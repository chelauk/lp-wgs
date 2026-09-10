#!/usr/bin/env python3
"""Collect ichorCNA SEG files across patients outside Nextflow."""

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

from prep_medicc_ichor import (
    compress_common_runs,
    read_seg,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Collect ichorCNA SEG files from files or directories and write a "
            "single composite segment table across patients."
        )
    )
    parser.add_argument(
        "inputs",
        nargs="*",
        help="SEG files and/or directories to scan. Not required with --manifest.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Directory for the composite SEG table and reports.",
    )
    parser.add_argument(
        "--composite-out",
        default="composite_ichor_segments.tsv",
        help="Composite SEG filename inside --outdir. Default: %(default)s",
    )
    parser.add_argument(
        "--pattern",
        default="*.seg*",
        help="Glob used when scanning directories. Default: %(default)s",
    )
    parser.add_argument(
        "--patient-regex",
        help=(
            "Regex applied to each SEG path. Must include a named (?P<patient>...) "
            "capture group. Overrides automatic patient inference."
        ),
    )
    parser.add_argument(
        "--copy-column",
        default="auto",
        help=(
            "Copy-number column suffix or full column name. Use 'auto' to prefer "
            "Corrected_Copy_Number, then copy.number. When Corrected_Copy_Number "
            "is selected and a row is missing, the matching copy.number value is "
            "used for that row."
        ),
    )
    parser.add_argument(
        "--coordinate-system",
        choices=("one-based-inclusive", "bed"),
        default="one-based-inclusive",
        help="Coordinate convention used by the ichorCNA SEG input.",
    )
    parser.add_argument(
        "--min-samples",
        type=int,
        default=2,
        help=(
            "Minimum SEG files required per patient when --write-medicc-per-patient "
            "is set. Default: %(default)s"
        ),
    )
    parser.add_argument(
        "--min-nondiploid-samples",
        type=int,
        default=2,
        help=(
            "Minimum non-diploid samples required per patient when "
            "--write-medicc-per-patient is set. Default: %(default)s"
        ),
    )
    parser.add_argument(
        "--skip-invalid-patients",
        action="store_true",
        help=(
            "Continue when a patient fails validation. By default the script exits "
            "on the first invalid patient."
        ),
    )
    parser.add_argument(
        "--manifest",
        help=(
            "Optional TSV with columns patient and seg_file. When provided, inputs "
            "are ignored for SEG discovery."
        ),
    )
    parser.add_argument(
        "--summary",
        default="composite_ichor_segments_summary.tsv",
        help="Batch summary filename inside --outdir. Default: %(default)s",
    )
    parser.add_argument(
        "--write-medicc-per-patient",
        action="store_true",
        help="Also write one MEDICC2 TSV and prep report per patient.",
    )
    return parser.parse_args()


def discover_seg_files(inputs, pattern):
    seg_files = []
    for raw_input in inputs:
        input_path = Path(raw_input)
        if input_path.is_dir():
            seg_files.extend(path for path in input_path.rglob(pattern) if path.is_file())
        elif input_path.is_file():
            seg_files.append(input_path)
        else:
            raise ValueError(f"{input_path}: input path does not exist")
    return sorted(set(path.resolve() for path in seg_files))


def read_manifest(path):
    grouped = defaultdict(list)
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path}: missing header")
        required = {"patient", "seg_file"}
        missing = sorted(required - set(reader.fieldnames))
        if missing:
            raise ValueError(f"{path}: missing required columns: {', '.join(missing)}")
        for line_number, row in enumerate(reader, start=2):
            patient = row["patient"].strip()
            seg_file = row["seg_file"].strip()
            if not patient:
                raise ValueError(f"{path}:{line_number}: patient is empty")
            if not seg_file:
                raise ValueError(f"{path}:{line_number}: seg_file is empty")
            seg_path = Path(seg_file).expanduser()
            if not seg_path.is_absolute():
                seg_path = Path(path).resolve().parent / seg_path
            if not seg_path.is_file():
                raise ValueError(f"{path}:{line_number}: {seg_path} does not exist")
            grouped[patient].append(seg_path.resolve())
    return grouped


def infer_patient_from_path(path):
    parts = path.parts
    if "low_pass_wgs" in parts:
        index = parts.index("low_pass_wgs")
        if index >= 2:
            return parts[index - 2]
    if len(path.parents) >= 2:
        return path.parents[1].name
    return path.parent.name


def patient_from_regex(path, regex):
    match = regex.search(str(path))
    if match is None:
        raise ValueError(f"{path}: does not match --patient-regex")
    try:
        patient = match.group("patient")
    except IndexError as error:
        raise ValueError("--patient-regex must define a named 'patient' group") from error
    if not patient:
        raise ValueError(f"{path}: --patient-regex produced an empty patient ID")
    return patient


def group_by_patient(seg_files, patient_regex):
    grouped = defaultdict(list)
    regex = re.compile(patient_regex) if patient_regex else None
    for seg_file in seg_files:
        patient = (
            patient_from_regex(seg_file, regex)
            if regex is not None
            else infer_patient_from_path(seg_file)
        )
        grouped[patient].append(seg_file)
    return grouped


def read_sample(seg_file, args):
    sample_id, copy_column, records = read_seg(
        seg_file,
        args.coordinate_system,
        args.copy_column,
        allow_corrected_copy_number_fallback=True,
    )
    return {
        "sample_id": sample_id,
        "copy_column": copy_column,
        "records": records,
        "seg_file": seg_file,
    }


def read_patient_samples(patient, seg_files, args):
    samples = []
    sample_ids = set()
    for seg_file in sorted(seg_files):
        sample = read_sample(seg_file, args)
        if sample["sample_id"] in sample_ids:
            raise ValueError(f"{patient}: duplicate sample_id '{sample['sample_id']}'")
        sample_ids.add(sample["sample_id"])
        samples.append(sample)
    return samples


def write_composite(grouped, outdir, args):
    output_path = outdir / args.composite_out
    summary_rows = []
    total_rows = 0

    with open(output_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "patient",
                "sample_id",
                "chrom",
                "start",
                "end",
                "copy_number",
                "diploid",
                "copy_column",
                "seg_file",
            ]
        )
        for patient, seg_files in sorted(grouped.items()):
            try:
                samples = read_patient_samples(patient, seg_files, args)
            except ValueError as error:
                if not args.skip_invalid_patients:
                    raise
                summary_rows.append(
                    {
                        "patient": patient,
                        "status": "skipped",
                        "samples": "",
                        "segments": "",
                        "composite_rows": "",
                        "medicc_tsv": "",
                        "report": "",
                        "message": str(error),
                    }
                )
                continue

            patient_rows = 0
            segment_counts = []
            for sample in samples:
                segment_counts.append(len(sample["records"]))
                for (chrom, start, end), copy_number in sample["records"]:
                    writer.writerow(
                        [
                            patient,
                            sample["sample_id"],
                            chrom,
                            start,
                            end,
                            copy_number,
                            2,
                            sample["copy_column"],
                            sample["seg_file"],
                        ]
                    )
                    patient_rows += 1

            total_rows += patient_rows
            summary_rows.append(
                {
                    "patient": patient,
                    "status": "written",
                    "samples": len(samples),
                    "segments": ",".join(str(count) for count in segment_counts),
                    "composite_rows": patient_rows,
                    "medicc_tsv": "",
                    "report": "",
                    "message": "",
                }
            )

    return output_path, summary_rows, total_rows


def write_patient(patient, seg_files, outdir, args):
    if len(seg_files) < args.min_samples:
        raise ValueError(
            f"{patient}: found {len(seg_files)} SEG file(s), need at least {args.min_samples}"
        )

    input_samples = read_patient_samples(patient, seg_files, args)
    samples = []
    interval_template = None
    copy_columns = {}
    for input_sample in input_samples:
        sample_id = input_sample["sample_id"]
        copy_column = input_sample["copy_column"]
        records = input_sample["records"]
        intervals = [interval for interval, _copy_number in records]
        if interval_template is None:
            interval_template = intervals
        elif intervals != interval_template:
            raise ValueError(
                f"{patient}: {input_sample['seg_file']}: intervals do not exactly match the first "
                "SEG file after coordinate normalization"
            )
        samples.append({"sample_id": sample_id, "records": records})
        copy_columns[sample_id] = copy_column

    nondiploid_samples = [
        sample["sample_id"]
        for sample in samples
        if any(copy_number != 2 for _interval, copy_number in sample["records"])
    ]
    if len(nondiploid_samples) < args.min_nondiploid_samples:
        raise ValueError(
            f"{patient}: found {len(nondiploid_samples)} non-diploid sample(s), "
            f"need at least {args.min_nondiploid_samples}"
        )

    original_segment_count = len(interval_template or [])
    samples, compressed_intervals = compress_common_runs(samples, interval_template)

    outdir.mkdir(parents=True, exist_ok=True)
    output_tsv = outdir / f"{patient}.tsv"
    report = outdir / f"{patient}.medicc2_ichor_prep.txt"

    with open(output_tsv, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "chrom", "start", "end", "Copies", "Diploid"])
        for sample in samples:
            for (chrom, start, end), copy_number in sample["records"]:
                writer.writerow([sample["sample_id"], chrom, start, end, copy_number, 2])

    with open(report, "w") as handle:
        handle.write(f"patient\t{patient}\n")
        handle.write(f"samples\t{len(samples)}\n")
        handle.write(f"segments_per_sample_in\t{original_segment_count}\n")
        handle.write(f"segments_per_sample_out\t{len(compressed_intervals)}\n")
        handle.write(f"coordinate_system_in\t{args.coordinate_system}\n")
        handle.write("coordinate_system_out\tbed\n")
        handle.write("seg_files\t")
        handle.write(",".join(str(path) for path in sorted(seg_files)))
        handle.write("\n")
        handle.write("copy_columns\t")
        handle.write(
            ",".join(
                f"{sample_id}:{copy_column}"
                for sample_id, copy_column in sorted(copy_columns.items())
            )
        )
        handle.write("\n")
        handle.write("nondiploid_samples\t" + ",".join(nondiploid_samples) + "\n")

    return {
        "patient": patient,
        "samples": len(samples),
        "segments_per_sample_in": original_segment_count,
        "segments_per_sample_out": len(compressed_intervals),
        "output_tsv": output_tsv,
        "report": report,
    }


def write_summary(path, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            lineterminator="\n",
            fieldnames=[
                "patient",
                "status",
                "samples",
                "segments",
                "composite_rows",
                "medicc_tsv",
                "report",
                "message",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    outdir = Path(args.outdir)

    if args.manifest:
        grouped = read_manifest(args.manifest)
    else:
        if not args.inputs:
            raise SystemExit("provide SEG files/directories or use --manifest")
        seg_files = discover_seg_files(args.inputs, args.pattern)
        if not seg_files:
            raise SystemExit("no SEG files found")
        grouped = group_by_patient(seg_files, args.patient_regex)

    outdir.mkdir(parents=True, exist_ok=True)
    composite_path, summary_rows, total_rows = write_composite(grouped, outdir, args)
    written = []

    if total_rows == 0:
        raise SystemExit("no composite segment rows were written")

    if args.write_medicc_per_patient:
        summary_by_patient = {row["patient"]: row for row in summary_rows}
        for patient, seg_files in sorted(grouped.items()):
            if summary_by_patient.get(patient, {}).get("status") == "skipped":
                continue
            try:
                result = write_patient(patient, seg_files, outdir, args)
                written.append(result)
                summary_by_patient[patient]["medicc_tsv"] = result["output_tsv"]
                summary_by_patient[patient]["report"] = result["report"]
            except ValueError as error:
                if not args.skip_invalid_patients:
                    raise
                summary_by_patient[patient]["status"] = "skipped"
                summary_by_patient[patient]["message"] = str(error)

    summary_path = outdir / args.summary
    write_summary(summary_path, summary_rows)

    print(f"wrote\t{composite_path}")
    for result in written:
        print(f"wrote\t{result['output_tsv']}\t{result['report']}")

    for row in summary_rows:
        if row["status"] == "skipped":
            print(f"skipped\t{row['message']}", file=sys.stderr)
    print(f"summary\t{summary_path}")


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        sys.exit(1)
    except ValueError as error:
        raise SystemExit(str(error)) from error
