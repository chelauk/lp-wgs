#!/usr/bin/env python3
"""Prepare ichorCNA SEG output for MEDICC2 total-copy-number input."""

import argparse
import csv
import math
import re
import sys
from pathlib import Path


COPY_COLUMN_CANDIDATES = ("Corrected_Copy_Number", "copy.number")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert patient-level ichorCNA SEG files to a MEDICC2 TSV."
    )
    parser.add_argument("--patient", required=True, help="Patient identifier.")
    parser.add_argument("--out", required=True, help="Output MEDICC2 TSV path.")
    parser.add_argument(
        "--report",
        default="medicc2_ichor_prep.txt",
        help="Validation report path.",
    )
    parser.add_argument(
        "--copy-column",
        default="auto",
        help=(
            "Copy-number column suffix or full column name. Use 'auto' to prefer "
            "Corrected_Copy_Number, then copy.number."
        ),
    )
    parser.add_argument(
        "--coordinate-system",
        choices=("one-based-inclusive", "bed"),
        default="one-based-inclusive",
        help="Coordinate convention used by the ichorCNA SEG input.",
    )
    parser.add_argument("seg_files", nargs="+", help="ichorCNA *.seg files.")
    return parser.parse_args()


def natural_key(chrom):
    value = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if value.isdigit():
        return (0, int(value))
    order = {"x": 23, "y": 24, "m": 25, "mt": 25}
    lower = value.lower()
    if lower in order:
        return (0, order[lower])
    return (1, lower)


def normalized_sample_id(path, header, copy_column):
    for suffix in COPY_COLUMN_CANDIDATES:
        suffix = "." + suffix
        if copy_column.endswith(suffix):
            return copy_column[: -len(suffix)]
    suffix = "." + copy_column
    for column in header:
        if column.endswith(suffix):
            return column[: -len(suffix)]
    stem = Path(path).name
    return re.sub(r"\.seg(?:\.txt)?$", "", stem)


def find_copy_column(header, requested):
    if requested != "auto":
        if requested in header:
            return requested
        matches = [column for column in header if column.endswith("." + requested)]
        if len(matches) == 1:
            return matches[0]
        if not matches:
            raise ValueError(f"copy-number column '{requested}' not found")
        raise ValueError(
            f"copy-number column suffix '{requested}' matched multiple columns: "
            + ", ".join(matches)
        )

    for suffix in COPY_COLUMN_CANDIDATES:
        matches = [column for column in header if column.endswith("." + suffix)]
        if len(matches) == 1:
            return matches[0]
        if suffix in header:
            return suffix
    raise ValueError(
        "could not find an ichorCNA copy-number column; tried suffixes: "
        + ", ".join(COPY_COLUMN_CANDIDATES)
    )


def parse_copy_number(value, path, line_number):
    try:
        parsed = float(value)
    except ValueError as error:
        raise ValueError(
            f"{path}:{line_number}: copy number '{value}' is not numeric"
        ) from error
    if not math.isfinite(parsed):
        raise ValueError(f"{path}:{line_number}: copy number '{value}' is not finite")
    rounded = round(parsed)
    if abs(parsed - rounded) > 1e-6:
        raise ValueError(
            f"{path}:{line_number}: copy number '{value}' is not an integer total CN"
        )
    if rounded < 0:
        raise ValueError(f"{path}:{line_number}: copy number '{value}' is negative")
    return int(rounded)


def read_seg(path, coordinate_system, requested_copy_column):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path}: missing header")
        required = {"chr", "start", "end"}
        missing = sorted(required - set(reader.fieldnames))
        if missing:
            raise ValueError(f"{path}: missing required columns: {', '.join(missing)}")

        copy_column = find_copy_column(reader.fieldnames, requested_copy_column)
        sample_id = normalized_sample_id(path, reader.fieldnames, copy_column)
        records = []
        seen = set()
        for line_number, row in enumerate(reader, start=2):
            chrom = row["chr"]
            try:
                start = int(row["start"])
                end = int(row["end"])
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: start/end must be integers"
                ) from error

            if coordinate_system == "one-based-inclusive":
                start -= 1

            if start < 0:
                raise ValueError(f"{path}:{line_number}: BED start is negative")
            if end <= start:
                raise ValueError(
                    f"{path}:{line_number}: end must be greater than start"
                )

            interval = (chrom, start, end)
            if interval in seen:
                raise ValueError(
                    f"{path}:{line_number}: duplicate interval {chrom}:{start}-{end}"
                )
            seen.add(interval)
            records.append((interval, parse_copy_number(row[copy_column], path, line_number)))

    records.sort(key=lambda item: (natural_key(item[0][0]), item[0][1], item[0][2]))
    return sample_id, copy_column, records


def main():
    args = parse_args()
    if len(args.seg_files) < 2:
        raise SystemExit("MEDICC2 preparation requires at least two ichorCNA SEG files")

    samples = []
    interval_template = None
    copy_columns = {}
    for seg_file in args.seg_files:
        sample_id, copy_column, records = read_seg(
            seg_file, args.coordinate_system, args.copy_column
        )
        if sample_id in {sample["sample_id"] for sample in samples}:
            raise SystemExit(f"duplicate sample_id '{sample_id}'")
        intervals = [interval for interval, _copy_number in records]
        if interval_template is None:
            interval_template = intervals
        elif intervals != interval_template:
            raise SystemExit(
                f"{seg_file}: intervals do not exactly match the first SEG file "
                "after coordinate normalization"
            )
        samples.append({"sample_id": sample_id, "records": records})
        copy_columns[sample_id] = copy_column

    nondiploid_samples = [
        sample["sample_id"]
        for sample in samples
        if any(copy_number != 2 for _interval, copy_number in sample["records"])
    ]
    if len(nondiploid_samples) < 2:
        raise SystemExit(
            "MEDICC2 requires at least two non-diploid samples; found "
            f"{len(nondiploid_samples)}"
        )

    with open(args.out, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "chrom", "start", "end", "Copies", "Diploid"])
        for sample in samples:
            for (chrom, start, end), copy_number in sample["records"]:
                writer.writerow([sample["sample_id"], chrom, start, end, copy_number, 2])

    with open(args.report, "w") as handle:
        handle.write(f"patient\t{args.patient}\n")
        handle.write(f"samples\t{len(samples)}\n")
        handle.write(f"segments_per_sample\t{len(interval_template or [])}\n")
        handle.write(f"coordinate_system_in\t{args.coordinate_system}\n")
        handle.write("coordinate_system_out\tbed\n")
        handle.write("copy_columns\t")
        handle.write(
            ",".join(
                f"{sample_id}:{copy_column}"
                for sample_id, copy_column in sorted(copy_columns.items())
            )
        )
        handle.write("\n")
        handle.write("nondiploid_samples\t" + ",".join(nondiploid_samples) + "\n")


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        sys.exit(1)
    except ValueError as error:
        raise SystemExit(str(error)) from error
