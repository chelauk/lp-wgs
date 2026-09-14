#!/usr/bin/env python3
"""Prepare ichorCNA SEG output for MEDICC2 total-copy-number input."""

import argparse
import csv
import math
import re
import sys
from pathlib import Path

SCRIPT_VERSION = "0.1.0"

COPY_COLUMN_CANDIDATES = ("Corrected_Copy_Number", "copy.number")
MISSING_COPY_NUMBER_VALUES = {"", ".", "NA", "N/A", "NaN", "nan", "NULL", "null"}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert patient-level ichorCNA SEG files to a MEDICC2 TSV."
    )
    parser.add_argument(
        "--version",
        action="version",
        version=SCRIPT_VERSION,
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


def find_fallback_copy_column(header, primary_column):
    corrected_suffix = ".Corrected_Copy_Number"
    if primary_column.endswith(corrected_suffix):
        sample_prefix = primary_column[: -len(corrected_suffix)]
        fallback_column = f"{sample_prefix}.copy.number"
        if fallback_column in header:
            return fallback_column
    if primary_column == "Corrected_Copy_Number" and "copy.number" in header:
        return "copy.number"
    return None


def copy_number_is_missing(value):
    return value is None or value.strip() in MISSING_COPY_NUMBER_VALUES


def parse_copy_number(value, path, line_number):
    value = value.strip()

    if value.upper() in {"NA", "NAN", "."}:
        return None

    try:
        parsed = float(value)
    except ValueError as error:
        raise ValueError(
            f"{path}:{line_number}: copy number '{value}' is not numeric"
        ) from error

    if not math.isfinite(parsed):
        return None

    rounded = round(parsed)
    if abs(parsed - rounded) > 1e-6:
        raise ValueError(
            f"{path}:{line_number}: copy number '{value}' is not an integer total CN"
        )
    if rounded < 0:
        raise ValueError(f"{path}:{line_number}: copy number '{value}' is negative")

    return int(rounded)


def read_seg(
    path,
    coordinate_system,
    requested_copy_column,
    allow_corrected_copy_number_fallback=False,
):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path}: missing header")
        required = {"chr", "start", "end"}
        missing = sorted(required - set(reader.fieldnames))
        if missing:
            raise ValueError(f"{path}: missing required columns: {', '.join(missing)}")

        copy_column = find_copy_column(reader.fieldnames, requested_copy_column)
        fallback_copy_column = None
        if allow_corrected_copy_number_fallback:
            fallback_copy_column = find_fallback_copy_column(
                reader.fieldnames, copy_column
            )
        sample_id = normalized_sample_id(path, reader.fieldnames, copy_column)
        records = []
        seen = set()
        used_fallback = False
        for line_number, row in enumerate(reader, start=2):
            chrom = normalize_autosome(row["chr"])
            if chrom is None:
                continue
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
            copy_number_value = row[copy_column]
            if (
                fallback_copy_column is not None
                and copy_number_is_missing(copy_number_value)
            ):
                copy_number_value = row[fallback_copy_column]
                used_fallback = True
            records.append(
                (interval, parse_copy_number(copy_number_value, path, line_number))
            )

    records.sort(key=lambda item: (natural_key(item[0][0]), item[0][1], item[0][2]))
    if used_fallback:
        copy_column = f"{copy_column};fallback={fallback_copy_column}"
    return sample_id, copy_column, records


def harmonize_samples(samples):
    """
    Retain genomic intervals present with non-missing CN in every sample.

    The ichorCNA inputs are expected to use the same fixed-bin coordinate grid.
    """
    interval_maps = []

    for sample in samples:
        interval_maps.append(
            {
                interval: copy_number
                for interval, copy_number in sample["records"]
            }
        )

    common_intervals = set(interval_maps[0])

    for interval_map in interval_maps[1:]:
        common_intervals &= set(interval_map)

    complete_intervals = [
        interval
        for interval in common_intervals
        if all(
            interval_map[interval] is not None
            for interval_map in interval_maps
        )
    ]

    complete_intervals.sort(
        key=lambda interval: (
            natural_key(interval[0]),
            interval[1],
            interval[2],
        )
    )

    if not complete_intervals:
        raise SystemExit(
            "no common intervals with complete copy-number calls remain"
        )

    total_intervals = len(
        set().union(*(set(interval_map) for interval_map in interval_maps))
    )
    removed_intervals = total_intervals - len(complete_intervals)

    harmonized_samples = []

    for sample, interval_map in zip(samples, interval_maps):
        harmonized_samples.append(
            {
                "sample_id": sample["sample_id"],
                "records": [
                    (interval, interval_map[interval])
                    for interval in complete_intervals
                ],
            }
        )

    return harmonized_samples, complete_intervals, removed_intervals


def compress_common_runs(samples, intervals):
    """Merge adjacent intervals with identical CN states across all samples."""
    copy_number_by_sample = [
        [copy_number for _interval, copy_number in sample["records"]]
        for sample in samples
    ]

    compressed_intervals = []
    compressed_copy_numbers = [[] for _sample in samples]
    current_chrom, current_start, current_end = intervals[0]
    current_state = tuple(
        sample_copy_numbers[0] for sample_copy_numbers in copy_number_by_sample
    )

    for index in range(1, len(intervals)):
        chrom, start, end = intervals[index]
        state = tuple(
            sample_copy_numbers[index] for sample_copy_numbers in copy_number_by_sample
        )
        if chrom == current_chrom and start == current_end and state == current_state:
            current_end = end
            continue

        compressed_intervals.append((current_chrom, current_start, current_end))
        for sample_index, copy_number in enumerate(current_state):
            compressed_copy_numbers[sample_index].append(copy_number)
        current_chrom, current_start, current_end = chrom, start, end
        current_state = state

    compressed_intervals.append((current_chrom, current_start, current_end))
    for sample_index, copy_number in enumerate(current_state):
        compressed_copy_numbers[sample_index].append(copy_number)

    compressed_samples = []
    for sample, sample_copy_numbers in zip(samples, compressed_copy_numbers):
        compressed_samples.append(
            {
                "sample_id": sample["sample_id"],
                "records": list(zip(compressed_intervals, sample_copy_numbers)),
            }
        )
    return compressed_samples, compressed_intervals


def main():
    args = parse_args()
    if len(args.seg_files) < 2:
        raise SystemExit("MEDICC2 preparation requires at least two ichorCNA SEG files")

    samples = []
    copy_columns = {}

    for seg_file in args.seg_files:
        sample_id, copy_column, records = read_seg(
            seg_file,
            args.coordinate_system,
            args.copy_column,
            allow_corrected_copy_number_fallback=True,
        )

        if sample_id in {sample["sample_id"] for sample in samples}:
            raise SystemExit(f"duplicate sample_id '{sample_id}'")

        samples.append(
            {
                "sample_id": sample_id,
                "records": records,
            }
        )
        copy_columns[sample_id] = copy_column

    if not samples:
        raise SystemExit("no ichorCNA SEG files were read")

    original_segment_count = sum(len(sample["records"]) for sample in samples)

    samples, interval_template, removed_interval_count = harmonize_samples(samples)

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

    samples, compressed_intervals = compress_common_runs(samples, interval_template)

    with open(args.out, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "chrom", "start", "end", "Copies", "Diploid"])
        for sample in samples:
            for (chrom, start, end), copy_number in sample["records"]:
                writer.writerow(
                    [sample["sample_id"], chrom, start, end, copy_number, 2]
                )

    with open(args.report, "w") as handle:
        handle.write(f"patient\t{args.patient}\n")
        handle.write(f"samples\t{len(samples)}\n")
        handle.write(f"segments_total_in\t{original_segment_count}\n")
        handle.write(f"harmonized_intervals\t{len(interval_template)}\n")
        handle.write(f"harmonized_intervals_removed\t{removed_interval_count}\n")
        handle.write(f"segments_per_sample_out\t{len(compressed_intervals)}\n")
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

