#!/usr/bin/env python3

# Copy this file to your directory of choice. 
# After running the monte_carlo_run.py script, run this script in command line within the same directory by typing: python monte_carlo_analysis.py
# This generates two output files with all the analysis results.
# PS - thank you ChatGPT for creating this script.

import csv
import math


def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def mean(values):
    return sum(values) / len(values)


def stdev(values):
    n = len(values)
    if n < 2:
        return 0.0
    mu = mean(values)
    return math.sqrt(sum((x - mu) ** 2 for x in values) / (n - 1))


def normal_pdf(x, mu, sigma):
    if sigma <= 0.0:
        return 0.0
    return (1.0 / (sigma * math.sqrt(2.0 * math.pi))) * math.exp(-0.5 * ((x - mu) / sigma) ** 2)


def read_csv_numeric_columns(csv_file):
    with open(csv_file, "r", newline="") as f:
        reader = csv.reader(f)
        rows = list(reader)

    if len(rows) < 2:
        raise ValueError("CSV file must contain a header and at least one data row.")

    headers = rows[0]
    data_rows = rows[1:]

    numeric_columns = {}
    numeric_headers = []

    for j, header in enumerate(headers):
        col_values = [row[j].strip() for row in data_rows if j < len(row)]
        nonblank = [v for v in col_values if v != ""]

        if len(nonblank) == 0:
            continue

        if all(is_float(v) for v in nonblank):
            numeric_headers.append(header)
            numeric_columns[header] = [float(v) for v in nonblank]

    return numeric_headers, numeric_columns


def make_histogram(values, n_bins=30, sigma_span=4.0):
    """
    Build histogram data and fitted normal PDF data.

    Returns:
        centers      : bin centers
        counts       : raw counts in each bin
        density      : normalized histogram density (area ~ 1)
        pdf_at_center: fitted normal PDF evaluated at bin centers
        bin_width    : bin width
    """
    n = len(values)
    mu = mean(values)
    sig = stdev(values)

    if n == 0:
        raise ValueError("Empty value list.")

    if sig == 0.0:
        centers = [mu]
        counts = [n]
        density = [1.0]
        pdf_at_center = [1.0]
        bin_width = 1.0
        return centers, counts, density, pdf_at_center, bin_width

    x_min = mu - sigma_span * sig
    x_max = mu + sigma_span * sig

    # Expand bounds if actual data exceed this range
    actual_min = min(values)
    actual_max = max(values)
    x_min = min(x_min, actual_min)
    x_max = max(x_max, actual_max)

    if x_max == x_min:
        x_max = x_min + 1.0

    bin_width = (x_max - x_min) / n_bins

    edges = [x_min + i * bin_width for i in range(n_bins + 1)]
    counts = [0] * n_bins

    for v in values:
        if v == x_max:
            idx = n_bins - 1
        else:
            idx = int((v - x_min) / bin_width)
            if idx < 0:
                idx = 0
            elif idx >= n_bins:
                idx = n_bins - 1
        counts[idx] += 1

    centers = [0.5 * (edges[i] + edges[i + 1]) for i in range(n_bins)]
    density = [c / (n * bin_width) for c in counts]
    pdf_at_center = [normal_pdf(x, mu, sig) for x in centers]

    return centers, counts, density, pdf_at_center, bin_width


def write_histogram_csv(output_file, headers, columns, n_bins=30, sigma_span=4.0):
    """
    Writes one wide CSV with histogram and fitted normal data for each column.

    For each variable h, output columns are:
        bin_center_h, count_h, density_h, normal_pdf_h
    """
    all_histograms = []
    max_len = 0

    for h in headers:
        centers, counts, density, pdf_at_center, bin_width = make_histogram(
            columns[h], n_bins=n_bins, sigma_span=sigma_span
        )
        all_histograms.append((h, centers, counts, density, pdf_at_center))
        max_len = max(max_len, len(centers))

    out_headers = []
    for h, _, _, _, _ in all_histograms:
        out_headers.extend([
            f"bin_center_{h}",
            f"count_{h}",
            f"density_{h}",
            f"normal_pdf_{h}",
        ])

    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(out_headers)

        for i in range(max_len):
            row = []
            for _, centers, counts, density, pdf_at_center in all_histograms:
                if i < len(centers):
                    row.extend([centers[i], counts[i], density[i], pdf_at_center[i]])
                else:
                    row.extend(["", "", "", ""])
            writer.writerow(row)


def write_summary_csv(output_file, headers, columns):
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["column", "mean", "stdev", "min", "max", "n"])

        for h in headers:
            vals = columns[h]
            writer.writerow([h, mean(vals), stdev(vals), min(vals), max(vals), len(vals)])


def main():
    input_file = "monte_carlo_final_states_severe.csv"
    summary_file = "monte_carlo_summary_severe.csv"
    histogram_file = "monte_carlo_histograms_severe.csv"

    headers, columns = read_csv_numeric_columns(input_file)

    # Skip run column if present
    headers = [h for h in headers if h.lower() != "run"]
    columns = {h: columns[h] for h in headers}

    write_summary_csv(summary_file, headers, columns)
    write_histogram_csv(histogram_file, headers, columns, n_bins=70, sigma_span=4.0)

    print(f"Wrote summary statistics to: {summary_file}")
    print(f"Wrote histogram data to: {histogram_file}")


if __name__ == "__main__":
    main()