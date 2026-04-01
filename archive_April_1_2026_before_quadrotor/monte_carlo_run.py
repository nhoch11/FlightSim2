#!/usr/bin/env python3

# To use this file, simply edit the LAST 5 LINES with your filenames and information.
# Copy this file to your directory of choice, and run in command line within that directory by typing: python monte_carlo_run.py
# This generates a file that is ready for the monte_carlo_analysis.py
# PS - thank you ChatGPT for creating this script.

import csv
import subprocess
from pathlib import Path


def get_last_csv_row(csv_path: Path):
    """
    Read a CSV file and return:
        header, last_data_row

    Assumes the first row is a header.
    """
    with csv_path.open("r", newline="") as f:
        reader = csv.reader(f)
        rows = list(reader)

    if len(rows) < 2:
        raise ValueError(f"{csv_path} does not contain a header and at least one data row.")

    header = rows[0]
    last_row = rows[-1]
    return header, last_row


def run_monte_carlo(
    executable,
    input_file,
    states_file,
    output_file,
    n_runs,
):
    executable = Path(executable)
    input_file = Path(input_file)
    states_file = Path(states_file)
    output_file = Path(output_file)

    if not executable.exists():
        raise FileNotFoundError(f"Executable not found: {executable}")
    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    output_header_written = False

    if output_file.exists():
        output_file.unlink()

    for run_idx in range(1, n_runs + 1):
        print(f"Run {run_idx}/{n_runs}")

        # Run the external program
        result = subprocess.run(
            [str(executable), str(input_file)],
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            print("Program failed.")
            print("STDOUT:")
            print(result.stdout)
            print("STDERR:")
            print(result.stderr)
            raise RuntimeError(f"Run {run_idx} failed with return code {result.returncode}")

        if not states_file.exists():
            raise FileNotFoundError(f"Expected output file not found after run {run_idx}: {states_file}")

        # Read the last row from states.csv
        header, last_row = get_last_csv_row(states_file)

        # Write or append to output file
        with output_file.open("a", newline="") as f:
            writer = csv.writer(f)

            if not output_header_written:
                writer.writerow(["run"] + header)
                output_header_written = True

            writer.writerow([run_idx] + last_row)

    print(f"\nDone. Final-state results written to: {output_file}")


if __name__ == "__main__":
    base_dir = Path(__file__).resolve().parent
    run_monte_carlo(
        executable=base_dir /"main.exe",
        input_file=base_dir /"9.8.1.json",
        states_file=base_dir /"9.8.1_F16_states.csv",
        output_file=base_dir /"monte_carlo_final_states_severe.csv",
        n_runs=10000,
    )