#!/usr/bin/env python3
import sys
from mpmath import mp, mpf

# Set precision (e.g., 128 decimal digits)
mp.dps = 128

def read_values(path):
    with open(path, 'r') as f:
        line = f.readline()
        values = list(map(mpf, line.strip().split()))
        if len(values) != 4:
            raise ValueError(f"File {path} does not contain exactly four values.")
        return values

def compare_files(file1, file2, tol_str="1e-9"):
    tol = mpf(tol_str)
    values1 = read_values(file1)
    values2 = read_values(file2)

    for i, (v1, v2) in enumerate(zip(values1, values2), start=1):
        diff = abs(v1 - v2)/v1
        if diff > tol:
            print(f"Mismatch at position {i}:")
            print(f"  File 1: {v1}")
            print(f"  File 2: {v2}")
            print(f"  |Difference| = {diff} > Tolerance = {tol}")
            return 1
    print("Files match within tolerance.")
    return 0

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: compare_txt_files.py <file1> <file2>")
        sys.exit(2)

    file1, file2 = sys.argv[1], sys.argv[2]
    sys.exit(compare_files(file1, file2))
