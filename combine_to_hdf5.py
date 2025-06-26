import os
import h5py

def collect_txt_to_hdf5(root_dir, output_hdf5_path):
    with h5py.File(output_hdf5_path, 'w') as h5file:
        for dirpath, _, filenames in os.walk(root_dir):
            for filename in filenames:
                if filename == "Z.txt":
                    full_path = os.path.join(dirpath, filename)

                    relative_path = os.path.relpath(full_path, root_dir)
                    parts = relative_path.split(os.sep)

                    # We expect: resultsGaussian/<prob>/<stddev>/<x>/<y>/<T_frac>/<prec>/<seed>/Z.txt
                    if len(parts) < 8 or parts[0] != "resultsGaussian":
                        print(f"Skipping unexpected path: {relative_path}")
                        continue

                    prob, stddev, x, y, T_frac, prec, seed = parts[1:8]
                    dataset_name = "Z"
                    group_path = f"{T_frac}/{prob}/{stddev}/{x}/{y}/{prec}/{seed}"

                    # Read Z.txt content (expecting 4 space-separated high-precision numbers)
                    with open(full_path, 'r') as f:
                        content = f.read().strip()
                        values = content.split()
                        if len(values) != 4:
                            print(f"Warning: {full_path} does not contain 4 values")
                            continue

                    # Store as strings for precision reason
                    dt = h5py.string_dtype(encoding='utf-8')
                    group = h5file.require_group(group_path)
                    group.create_dataset(dataset_name, data=values, dtype=dt)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Combine result txt files into HDF5.")
    parser.add_argument("result_dir", help="Root directory containing resultsGaussian")
    parser.add_argument("output_hdf5", help="Output HDF5 filename")
    args = parser.parse_args()

    collect_txt_to_hdf5(args.result_dir, args.output_hdf5)
