#!/usr/bin/env python3

import argparse
import os.path
import glob
import json


def dir_path(path):
    if os.path.isdir(path):
        return os.path.realpath(path)
    else:
        raise NotADirectoryError(path)


def main():
    # parse argument
    parser = argparse.ArgumentParser(
        description="Create inputs.json file for preprocessing."
    )
    parser.add_argument(
        "path", type=dir_path, help="sample or library directory"
    )
    args = parser.parse_args()
    path = os.path.realpath(args.path)
    library_id = os.path.basename(path)

    # list lanes
    lanes = glob.glob(path + "/**/*.fastq.gz", recursive=True)
    lanes.sort()

    # group lanes by sample
    sample_ids = []
    sample_dict = {}
    for lane in lanes:
        sample = os.path.basename(os.path.dirname(lane))
        if sample in sample_dict:
            sample_dict[sample].append(lane)
        else:
            sample_ids.append(sample)
            sample_dict[sample] = [lane]

    # write output in JSON format
    with open("inputs.json", "w") as output_file:
        json.dump(
            {
                "Main.library_id": library_id,
                "Main.sample_ids": sample_ids,
                "Main.samples": sample_dict,
            },
            output_file,
            indent=2,
        )

    # print output
    print("samples:")
    for i, sample in enumerate(sample_ids, start=1):
        print(i, sample)


if __name__ == "__main__":
    main()
