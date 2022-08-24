import argparse
import pandas as pd

from constants import COLUMNS


def split_annotation_csv(input_file, output_file, globus_id):
    df = pd.read_csv(input_file)
    df = df.loc[df[COLUMNS.GLOBUS_ID.value] == globus_id]
    df = df[[
        COLUMNS.CELL_ID.value,
        COLUMNS.ANNOTATION.value,
        COLUMNS.PREDICTION_SCORE.value
    ]]
    df.to_csv(output_file, index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i',
        '--input_file',
        type=str,
        required=True,
        help='Input CSV file'
    )
    parser.add_argument(
        '-o',
        '--output_file',
        type=str,
        required=True,
        help='Output CSV file'
    )
    parser.add_argument(
        '-gid',
        '--globus_id',
        type=str,
        required=True,
        help='Globus ID to use for filtering'
    )
    args = parser.parse_args()
    split_annotation_csv(
        args.input_file,
        args.output_file,
        args.globus_id
    )
