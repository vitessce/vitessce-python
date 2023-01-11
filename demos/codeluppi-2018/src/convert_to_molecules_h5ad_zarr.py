import argparse
import pandas as pd
import h5py
from anndata import AnnData
from vitessce.data_utils import (
    optimize_adata,
)


def convert_to_csv(args):
    with h5py.File(args.input, "r") as f:
        molecule_types = list(f.keys())

        molecules_df = pd.DataFrame(index=[], columns=["X", "Y", "Gene"])
        start_at = 0
        for molecule_type in molecule_types:
            molecule_arr = f[molecule_type][()]
            molecule_type_ids = [
                str(i) for i in range(start_at, start_at + len(molecule_arr))
            ]
            start_at += len(molecule_arr)
            molecule_type_df = pd.DataFrame(
                index=molecule_type_ids, columns=molecules_df.columns.values.tolist()
            )
            molecule_type_df["Gene"] = molecule_type
            molecule_type_df["X"] = [m[0] for m in molecule_arr]
            molecule_type_df["Y"] = [m[1] for m in molecule_arr]

            molecules_df = pd.concat([molecules_df, molecule_type_df])

    molecules_df.index = molecules_df.index.rename("molecule_id")

    # Molecules zarr
    m_obs_df = molecules_df[["Gene"]]
    m_var_df = pd.DataFrame(index=[], columns=[], data=[])

    m_adata = AnnData(X=None, obs=m_obs_df, var=m_var_df)
    m_adata.obsm["X_spatial"] = molecules_df[["X", "Y"]].values
    m_adata = optimize_adata(m_adata)

    # Write outputs
    m_adata.write_zarr(args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", type=str, required=True, help="Input hdf5 file"
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True, help="Output AnnData Zarr store"
    )
    args = parser.parse_args()
    convert_to_csv(args)
