## Processing scripts for Vitessce demo datasets

Previously, we have developed processing scripts for Vitessce demos in the `vitessce-data` repository. 

However, now we are standardizing the way that we process data stored in common single-cell file formats (.h5ad, .loom, etc.) for the Python package in this repository. Therefore, writing demo dataset processing code in this repository may allow us to iterate more quickly (i.e. a monorepo, which comes with the usual monorepo tradeoffs).

### Setup

Set up the `vitessce-python-demos` environment using conda.

```sh
cd demos
conda env create -f environment.yml
conda activate vitessce-python-demos
pip install -e "..[testing]"
```

### Run

```sh
snakemake -j 1
```

### Serve data locally

```sh
http-server --cors='*' --port 8000 .
```

### Run and deploy

```sh
snakemake -j 1 --config upload=true
```

### Add a new demo

To add a new demo, run

```sh
./create_demo.sh {new_demo_dirname}
```

Then,
- Add a subworkflow entry to the main demos workflow in `./Snakefile`
- Specify the output filenames in the `output` list in `./{new_demo_dir}/config.yml`
- Fill in the Snakefile in `./{new_demo_dirname}/Snakefile`
- Write a test Vitessce configuration in `./{new_demo_dirname}/vitessce.json`

See existing demo Snakefiles and scripts for examples.
