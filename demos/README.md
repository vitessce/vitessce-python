## Processing scripts for Vitessce demo datasets

Previously, we developed custom processing scripts for Vitessce demo data in the `vitessce-data` repository.

However, now that there are consensus single-cell file formats (e.g., h5ad) and an ecosystem of data processing packages in the community, we aim to leverage those. Writing demo dataset processing code using Scanpy and AnnData in small Snakemake workflows in this repository should allow us to iterate more quickly and share more demos.

### Setup

Set up the `vitessce-python-demos` environment using conda.

```sh
cd demos
conda env create -f environment.yml
conda activate vitessce-python-demos
pip install -e "..[dev]"
```

### Run

```sh
snakemake --cores all --rerun-triggers mtime
```

> NOTE: You will need to refresh the cellxgene urls as they are valid only for a week - please see snakemake files for instructions.

### Serve data locally

```sh
http-server --cors='*' --port 8000 .
```

### Add a new demo

To add a new demo, run

```sh
./create_demo.sh {new_demo_dirname}
```

Then,
- Add a subworkflow entry to the main demos workflow in `./Snakefile`
- Specify the output filenames in the `output` list in `./{new_demo_dir}/config.yml`
- Fill in the Snakefile in `./{new_demo_dir}/Snakefile`
- Write a test Vitessce configuration in `./{new_demo_dir}/vitessce.template.json`
- Test the demo by running the Vitessce frontend locally and navigating to `http://localhost:3000/?url=http://localhost:8000/{new_demo_dir}/vitessce.json`

See existing demo Snakefiles and scripts for examples.

Be sure to add comments (either in the demo-specific Snakefile or README.md file) about where raw file URLs were obtained.

Note that the name of `new_demo_dir` should match the intended key in `https://github.com/vitessce/vitessce/blob/main/src/demo/configs.js`, as this key can then be used in the documentation to point to the data processing scripts for each demo.


### Render all Vitessce config templates

To use the same configuration for both development (with localhost file URLs) and production/testing (with AWS/GCP file URLs), we write them as Jinja2 templates in the `vitessce.template.json` files where `{{ base_url }}` and `{{ base_url_gcp }}` are replaced at template render time.

We can render the template to `vitessce.local.json` and `vitessce.remote.json` using the following commands:

```sh
snakemake fill_templates -j 1
```

Be sure to re-run after making changes to `vitessce.template.json`.

#### Render an individual template

Local:

```sh
python fill_template.py -d codeluppi-2018 -t local > ./codeluppi-2018/vitessce.local.json
```

Remote:

```sh
python fill_template.py -d codeluppi-2018 -t remote -v 0.0.33 > ./codeluppi-2018/vitessce.remote.json
```

### Configure AWS and Google Cloud CLIs

Install `aws` CLI and add to your PATH ([reference](https://docs.aws.amazon.com/cli/latest/userguide/install-cliv2-linux.html)).

Install `gcloud` and `gsutil` and add to your PATH ([reference](https://cloud.google.com/storage/docs/gsutil_install#linux)).

Configure the AWS CLI by setting AWS environment variables ([reference](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-envvars.html)) or running `aws configure`  ([reference](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-files.html)).

Configure the Google Cloud CLI by running `gcloud auth login` ([reference](https://cloud.google.com/sdk/gcloud/reference/auth/login)).

### Run and deploy

```sh
snakemake --cores all --rerun-triggers mtime --config upload=true
```
