#!/usr/bin/env bash
set -o errexit
set -o pipefail

(( $# == 1 )) || die "Requires one argument, the demo name; Instead got $#."

SUBWORKFLOW_DIR=$1
SUBWORKFLOW_KEY=${SUBWORKFLOW_DIR//-/_}

mkdir "./$SUBWORKFLOW_DIR"
mkdir "./$SUBWORKFLOW_DIR/src"
touch "./$SUBWORKFLOW_DIR/Snakefile"
touch "./$SUBWORKFLOW_DIR/config.yml"
touch "./$SUBWORKFLOW_DIR/vitessce.template.json"

cat << EOF > "./$SUBWORKFLOW_DIR/Snakefile"
include: "../common.smk"
configfile: "config.yml"

rule all:
    input:
        [ (PROCESSED_DIR / f) for f in config['output'] ]
EOF

cat << EOF > "./$SUBWORKFLOW_DIR/config.yml"
output:
- $1.h5ad.zarr
EOF



cat << EOF
Add the following to the main demos Snakefile:

subworkflow $SUBWORKFLOW_KEY:
    workdir:
        "$SUBWORKFLOW_DIR"
EOF
