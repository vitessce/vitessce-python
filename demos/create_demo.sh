#!/usr/bin/env bash
set -o errexit
set -o pipefail

(( $# == 1 )) || die "Requires one argument, the demo name; Instead got $#."

mkdir "./$1"
mkdir "./$1/src"
touch "./$1/Snakefile"
touch "./$1/config.yml"
touch "./$1/vitessce.json"

cat << EOF > "./$1/Snakefile"
include: "../common.smk"
configfile: "config.yml"

rule all:
    input:
        [ (PROCESSED_DIR / f) for f in config['output'] ]
EOF

cat << EOF > "./$1/config.yml"
output:
- $1.h5ad.zarr
EOF

cat << EOF
Add the following to the main demos Snakefile:

subworkflow $1:
    workdir:
        "$1"
EOF