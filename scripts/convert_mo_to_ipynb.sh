cd docs/notebooks

for file in *.mo.py; do
    echo "Converting $file to ${file%.mo.py}.ipynb"
    uv run marimo export ipynb "$file" --output "__ipynb__/${file%.mo.py}.ipynb" —-sort top-down --force
done

# Uncomment the following lines if you want to use the list of existing .ipynb files in __ipynb__ directory for conversion
#for file in __ipynb__/*.ipynb; do
#    FILENAME=$(basename "${file%.ipynb}")
#    echo "Converting ${FILENAME}.mo.py to ${FILENAME}.ipynb"
#    uv run marimo export ipynb "${FILENAME}.mo.py" --output "__ipynb__/${FILENAME}.ipynb" --sort top-down --force
#done

cp example_configs.py "__ipynb__/example_configs.py"