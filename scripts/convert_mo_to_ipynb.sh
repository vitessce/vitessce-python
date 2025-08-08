cd docs/notebooks

for file in *.mo.py; do
    echo "Converting $file to ${file%.mo.py}.ipynb"
    uv run marimo export ipynb —-sort=top-down --force --output "__ipynb__/${file%.mo.py}.ipynb" "$file"
done

cp example_configs.py "__ipynb__/example_configs.py"