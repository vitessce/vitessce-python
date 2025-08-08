cd docs/notebooks

for file in *.mo.py; do
    echo "Converting $file to ${file%.mo.py}.ipynb"
    marimo export ipynb "$file" --output "__ipynb__/${file%.mo.py}.ipynb" --force
done

cp example_configs.py "__ipynb__/example_configs.py"