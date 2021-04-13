# analysis-pqc

Analysis for PQC measurements.

## Install

Install using pip in a virtual environment.

```bash
pip install git+https://github.com/hephy-dd/analysis-pqc.git@0.2.0
```

## Scripts

Create a virtual environment for using the scripts.

```bash
python3 -m venv env
. env/bin/activate  # . env/Scripts/activate on Windows
```

Install the package and its additional dependencies required 
by the analysis scripts.

```bash
pip install .
pip install -r scripts/requirements.txt
```

Run command `deactivate` to exit the virtual environment and
`. env/bin/activate` to actiate it.

Run `pip install .` again to apply local changes in the
`analysis_pqc` pacakge.

Alternatively run `python setup.py develop` once to install
the package as local development source.

### Analyze JSON
```bash
python scripts/pqc_analysis_json.py <path> <analysis>
```

### Analyze text
```bash
python scripts/pqc_analysis_txt.py -f <path/to/files>
```

### Convert text to JSON
```bash
python scripts/txt2json.py <input.txt> -o <output.json>
```
