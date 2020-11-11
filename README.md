# analysis-pqc

Analysis for PQC measurements.

## Install

Install using pip in a virtual environment.

```bash
pip install git+https://github.com/hephy-dd/analysis-pqc.git@0.1.1
```

## Local usage and development

Set `PYTHONPATH` environment variable to use the local package (if not using an virtual environment).

```bash
export PYTHONPATH=$PYTHONPATH:$(pwd)
```

**Note:** the package requires additional dependencies to be installed.

```bash
pip install -r requirements.txt
```

## Run scripts

Analyze JSON
```bash
python scripts/pqc_analysis_json.py <path> <analysis>
```

Analyze text
```bash
python scripts/pqc_analysis_txt.py -f <path/to/files>
```

Convert text to JSON
```bash
python scripts/txt2json.py <input.txt> -o <output.json>
```

**Note:** the above scripts require additional dependencies to be installed.

```bash
pip install -r scripts/requirements.txt
```
