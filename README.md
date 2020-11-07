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

## Run example

```bash
python pqc_analysis_scripts/pqc_analysis_json.py [path] [analysis]
python pqc_analysis_scripts/pqc_analysis_txt.py -f path/to/files
```
Note: the above example requires package `matplotlib` to be installed.
