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
python examples/analysis_example.py -f path/to/files
```
Note: the above example requires package `matplotlib` to be installed.
