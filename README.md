# analysis-pqc

Analysis for PQC measurements.

## creating tables from templates:

To enable one or more table output files, one must select the templates from scripts/templates and make a symlink to scripts/templates-enabled for each

templates that contain ```stdout``` will be sent to the stdout stream automatically, all others will be located in <outputdir>/tables/ 





## Install

Install using pip in a virtual environment.

```bash
not sure if this works...
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

Full line analysis (for whole batch)
```bash
python scripts/full_line.py [-h] [-o O] [-l] [-H] [-P] path
```

optional arguments:
  -h, --help  show this help message and exit
  -o PATH     override output directory location, a subdirectory will be created there: PATH/analysis_<batch-name>/
  -l          lazy evaluation: skip if the measurement folder is older than analysis folder
  -H          create histograms
  -P          create plots (for each single measurement used)


## Old scripts:

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
