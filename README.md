# analysis-pqc

Analysis for PQC measurements.

## creating tables from templates:

To enable one or more table output files, one must select the templates from ```scripts/templates``` and make a symlink to ```scripts/templates-enabled``` for each

the folder ```scripts/templates-enabled``` is in the gitignore and will not interfere with pulling new versions

## Install for full-line scripts

Clone the git repository

```bash
git clone https://github.com/hephy-dd/analysis-pqc.git -b <branch/tag>
cd analysis-pqc
```

Create a virtual environment

```bash
python3 -m venv env
. env/bin/activate # env/Scripts/activate.bat on Windows
python -m pip install -r requirements.txt
python -m pip install -r scripts/requirements.txt
python setup.py develop
```

Run full-line scripts

```bash
python3 scripts/full_line.py /home/measurements/PQC/Tracker/Production/Data/VPX35953/ -P -o ../test-pqc
```

## Install using pip

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
python scripts/full_line.py [-h] [-o DIR] [-l] [-H] [-P] path
```
```bash
required arguments:
  path        path to the folder of the batch - the "HPK_VPX12345_001_2-S_HM_WR" folders should be in this dir

optional arguments:
  -h, --help  show this help message and exit
  -o DIR      override output directory location, a sub-directory will be created at DIR/analysis_<batch-name>/
  -l          lazy evaluation: skip if the measurement folder is older than analysis folder
  -H          create histograms
  -P          create plots (for each single measurement used)
  -t EXPR     select templates to render (eg. -t*.tex -t*.html or -t* for all)
```

Templates that contain ```stdout``` will be sent to the stdout stream automatically, all others will be located in <outputdir>/tables/

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
