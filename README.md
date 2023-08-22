# analysis-pqc

Analysis functions for PQC measurements.

## Install using pip

Install using pip (>= 19.0) in a virtual environment.

```bash
pip install git+https://github.com/hephy-dd/analysis-pqc.git@0.8.1
```

**Note:** this does not install the bundled analysis scripts.

## Install analysis scripts

Clone the git repository

```bash
git clone https://github.com/hephy-dd/analysis-pqc.git -b <branch/tag>
cd analysis-pqc
```

Create a [virtual environment](https://docs.python.org/3/tutorial/venv.html).

```bash
python3 -m venv env
. env/bin/activate  # env/Scripts/activate.bat on Windows
pip install --upgrade pip
pip install --editable .
pip install -r scripts/requirements.txt
```

## Run full-line script

```bash
python scripts/full_line.py /PQC/Tracker/Production/Data/VPX35953/ -HP -t*.html -o ../test-pqc
```

Use flag `-t` to select templates to render output files by specifying a `glob`
expression with wildcard support. The flag can be set multiple times,
eg. `-t*.html` for all HTML reports, `-tresults.html` for a specific HTML
report, `-t*.tex -t*.html` for all LaTeX and HTML reports, or `-t*` for all.

XML files are generated with the -t*.xml flag. XML templates are located in scripts/templates.

For custom configurations clone and edit file `scripts/config/default.yaml` and
use flag `-c <name>` to load custom configuration from
`scripts/config/<name>.yaml`.

Run command `deactivate` to exit the virtual environment and
`. env/bin/activate` again to activate it.

Full line analysis (for whole batch)

```bash
python scripts/full_line.py [-h] [-o DIR] [-l] [-H] [-P] [-t EXPR] [-c NAME] path
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
  -c NAME     load custom configuration by name
```

Templates that contain ```stdout``` will be sent to the stdout stream automatically, all others will be located in DIR/analysis_<batch-name>/

## Other scripts

### Upload XML files to the database
Upload all .xml files in <path> (including subdirectories) to the database: `production` (cmsr) or `development` (int2r)
```bash
python3 scripts/upload_db.py -db <database> <path>
```

### Analyze JSON

```bash
python scripts/pqc_analysis_json.py <path> <analysis>
```

### Convert text to JSON
```bash
python scripts/txt2json.py <input.txt> -o <output.json>
```

**Note:** the scripts requires additional dependencies to be installed.

```bash
pip install -r scripts/requirements.txt
```
