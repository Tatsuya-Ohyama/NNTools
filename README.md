# NNTools

## Abstract
Group of programs to calculate nearest-neighbor parameters for nucleic acids

* `NNcalcTool.py`: Program to calculate various energies of nucleic acids from parameters
* `NNfitTool.py`: Program for fitting parameters from experimental values
* `NNcountTool.py`: Program to count NN pairs



## Usage
### NNcalcTool.py

```sh
$ NNcalcTool.py [-h] -s SEQUENCE_FILE.csv -p PARAMETER_FILE.csv -o OUTPUT_FILE.csv [-l LOG_FILE.log] [-O] [--make-template]
```

* `-h`, `--help`
	: show this help message and exit
* `-s SEQUENCE_FILE.csv`
	: sequence file of csv format

	* column: `Label` and `Sequence`
		* `Label`: experimental name, label, note and description.
		* `Sequence`: Sequence

* `-p PARAMETER_FILE.csv`
	: input file for parameters

	* column: `Parameter`, `dH`, `dS`, and `dG`
		* `Parameter`: `AA/TT`, `GC/CG`, or etc, and `re:` (regexp) and `reg:` (counted pattern by regexp)
		* `dH`: $\Delta H$ parameter
		* `dS`: $\Delta S$ parameter
		* `dG`: $\Delta G$ parameter

* `-o OUTPUT_FILE.csv`
	: output file
* `-l LOG_FILE.log`
	: log file (if not specify, a log file with the same name as -o option is generated)
* `-O`
	: overwrite forcibly
* `--make-template`
	: make template files (`template_ref_param.csv` and `template_sequence.csv`) and exit

	```sh
	$ NNcalcTool.py --make-template
	```

### NNfitTool.py

```sh
$ NNfitTool.py [-h] -x EXP.csv -r REF_PARAM.csv -o OUTPUT.csv [-O] [-d THRESHOLD] [-i INITIAL_INCREMENT] [-T TEMPERATURE] [-m EVALUATION_METHOD] [-S] [-I LOOP_COUNT] [--one-direction] [-e] [-t THREAD] [--verbose]
```

* Common
	* `-h`, `--help`
		: show this help message and exit

* Input:
	* `-x EXP.csv`
		: sequence and experimental value file

		* column: `Label`, `Sequence`, `dH`, `dH(error)`, `dS`, `dS(error)`, `dG`, and `dG(error)`
			* `Label`: experimental name, label, note and description.
			* `Sequence`: Sequence
			* `dH`: $\Delta H$ measured in the experiment.
			* `dH(error)`: error of $\Delta H$ measured in the experiment.
			* `dS`: $\Delta S$ measured in the experiment.
			* `dS(error)`: error of $\Delta S$ measured in the experiment.
			* `dG`: $\Delta G$ measured in the experiment.
			* `dG(error)`: error of $\Delta G$ measured in the experiment.

	* `-r REF_PARAM.csv`
		: referenced parameter values

		* column: `Parameter`, `dH`, `dS`, `dG`, space, `dH (change)`, `dS (change)`, `dG (change)`, space, `dH (Direction)`, `dS (Direction)`, and `dG (Direction)`
			* `Parameter`: `AA/TT`, `GC/CG`, or etc, and `re:` (regexp) and `reg:` (counted pattern by regexp)
				* e.g., `re:^A/^T` and `re:^T/^A` (initial parameter for A/T (both specification require))
				* e.g., `reg:.*?G.*?/.*?C.*?` (number of G/C pair parameter)
				* e.g., `reg:./.` (length parameter)
			* `dH`: initial $\Delta H$
			* `dS`: initial $\Delta S$
			* `dG`: initial $\Delta G$
			* `dH(change)` / `dS(change)` / `dG(change)`: whether the parameter changes is allowed during parameter fitting. (Default: `TRUE` (allow))
			* `dH (Direction)` / `dS (Direction)` / `dG (Direction)`: direction of the parameter change during parameter fitting (0: change to positive or negative / 1: change only to positive / 2: change only to negative)

* Output:
	* `-o OUTPUT.csv`
		: output file
	* `-O`
		: overwrite forcibly

* Config:
	* `-d THRESHOLD`
		: difference threshold of increment for searching (Default: 0.01)
	* `-i INITIAL_INCREMENT`
		: initial increment (Default: 1.0)
	* `-T TEMPERATURE`
		: temperature for experimental data (Default: 310.15)
	* `-m EVALUATION_METHOD`
		: evaluation method (r, r2, diff_mean, diff_std, diff_sum, diff_square) (Default: diff_square)
	* `-S`
		: Separately calculate dS (Default: OFF (dS is calculated by Gibbs free energy equation))
	* `-I LOOP_COUNT`
		: the number of looping optimize (Default: 1)
	* `-e`
		: consider with experimental value with error
	* `--one-direction`
		: Do not search for reverse order pattern (For example, this program searches AC/TG and reverse order pattern GT/CA as the same pattern. This option does not allow it.) (For mismatch or manual operation)

* Misc:
	* `-t THREAD`
		: number of threads for parallel calculation (Default: 1)(Efficient up to 3)
	* `--verbose`, `-v`
		: verbose (-v: display results / -vv: display calculation results)
	* `--make-template`
		: make template files (`template_ref_param.csv` and `template_exp.csv`) and exit

			```sh
			$ NNfitTool.py --make-template
			```


### NNcountTool.py

```sh
$ NNcountTool.py [-h] -f SEQUENCE_FILE.fasta [SEQUENCE_FILE.fasta ...] -p PARAMETER_FILE.csv -o OUTPUT_FILE.csv [-O] [--make-template]
```

* `-h, --help`
	: show this help message and exit
* `-f SEQUENCE_FILE.fasta [SEQUENCE_FILE.fasta ...]`
	: sequence file of FASTA format
* `-p PARAMETER_FILE.csv`
	: input file for parameters
* `-o OUTPUT_FILE.csv`
	: output file
* `-O`
	: overwrite forcibly
* `--make-template`
	: make template files (template_ref_param.csv) and exit



## Requirement
* Python 3
	* numpy
	* pandas
	* joblib


## License
The MIT License (MIT)

Copyright (c) 2019 Tatsuya Ohyama


## Authors
* Tatsuya Ohyama


## ChangeLog
### Ver. 7.2 (2026-02-03)
* Set the error value for overfitting to 0.001.
* Update template file.

### Ver. 7.1.1 (2026-02-03)
* Fix a bug that failed to reflect changes to the module path in NNcalcTools.py

### Ver. 7.1 (2024-06-23)
* Add `--one-direction` option.
* Change generation algorithm of pair from sequence.

### Ver. 7.0 (2024-06-02)
* Support for mismatch sequence
	* Add `template_exp_mismatch.csv`

### Ver. 6.21 (2023-10-23)
* Fix bug that do not use medium of parameter for prediction value at `<< Sequence >>` section

### Ver. 6.20 (2023-09-26)
* Small changes for NNfitTool.py

### Ver. 6.19 (2022-03-31)
* Add README.md
