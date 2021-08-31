# NNfitTool
Program to calculate Nearest Neighbor parameters fitted from experimental values

## Usage
```
usage: NNfitTool.py [-h] -x EXP.csv [-r REF_PARAM.csv] -o OUTPUT.csv [-O]
                    [-d THRESHOLD] [-i INITIAL_INCREMENT] [-T TEMPERATURE]
                    [-m EVALUATION_METHOD] [-S] [-I LOOP_COUNT] [-e | -es]
                    [-t THREAD] [--verbose] [--make-template]

NNfitTool.py

optional arguments:
  -h, --help            show this help message and exit

Input:
  -x EXP.csv            sequence and experimental value file
                        column: Label, Sequence, dH, dH(error), dS, dS(error), dG, and dG(error)
  -r REF_PARAM.csv      referenced parameter values
                        column: Parameter type, dH, dS, dG, space, dH(change flag), dS(change flag), and dG(change flag))
                        Parameter type: AA/TT, GC/CG, or etc, and re: (regexp) and reg: (counted pattern by regexp)
                        e.g., "re:^A/^T" and "re:^T/^A" (initial parameter for A/T (both specification require))
                        e.g., "reg:.*?G.*?/.*?C.*?" (number of G/C pair parameter)
                        e.g., "reg:./." (length parameter)

Output:
  -o OUTPUT.csv         output file
  -O                    overwrite forcibly

Config:
  -d THRESHOLD          difference threshold of increment for searching (Default: 0.01)
  -i INITIAL_INCREMENT  initial increment (Default: 1.0)
  -T TEMPERATURE        temperature for experimental data (Default: 310.15)
  -m EVALUATION_METHOD  evaluation method (r, r2, diff_mean, diff_std, diff_sum, diff_square) (Default: diff_square)
  -S                    Separately calculate dS (Default: OFF (dS is calculated by Gibbs free energy equation))
  -I LOOP_COUNT         the number of looping optimize (Default: 1)
  -e                    consider with experimental value with error
  -es                   strictly consider with experimental value with error

Misc:
  -t THREAD             number of threads for parallel calculation (Default: 1)
  --verbose, -v         verbose (-v: display results / -vv: display calculation results)
  --make-template       make template files (template_ref_param.csv and template_exp.csv) and exit
```

## Install
1. Build the Python operating environment.
1. Install requirement modules.
1. Copy this program.
    ```
    $ git clone https://github.com/Tatsuya-Ohyama/NNfitTool.git
    ```
1. Give execution rights to the program.
    ```
    $ cd NNfitTool/
    $ chmod +x NNfitTool.py
    ```
1. Check if it can be executed at a minimum.
    ```
    $ ./NNfitTool.py --help
    ```

## Requirement
* [numpy](https://pypi.org/project/numpy/)
* [joblib](https://pypi.org/project/joblib/)

## Licence
[MIT](https://github.com/tcnksm/tool/blob/master/LICENCE)

## Author
[Tatsuya Ohyama](https://github.com/Tatsuya-Ohyama)

