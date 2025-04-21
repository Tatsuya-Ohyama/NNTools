#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
NNfitTool.py -  Fitting program for parameters of Nearest Neighbor method
"""

import sys, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
import csv
import copy
import itertools
from decimal import Decimal, ROUND_HALF_UP, ROUND_HALF_EVEN
from joblib import Parallel, delayed
import pandas as pd

from mods.basicfunc import check_exist, check_overwrite
from mods.parameter import Parameter
from mods.sequence import Sequence
from mods.datagroup import DataGroup



# =============== variable =============== #
VERSION = "7.1"
iteration = 0



# =============== function =============== #
def calculation_worker(parameter, exp_data, mode, increment, threshold_increment, init_direction, verbose, error_sign=None):
	# calculate parameter
	# loop for energy type: dH, dS, and dG
	direction = copy.deepcopy(init_direction)
	if 1 <= verbose:
		print("_/" * 20)
		print("{0:^40}".format("Fitting {0}".format(exp_data.name)))
		print("_/" * 20)

	parameters_opt = [copy.deepcopy(parameter).set_name(parameter_type) for parameter_type in parameter_types]

	evaluation_val = [0.0] * len(parameter_types)
	evaluation_prev = [0.0] * len(parameter_types)
	cnt_i = 0
	while threshold_increment < increment:
		# loop while increment is larger than threshold
		cnt_i += 1
		for parameter_idx, parameter_type in enumerate(parameter_types):
			# loop for parameters
			evaluation_val_direction = []
			evaluation_val_tmp = []
			if direction[parameter_idx] == 0:
				# When not determine direction, determined direction
				# create 0,+,- changed parameter object
				parameter_plus = parameters_opt[parameter_idx].clone()
				parameter_plus.set_name(parameter_type + "_plus")
				parameter_plus.set_parameter(parameter_type, parameter_plus.get_parameter(parameter_type)[0] + increment)
				parameter_minus = parameters_opt[parameter_idx].clone()
				parameter_minus.set_name(parameter_type + "_minus")
				parameter_minus.set_parameter(parameter_type, parameter_minus.get_parameter(parameter_type)[0] - increment)

				# evaluation
				if error_sign is None:
					evaluation_prev[parameter_idx] = exp_data.get_stat(parameters_opt[parameter_idx], mode)
					evaluation_val_direction.append(evaluation_prev[parameter_idx])
					evaluation_val_direction.append(exp_data.get_stat(parameter_plus, mode))
					evaluation_val_direction.append(exp_data.get_stat(parameter_minus, mode))
				else:
					evaluation_prev[parameter_idx] = exp_data.get_stat(parameters_opt[parameter_idx], mode, error_sign=error_sign)
					evaluation_val_direction.append(evaluation_prev[parameter_idx])
					evaluation_val_direction.append(exp_data.get_stat(parameter_plus, mode, error_sign=error_sign))
					evaluation_val_direction.append(exp_data.get_stat(parameter_minus, mode, error_sign=error_sign))
				evaluation_val_tmp.append(evaluation_prev[parameter_idx])

				# determine direction
				min_val = min(evaluation_val_direction)
				min_val_idx = [i for i, x in enumerate(evaluation_val_direction) if min_val == x]
				if len(min_val_idx) != 1 or min_val_idx[0] == 0:
					# When all evaluation_val_tmp is the same even if parameter is changed
					# When evaluation_val_tmp value for base parameter is closest to 1, lock changing
					evaluation_val[parameter_idx] = evaluation_val_direction[0]
					evaluation_val_tmp.append(evaluation_val_direction[0])
					continue

				elif min_val_idx[0] == 1:
					# When evaluation for plus parameter is adopted
					direction[parameter_idx] = 1
					evaluation_val_tmp.append(evaluation_val_direction[1])

				elif min_val_idx[0] == 2:
					# When evaluation for minus parameter is adopted
					direction[parameter_idx] = -1
					evaluation_val_tmp.append(evaluation_val_direction[2])

				else:
					sys.stderr.write("ERROR: undefined condition at determination of direction. {0}\n".format(evaluation_val_tmp))
					sys.exit(1)

			# parepare increased parameter
			parameter_new = parameters_opt[parameter_idx].clone()
			parameter_new.set_name(parameter_type + "_new")
			parameter_new.set_parameter(parameter_type, Decimal(str(parameter_new.get_parameter(parameter_type)[0])) + Decimal(str(direction[parameter_idx] * increment)))

			# evaluation by mode
			if len(evaluation_val_tmp) == 0:
				# When direction determine step is skipped (ex. more than 2nd calculation for determined direction parameter)
				if error_sign is None:
					evaluation_val_tmp.append(exp_data.get_stat(parameters_opt[parameter_idx], mode))
					evaluation_val_tmp.append(exp_data.get_stat(parameter_new, mode))
				else:
					evaluation_val_tmp.append(exp_data.get_stat(parameters_opt[parameter_idx], mode, error_sign=error_sign))
					evaluation_val_tmp.append(exp_data.get_stat(parameter_new, mode, error_sign=error_sign))

			# choose parameter from statistics values (minimum value)
			min_val = min(evaluation_val_tmp)
			min_val_idx = [i for i, x in enumerate(evaluation_val_tmp) if min_val == x]
			if len(min_val_idx) != 1 or min_val_idx[0] == 0:
				# When statistic values for both parameter is the same or
				# prev parameter is closest to 1, only update evaluation_prev
				evaluation_val[parameter_idx] = evaluation_prev[parameter_idx]

			elif min_val_idx[0] == 1:
				# When statistic value for new parameter is closest to 1, update parameter
				parameter_new.set_name(parameter_type)
				parameters_opt[parameter_idx] = parameter_new
				evaluation_val[parameter_idx] = evaluation_val_tmp[1]

			else:
				sys.stderr.write("ERROR: undefined condition at evaluation. {0}\n".format(evaluation_val_tmp))
				sys.exit(1)


		if 2 <= verbose:
			print("-" * 62)
			print("{0} at {1} steps (dt = {2}) of {3} times  Mode: {4}".format(exp_label[exp_idx], cnt_i, increment, iteration, mode))

		# calculate diff statistic values between prev and present
		# and choose largest one (The value that is largely close to the experimental value)
		evaluation_diff = [abs(x - y) for x, y in zip(evaluation_prev, evaluation_val)]
		max_val = max(evaluation_diff)

		if max_val != 0.0:
			# max_val other than 0.0 (max_val = 0.0 means that diff value is converged)
			max_val_idx = [i for i, v in enumerate(evaluation_diff) if v == max_val][0]

			# update_parameter
			new_parameter_type = parameters_opt[max_val_idx].name
			new_parameter_val = parameters_opt[max_val_idx].get_parameter(new_parameter_type)
			for parameter in parameters_opt:
				parameter.set_parameter(new_parameter_type, copy.deepcopy(new_parameter_val))
			parameter = parameters_opt[max_val_idx]
			parameter.set_name(exp_data.name)

			if 2 <= verbose:
				print("{0:^8} {1:^10} {2:^3} {3:^13} {4:^13} {5:^4} {6:^5}".format("Type", "Parameter", "Chg", "E=sum(x'-x)^2", "Diff e (Prev)", "Sign", "Adopt"))
				print("{0:-^8} {1:-^10} {2:-^3} {3:-^13} {4:-^13} {5:-^4} {6:-^5}".format("", "", "", "", "", "", ""))
				print("{0:<8} {1:>10} {2:^3} {3:>13.3f}".format("(Prev)", "", "", evaluation_prev[0]))
				for i, (p, e1, e2, e_diff) in enumerate(zip(parameter_types, evaluation_prev, evaluation_val, evaluation_diff)):
					if i == max_val_idx:
						print("{0:<8} {1[0]:>10.3f} {2:^3} {3:>13.3f} {4:>13.3f}  {5:>2}  {6:^5}".format(p, parameter.get_parameter(p), str(parameter.is_change(p))[0], e2, e_diff, direction[i], "O"))
					else:
						print("{0:<8} {1[0]:>10.3f} {2:^3} {3:>13.3f} {4:>13.3f}  {5:>2}".format(p, parameter.get_parameter(p), str(parameter.is_change(p))[0], e2, e_diff, direction[i]))
				print("")
			evaluation_prev = [evaluation_val[max_val_idx] for parameter in parameter_types]

		else:
			# When all parameters were locked, unlock and change increment
			increment /= 2
			direction = copy.deepcopy(init_direction)

	if 1 <= verbose:
		print("")
		print("===== Last parameter =====")
		print("{0:^8} {1:^10} {2:^3} {3:^10} {4:^10}".format("Type", "Parameter", "Chg", "Error-", "Error+"))
		print("{0:-^8} {1:-^10} {2:-^3} {3:-^10} {4:-^10}".format("", "", "", "", ""))
		for p in parameter_types:
			print("{0:<8} {1[0]:>10.3f} {2:^3} {1[1]:>10.3f} {1[2]:>10.3f}".format(p, parameter.get_parameter(p), str(parameter.is_change(p))[0]))

		print("")
		print("===== Comparing experimental data =====")
		print("{0:^20} {1:^8} {2:^8} {3:^8}".format("Sequence", "Exp.", "Predict", "Diff"))
		print("{0:-^20} {1:-^8} {2:-^8} {3:-^8}".format("", "", "", ""))
		print(exp_data.get_energy(flag_sequence=True, obj_parameters=[parameters[2]]))
		print(exp_data.get_stat(parameters[2], "diff_abs"))
		for row, diff in zip(exp_data.get_energy(flag_sequence=True, parameters=[parameters[2]]), exp_data.get_stat(parameters[2], "diff_abs")):
			print("{0:<20} {1:>8.3f} {2:>8.3f} {3:>8.3f}".format(row[0], row[1], row[2], diff))

		print("")
		print("===== Curve fitting =====")
		print("Slope    :", exp_data.get_stat(parameter, "slope"))
		print("Intercept:", exp_data.get_stat(parameter, "intercept"))
		print("R   (1D) :", exp_data.get_stat(parameter, "r"))
		print("R^2 (1D) :", exp_data.get_stat(parameter, "r2"))
		print("E        :", exp_data.get_stat(parameter, "diff_square"))
		print("\n")

	return parameter


def read_parameters(input_file, parameters, directions):
	"""
	Function to read parameters

	Args:
		input_file (str): reference parameter file
		parameters (list): [obj_Parameter, ...]
		direcitons (list): [obj_Parameter, ...]

	Returns:
		list: [obj_Parameter, ...]
		list: [[direction(int), ...], ...]
		dict: base_pair
	"""
	flag_read = False
	base_pairs = {}
	has_problem_base_pair = False
	with open(input_file, "r") as obj_input:
		reader = csv.reader(obj_input)
		for row_val in reader:
			if row_val[0].lower() == "parameter":
				# read from "Parameter" at col 1
				flag_read = True
				continue

			if flag_read:
				if len(row_val) == 0 or row_val[0] == "":
					# read stop by empty row
					break

				parameter_types.append(row_val[0])
				row_val = [val for val in row_val if val != ""]
				if "/" in row_val[0] \
				and not row_val[0].startswith("init") \
				and not row_val[0].startswith("length") \
				and not row_val[0].startswith("symmetry") \
				and not row_val[0].startswith("re:") \
				and not row_val[0].startswith("reg:"):
					# lexical analysis for parameter label to base pair
					seq1, seq2 = row_val[0].split("/", 1)
					if len(seq1) != len(seq2):
						sys.stderr.write("ERROR: Length of pattern `{}` does not match.\n".format(row_val[0]))
						sys.exit(1)

					list_seq1 = list(seq1)
					list_seq2 = list(seq2)
					for base1, base2 in zip(list_seq1, list_seq2):
						if base1 in base_pairs.keys():
							if base_pairs[base1] != base2:
								sys.stderr.write("WARNING: Conflicting base pairs `{0}/{1}` and `{0}/{2}` are found in reading parameter `{3}`.\n".format(base1, base_pairs[base1], base2, row_val[0]))
								has_problem_base_pair = True

						else:
							base_pairs[base1] = base2

						if base2 in base_pairs.keys():
							if base_pairs[base2] != base1:
								sys.stderr.write("WARNING: Conflicting base pairs `{0}/{1}` and `{0}/{2}` are found in reading parameter `{3}`.\n".format(base2, base_pairs[base2], base1, row_val[0]))
								has_problem_base_pair = True

						else:
							base_pairs[base2] = base1

				parameters[0].append_parameter(row_val[0], float(row_val[1]))
				parameters[1].append_parameter(row_val[0], float(row_val[2]))
				parameters[2].append_parameter(row_val[0], float(row_val[3]))
				if 7 <= len(row_val):
					# change flag
					row_val[4:7] = [True if x.capitalize() == "True" else False for x in row_val[4:7]]
					parameters[0].set_change_stat(row_val[0], row_val[4])
					parameters[1].set_change_stat(row_val[0], row_val[5])
					parameters[2].set_change_stat(row_val[0], row_val[6])
				else:
					parameters[0].set_change_stat(row_val[0], True)
					parameters[1].set_change_stat(row_val[0], True)
					parameters[2].set_change_stat(row_val[0], True)

				if len(row_val) == 10:
					directions[0].append(float(row_val[7]))
					directions[1].append(float(row_val[8]))
					directions[2].append(float(row_val[9]))
				else:
					directions[0].append(0)
					directions[1].append(0)
					directions[2].append(0)

	if has_problem_base_pair:
		sys.stderr.write("WARNING: Because some conflicts were found in base pair information, sequence parameters are determined as mismatched sequences.\n")
		base_pairs = None

	return parameters, directions, base_pairs


# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="NNfitTool.py", formatter_class=argparse.RawTextHelpFormatter)
	input_group = parser.add_argument_group("Input")
	input_group.add_argument("-x", dest="EXPERIMENT_FILE", metavar="EXP.csv", required=True,
	help = """sequence and experimental value file
column: Label, Sequence, dH, dH(error), dS, dS(error), dG, and dG(error)
""")
	input_group.add_argument("-r", dest="REF_PARAM", metavar="REF_PARAM.csv", required=True,
	help = """referenced parameter values
column: `Parameter`, `dH`, `dS`, `dG`, space, `dH (change)`, `dS (change)`, `dG (change)`, space, `dH (Direction)`, `dS (Direction)`, and `dG (Direction)`
`Parameter`: AA/TT, GC/CG, or etc, and re: (regexp) and reg: (counted pattern by regexp)
e.g., "re:^A/^T" and "re:^T/^A" (initial parameter for A/T (both specification require))
e.g., "reg:.*?G.*?/.*?C.*?" (number of G/C pair parameter)
e.g., "reg:./." (length parameter)
""")

	output_group = parser.add_argument_group("Output")
	output_group.add_argument("-o", dest="OUTPUT_FILE", metavar="OUTPUT.csv", required=True, help="output file")
	output_group.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite forcibly")

	config_group = parser.add_argument_group("Config")
	config_group.add_argument("-d", dest="THRESHOLD_INCREMENT", metavar="THRESHOLD", type=float, default=0.01, help="difference threshold of increment for searching (Default: 0.01)")
	config_group.add_argument("-i", dest="INITIAL_INCREMENT", metavar="INITIAL_INCREMENT", type=float, default=1.0, help="initial increment (Default: 1.0)")
	config_group.add_argument("-T", dest="TEMPERATURE", metavar="TEMPERATURE", type=float, default=310.15, help="temperature for experimental data (Default: 310.15)")
	config_group.add_argument("-m", dest="MODE", metavar="EVALUATION_METHOD", default="diff_square", choices=["r", "r2", "diff_mean", "diff_std", "diff_sum", "diff_square"], help="evaluation method (r, r2, diff_mean, diff_std, diff_sum, diff_square) (Default: diff_square)")
	config_group.add_argument("-S", dest="FLAG_SEPARATE", action="store_true", default=False, help="Separately calculate dS (Default: OFF (dS is calculated by Gibbs free energy equation))")
	config_group.add_argument("-I", dest="OPTIMIZE_COUNT", metavar="LOOP_COUNT", type=int, default=1, help="the number of looping optimize (Default: 1)")
	config_group.add_argument("--one-direction", dest="FLAG_ONE_DIRECTION", action="store_true", default=False, help="Do not search for reverse order pattern (For example, this program searches AC/TG and reverse order pattern GT/CA as the same pattern. This option does not allow it.)")
	error_group = config_group.add_mutually_exclusive_group()
	error_group.add_argument("-e", dest="FLAG_ERROR", action="store_true", default=False, help="consider with experimental value with error")

	misc_group = parser.add_argument_group("Misc")
	misc_group.add_argument("-t", dest="THREAD", metavar="THREAD", type=int, default=1, help="number of threads for parallel calculation (Default: 1)(Efficient up to 3)")
	misc_group.add_argument("--verbose", "-v", dest="VERBOSE", action="count", default=0, help="verbose (-v: display results / -vv: display calculation results)")

	args = parser.parse_args()

	# calculate target
	target_list = [0, 2]	# without dS
	if args.FLAG_SEPARATE:
		# with dS
		target_list = [0, 1, 2]


	# initial parameter
	exp_label = ["dH", "dS", "dG"]
	parameters = [Parameter().set_name(label) for label in exp_label]
	parameters_init = []
	directions = [[] for idx in range(len(exp_label))]


	# loading reference parameter
	parameter_types = []
	check_exist(args.REF_PARAM, 2)
	parameters, directions, base_pairs = read_parameters(args.REF_PARAM, parameters, directions)
	for obj_parameter in parameters:
		obj_parameter.set_one_direction(args.FLAG_ONE_DIRECTION)

	# reading sequence and experimental data
	check_exist(args.EXPERIMENT_FILE, 2)
	exp_datas = [DataGroup(label) for label in exp_label]
	sys.stderr.write("Read (experimental file): {}".format(args.EXPERIMENT_FILE))
	sys.stderr.flush()
	df_exp = pd.read_csv(args.EXPERIMENT_FILE, header=0)
	df_exp = df_exp.astype({"dH": float, "dH (error)": float, "dS": float, "dS (error)": float, "dG": float, "dG (error)": float})
	df_exp = df_exp.fillna({"Comment": "", "dH": 0.0, "dH (error)": 0.0, "dS": 0.0, "dS (error)": 0.0, "dG": 0.0, "dG (error)": 0.0})
	format_type = 1
	if "Sequence1" in df_exp.columns and "Sequence2" in df_exp.columns:
		format_type = 2
		sys.stderr.write(" (complementary sequences are written by manual)\n".format(format_type))
		for row in df_exp.itertuples():
			obj_sequence = Sequence(row[1]).set_sequence(row[2], reversed(row[3]))
			obj_sequence.get_freq(parameters[0])	# make cache for freq
			exp_datas[0].append(obj_sequence, row[4], row[5])
			exp_datas[1].append(obj_sequence, row[6], row[7])
			exp_datas[2].append(obj_sequence, row[8], row[9])

	else:
		sys.stderr.write(" (no complementary sequences)\n".format(format_type))
		for row in df_exp.itertuples():
			if base_pairs is None:
				sys.stderr.write("ERROR: base pair information incomplete.\n")
				sys.exit(1)

			obj_sequence = Sequence(row[1]).set_sequence(row[2]).generate_complement(base_pairs)
			obj_sequence.get_freq(parameters[0])	# make cache for freq
			exp_datas[0].append(obj_sequence, row[3], row[4])
			exp_datas[1].append(obj_sequence, row[5], row[6])
			exp_datas[2].append(obj_sequence, row[7], row[8])


	# optimize parameter
	parameters_init = [copy.deepcopy(parameter) for parameter in parameters]
	for loop_idx in range(args.OPTIMIZE_COUNT):
		sys.stderr.write("Optimize parameters for {0} times.\n".format(loop_idx + 1))
		iteration = loop_idx + 1
		if 1 < args.THREAD:
			# multi-thread
			parameters_tmp = Parallel(n_jobs = args.THREAD)([
				delayed(calculation_worker)(
					parameter=parameters[exp_idx],
					exp_data=exp_datas[exp_idx],
					mode=args.MODE,
					increment=args.INITIAL_INCREMENT,
					threshold_increment=args.THRESHOLD_INCREMENT,
					init_direction=directions[exp_idx],
					verbose=0
				) for exp_idx in target_list])
			if args.FLAG_SEPARATE:
				parameters = parameters_tmp
			else:
				parameters[0] = parameters_tmp[0]
				parameters[2] = parameters_tmp[1]

		else:
			# single-thread
			for exp_idx in target_list:
				parameters[exp_idx] = calculation_worker(
					parameter=parameters[exp_idx],
					exp_data=exp_datas[exp_idx],
					mode=args.MODE,
					increment=args.INITIAL_INCREMENT,
					threshold_increment=args.THRESHOLD_INCREMENT,
					init_direction=directions[exp_idx],
					verbose=args.VERBOSE
				)

		# optimize for parameter by experimental values with error
		if args.FLAG_ERROR:
			sys.stderr.write("Optimize parameters with errors.\n")
			negative = [1 for x in range(len(exp_datas[0].get_sequence()))]
			positive = [-1 for x in range(len(exp_datas[0].get_sequence()))]
			new_parameters = []

			if args.THREAD is not None:
				# multi-thread
				new_parameters = Parallel(n_jobs=args.THREAD)([
					delayed(calculation_worker)(
						parameter=parameters[exp_idx],
						exp_data=exp_datas[exp_idx],
						mode=args.MODE,
						increment=args.INITIAL_INCREMENT,
						threshold_increment=args.THRESHOLD_INCREMENT,
						init_direction=directions[exp_idx],
						verbose=0,
						error_sign=sign
					) for exp_idx in target_list
						for sign in [negative, positive]
				])

			else:
				# single-thread
				for exp_idx in target_list:
					# exp values with errors
					for sign in [negative, positive]:
						new_parameters.append(calculation_worker(
							parameter=parameters[exp_idx],
							exp_data=exp_datas[exp_idx],
							mode=args.MODE,
							increment=args.INITIAL_INCREMENT,
							threshold_increment=args.THRESHOLD_INCREMENT,
							init_direction=directions[exp_idx],
							verbose=args.VERBOSE,
							error_sign=sign
						))

			for idx, exp_idx in enumerate(target_list):
				parameters[exp_idx].update_parameter_error("all", new_parameters[idx*2+0].get_parameter())
				parameters[exp_idx].update_parameter_error("all", new_parameters[idx*2+1].get_parameter())

	if not args.FLAG_SEPARATE and exp_datas[0].is_fitting and exp_datas[1].is_fitting:
		# calculate dS: (dH - dG) / T * 1000
		dS = {
			parameter_type: [
				(parameters[0].get_parameter(parameter_type, data_type="raw")[idx] - parameters[2].get_parameter(parameter_type, data_type = "raw")[idx]) / args.TEMPERATURE * 1000
				for idx in range(3)
			] for parameter_type in parameter_types}
		for parameter_type in parameter_types:
			parameters[1].set_parameter(parameter_type, dS[parameter_type])


	# output
	if args.FLAG_OVERWRITE == False:
		check_overwrite(args.OUTPUT_FILE)

	with open(args.OUTPUT_FILE, "w") as obj_output:
		writer = csv.writer(obj_output)
		writer.writerow(["<< Input >>"])
		writer.writerow(["Program version", VERSION])
		writer.writerow(["Experimental data", args.EXPERIMENT_FILE])
		writer.writerow(["Reference parameter", args.REF_PARAM])
		writer.writerow(["Initial iteration", args.INITIAL_INCREMENT])
		writer.writerow(["Increment threshold", args.THRESHOLD_INCREMENT])
		writer.writerow(["Iteration (whole)", args.OPTIMIZE_COUNT])
		writer.writerow(["Temperature", args.TEMPERATURE])
		writer.writerow(["Separate calculation", args.FLAG_SEPARATE])
		writer.writerow(["Error calculation", args.FLAG_ERROR])
		writer.writerow(["Evaluation mode", args.MODE])
		writer.writerow([""])

		writer.writerow(["Initial parameter", "dH", "dS", "dG", "", "Change (dH)", "Change (dS)", "Change (dG)", "", "Direction (dH)", "Direction (dS)", "Direction (dG)"])
		for idx, parameter_type in enumerate(parameter_types):
			parameter_dH = [Decimal(str(x)).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP) for x in parameters_init[0].get_parameter(data_type="fix")[parameter_type]]
			parameter_dS = [Decimal(str(x)).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP) for x in parameters_init[1].get_parameter(data_type="fix")[parameter_type]]
			parameter_dG = [Decimal(str(x)).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP) for x in parameters_init[2].get_parameter(data_type="fix")[parameter_type]]

			writer.writerow([
				parameter_type,
				parameter_dH[0],
				parameter_dS[0],
				parameter_dG[0],
				"",
				parameters_init[0].is_change(parameter_type),
				parameters_init[1].is_change(parameter_type),
				parameters_init[2].is_change(parameter_type),
				"",
				directions[0][idx],
				directions[1][idx],
				directions[2][idx]
				])
		writer.writerow([""])
		writer.writerow([""])

		writer.writerow(["<< Results >>"])
		if args.FLAG_SEPARATE:
			writer.writerow(["Parameter (optimized)", "dH", "dH (error)", "dS", "dS (error)", "dG", "dG (error)", "", "Change (dH)", "Change (dS)", "Change (dG)"])
		else:
			writer.writerow(["Parameter (optimized)", "dH", "dH (error)", "dS * 1000", "dS * 1000 (error)", "dG", "dG (error)", "", "Change (dH)", "Change (dS)", "Change (dG)"])
		for parameter_type in parameter_types:
			parameter_dH = [Decimal(str(x)).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP) for x in parameters[0].get_parameter(data_type="fix")[parameter_type]]
			parameter_dS = [Decimal(str(x)).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP) for x in parameters[1].get_parameter(data_type="fix")[parameter_type]]
			parameter_dG = [Decimal(str(x)).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP) for x in parameters[2].get_parameter(data_type="fix")[parameter_type]]

			writer.writerow([
				parameter_type,
				parameter_dH[0], parameter_dH[1],
				parameter_dS[0], parameter_dS[1],
				parameter_dG[0], parameter_dG[1],
				"",
				parameters[0].is_change(parameter_type),
				parameters[1].is_change(parameter_type),
				parameters[2].is_change(parameter_type)
				])
		writer.writerow([""])
		writer.writerow([""])

		writer.writerow(["<< Sequence >>"])
		writer.writerow([
			"Name",
			"Sequence A (5->3)",
			"Sequence B (3->5)",
			"Exp. (dH)",
			"Exp. (dS)",
			"Exp. (dG)",
			"Predict (dH)",
			"Predict (dS)",
			"Predict (dG)",
			"Diff. (dH)",
			"Diff. (dS)",
			"Diff. (dG)",
			""
			] + parameter_types)

		for sequence, dH, dS, dG in zip(
				exp_datas[0].get_sequence(),
				exp_datas[0].get_energy(flag_sequence=True, obj_parameters=[parameters[0]], data_type="fix"),
				exp_datas[1].get_energy(flag_sequence=True, obj_parameters=[parameters[1]], data_type="fix"),
				exp_datas[2].get_energy(flag_sequence=True, obj_parameters=[parameters[2]], data_type="fix")):
			name = sequence.name
			seqA = sequence.get_sequence("string")
			seqB = sequence.get_complement("string")
			exp_dH = dH[1]
			exp_dS = dS[1]
			exp_dG = dG[1]
			pred_dH = dH[3]
			pred_dS = dS[3]
			pred_dG = dG[3]
			diff_dH = pred_dH - exp_dH
			diff_dS = pred_dS - exp_dS
			diff_dG = pred_dG - exp_dG
			freq = sequence.get_freq(parameters[0])
			writer.writerow([name, seqA, seqB, exp_dH, exp_dS, exp_dG, pred_dH, pred_dS, pred_dG, diff_dH, diff_dS, diff_dG] + [""] + [freq[parameter_type] for parameter_type in parameter_types])
		writer.writerow([""])

		writer.writerow(["Correlation: exp. vs predict", "", "dH", "dS", "dG"])
		writer.writerow(["", "(Stat: R)"] + [exp_datas[idx].get_stat(parameters[idx], "r") for idx in range(3)])
		writer.writerow(["", "(Stat: R2)"] + [exp_datas[idx].get_stat(parameters[idx], "r2") for idx in range(3)])
		writer.writerow(["", "(Stat: Slope)"] + [exp_datas[idx].get_stat(parameters[idx], "slope") for idx in range(3)])
		writer.writerow(["", "(Stat: Intercept)"] + [exp_datas[idx].get_stat(parameters[idx], "intercept") for idx in range(3)])
		writer.writerow(["", "(Stat: Diff Mean)"] + [exp_datas[idx].get_stat(parameters[idx], "diff_mean") for idx in range(3)])
		writer.writerow(["", "(Stat: Diff Sum)"] + [exp_datas[idx].get_stat(parameters[idx], "diff_sum") for idx in range(3)])
		writer.writerow(["", "(Stat: Diff Sq)"] + [exp_datas[idx].get_stat(parameters[idx], "diff_square") for idx in range(3)])
