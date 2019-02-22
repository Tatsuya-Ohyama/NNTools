#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
NNfitTool.py -  Fitting program for parameters of Nearest Neighbor method
"""

import sys, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
import csv
import numpy as np
import copy
import itertools
from decimal import Decimal, ROUND_HALF_UP, ROUND_HALF_EVEN
from joblib import Parallel, delayed
from pprint import pprint

from basic_func import check_exist, check_overwrite
from classes.Parameter import Parameter
from classes.Sequence import Sequence
from classes.DataGroup import DataGroup


# =============== variable =============== #
parameter_types = ["AA/TT", "AT/TA", "TA/AT", "CA/GT", "GT/CA", "CT/GA", "GA/CT", "CG/GC", "GC/CG", "GG/CC", "init_GC", "init_AT", "symmetry", "5term_TA"]
VERSION = "4.3 (#91)"


# =============== function =============== #
def calculation_worker(parameter, exp_data, mode, increment, threshold_increment, verbose, error_sign = None):
	# calculate parameter
	# loop for energy type: dH, dS, and dG
	if 1 <= verbose:
		print("_/" * 20)
		print("{0:^40}".format("Fitting {0}".format(exp_data.get_name())))
		print("_/" * 20)

	parameters_opt = [copy.deepcopy(parameter).set_name(parameter_type) for parameter_type in parameter_types]
	direction = [0] * len(parameter_types)

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
					evaluation_prev[parameter_idx] = exp_data.get_stat(parameters_opt[parameter_idx], mode, error_sign = error_sign)
					evaluation_val_direction.append(evaluation_prev[parameter_idx])
					evaluation_val_tmp.append(evaluation_prev[parameter_idx])
					evaluation_val_direction.append(exp_data.get_stat(parameter_plus, mode, error_sign = error_sign))
					evaluation_val_direction.append(exp_data.get_stat(parameter_minus, mode, error_sign = error_sign))
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
					sys.stderr.write("ERROR: undefined condition.\n")
					sys.exit(1)

			# parepare increased parameter
			parameter_new = parameters_opt[parameter_idx].clone()
			parameter_new.set_name(parameter_type + "_new")
			parameter_new.set_parameter(parameter_type, parameter_new.get_parameter(parameter_type)[0] + direction[parameter_idx] * increment)

			# evaluation by mode
			if len(evaluation_val_tmp) == 0:
				# When direction determine step is skipped (ex. more than 2nd calculation for determined direction parameter)
				if error_sign is None:
					evaluation_val_tmp.append(exp_data.get_stat(parameters_opt[parameter_idx], mode))
					evaluation_val_tmp.append(exp_data.get_stat(parameter_new, mode))
				else:
					evaluation_val_tmp.append(exp_data.get_stat(parameters_opt[parameter_idx], mode, error_sign = error_sign))
					evaluation_val_tmp.append(exp_data.get_stat(parameter_new, mode, error_sign = error_sign))

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
				sys.stderr.write("ERROR: undefined condition. {0}\n".format(min_val_idx[0]))
				sys.exit(1)


		if 2 <= verbose:
			print("-" * 75)
			print("{0} at {1} steps (dt = {2})  Mode: {3}".format(exp_label[exp_idx], cnt_i, increment, mode))

		# calculate diff statistic values between prev and present
		# and choose largest one (The value that is largely close to the experimental value)
		evaluation_diff = [abs(x - y) for x, y in zip(evaluation_prev, evaluation_val)]
		max_val = max(evaluation_diff)

		if max_val != 0.0:
			# max_val other than 0.0 (max_val = 0.0 means that diff value is converged)
			max_val_idx = [i for i, v in enumerate(evaluation_diff) if v == max_val][0]

			# update_parameter
			new_parameter_type = parameters_opt[max_val_idx].get_name()
			new_parameter_val = parameters_opt[max_val_idx].get_parameter(new_parameter_type)
			for parameter in parameters_opt:
				parameter.set_parameter(new_parameter_type, copy.deepcopy(new_parameter_val))
			parameter = parameters_opt[max_val_idx]
			parameter.set_name(exp_data.get_name())

			if 2 <= verbose:
				print("{0:^8} {1:^10} {2:^3} {3:^13} {4:^13} {5:^5}".format("Type", "Parameter", "Chg", "E=sum(x'-x)^2", "Diff e (Prev)", "Adopt"))
				print("{0:-^8} {1:-^10} {2:-^3} {3:-^13} {4:-^13} {5:-^5}".format("", "", "", "", "", ""))
				print("{0:<8} {1:>10} {2:^3} {3:>13.3f}".format("(Prev)", "", "", evaluation_prev[0]))
				for i, (p, e1, e2, e_diff) in enumerate(zip(parameter_types, evaluation_prev, evaluation_val, evaluation_diff)):
					if i == max_val_idx:
						print("{0:<8} {1[0]:>10.3f} {2:^3} {3:>13.3f} {4:>13.3f} {5:^5}".format(p, parameter.get_parameter(p), str(parameter.is_change(p))[0], e2, e_diff, "O"))
					else:
						print("{0:<8} {1[0]:>10.3f} {2:^3} {3:>13.3f} {4:>13.3f}".format(p, parameter.get_parameter(p), str(parameter.is_change(p))[0], e2, e_diff))
				print("")
			evaluation_prev = [evaluation_val[max_val_idx] for parameter in parameter_types]

		else:
			# When all parameters were locked, unlock and change increment
			increment /= 2
			direction = [0 for x in parameter_types]

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
		for row, diff in zip(exp_data.get_energy(True, [parameters[2]]), exp_data.get_stat(parameters[2], "diff_abs")):
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



# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "NNfitTool.py", formatter_class=argparse.RawTextHelpFormatter)
	input_group = parser.add_argument_group("Input")
	input_group.add_argument("-x", dest = "experiment_file", metavar = "EXP.csv", required = True, help = "sequence and experimental value file")
	input_group.add_argument("-r", dest = "ref_param", metavar = "REF_PARAM.csv", help = "referenced parameter values")

	output_group = parser.add_argument_group("Output")
	output_group.add_argument("-o", dest = "output_file", metavar = "OUTPUT.csv", required = True, help = "output file")
	output_group.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly")

	config_group = parser.add_argument_group("Config")
	config_group.add_argument("-d", dest = "threshold_increment", metavar = "THRESHOLD", type = float, default = 0.01, help = "difference threshold of increment for searching (Default: 0.01)")
	config_group.add_argument("-i", dest = "initial_increment", metavar = "INITIAL_INCREMENT", type = float, default = 1.0, help = "initial increment (Default: 1.0)")
	config_group.add_argument("-T", dest = "temperature", metavar = "TEMPERATURE", type = float, default = 310.15, help = "temperature for experimental data (Default: 310.15)")
	config_group.add_argument("-m", dest = "mode", metavar = "EVALUATION_METHOD", default = "diff_square", choices = ["r", "r2", "diff_mean", "diff_std", "diff_sum", "diff_square"], help = "evaluation method (r, r2, diff_mean, diff_std, diff_sum, diff_square) (Default: diff_square)")
	config_group.add_argument("-S", dest = "flag_separate", action = "store_true", default = False, help = "Separately calculate dS (Default: OFF (dS is calculated by Gibbs free energy equation))")
	error_group = config_group.add_mutually_exclusive_group()
	error_group.add_argument("-e", dest = "flag_error", action = "store_true", default = False, help = "consider with experimental value with error")
	error_group.add_argument("-es", dest = "flag_error_strict", action = "store_true", default = False, help = "strictly consider with experimental value with error")

	misc_group = parser.add_argument_group("Misc")
	misc_group.add_argument("-t", dest = "thread", metavar = "THREAD", type = int, help = "number of threads for parallel calculation (Default: None)")
	misc_group.add_argument("--verbose", "-v", dest = "verbose", action = "count", default = 0, help = "verbose (-v: display results / -vv: display calculation results)")
	args = parser.parse_args()


	# calculate target
	target_list = [0, 2]	# without dS
	if args.flag_separate:
		# with dS
		target_list = [0, 1, 2]

	# initial parameter
	exp_label = ["dH", "dS", "dG"]
	parameters = [Parameter().set_name(label) for label in exp_label]
	parameters_init = []

	# loading reference parameter
	if args.ref_param is not None:
		check_exist(args.ref_param, 2)

		flag_read = False
		flag_init = False
		pos_sep = 0
		pos_offset = 1
		with open(args.ref_param, "r") as obj_input:
			reader = csv.reader(obj_input)
			for line_val in reader:
				if "Parameter" in line_val[0]:
					# read from "Parameter" at col 1
					flag_read = True
					flag_init = True
					continue
				if len(line_val) == 0 or line_val[0] == "":
					# read stop by empty row
					break

				if flag_read:
					if flag_init:
						pos_sep = line_val.index("")
						if pos_sep == 4:
							# without error
							pos_offset = 1

						elif pos_sep == 7:
							pos_offset = 2
						else:
							sys.stderr.write("ERROR: undefined format for reference parameter file.\n")
							sys.exit(1)
						flag_init = False

					parameters[0].set_parameter(line_val[0], float(line_val[pos_sep - pos_offset * 3]))
					parameters[1].set_parameter(line_val[0], float(line_val[pos_sep - pos_offset * 2]))
					parameters[2].set_parameter(line_val[0], float(line_val[pos_sep - pos_offset * 1]))
					if len(line_val) == 8:
						# format with error
						line_val[5:8] = [True if x.capitalize() == "True" else False for x in line_val[5:8]]
					elif len(line_val) == 11:
						# format without error
						line_val[8:11] = [True if x.capitalize() == "True" else False for x in line_val[8:11]]
					parameters[0].set_change_stat(line_val[0], line_val[pos_sep + 1])
					parameters[1].set_change_stat(line_val[0], line_val[pos_sep + 2])
					parameters[2].set_change_stat(line_val[0], line_val[pos_sep + 3])
	parameters_init = [copy.deepcopy(parameter) for parameter in parameters]

	# reading sequence and experimental data
	check_exist(args.experiment_file, 2)
	exp_datas = [DataGroup(label) for label in exp_label]
	with open(args.experiment_file, "r") as obj_input:
		reader = csv.reader(obj_input)

		# Ignore line number 1 (header) in CSV
		next(reader)

		for line_val in reader:
			exp_datas[0].append(
				Sequence(line_val[0], line_val[1]).set_energy_type(exp_label[0]),
				float(line_val[2]),
				float(line_val[3])
			)
			exp_datas[1].append(
				Sequence(line_val[0], line_val[1]).set_energy_type(
					exp_label[1]),
					float(line_val[4]),
					float(line_val[5])
			)
			exp_datas[2].append(
				Sequence(line_val[0], line_val[1]).set_energy_type(
					exp_label[2]),
					float(line_val[6]),
					float(line_val[7])
			)
	sys.stderr.write("Loading experimental values ... done.\n")

	# optimize parameter
	sys.stderr.write("Optimize parameters.\n")
	if args.thread is not None:
		# multi-thread
		parameters_tmp = Parallel(n_jobs = args.thread)([
			delayed(calculation_worker)(
				parameters[exp_idx],
				exp_datas[exp_idx],
				args.mode,
				args.initial_increment,
				args.threshold_increment,
				0
			) for exp_idx in target_list])
		if args.flag_separate:
			parameters = parameters_tmp
		else:
			parameters[0] = parameters_tmp[0]
			parameters[2] = parameters_tmp[1]

	else:
		# single-thread
		for exp_idx in target_list:
			parameters[exp_idx] = calculation_worker(
				parameters[exp_idx],
				exp_datas[exp_idx],
				args.mode,
				args.initial_increment,
				args.threshold_increment,
				args.verbose
			)
	sys.stderr.write("done.\n")

	# optimize for parameter by experimental values with error
	if args.flag_error:
		sys.stderr.write("Optimize parameters with errors.\n")
		negative = [1 for x in range(len(exp_datas[0].get_sequence()))]
		positive = [-1 for x in range(len(exp_datas[0].get_sequence()))]
		new_parameters = []
		if args.thread is not None:
			# multi-thread
			new_parameters = Parallel(n_jobs = args.thread)([
				delayed(calculation_worker)(
					parameters[exp_idx],
					exp_datas[exp_idx],
					args.mode,
					args.initial_increment,
					args.threshold_increment,
					0,
					sign
				) for exp_idx in target_list
					for sign in [negative, positive]
			])

		else:
			# single-thread
			for exp_idx in target_list:
				# exp values with errors
				for sign in [negative, positive]:
					new_parameters.append(calculation_worker(
						parameters[exp_idx],
						exp_datas[exp_idx],
						args.mode,
						args.initial_increment,
						args.threshold_increment,
						args.verbose,
						sign
					))

		for idx, exp_idx in enumerate(target_list):
			parameters[exp_idx].update_parameter_error("all", new_parameters[idx * 2 + 0].get_parameter())
			parameters[exp_idx].update_parameter_error("all", new_parameters[idx * 2 + 1].get_parameter())
		sys.stderr.write("done.\n")

	elif args.flag_error_strict:
		sys.stderr.write("Optimize parameters with errors by strict mode.\n")
		max_iter = len(list(itertools.product([-1, 1], repeat = len(exp_datas[0].get_sequence()))))
		if args.thread is not None:
			for exp_idx in target_list:
				sys.stderr.write("Calculalte for {0} with error: {1} steps\n".format(exp_label[exp_idx], max_iter))
				cnt = 0
				calc_set = []
				for job_idx, exp_error_pattern in enumerate(itertools.product([-1, 1], repeat = len(exp_datas[0].get_sequence()))):
					calc_set.append(exp_error_pattern)
					cnt += 1
					if 50 <= cnt:
						cnt = 0
						parameter_c = Parallel(n_jobs = args.thread)([
							delayed(calculation_worker)(
								parameters[exp_idx],
								exp_datas[exp_idx],
								args.mode,
								args.initial_increment,
								args.threshold_increment,
								0,
								exp_error_pattern
							) for error_pattern in calc_set
						])
						for new_parameter in parameter_c:
							parameters[exp_idx].update_parameter_error("all", new_parameter.get_parameter())
						calc_set = []
						sys.stderr.write("calculation of {0} with error for {1}-{2} ... done.\n".format(exp_label[exp_idx], job_idx + 1 - 50, job_idx + 1))

				if len(calc_set) != 0:
					parameter_c = Parallel(n_jobs = args.thread)([
						delayed(calculation_worker)(
							parameters[exp_idx],
							exp_datas[exp_idx],
							args.mode,
							args.initial_increment,
							args.threshold_increment,
							0,
							exp_error_pattern
						) for error_pattern in calc_set
					])
					for new_parameter in parameter_c:
						parameters[exp_idx].update_parameter_error("all", new_parameter.get_parameter())

		else:
			for exp_idx in target_list:
				# exp values with error
				for exp_error_pattern in itertools.product([-1, 1], repeat = len(exp_datas[0].get_sequence())):
					new_parameter = calculation_worker(
						parameters[exp_idx],
						exp_datas[exp_idx],
						args.mode,
						args.initial_increment,
						args.threshold_increment,
						args.verbose,
						exp_error_pattern
					)
					parameters[exp_idx].update_parameter_error("all", new_parameter.get_parameter())
		sys.stderr.write("done.\n")

	# calculate dS
	if not args.flag_separate:
		dS = {
			parameter_type: [
				(parameters[0].get_parameter(parameter_type, data_type = "raw")[idx] - parameters[2].get_parameter(parameter_type, data_type = "raw")[idx]) / args.temperature * 1000
				for idx in range(3)
			] for parameter_type in parameter_types}
		for parameter_type in parameter_types:
			parameters[1].set_parameter(parameter_type, dS[parameter_type])

	# output
	if args.flag_overwrite == False:
		check_overwrite(args.output_file)

	with open(args.output_file, "w") as obj_output:
		writer = csv.writer(obj_output)
		writer.writerow(["<< Input >>"])
		writer.writerow(["Program version", VERSION])
		writer.writerow(["Experimental data", args.experiment_file])
		writer.writerow(["Initial iteration", args.initial_increment])
		writer.writerow(["Increment threshold", args.threshold_increment])
		writer.writerow(["Temperature", args.temperature])
		writer.writerow(["Separate calculation", args.flag_separate])
		writer.writerow(["Evaluation mode", args.mode])
		writer.writerow(["Reference parameter", args.ref_param])
		writer.writerow([""])

		writer.writerow(["Initial parameter", "dH", "dS", "dG", "", "Change (dH)", "Change (dS)", "Change (dG)"])
		for parameter_type in parameter_types:
			parameter_dH = [Decimal(str(x)).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP) for x in parameters_init[0].get_parameter(data_type = "fix")[parameter_type]]
			parameter_dS = [Decimal(str(x)).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP) for x in parameters_init[1].get_parameter(data_type = "fix")[parameter_type]]
			parameter_dG = [Decimal(str(x)).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP) for x in parameters_init[2].get_parameter(data_type = "fix")[parameter_type]]

			writer.writerow([
				parameter_type,
				parameter_dH[0],
				parameter_dS[0],
				parameter_dG[0],
				"",
				parameters_init[0].is_change(parameter_type),
				parameters_init[1].is_change(parameter_type),
				parameters_init[2].is_change(parameter_type)
				])
		writer.writerow([""])
		writer.writerow([""])

		writer.writerow(["<< Results >>"])
		writer.writerow(["Parameter (optimized)", "dH", "dH (error)", "dS", "dS (error)", "dG", "dG (error)", "", "Change (dH)", "Change (dS)", "Change (dG)"])
		for parameter_type in parameter_types:
			parameter_dH = [Decimal(str(x)).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP) for x in parameters[0].get_parameter(data_type = "fix")[parameter_type]]
			parameter_dS = [Decimal(str(x)).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP) for x in parameters[1].get_parameter(data_type = "fix")[parameter_type]]
			parameter_dG = [Decimal(str(x)).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP) for x in parameters[2].get_parameter(data_type = "fix")[parameter_type]]

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
			"Sequence",
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

		for sequence, dH, dS, dG in zip(exp_datas[0].get_sequence(), exp_datas[0].get_energy(flag_sequence = True, obj_parameters = [parameters[0]]), exp_datas[1].get_energy(flag_sequence = True, obj_parameters = [parameters[1]]), exp_datas[2].get_energy(flag_sequence = True, obj_parameters = [parameters[2]])):
			name = sequence.get_name()
			seq = sequence.get_sequence("string")
			exp_dH = dH[1]
			exp_dS = dS[1]
			exp_dG = dG[1]
			pred_dH = dH[3]
			pred_dS = dS[3]
			pred_dG = dG[3]
			diff_dH = pred_dH - exp_dH
			diff_dS = pred_dS - exp_dS
			diff_dG = pred_dG - exp_dG
			freq = [sequence.get_freq()[parameter_type] for parameter_type in parameter_types]
			writer.writerow([name, seq, exp_dH, exp_dS, exp_dG, pred_dH, pred_dS, pred_dG, diff_dH, diff_dS, diff_dG] + [""] + freq)
		writer.writerow([""])

		writer.writerow(["Correlation: exp. vs predict", "", "dH", "dS", "dG"])
		writer.writerow(["", "(Stat: R)"] + [exp_datas[idx].get_stat(parameters[idx], "r") for idx in range(3)])
		writer.writerow(["", "(Stat: R2)"] + [exp_datas[idx].get_stat(parameters[idx], "r2") for idx in range(3)])
		writer.writerow(["", "(Stat: Slope)"] + [exp_datas[idx].get_stat(parameters[idx], "slope") for idx in range(3)])
		writer.writerow(["", "(Stat: Intercept)"] + [exp_datas[idx].get_stat(parameters[idx], "intercept") for idx in range(3)])
		writer.writerow(["", "(Stat: Diff Mean)"] + [exp_datas[idx].get_stat(parameters[idx], "diff_mean") for idx in range(3)])
		writer.writerow(["", "(Stat: Diff Sum)"] + [exp_datas[idx].get_stat(parameters[idx], "diff_sum") for idx in range(3)])
		writer.writerow(["", "(Stat: Diff Sq)"] + [exp_datas[idx].get_stat(parameters[idx], "diff_square") for idx in range(3)])
