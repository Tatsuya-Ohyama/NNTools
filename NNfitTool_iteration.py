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
from decimal import Decimal, ROUND_HALF_UP, ROUND_HALF_EVEN
from joblib import Parallel, delayed

from basic_func import check_exist, check_overwrite
from classes.Parameter import Parameter
from classes.Sequence import Sequence
from classes.DataGroup import DataGroup


# =============== variable =============== #
parameter_types = parameter_list = ["AA/TT", "AT/TA", "TA/AT", "CA/GT", "GT/CA", "CT/GA", "GA/CT", "CG/GC", "GC/CG", "GG/CC", "init_GC", "init_AT", "symmetry", "5term_TA"]


# =============== function =============== #
def calculation_worker(parameter, exp_data, mode, increment, threshold_increment, verbose, flag_thread):
	# calculate parameter
	# loop for energy type: dH, dS, and dG
	if not flag_thread and 1 <= verbose:
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
			if direction[parameter_idx] == False:
				# When not determine direction, determined direction
				# create 0,+,- changed parameter object
				parameter_plus = parameters_opt[parameter_idx].clone()
				parameter_plus.set_name(parameter_type + "_plus")
				parameter_plus.set_parameter(parameter_type, parameter_plus.get_parameter(parameter_type) + increment)
				parameter_minus = parameters_opt[parameter_idx].clone()
				parameter_minus.set_name(parameter_type + "_minus")
				parameter_minus.set_parameter(parameter_type, parameter_minus.get_parameter(parameter_type) - increment)

				# evaluation
				evaluation_val_tmp = []
				evaluation_val_tmp.append(exp_data.get_stat(parameters_opt[parameter_idx], mode))
				evaluation_val_tmp.append(exp_data.get_stat(parameter_plus, mode))
				evaluation_val_tmp.append(exp_data.get_stat(parameter_minus, mode))
				evaluation_prev[parameter_idx] = exp_data.get_stat(parameters_opt[parameter_idx], mode)

				# determine direction
				min_val = min(evaluation_val_tmp)
				min_val_idx = [i for i, x in enumerate(evaluation_val_tmp) if min_val == x]
				if len(min_val_idx) != 1 or min_val_idx[0] == 0:
					# When all evaluation_val_tmp is the same even if parameter is changed
					# When evaluation_val_tmp value for base parameter is closest to 1, lock changing
					evaluation_val[parameter_idx] = exp_data.get_stat(parameters_opt[parameter_idx], mode)

				elif min_val_idx[0] == 1:
					# When evaluation for plus parameter is adopted
					direction[parameter_idx] = 1

				elif min_val_idx[0] == 2:
					# When evaluation for minus parameter is adopted
					direction[parameter_idx] = -1

				else:
					sys.stderr.write("ERROR: undefined condition.\n")
					sys.exit(1)

			# parepare increased parameter
			parameter_new = parameters_opt[parameter_idx].clone()
			parameter_new.set_name(parameter_type + "_new")
			parameter_new.set_parameter(parameter_type, parameter_new.get_parameter(parameter_type) + direction[parameter_idx] * increment)

			# evaluation by mode
			evaluation_val_tmp = []
			evaluation_val_tmp.append(exp_data.get_stat(parameters_opt[parameter_idx], mode))
			evaluation_val_tmp.append(exp_data.get_stat(parameter_new, mode))

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
				sys.stderr.write("ERROR: undefined condition.\n")
				sys.exit(1)

		if not flag_thread and 2 <= verbose:
			print("-" * 53)
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

			if not flag_thread and 2 <= verbose:
				print("{0:^8} {1:^10} {2:^13} {3:^13} {4:^5}".format("Type", "Parameter", "e=sum(x'-x)^2", "Diff e (Prev)", "Adopt"))
				print("{0:-^8} {1:-^10} {2:-^13} {3:-^13} {4:-^5}".format("", "", "", "", "", ""))
				print("{0:<8} {1:>10} {2:>13.3f}".format("(Prev)", "", evaluation_prev[0]))
				for i, (p, e1, e2, e_diff) in enumerate(zip(parameter_types, evaluation_prev, evaluation_val, evaluation_diff)):
					if i == max_val_idx:
						print("{0:<8} {1:>10.3f} {2:>13.3f} {3:>13.3f} {4:^5}".format(p, parameter.get_parameter(p), e2, e_diff, "O"))
					else:
						print("{0:<8} {1:>10.3f} {2:>13.3f} {3:>13.3f}".format(p, parameter.get_parameter(p), e2, e_diff))
				print("")
			evaluation_prev = [evaluation_val[max_val_idx] for parameter in parameter_types]

		else:
			# When all parameters were locked, unlock and change increment
			increment /= 2
			direction = [0 for x in parameter_types]

	if not flag_thread and 1 <= verbose:
		print("")
		print("===== Last parameter =====")
		print("{0:^8} {1:^10}".format("Type", "Parameter"))
		print("{0:-^8} {1:-^10}".format("", ""))
		for p in parameter_types:
			print("{0:<8} {1:>10.3f}".format(p, parameter.get_parameter(p)))

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
	parser.add_argument("-i", dest = "input_file", metavar = "INPUT.csv", required = True, help = "sequence and experimental value file")
	parser.add_argument("-r", dest = "ref_param", metavar = "REF_PARAM.csv", help = "referenced parameter values")
	parser.add_argument("-o", dest = "output_file", metavar = "OUTPUT.csv", required = True, help = "output file")
	parser.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly")
	parser.add_argument("-d", dest = "threshold_increment", metavar = "THRESHOLD", type = float, default = 0.00001, help = "difference threshold of increment for searching (Default: 0.00001)")
	parser.add_argument("-ii", dest = "initial_increment", metavar = "INITIAL_INCREMENT", type = float, default = 0.01, help = "initial increment (Default: 1.0)")
	parser.add_argument("-T", dest = "temperature", metavar = "TEMPERATURE", type = float, default = 310.0, help = "temperature for experimental data (Default: 310.0)")
	parser.add_argument("-t", dest = "flag_thread", action = "store_true", default = False, help = "parallel calculation (Default: False)")
	parser.add_argument("-m", dest = "mode", metavar = "EVALUATION_METHOD", default = "diff_square", choices = ["r", "r2", "diff_mean", "diff_std", "diff_sum", "diff_square"], help = "evaluation method (r, r2, diff_mean, diff_std, diff_sum, diff_square) (Default: diff_square)")
	parser.add_argument("--verbose", "-v", dest = "verbose", action = "count", default = 0, help = "verbose (-v: display results / -vv: display calculation results)")
	args = parser.parse_args()

	# initial parameter
	exp_label = ["dH", "dS", "dG"]
	parameters = [Parameter(label) for label in exp_label]

	# loading reference parameter
	if args.ref_param is not None and check_exist(args.ref_param, 2, False):
		with open(args.ref_param, "r") as obj_input:
			reader = csv.reader(obj_input)
			for line_idx, line_val in enumerate(reader):
				if line_idx != 0:
					parameters[0].set_parameter(line_val[0], float(line_val[1]))
					parameters[1].set_parameter(line_val[0], float(line_val[2]))
					parameters[2].set_parameter(line_val[0], float(line_val[3]))
					if len(line_val) == 7:
						line_val[4:7] = [True if x.capitalize() == "True" else False for x in line_val[4:7]]
						parameters[0].set_change_stat(line_val[0], line_val[4])
						parameters[1].set_change_stat(line_val[0], line_val[5])
						parameters[2].set_change_stat(line_val[0], line_val[6])

	# reading sequence and experimental data
	check_exist(args.input_file, 2)
	exp_datas = [DataGroup(label) for label in exp_label]
	with open(args.input_file, "r") as obj_input:
		reader = csv.reader(obj_input)

		# Ignore line number 1 (header) in CSV
		next(reader)

		for line_val in reader:
			exp_datas[0].append(
				Sequence(line_val[0], line_val[1]).set_energy_type(exp_label[0]),
				float(line_val[2])
			)
			exp_datas[1].append(
				Sequence(line_val[0], line_val[1]).set_energy_type(exp_label[1]),
				float(line_val[3])
			)
			exp_datas[2].append(
				Sequence(line_val[0], line_val[1]).set_energy_type(exp_label[2]),
				float(line_val[4])
			)


	# optimize parameter for dH and dG
	if args.flag_thread:
		parameters = Parallel(n_jobs = 2)([
			delayed(calculation_worker)(
				parameters[exp_idx],
				exp_datas[exp_idx],
				args.mode,
				args.initial_increment,
				args.threshold_increment,
				args.verbose, args.flag_thread
			) for exp_idx in [0, 2]])
	else:
		for exp_idx in [0, 2]:
			parameters[exp_idx] = calculation_worker(
				parameters[exp_idx],
				exp_datas[exp_idx],
				args.mode,
				args.initial_increment,
				args.threshold_increment,
				args.verbose, args.flag_thread
			)

	# calculate dS
	for parameter_type in parameter_types:
		dS = (parameters[0].get_parameter(parameter_type) - parameters[2].get_parameter(parameter_type)) / args.temperature * 1000
		parameters[1].set_parameter(parameter_type, dS)

	# output
	if args.flag_overwrite == False:
		check_overwrite(args.output_file)

	with open(args.output_file, "w") as obj_output:
		writer = csv.writer(obj_output)
		writer.writerow(["<< Input >>"])
		writer.writerow(["Input", args.input_file])
		writer.writerow(["Initial iteration", args.initial_increment])
		writer.writerow(["Increment threshold", args.threshold_increment])
		writer.writerow([""])
		writer.writerow([""])

		writer.writerow(["<< Parameter >>"])
		writer.writerow(["", "dH", "dS", "dG", "", "Change (dH)", "Change (dS)", "Change (dG)"])
		for parameter_type in parameter_types:
			writer.writerow([
				parameter_type,
				Decimal(str(parameters[0].get_parameter()[parameter_type])).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP),
				Decimal(str(parameters[1].get_parameter()[parameter_type])).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP),
				Decimal(str(parameters[2].get_parameter()[parameter_type])).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP),
				"",
				parameters[0].is_change(parameter_type),
				parameters[1].is_change(parameter_type),
				parameters[2].is_change(parameter_type)
				])
		writer.writerow([""])
		writer.writerow(["(Stat: R)"] + [exp_datas[idx].get_stat(parameters[idx], "r") for idx in range(3)])
		writer.writerow(["(Stat: R2)"] + [exp_datas[idx].get_stat(parameters[idx], "r2") for idx in range(3)])
		writer.writerow(["(Stat: Slope)"] + [exp_datas[idx].get_stat(parameters[idx], "slope") for idx in range(3)])
		writer.writerow(["(Stat: Intercept)"] + [exp_datas[idx].get_stat(parameters[idx], "intercept") for idx in range(3)])
		writer.writerow(["(Stat: Diff Mean)"] + [exp_datas[idx].get_stat(parameters[idx], "diff_mean") for idx in range(3)])
		writer.writerow(["(Stat: Diff Sum)"] + [exp_datas[idx].get_stat(parameters[idx], "diff_sum") for idx in range(3)])
		writer.writerow(["(Stat: Diff Sq)"] + [exp_datas[idx].get_stat(parameters[idx], "diff_square") for idx in range(3)])
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
		datas = []
		datas.append(exp_datas[0].get_energy(flag_sequence = True))
		datas.append(exp_datas[1].get_energy(flag_sequence = True))
		datas.append(exp_datas[2].get_energy(flag_sequence = True))
		datas.append(exp_datas[0].get_energy(flag_sequence = True, obj_parameters = [parameters[0]]))
		datas.append(exp_datas[1].get_energy(flag_sequence = True, obj_parameters = [parameters[1]]))
		datas.append(exp_datas[2].get_energy(flag_sequence = True, obj_parameters = [parameters[2]]))

		for sequence in exp_datas[0].get_sequence():
			name = sequence.get_name()
			seq = sequence.get_sequence("string")
			exp_dH = [x[1] for x in datas[0] if x[0] == seq][0]
			exp_dS = [x[1] for x in datas[1] if x[0] == seq][0]
			exp_dG = [x[1] for x in datas[2] if x[0] == seq][0]
			dH = [x[2] for x in datas[3] if x[0] == seq][0]
			dS = [x[2] for x in datas[4] if x[0] == seq][0]
			dG = [x[2] for x in datas[5] if x[0] == seq][0]
			diff_dH = dH - exp_dH
			diff_dS = dS - exp_dS
			diff_dG = dG - exp_dG
			freq = [sequence.get_freq()[parameter_type] for parameter_type in parameter_types]
			writer.writerow([name, seq, exp_dH, exp_dS, exp_dG, dH, dS, dG, diff_dH, diff_dS, diff_dG] + [""] + freq)
