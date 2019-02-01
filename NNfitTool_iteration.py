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

from basic_func import check_exist, check_overwrite
from classes.Parameter import Parameter
from classes.Sequence import Sequence
from classes.DataGroup import DataGroup2


# =============== variable =============== #


# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "NNfitTool.py", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-i", dest = "input_file", metavar = "INPUT.csv", required = True, help = "sequence and experimental value file")
	parser.add_argument("-o", dest = "output_file", metavar = "OUTPUT.csv", required = True, help = "output file")
	parser.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly")
	parser.add_argument("-d", dest = "threshold_increment", metavar = "THRESHOLD", type = float, default = 0.00001, help = "difference threshold of increment for searching (Default: 0.00001)")
	parser.add_argument("-ii", dest = "initial_increment", metavar = "INITIAL_INCREMENT", type = float, default = 1.0, help = "initial increment (Default: 1.0)")
	args = parser.parse_args()

	check_exist(args.input_file, 2)

	# initial parameter
	exp_label = ["dH", "dS", "dG"]
	parameters = [Parameter(label) for label in exp_label]

	# reading sequence and experimental data
	exp_datas = [DataGroup2(label) for label in exp_label]
	with open(args.input_file, "r") as obj_input:
		reader = csv.reader(obj_input)

		# Ignore line number 1 (header) in CSV
		next(reader)

		for line_val in reader:
			sequence = Sequence(line_val[0], line_val[1])
			exp_datas[0].append(sequence, float(line_val[2]))
			exp_datas[1].append(sequence, float(line_val[3]))
			exp_datas[2].append(sequence, float(line_val[4]))


	# optimize parameter
	for exp_idx in range(3):
		# loop for energy type: dH, dS, and dG
		print("=============== Fitting {0} ===============".format(exp_label[exp_idx]))
		point_data = exp_datas[exp_idx]

		parameter_types = [parameter_type for parameter_type in parameters[exp_idx].get_parameter().keys()]
		parameters_opt = [Parameter(parameter_type) for parameter_type in parameter_types]
		increment = args.initial_increment
		direction = [0 for x in parameter_types]

		evaluation_val = [0.0 for x in parameter_types]
		evaluation_prev = [0.0 for x in parameter_types]
		cnt_i = 0
		while args.threshold_increment < increment:
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
					evaluation_val_tmp.append(exp_datas[exp_idx].get_stat(parameters_opt[parameter_idx], "diff_square"))
					evaluation_val_tmp.append(exp_datas[exp_idx].get_stat(parameter_plus, "diff_square"))
					evaluation_val_tmp.append(exp_datas[exp_idx].get_stat(parameter_minus, "diff_square"))
					evaluation_prev[parameter_idx] = exp_datas[exp_idx].get_stat(parameters_opt[parameter_idx], "diff_square")

					# determine direction
					min_val = min(evaluation_val_tmp)
					min_val_idx = [i for i, x in enumerate(evaluation_val_tmp) if min_val == x]
					if len(min_val_idx) != 1 or min_val_idx[0] == 0:
						# When all evaluation_val_tmp is the same even if parameter is changed
						# When evaluation_val_tmp value for base parameter is closest to 1, lock changing
						evaluation_val[parameter_idx] = exp_datas[exp_idx].get_stat(parameters_opt[parameter_idx], "diff_square")

					elif min_val_idx[0] == 1:
						# When evaluation for plus parameter is adopted
						direction[parameter_idx] = 1

					elif min_val_idx[0] == 2:
						# When evaluation for minus parameter is adopted
						direction[parameter_idx] = -1

					else:
						sys.stderr.write("ERROR: undefined condition.\n")
						sys.exit(1)

				parameter_new = parameters_opt[parameter_idx].clone()
				parameter_new.set_name(parameter_type + "_new")
				parameter_new.set_parameter("all", parameters_opt[parameter_idx].get_parameter())
				parameter_new.set_parameter(parameter_type, parameter_new.get_parameter(parameter_type) + direction[parameter_idx] * increment)

				# evaluation
				evaluation_val_tmp = []
				evaluation_val_tmp.append(exp_datas[exp_idx].get_stat(parameters_opt[parameter_idx], "diff_square"))
				evaluation_val_tmp.append(exp_datas[exp_idx].get_stat(parameter_new, "diff_square"))

				# calculate evaluation_val_tmp
				min_val = min(evaluation_val_tmp)
				min_val_idx = [i for i, x in enumerate(evaluation_val_tmp) if min_val == x]
				if len(min_val_idx) != 1 or min_val_idx[0] == 0:
					# When all factors of evaluation_val_tmp take the same value, lock changing
					# When evaluation_val_tmp for prev parameter is closest to 1, lock changing
					evaluation_val[parameter_idx] = evaluation_prev[parameter_idx]

				elif min_val_idx[0] == 1:
					# When r value for new parameter is closest to 1, update base parameter to low parameter
					parameter_new.set_name(parameter_type)
					parameters_opt[parameter_idx] = parameter_new
					evaluation_val[parameter_idx] = evaluation_val_tmp[1]

				else:
					sys.stderr.write("ERROR: undefined condition.\n")
					sys.exit(1)

			print("-" * 58)
			print("{0}     Iteration: {1} (dt = {2})".format(exp_label[exp_idx], cnt_i, increment))

			# evaluation for error
			evaluation_diff = [abs(x - y) for x, y in zip(evaluation_prev, evaluation_val)]
			max_val = max(evaluation_diff)

			if max_val != 0.0:
				max_val_idx = [i for i, v in enumerate(evaluation_diff) if v == max_val][0]

				# update_parameter
				new_parameter_type = parameters_opt[max_val_idx].get_name()
				new_parameter_val = parameters_opt[max_val_idx].get_parameter(new_parameter_type)
				for parameter in parameters_opt:
					parameter.set_parameter(new_parameter_type, copy.deepcopy(new_parameter_val))
				parameters[exp_idx] = parameters_opt[max_val_idx]

				print("{0:^8} {1:^10} {2:^10} {3:^10} {4:^10} {5:^5}".format("Type", "Parameter", "e1", "e2", "e_diff", "Adopt"))
				print("{0:-^8} {1:-^10} {2:-^10} {3:-^10} {4:-^10} {5:-^5}".format("", "", "", "", "", ""))
				for i, (p, e1, e2, e_diff) in enumerate(zip(parameter_types, evaluation_prev, evaluation_val, evaluation_diff)):
					if i == max_val_idx:
						print("{0:<8} {1:>10.3f} {2:>10.3f} {3:>10.3f} {4:>10.3f} {5:^5}".format(p, parameters[exp_idx].get_parameter(p), e1, e2, e_diff, "O"))
					else:
						print("{0:<8} {1:>10.3f} {2:>10.3f} {3:>10.3f} {4:>10.3f}".format(p, parameters[exp_idx].get_parameter(p), e1, e2, e_diff))
				print("")
				evaluation_prev = [evaluation_val[max_val_idx] for parameter in parameter_types]

			else:
				# When all parameters were locked, unlock and change increment
				increment /= 2

	print("")
	print("===== Last parameter =====")
	print("{0:^8} {1:^10}".format("Type", "Parameter"))
	print("{0:-^8} {1:-^10}".format("", ""))
	for p in parameter_types:
		print("{0:<8} {1:>10.3f}".format(p, parameters[exp_idx].get_parameter(p)))

	print("")
	print("===== Comparing experimental data =====")
	print("{0:^20} {1:^8} {2:^8} {3:^8}".format("Sequence", "Exp.", "Predict", "Diff"))
	print("{0:-^20} {1:-^8} {2:-^8} {3:-^8}".format("", "", "", ""))
	for row, diff in zip(exp_datas[2].get_energy(True, [parameters[2]]), exp_datas[2].get_stat(parameters[2], "diff_abs")):
		print("{0:<20} {1:>8.3f} {2:>8.3f} {3:>8.3f}".format(row[0], row[1], row[2], diff))

	print("")
	print("===== Curve fitting =====")
	print("Slope    :", exp_datas[exp_idx].get_stat(parameters[exp_idx], "slope"))
	print("Intercept:", exp_datas[exp_idx].get_stat(parameters[exp_idx], "intercept"))
	print("R   (1D) :", exp_datas[exp_idx].get_stat(parameters[exp_idx], "r"))
	print("R^2 (1D) :", exp_datas[exp_idx].get_stat(parameters[exp_idx], "r2"))
	print("E        :", exp_datas[exp_idx].get_stat(parameters[exp_idx], "diff_square"))

	# output
	if args.flag_overwrite == False:
		check_overwrite(args.output_file)

	with open(args.output_file, "w") as obj_output:
		writer = csv.writer(obj_output)
		writer.writerow(["", "dH", "dS", "dG"])
		parameters = [x.get_parameter() for x in parameters]
		for parameter_type in parameters[0].keys():
			writer.writerow([parameter_type, Decimal(str(parameters[0][parameter_type])).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP), Decimal(str(parameters[1][parameter_type])).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP), Decimal(str(parameters[1][parameter_type])).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP), Decimal(str(parameters[2][parameter_type])).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP)])
