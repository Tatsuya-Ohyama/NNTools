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
from classes.DataGroup import DataGroup


# =============== variable =============== #


# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "NNfitTool.py", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-i", dest = "input_file", metavar = "INPUT.csv", required = True, help = "sequence and experimental value file")
	parser.add_argument("-o", dest = "output_file", metavar = "OUTPUT.csv", required = True, help = "output file")
	parser.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly")
	parser.add_argument("-d", dest = "threshold_increment", metavar = "THRESHOLD", type = float, default = 0.00001, help = "difference threshold of increment for searching (Default: 0.00001)")
	parser.add_argument("-ii", dest = "initial_increment", metavar = "INITIAL_INCREMENT", type = float, default = 0.01, help = "initial increment (Default: 0.01)")
	args = parser.parse_args()

	check_exist(args.input_file, 2)

	# initial parameter
	parameters = [Parameter("dH"), Parameter("dS"), Parameter("dG"), Parameter("dH_opt"), Parameter("dS_opt"), Parameter("dG_opt")]

	# reading sequence and experimental data
	sequences = []
	point_datas = [DataGroup(), DataGroup(), DataGroup()]
	with open(args.input_file, "r") as obj_input:
		reader = csv.reader(obj_input)

		# Ignore line number 1 (header) in CSV
		next(reader)

		for line_val in reader:
			sequence = Sequence(line_val[2], line_val[0])
			sequences.append(sequence)
			point_datas[0].add_data("row", [sequence.get_name()], ["".join(sequence.get_sequence()), line_val[1]])

	point_datas[0].set_label("column", ["Sequence", "Exp."])
	point_datas[0].set_dtype(np.float, "label", "Exp.")
	point_datas[1] = copy.deepcopy(point_datas[0])
	point_datas[2] = copy.deepcopy(point_datas[0])


	# optimize parameter
	label_base = ["Sequence", "Exp."]
	for energy_idx in range(3):
		# loop for energy type: dH, dS, and dG
		point_data = point_datas[energy_idx]

		parameter_types = [parameter_type for parameter_type in parameters[energy_idx + 3].get_parameter().keys()]
		candidate_parameters = [Parameter(parameter_type) for parameter_type in parameter_types]
		flag_change = [True for x in parameter_types]
		increment = args.initial_increment
		direction = [False for x in parameter_types]

		while :
			# loop while r2 or increment value .....
			r2_eval = []
			for parameter_idx, parameter_type in enumerate(parameter_types):
				# loop for parameters

				if flag_change[parameter_idx]:
					# for parameter changed with current increment

					# change parameter
					if direction[parameter_idx] == False:
						# not determine direction
						parameter_base = Parameter(parameter_type + "_base")
						parameter_base.set_parameter({k: v if k == parameter.get_name() else v for k, v in parameter_base.get_parameter().items()})
						parameter_plus = Parameter(parameter_type + "_plus")
						parameter_plus.set_parameter({k: v + increment if k == parameter.get_name() else v for k, v in parameter_base.get_parameter().items()})
						parameter_minus = Parameter(parameter_type + "_minus")
						parameter_minus.set_parameter({k: v - increment if k == parameter.get_name() else v for k, v in parameter_base.get_parameter().items()})

						# add energy
						point_data.add_data("column", [parameter_type + "_base"], [sequence.set_parameter(parameter_plus).get_energy() for sequence in sequences], np.float)
						point_data.add_data("column", [parameter_type + "_plus"], [sequence.set_parameter(parameter_plus).get_energy() for sequence in sequences], np.float)
						point_data.add_data("column", [parameter_type + "_minus"], [sequence.set_parameter(parameter_minus).get_energy() for sequence in sequences], np.float)

						# calculate r2 & determine direction
						diff_r2 = []
						diff_r2.append(abs(1 - point_data.get_factor("Exp.", parameter_type + "_base")[3]))
						diff_r2.append(abs(1 - point_data.get_factor("Exp.", parameter_type + "_plus")[3]))
						diff_r2.append(abs(1 - point_data.get_factor("Exp.", parameter_type + "_minus")[3]))
						min_idx = [i for i, x in enumerate(diff_r2) if min(diff_r2) == x]
						if len(min_idx) != 1 or min_idx[0] == 0:
							# When all diff_r2 is the same even if parameter is changed
							# When diff_r2 value for base parameter is closest to 1, change increment
							flag_change[parameter_idx] = False
							r2_eval.append(diff_r2[0])


						elif min_idx[0] == 1:
							# When diff_r2 value for high parameter is closest to 1, update base parameter to high parameter
							increment[parameter_idx] *= 1
							candidate_parameters[parameter_idx] = parameter_high
							point_data.remove_data(parameter_type + "_minus")
							point_data.set_label("column", [parameter_type for x in point_data.get_label("column") if x == parameter_type + "_plus"])
							direction[parameter_idx] = True
							r2_eval.append(diff_r2[1])

						elif min_idx[0] == 2:
							# When diff_r2 value for low parameter is closest to 1, update base parameter to low parameter
							increment[parameter_idx] *= -1
							candidate_parameters[parameter_idx] = parameter_low
							point_data.remove_data(parameter_type + "_plus")
							point_data.set_label("column", [parameter_type for x in point_data.get_label("column") if x == parameter_type + "_minus"])
							direction[parameter_idx] = True
							r2_eval.append(diff_r2[2])

						else:
							sys.stderr.write("ERROR: undefined condition.\n")
							sys.exit(1)

					else:
						"""
						ここの実装
						"""
						# determined direction
						parameter = Parameter(parameter_type + "_new")
						parameter.set_parameter({k: v if k == parameter.get_name() else v for k, v in parameters[parameter_idx].get_parameter().items()})

						# add energy
						point_data.add_data("column", [parameter_type + "_base"], [sequence.set_parameter(parameter_plus).get_energy() for sequence in sequences], np.float)

						# calculate r2 & determine direction
						diff_r2 = []
						diff_r2.append(abs(1 - point_data.get_factor("Exp.", parameter_type + "_base")[3]))
						diff_r2.append(abs(1 - point_data.get_factor("Exp.", parameter_type + "_plus")[3]))
						diff_r2.append(abs(1 - point_data.get_factor("Exp.", parameter_type + "_minus")[3]))
						min_idx = [i for i, x in enumerate(diff_r2) if min(diff_r2) == x]
						if len(min_idx) != 1 or min_idx[0] == 0:
							# When all diff_r2 is the same even if parameter is changed
							# When diff_r2 value for base parameter is closest to 1, change increment
							flag_change[parameter_idx] = False


						elif min_idx[0] == 1:
							# When diff_r2 value for high parameter is closest to 1, update base parameter to high parameter
							increment[parameter_idx] *= 1
							candidate_parameters[parameter_idx] = parameter_high
							point_data.remove_data(parameter_type + "_minus")
							point_data.set_label("column", [parameter_type for x in point_data.get_label("column") if x == parameter_type + "_plus"])
							direction[parameter_idx] = True

						elif min_idx[0] == 2:
							# When diff_r2 value for low parameter is closest to 1, update base parameter to low parameter
							increment[parameter_idx] *= -1
							candidate_parameters[parameter_idx] = parameter_low
							point_data.remove_data(parameter_type + "_plus")
							point_data.set_label("column", [parameter_type for x in point_data.get_label("column") if x == parameter_type + "_minus"])
							direction[parameter_idx] = True

						else:
							sys.stderr.write("ERROR: undefined condition.\n")
							sys.exit(1)


			if len([True for x in flag_change if True]) == 0:
				# When all parameter did not changed, increment change
				increment /= 2
				flag_change = [True for x in flag_change]







	if args.flag_overwrite == False:
		check_overwrite(args.output_file)

	with open(args.output_file, "w") as obj_output:
		writer = csv.writer(obj_output)
		writer.writerow(["", "dH (ref.)", "dH (opt.)", "dS (ref.)", "dS (opt.)", "dG (ref.)", "dG (opt.)"])
		parameters = [x.get_parameter() for x in parameters]
		for parameter_type in parameters[0].keys():
			writer.writerow([parameter_type, Decimal(str(parameters[0][parameter_type])).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP), Decimal(str(parameters[3][parameter_type])).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP), Decimal(str(parameters[1][parameter_type])).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP), Decimal(str(parameters[4][parameter_type])).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP), Decimal(str(parameters[2][parameter_type])).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP), Decimal(str(parameters[5][parameter_type])).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP)])
