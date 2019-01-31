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
	parser.add_argument("-ii", dest = "initial_increment", metavar = "INITIAL_INCREMENT", type = float, default = 1.0, help = "initial increment (Default: 1.0)")
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

		for line_eval in reader:
			sequence = Sequence(line_eval[2], line_eval[0])
			sequences.append(sequence)
			point_datas[0].add_data("row", [sequence.get_name()], ["".join(sequence.get_sequence()), line_eval[1]])

	point_datas[0].set_label("column", ["Sequence", "Exp."])
	point_datas[0].set_dtype(np.float, "label", "Exp.")
	point_datas[1] = copy.deepcopy(point_datas[0])
	point_datas[2] = copy.deepcopy(point_datas[0])


	# optimize parameter
	energy_idx = 2
	if energy_idx == 2:
	# for energy_idx in range(3):
		# loop for energy type: dH, dS, and dG
		point_data = point_datas[energy_idx]

		parameter_types = [parameter_type for parameter_type in parameters[energy_idx + 3].get_parameter().keys()]
		candidate_parameters = [Parameter(parameter_type) for parameter_type in parameter_types]
		flag_change = [True for x in parameter_types]
		increment = args.initial_increment
		direction = [0 for x in parameter_types]
		flag_first = True

		e_eval = []
		e_prev = []
		while args.threshold_increment < increment:
			# loop while increment is larger than threshold

			# add initial parameter to point data

			for parameter_idx, parameter_type in enumerate(parameter_types):
				# loop for parameters

				if flag_first:
					point_data.add_data("column", [parameter_type], [sequence.set_parameter(candidate_parameters[parameter_idx]).get_energy() for sequence in sequences], np.float)
					e_eval = [0.0 for x in parameter_types]

				if flag_change[parameter_idx]:
					# for parameter changed with current increment

					# change parameter
					if direction[parameter_idx] == False:
						# not determine direction

						# create +/- increment parameter object
						parameter_plus = Parameter(parameter_type + "_plus")
						parameter_plus.set_parameter("all", candidate_parameters[parameter_idx].get_parameter())
						parameter_plus.set_parameter(parameter_type, parameter_plus.get_parameter(parameter_type) + increment)
						parameter_minus = Parameter(parameter_type + "_minus")
						parameter_minus.set_parameter("all", candidate_parameters[parameter_idx].get_parameter())
						parameter_minus.set_parameter(parameter_type, parameter_minus.get_parameter(parameter_type) - increment)

						# add energy
						point_data.add_data("column", [parameter_type + "_plus"], [sequence.set_parameter(parameter_plus).get_energy() for sequence in sequences], np.float)
						point_data.add_data("column", [parameter_type + "_minus"], [sequence.set_parameter(parameter_minus).get_energy() for sequence in sequences], np.float)

						# calculate error
						e = [point_data.get_diff("label", "Exp.", parameter_type)[1], point_data.get_diff("label", "Exp.", parameter_type + "_plus")[1], point_data.get_diff("label", "Exp.", parameter_type + "_minus")[1]]
						print(100,e)
						e_prev = e
						min_e = min(e)
						min_e_idx = [i for i, x in enumerate(e) if min_e == x]
						if len(min_e_idx) != 1 or min_e_idx[0] == 0:
							# When all e is the same even if parameter is changed
							# When e value for base parameter is closest to 1, lock changing
							flag_change[parameter_idx] = False
							e_eval[parameter_idx] = e[0]
							point_data.remove_data("column", "label", parameter_type + "_plus")
							point_data.remove_data("column", "label", parameter_type + "_minus")

						elif min_e_idx[0] == 1:
							# When r value for high parameter is closest to 1, update base parameter to high parameter
							parameter_plus.set_name(parameter_type)
							candidate_parameters[parameter_idx] = parameter_plus
							point_data.remove_data("column", "label", parameter_type + "_minus")
							point_data.remove_data("column", "label", parameter_type)
							point_data.rename_label("column", parameter_type + "_plus", parameter_type)
							direction[parameter_idx] = 1
							e_eval[parameter_idx] = e[1]
							print("shift plus")

						elif min_e_idx[0] == 2:
							# When r value for low parameter is closest to 1, update base parameter to low parameter
							parameter_minus.set_name(parameter_type)
							candidate_parameters[parameter_idx] = parameter_minus
							point_data.remove_data("column", "label", parameter_type + "_plus")
							point_data.remove_data("column", "label", parameter_type)
							point_data.rename_label("column", parameter_type + "_minus", parameter_type)
							direction[parameter_idx] = -1
							e_eval[parameter_idx] = e[2]
							print("shift minus")

						else:
							sys.stderr.write("ERROR: undefined condition.\n")
							sys.exit(1)

					else:
						# determined direction
						# change parameter
						parameter_new = Parameter(parameter_type + "_new")
						parameter_new.set_parameter("all", {k: v + direction[parameter_idx] * increment if k == candidate_parameters[parameter_idx].get_name() else v for k, v in candidate_parameters[parameter_idx].get_parameter().items()})

						# add energy
						point_data.add_data("column", [parameter_type + "_new"], [sequence.set_parameter(parameter_new).get_energy() for sequence in sequences], np.float)

						# calculate error
						e = [point_data.get_diff("label", "Exp.", parameter_type)[1], point_data.get_diff("label", "Exp.", parameter_type + "_new")[1]]
						print(101, e)
						e_prev = e
						min_e = min(e)
						min_e_idx = [i for i, x in enumerate(e) if min_e == x]
						if len(min_e_idx) != 1 or min_e_idx[0] == 0:
							# When r value for base parameter is closest to 1, change increment
							# When r value for prev parameter is closest to 1, lock changing
							flag_change[parameter_idx] = False
							point_data.remove_data("column", "label", parameter_type + "_new")
							print("fail")

						elif min_e_idx[0] == 1:
							# When r value for new parameter is closest to 1, update base parameter to low parameter
							print("old", parameter.get_parameter())
							parameter_new.set_name(parameter_type)
							candidate_parameters[parameter_idx] = parameter_new
							point_data.remove_data("column", "label", parameter_type)
							point_data.rename_label("column", parameter_type + "_new", parameter_type)
							e_eval[parameter_idx] = e[1]
							print("new", parameter_new.get_parameter())
							print("update")

						else:
							sys.stderr.write("ERROR: undefined condition.\n")
							sys.exit(1)
				# print(point_data.get_data())
				# print(e_prev, e)
				# print(parameter_type)
				sys.stdin.readline()

			# evaluation for error
			print(e_prev, e)
			min_e = max(e_eval)
			min_e_idx = [i for i, v in enumerate(e_eval) if v == min_e][0]

			# update_parameter
			new_parameter_type = candidate_parameters[min_e_idx].get_name()
			new_parameter_val = candidate_parameters[min_e_idx].get_parameter(new_parameter_type)
			for parameter in candidate_parameters:
				parameter.set_parameter(new_parameter_type, new_parameter_val)
			parameters[energy_idx + 3] = candidate_parameters[min_e_idx]
			flag_first = False

			print(parameters[energy_idx + 3].get_parameter())
			print("e", e_eval[min_e_idx])
			for i, (p, r) in enumerate(zip(parameter_types, e_eval)):
				if i == min_e_idx:
					print(p, r, "O")
				else:
					print(p, r)

			if flag_change.count(False) == len(parameter_types):
				# When all parameters were locked, unlock and change increment
				increment /= 2
				flag_change = [True for x in flag_change]
				print("next_step", increment)
				print("=" * 50)
				sys.stdin.readline()

	point_datas[2].save_csv("point.csv")




	if args.flag_overwrite == False:
		check_overwrite(args.output_file)

	with open(args.output_file, "w") as obj_output:
		writer = csv.writer(obj_output)
		writer.writerow(["", "dH (ref.)", "dH (opt.)", "dS (ref.)", "dS (opt.)", "dG (ref.)", "dG (opt.)"])
		parameters = [x.get_parameter() for x in parameters]
		for parameter_type in parameters[0].keys():
			writer.writerow([parameter_type, Decimal(str(parameters[0][parameter_type])).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP), Decimal(str(parameters[3][parameter_type])).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP), Decimal(str(parameters[1][parameter_type])).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP), Decimal(str(parameters[4][parameter_type])).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP), Decimal(str(parameters[2][parameter_type])).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP), Decimal(str(parameters[5][parameter_type])).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP)])
