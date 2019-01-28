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
	parser.add_argument("-r", dest = "parameter_file", metavar = "PARAMETER.csv", required = True, help = "initial parameter file")
	parser.add_argument("-o", dest = "output_file", metavar = "OUTPUT.csv", required = True, help = "output file")
	parser.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly")
	parser.add_argument("-d", dest = "threshold_increment", metavar = "THRESHOLD", type = float, default = 0.00001, help = "difference threshold of increment for searching (Default: 0.00001)")
	parser.add_argument("-ii", dest = "initial_increment", metavar = "INITIAL_INCREMENT", type = float, default = 0.01, help = "initial increment (Default: 0.01)")
	args = parser.parse_args()

	check_exist(args.input_file, 2)
	check_exist(args.parameter_file, 2)

	# reading parameter
	parameters = [Parameter("dH"), Parameter("dS"), Parameter("dG"), Parameter("dH_opt"), Parameter("dS_opt"), Parameter("dG_opt")]
	with open(args.parameter_file, "r") as obj_input:
		reader = csv.reader(obj_input)

		# Ignore line number 1 (header) in CSV
		next(reader)

		for line_val in reader:
			parameters[0].set_parameter(line_val[0], line_val[1])
			parameters[1].set_parameter(line_val[0], line_val[2])
			parameters[2].set_parameter(line_val[0], line_val[3])
			parameters[3].set_parameter(line_val[0], line_val[1])
			parameters[4].set_parameter(line_val[0], line_val[2])
			parameters[5].set_parameter(line_val[0], line_val[3])

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


	# import random
	# parameter_types = ["AA/TT", "AT/TA", "TA/AT", "CA/GT", "GT/CA", "CT/GA", "GA/CT", "CG/GC", "GC/CG", "GG/CC", "init_GC", "init_AT", "symmetry", "5term_TA"]
	# random.shuffle(parameter_types)
	label_base = ["Sequence", "Exp."]
	for idx in range(3):
		flag_first = True
		diff_r2 = [args.threshold_increment + 1 for i in range(3)]
		parameter_base = parameters[idx + 3]
		parameter_high = None
		parameter_low = None
		point_data = point_datas[idx]

		for parameter_type in parameter_base.get_parameter().keys():
		# for parameter_type in parameter_types:
			# optimize for each parameter
			increment = args.initial_increment
			diff_r2_prev = 1.0
			print("=" * 30, "PARAMETER", parameter_type, "=" * 30)
			direction = 0
			while args.threshold_increment < increment:
				# loop while r2 is larger than threshold

				if direction == 0:
					# calculate base parameter
					point_data.add_data("column", ["Base"], [sequence.set_parameter(parameter_base).get_energy() for sequence in sequences], np.float)

				if 0 <= direction <= 1:
					# calculate high parameter
					parameter_high = Parameter("High").set_parameter("all", {k: v + increment if k == parameter_type else v for k, v in parameter_base.get_parameter().items()})
					point_data.add_data("column", ["High"], [sequence.set_parameter(parameter_high).get_energy() for sequence in sequences], np.float)

				if -1 <= direction <= 0:
					# calculate low parameter
					parameter_low = Parameter("Low").set_parameter("all", {k: v - increment if k == parameter_type else v for k, v in parameter_base.get_parameter().items()})
					point_data.add_data("column", ["Low"], [sequence.set_parameter(parameter_low).get_energy() for sequence in sequences], np.float)

				r2 = [point_data.get_factor("label", "Exp.", "Base")[3], point_data.get_factor("label", "Exp.", "High")[3], point_data.get_factor("label", "Exp.", "Low")[3]]
				diff_r2 = [abs(x - 1.0) for x in r2]
				print("parameter", parameter_base.get_parameter()[parameter_type], parameter_high.get_parameter()[parameter_type], parameter_low.get_parameter()[parameter_type])

				print("diff_r2", diff_r2)
				if min(diff_r2) == diff_r2[0]:
					# When diff_r2 value for base parameter is closest to 1, change increment
					increment /= 2
					direction = 0
					point_data.remove_data("column", "label", "Base")
					point_data.remove_data("column", "label", "High")
					point_data.remove_data("column", "label", "Low")
					continue
					print("even")
				elif min(diff_r2) == diff_r2[1]:
					# When diff_r2 value for high parameter is closest to 1, update base parameter to high parameter
					parameter_base = parameter_high
					diff_r2[2] = diff_r2[0]
					diff_r2[0] = diff_r2[1]
					direction = 1
					point_data.remove_data("column", "label", "Low")
					point_data.set_label("column", ["Sequence", "Exp.", "Low", "Base"])
					print("high")
				elif min(diff_r2) == diff_r2[2]:
					# When diff_r2 value for low parameter is closest to 1, update base parameter to low parameter
					parameter_base = parameter_low
					diff_r2[1] = diff_r2[0]
					diff_r2[0] = diff_r2[2]
					point_data.remove_data("column", "label", "High")
					point_data.set_label("column", ["Sequence", "Exp.", "High", "Base"])
					direction = -1
					print("low")
				else:
					# When all diff_r2 is the same even if parameter is changed
					sys.stderr.write("R2 value converged.\n")
					print("r2", point_data.get_factor("label", "Exp.", "Base"))
					point_data.remove_data("column", "label", "Base")
					point_data.remove_data("column", "label", "High")
					point_data.remove_data("column", "label", "Low")
					direction = 0
					break
				point_data.save_csv("tmp.csv")

				diff_r2_prev = min(diff_r2)

		parameters[idx + 3] = parameter_base



	if args.flag_overwrite == False:
		check_overwrite(args.output_file)

	with open(args.output_file, "w") as obj_output:
		writer = csv.writer(obj_output)
		writer.writerow(["", "dH (ref.)", "dH (opt.)", "dS (ref.)", "dS (opt.)", "dG (ref.)", "dG (opt.)"])
		parameters = [x.get_parameter() for x in parameters]
		for parameter_type in parameters[0].keys():
			writer.writerow([parameter_type, Decimal(str(parameters[0][parameter_type])).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP), Decimal(str(parameters[3][parameter_type])).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP), Decimal(str(parameters[1][parameter_type])).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP), Decimal(str(parameters[4][parameter_type])).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP), Decimal(str(parameters[2][parameter_type])).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP), Decimal(str(parameters[5][parameter_type])).quantize(Decimal('0.001'), rounding = ROUND_HALF_UP)])
