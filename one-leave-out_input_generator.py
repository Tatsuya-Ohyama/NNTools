#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
one-leave-out_input_generator.py
Program to generate input file for one leave out method
"""

import sys, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
import os
import csv

from basic_func import check_exist, check_overwrite


# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "Program to generate input file for one leave out method", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-i", dest = "INPUT_FILE", metavar = "INPUT.csv", required = True, help = "original input file for NNfitTool.py")
	parser.add_argument("-O", dest = "FLAG_OVERWRITE", action = "store_true", default = False, help = "overwrite forcibly")
	args = parser.parse_args()

	check_exist(args.INPUT_FILE, 2)

	label = []
	values = []
	with open(args.INPUT_FILE, "r") as obj_input:
		reader = csv.reader(obj_input)
		for row_idx, row_val in enumerate(reader, 1):
			if row_idx == 1:
				label.extend(row_val)
			else:
				tmp_values = row_val[:2] + [float(v) for v in row_val[2:] if len(v) != 0]
				values.append(tmp_values)

	digit = len(str(len(values)))
	output_prefix = os.path.splitext(args.INPUT_FILE)[0]
	output_format = output_prefix + "_wo_{0:0" + str(digit) + "}.csv"
	for i in range(len(values)):
		output_file = output_format.format(i + 1)

		if args.FLAG_OVERWRITE == False:
			check_overwrite(output_file)

		with open(output_file, "w") as obj_output:
			writer = csv.writer(obj_output)
			writer.writerow(label)
			writer.writerows(values[0 : i] + values[i + 1:])

		sys.stderr.write("Create: {0}\n".format(output_file))
