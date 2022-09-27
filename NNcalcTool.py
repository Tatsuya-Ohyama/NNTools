#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
NNcalcTool.py - Program for calculation of stability by nearest neighbor method
"""

import sys, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
import csv
import os
import datetime

from classes.basicfunc import check_exist, check_overwrite
from classes.parameter import Parameter
from classes.sequence import Sequence



# =============== variable =============== #
VERSION = "1.0"
TEMPLATE_PARAM = "template_ref_param.csv"
TEMPLATE_SEQUENCE = "template_sequence.csv"



# =============== function =============== #
def make_template(flag_overwrite):
	"""
	create template files for ref_param.csv and ref_exp.csv
	"""
	if flag_overwrite == False:
		check_overwrite(TEMPLATE_PARAM)
	with open(TEMPLATE_PARAM, "w") as obj_output:
		writer = csv.writer(obj_output)
		writer.writerow(["Parameter", "dH", "dS", "dG"])
	sys.stderr.write("{0} is created.\n".format(TEMPLATE_PARAM))

	if flag_overwrite == False:
		check_overwrite(TEMPLATE_SEQUENCE)
	with open(TEMPLATE_SEQUENCE, "w") as obj_output:
		writer = csv.writer(obj_output)
		writer.writerow(["Label", "Sequence"])
	sys.stderr.write("{0} is created.\n".format(TEMPLATE_SEQUENCE))



# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="neighbor method", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-s", dest="FILE_SEQUENCE", metavar="SEQUENCE_FILE.csv", required="--make-template" not in sys.argv, help="sequence file of csv format")
	parser.add_argument("-p", dest="FILE_PARAMETER", metavar="PARAMETER_FILE.csv", required="--make-template" not in sys.argv, help="input file for parameters")
	parser.add_argument("-o", dest="FILE_OUTPUT", metavar="OUTPUT_FILE.csv", required="--make-template" not in sys.argv, help="output file")
	parser.add_argument("-l", dest="FILE_LOG", metavar="LOG_FILE.log", help="log file (if not specify, a log file with the same name as -o option is generated)")
	parser.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite forcibly")
	parser.add_argument("--make-template", dest="FLAG_MAKE_TEMPLATE", action="store_true", default=False, help="make template files ({0} and {1}) and exit".format(TEMPLATE_PARAM, TEMPLATE_SEQUENCE))
	args = parser.parse_args()

	if args.FLAG_MAKE_TEMPLATE:
		make_template(args.FLAG_OVERWRITE)
		sys.exit(0)

	check_exist(args.FILE_PARAMETER, 2)
	check_exist(args.FILE_SEQUENCE, 2)

	# reading parameters
	parameters = []
	base_pair = {}
	with open(args.FILE_PARAMETER, "r") as obj_input:
		flag_read = False
		flag_init = False
		reader = csv.reader(obj_input)
		for line_val in reader:
			if "Parameter" in line_val[0]:
				flag_read = True
				parameters = [Parameter().set_name(line_val[i]) for i in range(1, len(line_val))]
				continue

			if len(line_val[0]) == 0:
				break

			if flag_read:
				if "/" in line_val[0] and not line_val[0].startswith("init") and not line_val[0].startswith("length") and not line_val[0].startswith("symmetry") and not line_val[0].startswith("re:") and not line_val[0].startswith("reg:"):
					# generate base_pair
					base = line_val[0].split("/", 2)
					tmp_base_pair = {}
					tmp_base_pair[base[0][0:1]] = base[1][0:1]
					tmp_base_pair[base[0][1:2]] = base[1][1:2]
					for k, v in tmp_base_pair.items():
						if k in base_pair.keys():
							if base_pair[k] != v:
								sys.stderr.write("ERROR: base pair is duplicated: {0}-{1} vs {0}-{2}.\n".format(k, base_pair[k], tmp_base_pair[k]))
								sys.exit(1)
						else:
							base_pair[k] = v

				for param_idx in range(len(parameters)):
					parameters[param_idx].append_parameter(line_val[0], float(line_val[param_idx + 1].strip()))

	new_base_pair = {}
	for k, v in base_pair.items():
		new_base_pair[k] = v
		if v not in base_pair.keys():
			new_base_pair[v] = k
	base_pair = new_base_pair


	# reading sequence
	sequences = []
	with open(args.FILE_SEQUENCE, "r") as obj_input:
		reader = csv.reader(obj_input)
		flag_read = False
		for line_val in reader:
			if "Sequence" in line_val or "Sequences" in line_val:
				flag_read = True
				continue

			if flag_read:
				sequences.append(Sequence(line_val[0]).set_sequence(line_val[1], base_pair))

	# output
	if args.FILE_LOG is not None:
		path = os.path.splitext(args.FILE_OUTPUT)
		args.FILE_LOG = path[0] + ".log"

		if args.FLAG_OVERWRITE == False:
			check_overwrite(args.FILE_OUTPUT)
			check_overwrite(args.FILE_LOG)

		with open(args.FILE_LOG, "w") as obj_output:
			obj_output.write("This log file was generated at " + datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "\n\n")
			obj_output.write("<< Input >>\n")
			obj_output.write("Program version: {0}\n".format(VERSION))
			obj_output.write("Sequence file: {0} ({1} sequences)\n".format(args.FILE_SEQUENCE, len(sequences)))
			obj_output.write("Parameter file: {0}\n".format(args.FILE_PARAMETER))
			obj_output.write("|{0:^15}|{1[0]:^8}|{1[1]:^8}|{1[2]:^8}|\n".format("", [param.get_name() for param in parameters]))
			obj_output.write("|{0:-^15}|{0:-^8}|{0:-^8}|{0:-^8}|\n".format(""))
			for parameter_type in parameters[0].get_parameter(data_type="name"):
				obj_output.write("|{0:^15}|{1[0]:>8.3f}|{1[1]:>8.3f}|{1[2]:>8.3f}|\n".format(parameter_type, [param.get_parameter(parameter_type)[0] for param in parameters]))
			obj_output.write("\n")
			obj_output.write("Output file: {0}\n".format(args.FILE_OUTPUT))
			obj_output.write("\n")

	with open(args.FILE_OUTPUT, "w") as obj_output:
		writer = csv.writer(obj_output)
		parameter_types = list(parameters[0].get_parameter(data_type="name"))
		writer.writerow(["Comment", "Sequence"] + [param.name for param in parameters] + [""] + parameter_types)
		for sequence in sequences:
			freq = sequence.get_freq(parameters[0], base_pair)
			energy = [sequence.get_energy(param, base_pair) for param in parameters]
			writer.writerow([sequence.name, sequence.get_sequence("string")] + energy + [""] + [freq[param] for param in parameter_types])
