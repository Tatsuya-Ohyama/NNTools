#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
NNscanTool.py - Program for searching stable area by nearest neighbor method
"""

import sys, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
import csv
import os
import datetime
import matplotlib.pyplot as plt
import numpy as np

from mods.basicfunc import check_exist, check_overwrite
from mods.parameter import Parameter
from mods.sequence import Sequence
from mods.graph import Graph2D, matplot_init



# =============== variable =============== #
VERSION = "1.4"
TEMPLATE_PARAM = "template_ref_param.csv"
TEMPLATE_SEQUENCE = "template_sequence.csv"



# =============== function =============== #
def make_template(flag_overwrite):
	"""
	Function to create template files for ref_param.csv and ref_exp.csv

	Args:
		flag_overwrite (bool): overwrite forcibly
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


def read_sequence_txt(input_file):
	"""
	Function to read sequence from txt format file

	Args:
		input_file (str): input file of txt format

	Returns:
		sequence (str)
	"""
	sequence = ""
	with open(input_file, "r") as obj_input:
		sequence = "".join([line_val.strip() for line_val in obj_input.readlines()])
	return sequence.upper()


def read_sequence_genbank(input_file):
	"""
	Function to read sequence from genbank format file

	Args:
		input_file (str): input file of genbank format file

	Returns:
		sequence (str)
	"""
	flag_read = False
	sequence = ""
	with open(input_file, "r") as obj_input:
		for line_val in obj_input:
			if line_val.startswith("ORIGIN"):
				flag_read = True
			elif line_val.startswith("//"):
				break
			elif flag_read:
				sequence += line_val[10:75].strip().replace(" ", "")
	return sequence.upper()


def read_parameters(input_file):
	"""
	Function to read parameters from csv file

	Args:
		input_file (str): parameter file
	"""
	parameters = []
	base_pair = {}
	with open(input_file, "r") as obj_input:
		flag_read = False
		flag_init = False
		reader = csv.reader(obj_input)
		for line_val in reader:
			if "Parameter" in line_val[0]:
				flag_read = True
				parameters = [Parameter().set_name(line_val[i]) for i in range(1, len(line_val))]
				continue

			if flag_read:
				if "/" in line_val[0] and not line_val[0].startswith("init") and not line_val[0].startswith("length") and not line_val[0].startswith("symmetry") and not line_val[0].startswith("re:"):
					# generate base_pair
					base = line_val[0].split("/", 2)
					tmp_base_pair = {}
					tmp_base_pair[base[0][0:1]] = base[1][0:1]
					tmp_base_pair[base[0][1:2]] = base[1][1:2]
					for k, v in tmp_base_pair.items():
						if k in base_pair.keys():
							if base_pair[k] != v:
								sys.stderr.write("ERROR: base pair are duplicated: {0}-{1} vs {0}-{2}.\n".format(k, base_pair[k], tmp_base_pair[k]))
								sys.exit(1)
						else:
							base_pair[k] = v

				for param_idx in range(len(parameters)):
					# register parameter
					parameters[param_idx].append_parameter(line_val[0], float(line_val[param_idx+1].strip()))

	new_base_pair = {}
	for k, v in base_pair.items():
		new_base_pair[k] = v
		if v not in base_pair.keys():
			new_base_pair[v] = k
	base_pair = new_base_pair

	return parameters, base_pair


# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="NNscanTool.py - Program for searching stable area by nearest neighbor method", formatter_class=argparse.RawTextHelpFormatter)
	sequence_group = parser.add_mutually_exclusive_group(required="--make-template" not in sys.argv)
	sequence_group.add_argument("-s", dest="FILE_SEQUENCE_TXT", metavar="SEQUENCE_FILE.txt", help="sequence file of txt format")
	sequence_group.add_argument("-g", dest="FILE_SEQUENCE_GEN", metavar="SEQUENCE_FILE.gb", help="sequence file of genbank format")
	parser.add_argument("-p", dest="FILE_PARAMETER", metavar="PARAMETER_FILE.csv", required="--make-template" not in sys.argv, help="input file for parameters")
	parser.add_argument("-b", dest="BLOCK_LENGTH", metavar="SEQUENCE_RANGE", type=int, help="block length for sequence (Default: full length)")
	parser.add_argument("-o", dest="FILE_OUTPUT", metavar="OUTPUT_FILE.csv", required="--make-template" not in sys.argv, help="output file")
	parser.add_argument("-l", dest="FILE_LOG", metavar="LOG_FILE.log", help="log file (if not specify, a log file with the same name as -o option is generated)")
	parser.add_argument("--graph", dest="FILE_GRAPH", metavar="GRAPH_BASE_NAME", help="base name of graph (.png)")
	parser.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite forcibly")
	parser.add_argument("--make-template", dest="FLAG_MAKE_TEMPLATE", action="store_true", default=False, help="make template files ({0} and {1}) and exit".format(TEMPLATE_PARAM, TEMPLATE_SEQUENCE))
	parser.add_argument("--version", action="version", version=VERSION)
	args = parser.parse_args()

	if args.FLAG_MAKE_TEMPLATE:
		make_template(args.FLAG_OVERWRITE)
		sys.exit(0)

	check_exist(args.FILE_PARAMETER, 2)

	# read sequence
	str_sequence = ""
	if args.FILE_SEQUENCE_TXT is not None:
		check_exist(args.FILE_SEQUENCE_TXT, 2)
		str_sequence = read_sequence_txt(args.FILE_SEQUENCE_TXT)
	elif args.FILE_SEQUENCE_GEN is not None:
		check_exist(args.FILE_SEQUENCE_GEN, 2)
		str_sequence = read_sequence_genbank(args.FILE_SEQUENCE_GEN)

	# read parameter
	parameters, base_pair = read_parameters(args.FILE_PARAMETER)

	# generate block sequence and calculate stability
	BLOCK_LENGTH = args.BLOCK_LENGTH
	if args.BLOCK_LENGTH is None:
		BLOCK_LENGTH = len(str_sequence)


	# output
	FILE_LOG = args.FILE_LOG
	if FILE_LOG is None:
		path = os.path.splitext(args.FILE_OUTPUT)
		FILE_LOG = "{}.log".format(path[0])
	if args.FLAG_OVERWRITE == False:
		check_overwrite(args.FILE_OUTPUT)
		check_overwrite(args.FILE_LOG)

	with open(FILE_LOG, "w") as obj_output:
		obj_output.write("This log file was generated at " + datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "\n\n")
		obj_output.write("<< Input >>\n")
		obj_output.write("Program version: {0}\n".format(VERSION))
		obj_output.write("\n")

		if args.FILE_SEQUENCE_GEN is not None:
			obj_output.write("Sequence file: {0} ({1})\n".format(args.FILE_SEQUENCE_GEN, "GenBank"))
		else:
			obj_output.write("Sequence file: {0} ({1})\n".format(args.FILE_SEQUENCE_TXT, "TXT"))
		for idx in range(0, len(str_sequence), 50):
			if len(str_sequence[idx :]) < 50:
				obj_output.write("> {0}\n".format(str_sequence[idx :]))
			else:
				obj_output.write("> {0}\n".format(str_sequence[idx : idx + 50]))
		obj_output.write("{0} bp.\n".format(len(str_sequence)))
		obj_output.write("\n")

		obj_output.write("Block length: {0}\n".format(args.BLOCK_LENGTH))

		obj_output.write("Block number: {0}\n".format(range(0, len(str_sequence) - BLOCK_LENGTH + 1)))
		obj_output.write("\n")

		obj_output.write("Parameter file: {0}\n".format(args.FILE_PARAMETER))
		obj_output.write("|{0:^15}|{1[0]:^8}|{1[1]:^8}|{1[2]:^8}|\n".format("", [param.name for param in parameters]))
		obj_output.write("|{0:-^15}|{0:-^8}|{0:-^8}|{0:-^8}|\n".format(""))
		for parameter_type in parameters[0].get_parameter(data_type="name"):
			obj_output.write("|{0:^15}|{1[0]:>8.3f}|{1[1]:>8.3f}|{1[2]:>8.3f}|\n".format(parameter_type, [param.get_parameter(parameter_type)[0] for param in parameters]))
		obj_output.write("\n")
		obj_output.write("Output file: {0}\n".format(args.FILE_OUTPUT))
		obj_output.write("\n")

	x = []
	list_y = []
	with open(args.FILE_OUTPUT, "w") as obj_output:
		writer = csv.writer(obj_output)
		parameter_types = list(parameters[0].get_parameter(data_type="name"))
		writer.writerow(["Block No.", "Start (Sequence)", "End (Sequence)", "Sequence"] + [param.name for param in parameters] + [""] + parameter_types)
		for i in range(0, len(str_sequence)-BLOCK_LENGTH+1):
			sequence = Sequence(i)
			sequence.set_sequence(str_sequence[i : i+BLOCK_LENGTH], base_pair)
			freq = sequence.get_freq(parameters[0], base_pair)
			energy = [sequence.get_energy(param, base_pair) for param in parameters]
			writer.writerow([i+1, i+1, i+BLOCK_LENGTH] + [sequence.get_sequence("string")] + energy + [""] + [freq[param] for param in parameter_types])
			x.append(i+1)
			list_y.append(energy)

	# graph
	if args.FILE_GRAPH is not None:
		matplot_init()

		for param_idx in range(len(parameters)):
			output_graph_file = "{}_{}.png".format(args.FILE_GRAPH, parameters[param_idx].name)
			if args.FLAG_OVERWRITE == False:
				check_overwrite(output_graph_file)

			name = parameters[param_idx].name

			tic_x = [1, len(x), len(x)/10]
			y = [v[param_idx] for v in list_y]

			obj_graph = Graph2D(600, 300)
			obj_graph.append_ax(1, 1, 1)

			obj_graph.set_label("x", "Block number")
			obj_graph.set_range("x", tic_x[0], tic_x[1])

			obj_graph.set_label("y", name)

			obj_graph.add_grid("y")
			obj_graph.add_zeroaxis_bar("y", 0)

			obj_graph.ax.plot(
				x,
				y,
				color="#000000",
				linewidth=2.0,
				label=name
			)

			obj_graph.make_graph("CUI", output_file=output_graph_file)
