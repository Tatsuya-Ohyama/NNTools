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
import time

from basic_func import check_exist, check_overwrite
from classes.Parameter import Parameter
from classes.Sequence import Sequence


# =============== variable =============== #
VERSION = "1.2"


# =============== function =============== #
def read_sequence_txt(input_file):
	"""
	read sequence from txt format file
	@param input_file: input file of txt format
	@return sequence (str)
	"""
	sequence = ""
	with open(args.file_sequence, "r") as obj_input:
		sequence = "".join([line_val.strip() for line_val in obj_input.readlines()])
	return sequence.upper()


def read_sequence_genbank(input_file):
	"""
	read sequence from genbank format file
	@param input_file: input file of genbank format file
	@return sequence (str)
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


# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "NNscanTool.py - Program for searching stable area by nearest neighbor method", formatter_class=argparse.RawTextHelpFormatter)
	sequence_group = parser.add_mutually_exclusive_group(required = True)
	sequence_group.add_argument("-s", dest = "file_sequence_txt", metavar = "SEQUENCE_FILE.txt", help = "sequence file of txt format")
	sequence_group.add_argument("-g", dest = "file_sequence_gen", metavar = "SEQUENCE_FILE.gb", help = "sequence file of genbank format")
	parser.add_argument("-p", dest = "file_parameter", metavar = "PARAMETER_FILE.csv", required = True, help = "input file for parameters")
	parser.add_argument("-b", dest = "block_length", metavar = "SEQUENCE_RANGE", type = int, help = "block length for sequence (Default: full length)")
	parser.add_argument("-o", dest = "file_output", metavar = "OUTPUT_FILE.csv", required = True, help = "output file")
	parser.add_argument("-l", dest = "file_log", metavar = "LOG_FILE.log", help = "log file (if not specify, a log file with the same name as -o option is generated)")
	parser.add_argument("--graph", dest = "file_graph", metavar = "GRAPH_BASE_NAME", help = "base name of graph (.png)")
	parser.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly")
	args = parser.parse_args()

	check_exist(args.file_parameter, 2)

	str_sequence = ""
	if args.file_sequence_txt is not None:
		check_exist(args.file_sequence_txt, 2)
		str_sequence = read_sequence_txt(args.file_sequence_txt)
	elif args.file_sequence_gen is not None:
		check_exist(args.file_sequence_gen, 2)
		str_sequence = read_sequence_genbank(args.file_sequence_gen)


	# reading parameters
	parameters = []
	base_pair = {}
	with open(args.file_parameter, "r") as obj_input:
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
					parameters[param_idx].append_parameter(line_val[0], float(line_val[param_idx + 1].strip()))

	new_base_pair = {}
	for k, v in base_pair.items():
		new_base_pair[k] = v
		if v not in base_pair.keys():
			new_base_pair[v] = k
	base_pair = new_base_pair


	# generate block sequence and calculate stability
	block_length = args.block_length
	if args.block_length is None:
		block_length = len(str_sequence)


	# output
	file_log = args.file_log
	if file_log is None:
		path = os.path.splitext(args.file_output)
		file_log = path[0] + ".log"
	if args.flag_overwrite == False:
		check_overwrite(args.file_output)
		check_overwrite(args.file_log)

	with open(file_log, "w") as obj_output:
		obj_output.write("This log file was generated at " + datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "\n\n")
		obj_output.write("<< Input >>\n")
		obj_output.write("Program version: {0}\n".format(VERSION))
		obj_output.write("\n")

		if args.file_sequence_gen is not None:
			obj_output.write("Sequence file: {0} ({1})\n".format(args.file_sequence_gen, "GenBank"))
		else:
			obj_output.write("Sequence file: {0} ({1})\n".format(args.file_sequence_txt, "TXT"))
		for idx in range(0, len(str_sequence), 50):
			if len(str_sequence[idx :]) < 50:
				obj_output.write("> {0}\n".format(str_sequence[idx :]))
			else:
				obj_output.write("> {0}\n".format(str_sequence[idx : idx + 50]))
		obj_output.write("{0} bp.\n".format(len(str_sequence)))
		obj_output.write("\n")

		obj_output.write("Block length: {0}\n".format(args.block_length))

		obj_output.write("Block number: {0}\n".format(range(0, len(str_sequence) - block_length + 1)))
		obj_output.write("\n")

		obj_output.write("Parameter file: {0}\n".format(args.file_parameter))
		obj_output.write("|{0:^15}|{1[0]:^8}|{1[1]:^8}|{1[2]:^8}|\n".format("", [param.get_name() for param in parameters]))
		obj_output.write("|{0:-^15}|{0:-^8}|{0:-^8}|{0:-^8}|\n".format(""))
		for parameter_type in parameters[0].get_parameter(data_type = "name"):
			obj_output.write("|{0:^15}|{1[0]:>8.3f}|{1[1]:>8.3f}|{1[2]:>8.3f}|\n".format(parameter_type, [param.get_parameter(parameter_type)[0] for param in parameters]))
		obj_output.write("\n")
		obj_output.write("Output file: {0}\n".format(args.file_output))
		obj_output.write("\n")

	x = []
	list_y = []
	with open(args.file_output, "w") as obj_output:
		writer = csv.writer(obj_output)
		parameter_types = list(parameters[0].get_parameter(data_type = "name"))
		writer.writerow(["Block No.", "Start (Sequence)", "End (Sequence)", "Sequence"] + [param.get_name() for param in parameters] + [""] + parameter_types)
		for i in range(0, len(str_sequence) - block_length + 1):
			sequence = Sequence(i)
			sequence.set_sequence(str_sequence[i : i + block_length], base_pair)
			freq = sequence.get_freq(parameters[0], base_pair)
			energy = [sequence.get_energy(param, base_pair) for param in parameters]
			writer.writerow([i + 1, i + 1, i + block_length] + [sequence.get_sequence("string")] + energy + [""] + [freq[param] for param in parameter_types])
			x.append(i + 1)
			list_y.append(energy)

	# graph
	if args.file_graph is not None:
		files_graph = []
		for param in parameters:
			if args.file_graph is None:
				path = os.path.splitext(args.file_output)
				path = path[0] + "_" + param.get_name() + ".png"
				if args.flag_overwrite == False:
					check_overwrite(path)
				files_graph.append(path)

		for param_idx in range(len(parameters)):
			font_size = 22
			y = [v[param_idx] for v in list_y]
			xtic = [1, len(x), len(x) / 10]
			ystep = (max(y) - min(y)) / 10
			ytic = [min(y) - ystep, max(y) + ystep, ystep]
			name = parameters[param_idx].get_name()

			fig = plt.figure(1, figsize=(8.345, 6.95))
			plt.rcParams['axes.linewidth'] = 2.0
			plt.rcParams['axes.axisbelow'] = True
			plt.rcParams['font.family'] = "Segoe UI"
			plt.rcParams['xtick.direction'] = "in"
			plt.rcParams.update({"mathtext.default": "regular"})

			ax = fig.add_subplot(111)
			ax.grid(True, axis = "y", color = "#000000", linewidth = 1, dashes = (5,2.5))
			ax.tick_params(axis = "x", which = "major", labelsize = font_size, pad = 5, rotation = 90)
			ax.tick_params(axis = "y", which = "major", labelsize = font_size, pad = 5)

			ax.set_xlabel("Block number", fontsize = font_size)
			ax.set_xlim(xtic[0], xtic[1])
			ax.set_xticks(np.arange(xtic[0], xtic[1] + xtic[2], step = xtic[2]))

			ax.set_ylabel(name, fontsize = font_size)
			ax.set_ylim(ytic[0], ytic[1])
			ax.set_yticks(np.arange(ytic[0], ytic[1] + ytic[2], step = ytic[2]))
			ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))

			ax.axhline(y = 0, color = "#000000")
			ax.plot(x,
					y,
					color = "#000000",
					linewidth = 2.0,
					label = name
				)

			fig.tight_layout()
			plt.savefig(files_graph[param_idx])
			plt.close()
