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
import re

from classes.basicfunc import check_exist, check_overwrite
from classes.parameter import Parameter
from classes.sequence import Sequence



# =============== variable =============== #
VERSION = "1.1"
TEMPLATE_PARAM = "template_ref_param.csv"
TEMPLATE_SEQUENCE = "template_sequence.csv"
LIMIT_LEN_SEQUENCE = 100
STANDARD_BASES = ["A", "C", "G", "T", "U"]
RE_STANDARD_SEQUENCE = re.compile(r"([ACGTU]+)")



# =============== function =============== #
def make_template(flag_overwrite):
	"""
	Function to create template files for ref_param.csv and ref_exp.csv
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



def read_sequence_csv(input_file, base_pair):
	"""
	Function to read sequences from .csv file

	Args:
		input_file (str): csv file path
		base_pair (dict): base pair

	Returns:
		list: [SequenceObject, ...]
	"""
	sequences = []
	with open(input_file, "r") as obj_input:
		reader = csv.reader(obj_input)
		flag_read = False
		for line_val in reader:
			if "Sequence" in line_val or "Sequences" in line_val:
				flag_read = True
				continue

			if flag_read:
				sequences.append(Sequence(line_val[0]).set_sequence(line_val[1].upper(), base_pair))
	return sequences



def read_sequence_fasta(input_file, base_pair, flag_skip):
	"""
	Function to read sequences from .fasta or .fna file

	Args:
		input_file (str): fasta file path
		base_pair (dict): base pair
		flag_skip (bool): if True, skip ambiguous base and separate sequence by the base

	Returns:
		list: [SequenceObject, ...]
	"""
	sequences = []
	with open(input_file, "r") as obj_input:
		sequence = ""
		for line_val in obj_input:
			if line_val.startswith(">"):
				if sequence != "":
					# when previous sequence remained, register sequence
					if flag_skip:
						standard_sequences = RE_STANDARD_SEQUENCE.findall(sequence)
						if len(standard_sequences) != 1:
							# when multiple sequences were found
							num = 1
							base_name = sequences[-1].name
							sequences[-1].set_name("{0} ({1})".format(base_name, num))
							sequences[-1].set_sequence(standard_sequences[0], base_pair)
							for seq in standard_sequences[1:]:
								if len(seq) == 1:
									continue

								num += 1
								sequences.append(Sequence("{0} ({1})".format(base_name, num)).set_sequence(seq, base_pair))
							sequence = ""
							continue

					sequences[-1].set_sequence(sequence, base_pair)
					sequence = ""

				# create new sequence object
				sequences.append(Sequence(line_val.strip().replace(">", "", 1)))

			else:
				# save sequence or connect sequence with previous line
				sequence += line_val.strip().upper()

		if sequence != "":
			sequences[-1].set_sequence(sequence, base_pair)

	return sequences



# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="neighbor method", formatter_class=argparse.RawTextHelpFormatter)
	parser_sequence = parser.add_mutually_exclusive_group(required="--make-template" not in sys.argv)
	parser_sequence.add_argument("-s", dest="FILE_SEQUENCE", metavar="SEQUENCE_FILE.csv", nargs="+", help="sequence file of csv format")
	parser_sequence.add_argument("-f", dest="FILE_FASTA", metavar="SEQUENCE_FILE.fasta", nargs="+", help="sequence file of FASTA format")
	parser.add_argument("-p", dest="FILE_PARAMETER", metavar="PARAMETER_FILE.csv", required="--make-template" not in sys.argv, help="input file for parameters")
	parser.add_argument("-o", dest="FILE_OUTPUT", metavar="OUTPUT_FILE.csv", required="--make-template" not in sys.argv, help="output file")
	parser.add_argument("-l", dest="FILE_LOG", metavar="LOG_FILE.log", help="log file (if not specify, a log file with the same name as -o option is generated)")
	parser.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite forcibly")
	parser.add_argument("--make-template", dest="FLAG_MAKE_TEMPLATE", action="store_true", default=False, help="make template files ({0} and {1}) and exit".format(TEMPLATE_PARAM, TEMPLATE_SEQUENCE))
	parser.add_argument("--skip-ambiguous", dest="FLAG_SKIP", action="store_true", default=False, help="skip ambiguous base (separate sequence by ambiguous base)")
	args = parser.parse_args()

	if args.FLAG_MAKE_TEMPLATE:
		make_template(args.FLAG_OVERWRITE)
		sys.exit(0)

	check_exist(args.FILE_PARAMETER, 2)

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
					tmp_base_pair[base[0][0:1]] = base[1][0:1].upper()
					tmp_base_pair[base[0][1:2]] = base[1][1:2].upper()
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
	if args.FILE_SEQUENCE is not None:
		for file_sequence in args.FILE_SEQUENCE:
			check_exist(file_sequence, 2)
			sequences.extend(read_sequence_csv(file_sequence, base_pair))

	elif args.FILE_FASTA is not None:
		for file_fasta in args.FILE_FASTA:
			check_exist(file_fasta, 2)
			sequences.extend(read_sequence_fasta(file_fasta, base_pair, args.FLAG_SKIP))


	# output
	if args.FILE_LOG is not None:
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
			seq = sequence.get_sequence("string")
			if len(seq) > LIMIT_LEN_SEQUENCE:
				seq = seq[:LIMIT_LEN_SEQUENCE] + "..."
			writer.writerow([sequence.name, seq] + energy + [""] + [freq[param] for param in parameter_types])
