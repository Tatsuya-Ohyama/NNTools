#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
NNcountTool.py - Program to count NN pairs
"""

import sys, signal
sys.dont_write_bytecode = True
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
import csv
import re

from classes.basicfunc import check_exist, check_overwrite



# =============== constant =============== #
VERSION = 1.0
TEMPLATE_PARAM = "template_ref_param.csv"
VALID_BASES = ["A", "C", "G", "T", "U"]
INVALID_BASES = ["R", "Y", "K", "M", "S", "W", "B", "D", "H", "V", "N", "-"]



# =============== function =============== #
def make_template(flag_overwrite):
	"""
	Function to create template files for ref_param.csv
	"""
	if flag_overwrite == False:
		check_overwrite(TEMPLATE_PARAM)
	with open(TEMPLATE_PARAM, "w") as obj_output:
		writer = csv.writer(obj_output)
		writer.writerow(["Parameter", "dH", "dS", "dG"])
	sys.stderr.write("{0} is created.\n".format(TEMPLATE_PARAM))


def read_parameter_file(parameter_file):
	"""
	Function to read parameter file

	Args:
		parameter_file (str): parameter file path

	Returns:
		dict: {"AA/TT": ["AA", "TT"], ...}
	"""
	parameters = {}
	with open(parameter_file) as obj_input:
		reader = csv.reader(obj_input)
		flag_read = False
		for row_val in reader:
			if row_val[0] == "Parameter":
				flag_read = True
				continue

			if flag_read and len(row_val) != 0:
				parameter_name = row_val[0]
				if "/" in row_val[0] and len(row_val[0]) == 5:
					elems = row_val[0].upper().split("/")
					parameters[parameter_name] = list(set([elems[0], elems[1][::-1]]))
	return parameters


def read_fasta(fasta_file, parameters):
	sequences = []
	with open(fasta_file, "r") as obj_input:
		sequence = ""
		for line_val in obj_input:
			if line_val.startswith(">"):
				if len(sequence) != 0:
					sequences[-1].set_parameters(analysis_sequence(sequence, parameters))
					sequences[-1].set_length(len(sequence))
					sequences[-1].set_invalid_bases(count_invalid_bases(sequence))
					sequence = ""
					sys.stderr.write("done.\n")

				sequences.append(Sequence(line_val.strip().replace(">", "", 1), parameters.keys()))
				sys.stderr.write("{0} ... ".format(sequences[-1].name))
				sys.stderr.flush()
				continue

			sequence += line_val.upper().strip()

		if len(sequence) != 0:
			sequences[-1].set_parameters(analysis_sequence(sequence, parameters))
			sequences[-1].set_length(len(sequence))
			sequences[-1].set_invalid_bases(count_invalid_bases(sequence))
			sys.stderr.write("done.\n")

	return sequences


def analysis_sequence(sequence, parameters):
	"""
	Function count NN pairs

	Args:
		sequence (str): sequence
		parameters (dict): {"AA/TT": [AA, TT], ...}

	Returns:
		dict: {"AA/TT": count(int), ...}
	"""
	parameter_count = {}
	for name, patterns in parameters.items():
		parameter_count[name] = 0
		for pattern in patterns:
			parameter_count[name] += sequence.count(pattern)

	return parameter_count


def count_invalid_bases(sequence):
	"""
	Function to count invalid bases

	Args:
		sequence (str): sequence

	Returns:
		int: number of invalid bases
	"""
	return {invalid_base: sequence.count(invalid_base) for invalid_base in INVALID_BASES}


def output_csv(output_file, sequences, parameters):
	"""
	Function to output statistics for sequence

	Args:
		output_file (str): output file
		sequences (list): [SequenceObject, ...]
		parameters (dict): {"AA/TT": [AA, TT], ...}
	"""
	parameter_names = list(parameters.keys())
	with open(output_file, "w") as obj_output:
		writer = csv.writer(obj_output)
		writer.writerow(["Name", "Length", "Invalid"] + parameter_names + INVALID_BASES)
		for obj_sequence in sequences:
			writer.writerow(
				[
					obj_sequence.name,
					obj_sequence.length,
					sum(obj_sequence.invalid_bases.values()),
				] \
				+ [obj_sequence.parameters[parameter_name] for parameter_name in parameter_names] \
				+ [obj_sequence.invalid_bases[invalid_base] for invalid_base in INVALID_BASES]
			)



# =============== class =============== #
class Sequence:
	def __init__(self, name, parameter_names):
		self._name = None
		self._parameters = {}
		self._length = 0
		self._invalid_bases = {}

		self.set_name(name)
		self._parameters = {parameter_name: 0 for parameter_name in parameter_names}
		self._invalid_bases = {base: 0 for base in INVALID_BASES}


	@property
	def name(self):
		return self._name

	@property
	def parameters(self):
		return self._parameters

	@property
	def length(self):
		return self._length

	@property
	def invalid_bases(self):
		return self._invalid_bases


	def set_name(self, name):
		"""
		Method to set sequence name

		Args:
			name (str): name

		Returns:
			self
		"""
		self._name = name
		return self


	def set_length(self, length):
		"""
		Method to set sequence length

		Args:
			length (int): sequence length

		Returns:
			self
		"""
		self._length = length
		return self


	def set_invalid_bases(self, invalid_bases):
		"""
		Method to set number of invalid bases

		Args:
			invalid_bases (dict): {INVALID_BASE(str): count(int), ...}

		Returns:
			self
		"""
		self._invalid_bases = invalid_bases
		return self


	def set_invalid_base(self, name, invalid_base):
		"""
		Method to set number of a invalid bases

		Args:
			name (str): invalid base name
			invalid_base (int): number of invalid bases

		Returns:
			self
		"""
		self._invalid_bases[name] = invalid_base
		return self


	def set_parameters(self, parameters):
		"""
		Method to set parameters

		Args:
			parameters (list): ["AA/TT", "TA/AT", ...]

		Returns:
			self
		"""
		self._parameters = parameters
		return self


	def update_parameter(self, parameter_name, count):
		"""
		Method to update a parameter

		Args:
			parameter_name (str): parameter name
			count (int): parameter count

		Returns:
			self
		"""
		self._parameters[parameter_name] = count
		return self


	def increase_parameter(self, parameter_name):
		"""
		Method to increase parameter count

		Args:
			parameter_name (str): parameter name

		Returns:
			self
		"""
		self._parameters[parameter_name] += 1
		return self



# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Program to count NN pairs", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-f", dest="FASTA_FILES", metavar="SEQUENCE_FILE.fasta", nargs="+", required="--make-template" not in sys.argv, help="sequence file of FASTA format")
	parser.add_argument("-p", dest="PARAMETER_FILE", metavar="PARAMETER_FILE.csv", required="--make-template" not in sys.argv, help="input file for parameters")
	parser.add_argument("-o", dest="OUTPUT_FILE", metavar="OUTPUT_FILE.csv", required="--make-template" not in sys.argv, help="output file")
	parser.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite forcibly")
	parser.add_argument("--make-template", dest="FLAG_MAKE_TEMPLATE", action="store_true", default=False, help="make template files ({0}) and exit".format(TEMPLATE_PARAM))
	args = parser.parse_args()

	if args.FLAG_MAKE_TEMPLATE:
		make_template(args.FLAG_OVERWRITE)
		sys.exit(0)

	check_exist(args.PARAMETER_FILE, 2)

	# reading parameters
	parameters = read_parameter_file(args.PARAMETER_FILE)

	# reading sequence
	sequences = []
	for fasta_file in args.FASTA_FILES:
		sequences.extend(read_fasta(fasta_file, parameters))

	if args.FLAG_OVERWRITE == False:
		check_overwrite(args.OUTPUT_FILE)

	output_csv(args.OUTPUT_FILE, sequences, parameters)
