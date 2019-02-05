#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pickle

from classes.Parameter import Parameter

# =============== variable =============== #
BASE_PAIR = {"A": "T", "T": "A", "G": "C", "C": "G"}
parameter_list = ["AA/TT", "AT/TA", "TA/AT", "CA/GT", "GT/CA", "CT/GA", "GA/CT", "CG/GC", "GC/CG", "GG/CC", "init_GC", "init_AT", "symmetry", "5term_TA"]


# =============== class =============== #
class Sequence:
	""" Sequence class """
	def __init__(self, name, sequence = None):
		# member variables
		self._name = ""
		self._sequence = ""
		self._energy_type = ""
		self._nucleic_type = "DNA"
		self._is_self_complement = False

		# initiation
		self.set_name(name)
		self.set_sequence(sequence)


	def save_pickle(self, output_file):
		"""
		Pickle ファイルに保存するメソッド
		@param output_file: 出力する pickle ファイルのパス
		@return: 自身を返す (チェーンメソッドのため)"""
		with open(output_file, "wb") as obj_output:
			pickle.dump(self, obj_output)
			sys.stderr.write("INFO: save pickle file to '{0}'\n".format(output_file))
		return self


	def restore_pickle(self, input_file):
		"""
		Pickle ファイルから復元するメソッド
		@param input_file: pickle ファイルのパス
		@return: 自身を返す (チェーンメソッドのため)
		"""
		with open(input_file, "rb") as obj_input:
			self = pickle.load(obj_input)
			sys.stderr.write("INFO: restore object from pickle file '{0}'\n".format(input_file))
		return self


	def set_name(self, name):
		"""
		set name
		@param name: name
		@return self
		"""
		self._name = name
		return self


	def set_energy_type(self, energy_type):
		"""
		set energy_type
		@param energy_type: dH, dS, dG
		@return self
		"""
		self._energy_type = energy_type
		return self


	def set_sequence(self, sequence, nucleic_type = "DNA"):
		"""
		set sequence method
		@param sequence: sequence
		@param nucleic_type: DNA or RNA
		@return: self
		"""
		if sequence is not None:
			self._sequence = list(sequence)
			self._nucleic_type = nucleic_type
			complement = [BASE_PAIR[base] for base in self._sequence]
			self._is_self_complement = self._sequence == list(reversed(complement))
		return self


	def get_name(self):
		"""
		return name
		@return name
		"""
		return self._name


	def get_sequence(self, sequence_type = "list"):
		"""
		return sequence
		@param sequence_type: "list" or "string" (Default: list)
		@return sequence
		"""
		if sequence_type == "string":
			return "".join(self._sequence)
		elif sequence_type == "list":
			return self._sequence
		else:
			sys.stderr.write("ERROR: undefined sequence_type at get_sequence() in Sequence class.\n")
			sys.exit(1)


	def is_complement(self):
		"""
		return self complement or not
		@return True or False
		"""
		return self._is_self_complement


	def get_energy_type(self):
		"""
		return energy_type
		@return energy_type
		"""
		return self._energy_type


	def get_freq(self):
		"""
		return pair frequency
		@return pair frequency list
		"""
		freq = {
			"AA/TT": 0,
			"AT/TA": 0,
			"TA/AT": 0,
			"CA/GT": 0,
			"GT/CA": 0,
			"CT/GA": 0,
			"GA/CT": 0,
			"CG/GC": 0,
			"GC/CG": 0,
			"GG/CC": 0,
			"init_GC": 0,
			"init_AT": 0,
			"symmetry": 0,
			"5term_TA": 0
		}
		if self._sequence[0] in ["G", "C"]:
			freq["init_GC"] += 1
		elif self._sequence[0] in ["A", "T"]:
			freq["init_AT"] += 1
			if self._sequence[0] == "T":
				freq["5term_TA"] += 1

		if self._is_self_complement:
			freq["symmetry"] += 1

		for base_idx in range(len(self._sequence) - 1):
			pair_forward = [self._sequence[base_idx], self._sequence[base_idx + 1]]
			pair_reverse = [BASE_PAIR[self._sequence[base_idx]], BASE_PAIR[self._sequence[base_idx + 1]]]
			pair_type = "/".join(["".join(pair_forward), "".join(pair_reverse)])
			if pair_type not in parameter_list:
				pair_type = "/".join(["".join(reversed(pair_reverse)), "".join(reversed(pair_forward))])
			freq[pair_type] += 1

		return freq


	def get_energy(self, obj_parameter):
		"""
		return energy value
		@param obj_parameter: Parameter object
		@return energy_value
		"""
		# エネルギーを計算させる
		energy = 0.0
		if self._sequence[0] in ["G", "C"]:
			# init_GC
			energy += obj_parameter.get_parameter("init_GC")
		elif self._sequence[0] in ["A", "T"]:
			# init_AT
			energy += obj_parameter.get_parameter("init_AT")
			if self._sequence[0] == "T":
				# 5term_TA
				energy += obj_parameter.get_parameter("5term_TA")
		else:
			sys.stderr.write("ERROR: undefined initiation base pair.\n")
			sys.exit(1)


		for base_idx in range(len(self._sequence) - 1):
			pair_forward = [self._sequence[base_idx], self._sequence[base_idx + 1]]
			pair_reverse = [BASE_PAIR[self._sequence[base_idx]], BASE_PAIR[self._sequence[base_idx + 1]]]
			pair_type = "/".join(["".join(pair_forward), "".join(pair_reverse)])
			if pair_type not in parameter_list:
				pair_type = "/".join(["".join(reversed(pair_reverse)), "".join(reversed(pair_forward))])

			energy += obj_parameter.get_parameter(pair_type)

		if self._is_self_complement:
			energy += obj_parameter.get_parameter("symmetry")

		return energy


# =============== main =============== #
# if __name__ == '__main__':
# 	main()
