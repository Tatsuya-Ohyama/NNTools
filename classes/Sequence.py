#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pickle

from classes.Parameter import Parameter

# =============== variable =============== #
BASE_PAIR = {"A": "T", "T": "A", "G": "C", "C": "G", }


# =============== class =============== #
class Sequence:
	""" Sequence class """
	def __init__(self, name, sequence = None, parameter = None):
		# member variables
		self._name = ""
		self._sequence = ""
		self._parameter = None
		self._nucleic_type = "DNA"
		self._flag_self_complementary = False

		# initiation
		self.set_name(name)
		self.set_sequence(sequence)
		self.set_parameter(parameter)


	def save_pickle(self, output_file):
		"""
		Pickle ファイルに保存するメソッド
		@param output_file: 出力する pickle ファイルのパス
		@return: 自身を返す (チェーンメソッドのため)"""
		with open(output_file, "wb") as obj_output:
			pickle.dump(self, obj_output)
		return self


	def restore_pickle(self, input_file):
		"""
		Pickle ファイルから復元するメソッド
		@param input_file: pickle ファイルのパス
		@return: 自身を返す (チェーンメソッドのため)
		"""
		with open(input_file, "rb") as obj_input:
			self = pickle.load(obj_input)
		return self


	def set_name(self, name):
		"""
		set name
		@param name: name
		@return self
		"""
		self._name = name
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
		return self


	def set_self_complementary(self, flag_self_complementary = False):
		"""
		set flag_self_complementary
		@param flag_self_complementary: True or False (Default: False)
		@return self
		"""
		self._flag_self_complementary = flag_self_complementary
		return self


	def set_parameter(self, obj_parameter = None):
		"""
		set parameter method
		@param obj_parameter: ParameterValue object
		@return self
		"""
		if obj_parameter is not None:
			self._parameter = obj_parameter
		return self


	def get_name(self):
		"""
		return name
		@return name
		"""
		return self._name


	def get_sequence(self):
		"""
		return sequence
		@return sequence
		"""
		return self._sequence


	def get_energy(self):
		"""
		return energy value
		@return energy_value
		"""
		# エネルギーを計算させる
		energy = 0.0
		if self._sequence[0] in ["G", "C"]:
			# init_GC
			energy += self._parameter.get_parameter("init_GC")
		elif self._sequence[0] in ["A", "T"]:
			# init_AT
			energy += self._parameter.get_parameter("init_AT")
			if self._sequence[0] == "T":
				# 5term_AT
				energy += self._parameter.get_parameter("5term_AT")
				pass
		else:
			sys.stderr.write("ERROR: undefined initiation base pair.\n")
			sys.exit(1)


		for base_idx in range(len(self._sequence) - 1):
			pair_forward = self._sequence[base_idx] + self._sequence[base_idx + 1]
			pair_reverse = BASE_PAIR[self._sequence[base_idx]] + BASE_PAIR[self._sequence[base_idx + 1]]
			pair_type = "/".join([pair_forward, pair_reverse])
			energy += self._parameter.get_parameter(pair_type)

		if self._flag_self_complementary:
			energy = [x + y for x, y in zip(energy, self._parameter.get_parameter("all", "symmetry"))]

		return energy


# =============== main =============== #
# if __name__ == '__main__':
# 	main()
