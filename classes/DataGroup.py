#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pickle
import numpy as np
import statistics


# =============== class =============== #
class DataGroup:
	""" DataGroup class """
	def __init__(self, name):
		"""
		@param name: this object name
		"""
		# member variables
		self._sequences = []
		self._energy = []


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


	def append(self, obj_sequence, exp_value):
		"""
		append sequence object and experimental data
		@param obj_sequence: Sequence object
		@param exp_value: experimental value
		@return: self
		"""
		self._sequences.append(obj_sequence)
		self._energy.append(exp_value)
		return self


	def get_sequence(self, data_type = None):
		"""
		return Sequence object
		@param data_type: None or "sequence"
		@return sequence object list for None or sequence string list for "sequence"
		"""
		if data_type is None:
			return self._sequences
		elif data_type == "sequence":
			return [x.get_sequence() for x in self._sequences]


	def get_energy(self, flag_sequence = False, obj_parameters = []):
		"""
		return experimental value list
		@param flag_sequence: return energy with sequence (Default: False)
		@param obj_parameters: Parameter object list (Default: [])
		@return energy value list
		"""
		energy = []
		if flag_sequence:
			energy = [[sequence.get_sequence("string") for sequence in self._sequences]]
		energy += [self._energy]
		energy += [[sequence.set_parameter(parameter).get_energy() for sequence in self._sequences] for parameter in obj_parameters]
		energy = [[energy[row_idx][col_idx] for row_idx in range(len(energy))] for col_idx in range(len(energy[0]))]
		return energy


	def get_stat(self, obj_parameter, data_type = None, deg = 1):
		"""
		return statistics
		@param data_type: None, "r", "r2", "slope", "intercept", "diff_abs", "diff_mean", "diff_sum, diff_square" (Default: None)
		@param obj_parameter: degree of the fitting polynomial (Default: 1)
		@param deg:
		@return statistics value or return [r, r2, slope, intercept, diff_abs, diff_mean, diff_sum, diff_square] list when data_type is None
		"""
		x = np.array(self._energy)
		y = np.array([sequence.set_parameter(obj_parameter).get_energy() for sequence in self._sequences])
		result = [float(x) for x in np.polyfit(x, y, deg).tolist()]
		result.append(np.corrcoef(x, y)[0, 1])
		result.append(result[-1] ** 2)
		result.append(np.abs(x - y))
		result.append(np.mean(x - y))
		result.append(np.sum(np.abs(x - y)))
		result.append(np.sum((x - y) ** 2))

		if data_type is None:
			return result
		elif data_type == "slope":
			return result[0]
		elif data_type == "intercept":
			return result[1]
		elif data_type == "r":
			return result[2]
		elif data_type == "r2":
			return result[3]
		elif data_type == "diff_abs":
			return result[4]
		elif data_type == "diff_mean":
			return result[5]
		elif data_type == "diff_sum":
			return result[6]
		elif data_type == "diff_square":
			return result[7]
		else:
			sys.stderr.write("ERROR: undefined data_type at get_stat() in DataGroup class.\n")
			sys.exit(1)


# =============== main =============== #
# if __name__ == '__main__':
# 	main()
