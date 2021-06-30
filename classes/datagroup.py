#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
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
		self._name = ""
		self._sequences = []
		self._energy = []
		self._error = []
		self._error_sign = []
		self._base_pairs = {}
		self._is_fitting = None

		self.set_name(name)


	@property
	def name(self):
		return self._name

	@property
	def sequences(self):
		return self._sequences

	@property
	def energy(self):
		return self._energy

	@property
	def error(self):
		return self._error

	@property
	def error_sign(self):
		return self._error_sign

	@property
	def base_pairs(self):
		return self._base_pairs

	@property
	def is_fitting(self):
		return self._is_fitting


	def set_name(self, name):
		"""
		Method to set name

		Args:
			name (str): data name

		Returns:
			self
		"""
		self._name = name
		return self


	def set_base_pair(self, base_pairs):
		"""
		Method to set base pair information

		Args:
			base_pairs (dict): base pairs

		Returns:
			self
		"""
		self._base_pairs = base_pairs
		return self


	def append(self, obj_sequence, exp_value, exp_value_e=0.0):
		"""
		Method to append sequence object and experimental data

		Args:
			obj_sequence (objSequence): Sequence Object
			exp_value (float): experimental value
			exp_value_e (float, optional): error value of experimental value (Default: 0.0)

		Returns:
			self
		"""
		self._sequences.append(obj_sequence)
		self._energy.append(exp_value)
		self._error.append(exp_value_e)
		self._error_sign.append(0)
		if len([v for v in self._energy if v == 0.0]) == len(self._energy):
			# If all the values are 0.0, it can not be fitted and will not be calculated later
			self._is_fitting = False
		else:
			self._is_fitting = True

		return self


	def get_sequence(self, data_type=None):
		"""
		Method to return Sequence object

		Args:
			data_type (None or str, optional): None or "sequence" (Default: None)

		Returns:
			list of objSequence for None or sequence string list for "sequence"
		"""
		if data_type is None:
			return self._sequences
		elif data_type == "sequence":
			return [x.get_sequence() for x in self._sequences]


	def get_energy(self, flag_sequence=False, obj_parameters=[]):
		"""
		Method to return experimental value list

		Args:
			flag_sequence (bool, optional): return energy with sequence (Default: False)
			obj_parameters (list, optional): Paramete robject list (Default: [])

		Returns:
			list: energy value
		"""
		energy_data = []
		for idx in range(len(self._sequences)):
			data = []
			if flag_sequence:
				data.append(self._sequences[idx].get_sequence("string"))
			data.append(self._energy[idx])
			data.append(self._error[idx])
			for parameter in obj_parameters:
				data.append(self._sequences[idx].get_energy(parameter, self._base_pairs))
			energy_data.append(data)
		return energy_data


	def get_stat(self, obj_parameter, mode=None, deg=1, error_sign=None):
		"""
		Method to return statistics

		Args:
			obj_parameter (str): parameter name
			mode (None or str, optional): None, "r", "r2", "slope", "intercept", "diff_abs", "diff_mean", "diff_std", "diff_sum, diff_square" (Default: None)
			deg (int, optional): dimension for curve fitting (Default: 1)
			error_sign (list, optional): sign for error (Default: None)

		Returns:
			list: [r, r2, slope, intercept, diff_abs, diff_mean, diff_std, diff_sum, diff_square] list when data_type is None
		"""
		if self._is_fitting:
			if type(error_sign) != list:
				error_sign = self._error_sign
			x = np.array([self._energy[idx] + self._error[idx] * error_sign[idx] for idx in range(len(self._energy))])
			y = np.array([sequence.get_energy(obj_parameter, self._base_pairs) for sequence in self._sequences])
			result = [float(x) for x in np.polyfit(x, y, deg).tolist()]
			if y[y == 0.0].shape[0] == len(self._sequences) or np.std(x) == 0.0 or np.std(y) == 0.0:
				result.append(0.0)
			else:
				result.append(np.corrcoef(x, y)[0, 1])

			result.append(result[-1] ** 2)
			result.append(np.abs(x - y))
			result.append(np.mean(x - y))
			result.append(np.std(x - y))
			result.append(np.sum(np.abs(x - y)))
			result.append(np.sum((x - y) ** 2))

		else:
			result = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

		if mode is None:
			return result
		elif mode == "slope":
			return result[0]
		elif mode == "intercept":
			return result[1]
		elif mode == "r":
			return result[2]
		elif mode == "r2":
			return result[3]
		elif mode == "diff_abs":
			return result[4]
		elif mode == "diff_mean":
			return result[5]
		elif mode == "diff_std":
			return result[6]
		elif mode == "diff_sum":
			return result[7]
		elif mode == "diff_square":
			return result[8]
		else:
			sys.stderr.write("ERROR: undefined mode at get_stat() in DataGroup class.\n")
			sys.exit(1)
