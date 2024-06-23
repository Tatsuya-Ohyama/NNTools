#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import copy



# =============== class =============== #
class Parameter:
	""" Parameter class for one energy type """
	def __init__(self):
		# member variables
		self._name = ""
		self._parameters = {}	# key: parameter_type, value: [raw, +error, -error]
		self._change = {}
		self._one_direction = False


	@property
	def name(self):
		return self._name

	@property
	def parameters(self):
		return self._parameters

	@property
	def change(self):
		return self._change

	@property
	def one_direction(self):
		return self._one_direction


	def set_name(self, name):
		"""
		Method to set parameter name

		Args:
			name (str): parameter name

		Returns:
			self
		"""
		self._name = name
		return self


	def append_parameter(self, parameter_type, parameter_val):
		"""
		Method to append parameter with error

		Args:
			parameter_type (str): parameter name ("AA/TT", "AT/TA", ..., "init_XX", "symmetry", or "5term_TA (by regexp)")
			parameter_val (list): [param(float), param_min(float), param_max(float)]

		Returns:
			self
		"""
		self._parameters[parameter_type] = [parameter_val for x in range(3)]
		self._change[parameter_type] = True
		return self


	def set_parameter(self, parameter_type, parameter_val=None):
		"""
		Method to set parameter with error

		Args:
			parameter_type (str): parameter name (all, "AA/TT", "AT/TA", "TA/AT", "CA/GT", "GT/CA", "CT/GA", "GA/CT", "CG/GC", "GC/CG", "GG/CC", "init_GC", "init_AT", "symmetry", or "5term_TA (by regexp)")
			parameter_val (float, optional): [param(float), param_min(float), param_max(float)] (Default: None)

		Returns:
			self
		"""
		if parameter_type == "all":
			# All data are received without modification
			self._parameters = {x: [float(y[idx]) for idx in range(3)] if self._change[x] else self._parameters[x] for x, y in parameter_val.items()}
			return self

		if parameter_type not in self._parameters.keys():
			# parameter include by complementary sequence
			parameter_type = "".join(reversed(parameter_type))

		if parameter_type in self._parameters.keys():
			# include parameter name
			if self._change[parameter_type]:
				if type(parameter_val) == list:
					self._parameters[parameter_type] = [float(x) for x in parameter_val]
				else:
					self._parameters[parameter_type] = [float(parameter_val)] * 3
		else:
			sys.stderr.write("ERROR: Undefined parameter_type in set_parameter() of ParameterData class ({0}).\n".format(parameter_type))
			sys.exit(1)

		return self


	def set_one_direction(self, one_direction):
		"""
		Method to set one_direction

		Args:
			one_direction (bool): _description_

		Returns:
			self
		"""
		self._one_direction = one_direction
		return self


	def update_parameter_error(self, parameter_type, parameter_value):
		"""
		Method to calculate and set error

		Args:
			parameter_type (str): parameter name or "all"
			parameter_value (float or list): parameter value(float) for one parameter type or [param(float), param_min(float), param_max(float)] for "all"

		Returns:
			self
		"""
		parameter_types = []
		parameter_values = []
		if parameter_type == "all":
			parameter_types = self._parameters.keys()
			parameter_values = parameter_value
		else:
			if parameter_type not in self._parameters.keys():
				sys.stderr("ERROR: Undefined parameter_type at set_parameter_error() in Parameter class ({0}).\n".format(parameter_type))
				sys.exit(1)
			parameter_types = [parameter_type]
			parameter_values = [parameter_value]

		for parameter_type in parameter_types:
			# loop for parameter types
			if self._change[parameter_type]:
				# when chang flag is True
				if parameter_values[parameter_type][0] < self._parameters[parameter_type][1]:
					# update minimum value
					self._parameters[parameter_type][1] = parameter_values[parameter_type][0]
				elif self._parameters[parameter_type][2] < parameter_values[parameter_type][0]:
					# update maximum value
					self._parameters[parameter_type][2] = parameter_values[parameter_type][0]
		return self


	def set_change_stat(self, parameter_type, state):
		"""
		Method to change state for changing parameter

		Args:
			parameter_type (str): parameter type or "all"
			state (bool): True or False

		Returns:
			self
		"""
		if parameter_type == "all":
			self._change = {k: state for k in self._change.keys()}
		else:
			self._change[parameter_type] = state
		return self


	def remove_parameter(self, parameter_name):
		"""
		Method to remove parameter

		Args:
			parameter_name (str): parameter name

		Returns:
			self
		"""
		if parameter_name in self._parameters.keys():
			del(self._parameters[parameter_name])
		else:
			sys.stderr.write("ERROR: Not found `{0}`.\n".format(parameter_name))
			sys.exit(1)
		return self


	def clone(self):
		"""
		Method to return clone parameter set

		Returns:
			obj_Parameter
		"""
		return copy.deepcopy(self)


	def is_change(self, parameter_type=None):
		"""
		Method to return change state

		Args:
			parameter_type (str, optional): None or parameter type (Default: None)

		Returns:
			list for None or bool for parameter type
		"""
		if parameter_type is None:
			return self._change
		elif parameter_type in self._change.keys():
			return self._change[parameter_type]
		else:
			sys.stderr.write("ERROR: Undefined parameter type at is_chage() in Parameter class ({0}).\n".format(parameter_type))
			sys.exit(1)


	def get_parameter(self, parameter_type=None, data_type="raw"):
		"""
		Method to return parameter

		Args:
			parameter_type (str, optional): parameter type (Default: None)
			data_type (str, optional): raw ([value, minimum value with error, maximum value with error]) or fix (value, error+/-) (Default: "raw")

		Returns:
			list for None(all) or float value for each parameter
		"""
		values = {}
		if data_type == "raw":
			values = self._parameters
		elif data_type == "fix":
			error = {k: abs(v[1] - v[2]) / 2 for k, v in self._parameters.items()}
			min_val = [min(v[1:]) for v in self._parameters.values()]
			values = {k: [m + error[k], error[k]] for m, (k, v) in zip(min_val, self._parameters.items())}
		elif data_type == "name":
			return self._parameters.keys()
		else:
			sys.stderr.write("ERROR: Undefined data_type at get_parameter() in Parameter class.\n".format(data_type))
			sys.exit(1)

		if parameter_type is None:
			return values
		elif parameter_type in self._parameters.keys():
			# exist parameter_type
			return values[parameter_type]
		else:
			parameter_key = "".join(reversed(parameter_type))
			if parameter_key in self._parameters.keys():
				return values[parameter_key]
			else:
				# all parameters are specified
				sys.stderr.write("ERROR: Undefined parameter_type in get_parameter() in Parameter class ({0}).\n".format(parameter_type))
				sys.exit(1)
