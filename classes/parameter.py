#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
<<<<<<< HEAD
import copy



=======
import pickle
import copy


>>>>>>> origin/master
# =============== class =============== #
class Parameter:
	""" Parameter class for one energy type """
	def __init__(self):
		# member variables
		self._name = ""
		self._parameters = {}	# key: parameter_type, value: [raw, +error, -error]
		self._change = {}


<<<<<<< HEAD
	@property
	def name(self):
		return self._name

	@property
	def parameters(self):
		return self._parameters

	@property
	def change(self):
		return self._change
=======
	def save_pickle(self, output_file):
		"""
		save to pickle
		@param output_file: output pickle file path
		@return self (for chain method)
		"""
		with open(output_file, "wb") as obj_output:
			pickle.dump(self, obj_output)
			sys.stderr.write("INFO: save pickle file to '{0}'\n".format(output_file))
		return self


	def restore_pickle(self, input_file):
		"""
		restore from pickle
		@param input_file: pickle file path
		@return self (for chain method)
		"""
		with open(input_file, "rb") as obj_input:
			self = pickle.load(obj_input)
			sys.stderr.write("INFO: restore object from pickle file '{0}'\n".format(input_file))
		return self
>>>>>>> origin/master


	def set_name(self, name):
		"""
<<<<<<< HEAD
		Method to set parameter name

		Args:
			name (str): parameter name

		Returns:
			self
=======
		set name
		@param name: name
		@return self
>>>>>>> origin/master
		"""
		self._name = name
		return self


	def append_parameter(self, parameter_type, parameter_val):
		"""
<<<<<<< HEAD
		Method to append parameter with error

		Args:
			parameter_type (str): parameter name ("AA/TT", "AT/TA", ..., "init_XX", "symmetry", or "5term_TA (by regexp)")
			parameter_val (list): [param(float), param_min(float), param_max(float)]

		Returns:
			self
=======
		append parameter with error method
		@param parameter_type: parameter name ("AA/TT", "AT/TA", ..., "init_XX", "symmetry", or "5term_TA (by regexp)")
		@param parameter_val: parameter value list ([parameter, minimum parameter with error, maximum value with error] or parameter)
		@return self
>>>>>>> origin/master
		"""
		self._parameters[parameter_type] = [parameter_val for x in range(3)]
		self._change[parameter_type] = True
		return self


<<<<<<< HEAD
	def set_parameter(self, parameter_type, parameter_val=None):
		"""
		Method to set parameter with error

		Args:
			parameter_type (str): parameter name (all, "AA/TT", "AT/TA", "TA/AT", "CA/GT", "GT/CA", "CT/GA", "GA/CT", "CG/GC", "GC/CG", "GG/CC", "init_GC", "init_AT", "symmetry", or "5term_TA (by regexp)")
			parameter_val (float, optional): [param(float), param_min(float), param_max(float)] (Default: None)

		Returns:
			self
=======
	def set_parameter(self, parameter_type, parameter_val = None):
		"""
		set parameter with error method
		@param parameter_type: parameter name (all, "AA/TT", "AT/TA", "TA/AT", "CA/GT", "GT/CA", "CT/GA", "GA/CT", "CG/GC", "GC/CG", "GG/CC", "init_GC", "init_AT", "symmetry", or "5term_TA (by regexp)")
		@param parameter_val: parameter value ([parameter, minimum parameter with error, maximum value with error] or parameter)
		@return self
>>>>>>> origin/master
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
			sys.stderr.write("ERROR: undefined parameter_type in set_parameter() of ParameterData class ({0}).\n".format(parameter_type))
			sys.exit(1)

		return self


	def update_parameter_error(self, parameter_type, parameter_value):
		"""
<<<<<<< HEAD
		Method to calculate and set error

		Args:
			parameter_type (str): parameter name or "all"
			parameter_value (float or list): parameter value(float) for one parameter type or [param(float), param_min(float), param_max(float)] for "all"

		Returns:
			self
=======
		calculate and set error method
		@param parameter_type: parameter type or "all"
		@param parameter_value: parameter value for one parameter type or parameter value list for "all"
		@return self
>>>>>>> origin/master
		"""
		parameter_types = []
		parameter_values = []
		if parameter_type == "all":
			parameter_types = self._parameters.keys()
			parameter_values = parameter_value
		else:
			if parameter_type not in self._parameters.keys():
				sys.stderr("ERROR: undefined parameter_type at set_parameter_error() in Parameter class ({0}).\n".format(parameter_type))
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
<<<<<<< HEAD
		Method to change state for changing parameter

		Args:
			parameter_type (str): parameter type or "all"
			state (bool): True or False

		Returns:
			self
=======
		change state for changing parameter
		@param parameter_type: parameter type or "all"
		@param state: True or False
		@return self
>>>>>>> origin/master
		"""
		if parameter_type == "all":
			self._change = {k: state for k in self._change.keys()}
		else:
			self._change[parameter_type] = state
		return self


	def remove_parameter(self, parameter_name):
		"""
<<<<<<< HEAD
		Method to remove parameter

		Args:
			parameter_name (str): parameter name

		Returns:
			self
=======
		remove parameter
		@param parameter_name: parameter name
		@return self
>>>>>>> origin/master
		"""
		if parameter_name in self._parameters.keys():
			del(self._parameters[parameter_name])
		else:
			sys.stderr.write("ERROR: {0} does not found.\n".format(parameter_name))
			sys.exit(1)
		return self


	def clone(self):
		"""
<<<<<<< HEAD
		Method to return clone parameter set

		Returns:
			obj_Parameter
=======
		return clone(self)
		@return self
>>>>>>> origin/master
		"""
		return copy.deepcopy(self)


<<<<<<< HEAD
	def is_change(self, parameter_type=None):
		"""
		Method to return change state

		Args:
			parameter_type (str, optional): None or parameter type (Default: None)

		Returns:
			list for None or bool for parameter type
=======
	def get_name(self):
		"""
		return name
		@return name
		"""
		return self._name


	def is_change(self, parameter_type = None):
		"""
		return change state
		@param parameter_type: None or parameter type
		@return change state list for None, or boolean type value for parameter type
>>>>>>> origin/master
		"""
		if parameter_type is None:
			return self._change
		elif parameter_type in self._change.keys():
			return self._change[parameter_type]
		else:
			sys.stderr.write("ERROR: undefined parameter type at is_chage() in Parameter class ({0}).\n".format(parameter_type))
			sys.exit(1)


<<<<<<< HEAD
=======
	def get_parameter(self, parameter_type = None, data_type = "raw"):
>>>>>>> origin/master
		"""
		return parameter
		@param parameter_type: parameter type
		@param data_type: raw ([value, minimum value with error, maximum value with error]) or fix (value, error+/-)
		@return parameter (list for None(all) or float value for each parameter)
		"""
<<<<<<< HEAD
	def get_parameter(self, parameter_type=None, data_type="raw"):
		"""
		Method to return parameter

		Args:
			parameter_type (str, optional): parameter type (Default: None)
			data_type (str, optional): raw ([value, minimum value with error, maximum value with error]) or fix (value, error+/-) (Default: "raw")

		Returns:
			list for None(all) or float value for each parameter
		"""
=======
>>>>>>> origin/master
		values = {}
		if data_type == "raw":
			values = self._parameters
		elif data_type == "fix":
			error = {k: abs(v[1] - v[2]) / 2 for k, v in self._parameters.items()}
<<<<<<< HEAD
			min_val = [min(v[1:]) for v in self._parameters.values()]
			values = {k: [m + error[k], error[k]] for m, (k, v) in zip(min_val, self._parameters.items())}
=======
			values = {k: [v[1] + error[k], error[k]] for k, v in self._parameters.items()}
>>>>>>> origin/master
		elif data_type == "name":
			return self._parameters.keys()
		else:
			sys.stderr.write("ERROR: undefined data_type at get_parameter() in Parameter class.\n".format(data_type))
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
				sys.stderr.write("ERROR: undefined parameter_type in get_parameter() in Parameter class ({0}).\n".format(parameter_type))
				sys.exit(1)
<<<<<<< HEAD
=======



# =============== main =============== #
# if __name__ == '__main__':
# 	main()
>>>>>>> origin/master
