#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pickle
import copy

# =============== variables =============== #
parameter_list = ["AA/TT", "AT/TA", "TA/AT", "CA/GT", "GT/CA", "CT/GA", "GA/CT", "CG/GC", "GC/CG", "GG/CC", "init_GC", "init_AT", "symmetry", "5term_TA"]
complement_list = {
	"AA/TT": "TT/AA",
	"AT/TA": "AT/TA",
	"TA/AT": "TA/AT",
	"CA/GT": "TG/AC",
	"GT/CA": "AC/TG",
	"CT/GA": "AG/TC",
	"GA/CT": "TC/AG",
	"CG/GC": "CG/GC",
	"GC/CG": "GC/CG",
	"GG/CC": "CC/GG"
}


# =============== class =============== #
class Parameter:
	""" Parameter class for one energy type """
	def __init__(self):
		# member variables
		self._name = ""
		self._parameters = {
			"AA/TT":    [0.0, 0.0],
			"AT/TA":    [0.0, 0.0],
			"TA/AT":    [0.0, 0.0],
			"CA/GT":    [0.0, 0.0],
			"GT/CA":    [0.0, 0.0],
			"CT/GA":    [0.0, 0.0],
			"GA/CT":    [0.0, 0.0],
			"CG/GC":    [0.0, 0.0],
			"GC/CG":    [0.0, 0.0],
			"GG/CC":    [0.0, 0.0],
			"init_GC":  [0.0, 0.0],
			"init_AT":  [0.0, 0.0],
			"symmetry": [0.0, 0.0],
			"5term_TA": [0.0, 0.0]
		}
		self._change = {k: True for k in self._parameters.keys()}


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


	def set_parameter(self, parameter_type, parameter_val = None):
		"""
		set parameter method
		@param parameter_type: parameter name (all, "AA/TT", "AT/TA", "TA/AT", "CA/GT", "GT/CA", "CT/GA", "GA/CT", "CG/GC", "GC/CG", "GG/CC", "init_GC", "init_AT", "symmetry", or "5term_TA")
		@param parameter_val: parameter value ([parameter, +error, -error] or parameter)
		@return self
		"""
		if parameter_type == "all":
			# すべての場合、そのまま引き受ける
			self._parameters = {x: [float(y[idx]) for idx in range(2)] if self._change[x] else self._parameters[x] for x, y in parameter_val.items()}
			return self

		if parameter_type not in parameter_list:
			# パラメータが相補鎖の形式で含まれる場合
			parameter_type = [k for k, v in complement_list.items()][0]

		if parameter_type in parameter_list:
			# パラメータ名が含まれる場合
			if self._change[parameter_type]:
				if type(parameter_val) == list:
					self._parameters[parameter_type] = [float(x) for x in parameter_val]
				else:
					self._parameters[parameter_type] = [float(parameter_val), 0.0, 0.0]
		else:
			sys.stderr.write("ERROR: undefined parameter_type in set_parameter() of ParameterData class ({0}).\n".format(parameter_type))
			sys.exit(1)

		return self

	def set_parameter_error(self, parameter_type, parameter_value):
		"""
		set parameter error
		@param parameter_type: parameter type or "all"
		@param parameter_value: parameter value for one parameter type or parameter value list for "all"
		@return self
		"""
		parameter_types = []
		parameter_values = []
		if parameter_type == "all":
			parameter_types = parameter_list
			parameter_values = parameter_value
		else:
			if parameter_type not in parameter_list:
				sys.stderr("ERROR: undefined parameter_type at set_parameter_error() in Parameter class ({0}).\n".format(parameter_type))
				sys.exit(1)
			parameter_types = [parameter_type]
			parameter_values = [parameter_value]

		for idx in range(len(parameter_types)):
			parameter_values[parameter_types[idx]][0] -= self._parameters[parameter_types[idx]][0]
			if self._parameters[parameter_types[idx]][1] < abs(parameter_values[parameter_types[idx]][0]):
				# update negative error
				self._parameters[parameter_types[idx]][1] = parameter_values[parameter_types[idx]][0]
		return self


	def set_change_stat(self, parameter_type, state):
		"""
		change state for changing parameter
		@param parameter_type: parameter type
		@param state: True or False
		@return self
		"""
		self._change[parameter_type] = state
		return self


	def remove_parameter(self, parameter_name):
		"""
		remove parameter
		@param parameter_name: parameter name
		@return self
		"""
		if parameter_name in self._parameters.keys():
			del(self._parameters[parameter_name])
		else:
			sys.stderr.write("ERROR: {0} does not found.\n".format(parameter_name))
			sys.exit(1)
		return self


	def clone(self):
		"""
		return clone(self)
		@return self
		"""
		return copy.deepcopy(self)


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
		"""
		if parameter_type is None:
			return self._change
		elif parameter_type in self._change.keys():
			return self._change[parameter_type]
		else:
			sys.stderr.write("ERROR: undefined parameter type.\n")
			sys.exit(1)


	def get_parameter(self, parameter_type = None):
		"""
		return parameter
		@param parameter_type: parameter type
		@return parameter (list for None(all) or float value for each parameter)
		"""
		if parameter_type in parameter_list:
			# パラメータタイプが存在する場合
			return self._parameters[parameter_type]
		elif parameter_type in complement_list.values():
			parameter_key = [k for k, v in complement_list.items()][0]
			return self._parameters[parameter_key]
		else:
			# すべてのパラメータが指定された場合
			return self._parameters



# =============== main =============== #
# if __name__ == '__main__':
# 	main()
