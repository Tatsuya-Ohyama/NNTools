#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pickle


# =============== class =============== #
class ParameterValue:
	""" Sequence class """
	def __init__(self, name):
		# member variables
		self._name = ""
		self._parameters = {
			"AA/TT": 0.0,
			"AT/TA": 0.0,
			"TA/AT": 0.0,
			"CA/GT": 0.0,
			"GT/CA": 0.0,
			"CT/GA": 0.0,
			"GA/CT": 0.0,
			"CG/GC": 0.0,
			"GC/CG": 0.0,
			"GG/CC": 0.0,
			"init_GC": 0.0,
			"init_AT": 0.0,
			"symmetry": 0.0,
			"5term_TA": 0.0
		}


		# initiation
		self.set_name(name)


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


	def set_parameter(self, parameter_type, parameter_val = None):
		"""
		set parameter method
		@param parameter_type: parameter name (all, "AA/TT", "AT/TA", "TA/AT", "CA/GT", "GT/CA", "CT/GA", "GA/CT", "CG/GC", "GC/CG", "GG/CC", "init_GC", "init_AT", "symmetry", or "5term_TA")
		@param parameter_val: parameter value (parameter list [dH, dS, dG])
		@return self
		"""
		if parameter_type in self._parameters.keys:
			# パラメータ名が含まれる場合
			self._parameters[parameter_type] = [float(x) for x in parameter_val]
		elif parameter_type == "all":
			# すべての場合、そのまま引き受ける
			self._parameters = {x: [float(z) for z in y] for x, y in parameter_val.items()}
		else:
			sys.stderr.write("ERROR: undefined parameter_type in set_parameter() of ParameterData class.\n")
			sys.exit(1)

		return self


	def get_name(self):
		"""
		return name
		@return name
		"""
		return self._name


	def get_parameter(self, parameter_type = None):
		"""
		return parameter
		@param parameter_type: parameter type (None (all), dH, dS, or dG)
		@return parameter (None(all): [dH, dS, dG], or float value of dH, dS, and dG)
		"""
		if parameter_type in self._parameters.keys():
			# パラメータタイプが存在する場合
			return self._parameters[parameter_type]
		elif parameter_type == "all":
			# すべてのパラメータが指定された場合
			return self._parameters
		else:
			sys.stderr.write("ERROR: undefined parameter_type in get_parameter() of ParameterData class.\n")
			sys.exit(1)



# =============== main =============== #
# if __name__ == '__main__':
# 	main()
