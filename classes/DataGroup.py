#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pickle
import numpy as np
import statistics


# =============== class =============== #
class DataGroup:
	""" LinearFunction class """
	def __init__(self, data = None):
		# member variables
		self._datas = []
		self._slope = 0.0
		self._intercept = 0.0
		self._r2 = 0.0

		# initiation
		self.set_data(data)


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


	def set_data(self, data):
		"""
		set data method
		@param data: new data
		@return: self
		"""
		if data is not None:
			self._datas = data
		return self


	def add_data(self, data):
		"""
		add data method
		@param data: additional data
		@return self
		"""
		if data is not None:
			self._datas.append(data)
		return self


	def get_data(self):
		"""
		return data
		@return data list
		"""
		return self._datas


	def get_factor(self,deg = 1):
		"""
		return list of optimized slope factor, intercept, and R2 value
		@param deg: Degree of the fitting polynomial
		@return self: [polynomial_coefficients, residuals, r2]
		"""
		x = np.array([x[0] for x in self._datas], dtype = np.float)
		y = np.array([x[1] for x in self._datas], dtype = np.float)
		coeff = np.polyfit(x, y, deg).tolist()
		r2 = np.corrcoef(x, y)[0,1]
		return coeff + [r2]


# =============== main =============== #
# if __name__ == '__main__':
# 	main()
