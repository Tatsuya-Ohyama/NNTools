#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pickle
import numpy as np
import statistics
import pandas as pd


# =============== class =============== #
class DataGroup:
	""" DataGroup class """
	def __init__(self, data = None):
		# member variables
		self._datas = pd.DataFrame()

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
			self._datas = pd.DataFrame(data)
		return self


	def add_data(self, data, data_type = None):
		"""
		add data method
		@param data: additional data
		@return self
		"""
		if data is not None:
			if data_type is None:
				data_type = "object"
			self._datas = self._datas.append(pd.Series(data, dtype = data_type), ignore_index = True)
		return self


	def set_label(self, label_type, label_list):
		"""
		set label for row or column
		@param label_type: "row" or "column"
		@param label_list: label list
		@return self
		"""
		if label_type == "row":
			self._datas.index = label_list
		elif label_type == "column":
			self._datas.columns = label_list
		else:
			sys.stderr.write("ERROR: undefined label_type at set_label() in DataGroup class.\n")
			sys.exit(1)


	def set_dtype(self, data_type, label_type = None, label_val = None):
		"""
		set dtype
		@param data_type: data type
		@param label_type: label value type for column ("label" or "index")
		@param label_val: label value for changing type
		@return self
		"""
		if label_type is not None:
			if label_type == "label":
				if label_val is not None:
					self._datas.at[label_val].astype(data_type)
				else:
					sys.stderr.write("ERROR: label_val is not defined at set_dtype() in DataGroup class.\n")
					sys.exit(1)
			elif label_type == "index":
				if label_val is not None:
					self._datas[label_val].astype(data_type)
				else:
					sys.stderr.write("ERROR: label_val is not defined at set_dtype() in DataGroup class.\n")
					sys.exit(1)
			else:
				sys.stderr.write("ERROR: undefined label_type at set_dtype() in DataGroup class.\n")
				sys.exit(1)
		else:
			self._datas = self._datas.astype(data_type)
		return self


	def update_data(self, label_type, index, column, new_val):
		"""
		update value at specified x
		@param label_type: "label" or "index"
		@param index: row position
		@param column: column position
		@param new_val: value
		@return self
		"""
		if labe_type == "label":
			self._datas.at[index, column] = new_val
		elif label_type == "index":
			self._datas.iat[index, column] = new_val
		else:
			sys.stderr.write("ERROR: undefined label_type at update_data() in DataGroup class.\n")
			sys.exit(1)
		return self


	def get_label(self, label_type):
		"""
		return label list
		@param label_type: "row" or "column"
		@return label list
		"""
		if label_type == "row":
			return self._datas.index
		elif label_type == "column":
			return self._datas.columns
		else:
			sys.stderr.write("ERROR: undefined label_type at get_label() in DataGroup class.\n")
			sys.exit(1)


	def get_dtype(self):
		"""
		return dtype
		@return dtype
		"""
		return self._datas.dtypes


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
		x = self._datas[self._datas.columns[0]]
		y = self._datas[self._datas.columns[1]]
		coeff = np.polyfit(x, y, deg)
		r2 = np.corrcoef(x, y)[0,1]
		return coeff.tolist() + [r2]


# =============== main =============== #
# if __name__ == '__main__':
# 	main()
