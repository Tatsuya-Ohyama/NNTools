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


	def add_data(self, direction, name, data, data_type = None):
		"""
		add data method
		@param direction: "row" or "column"
		@param name: column name
		@param data: additional data
		@param data_type: dtype (Default: "object")
		@return self
		"""
		if data is not None:
			if data_type is None:
				data_type = "object"
			if direction == "row":
				# add row
				if self._datas.shape == (0,0):
					# empty DataFrame
					self._datas = pd.DataFrame([data], index = name, dtype = data_type)
				else:
					# add
					new_row = pd.DataFrame([data], index = name, dtype = data_type)
					new_row.columns = self._datas.columns
					self._datas = pd.concat([self._datas, new_row], axis = 0)
			elif direction == "column":
				# add column
				if self._datas.shape == (0,0):
					# empty DataFrame
					self._datas = pd.DataFrame(data, columns = name, dtype = data_type)
				else:
					# add
					new_column = pd.DataFrame(data, columns = name, dtype = data_type)
					new_column.index = self._datas.index
					self._datas = pd.concat([self._datas, new_column], axis = 1)
		return self


	def set_label(self, label_direction, label_list):
		"""
		set label for row or column
		@param label_direction: "row" or "column"
		@param label_list: label list
		@return self
		"""
		if label_direction == "row":
			self._datas.index = label_list
		elif label_direction == "column":
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
					self._datas[label_val] = self._datas[label_val].astype(data_type)
				else:
					sys.stderr.write("ERROR: label_val is not defined at set_dtype() in DataGroup class.\n")
					sys.exit(1)
			elif label_type == "index":
				if label_val is not None:
					self._datas[self._datas.columns[label_val]] = self._datas[self._datas.columns[label_val]].astype(data_type)
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
		if label_type == "label":
			if index in self._datas.index and column in self._datas.columns:
				self._datas.at[index, column] = new_val
			else:
				sys.stderr.write("ERROR: index '{0}', column '{1}' does not found.\n".format(index, column))
				sys.exit(1)
		elif label_type == "index":
			self._datas.iat[index, column] = new_val
		else:
			sys.stderr.write("ERROR: undefined label_type at update_data() in DataGroup class.\n")
			sys.exit(1)
		return self


	def remove_data(self, label_direction, label_type, label_name):
		"""
		remove data
		@param label_direction: "row" or "column"
		@param label_type: "index" or "label"
		@param label_name: label name
		@return self
		"""
		if label_direction == "row":
			# 行の削除
			if label_type == "index":
				# インデックス指定
				self._datas = self._datas.drop(self._datas.index[label_name], axis = 0)
			elif label_type == "label":
				# ラベル指定
				self._datas = self._datas.drop(label_name, axis = 0)
			else:
				sys.stderr.write("ERROR: undefined label_type at remove_data() in DataGroup class.\n")
				sys.exit(1)
		elif label_direction == "column":
			# 列の削除
			if label_type == "index":
				# インデックス指定
				self._datas = self._datas.drop(self._datas.columns[label_name], axis = 1)
			elif label_type == "label":
				# ラベル指定
				self._datas = self._datas.drop(label_name, axis = 1)
			else:
				sys.stderr.write("ERROR: undefined label_type at remove_data() in DataGroup class.\n")
				sys.exit(1)
		else:
			sys.stderr.write("ERROR: undefined label_direction at remove_data() in DataGroup class.\n")
			sys.exit(1)
		return self


	def get_label(self, label_direction):
		"""
		return label list
		@param label_direction: "row" or "column"
		@return label list
		"""
		if label_direction == "row":
			return self._datas.index
		elif label_direction == "column":
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


	def get_factor(self, label_type, label_x, label_y, deg = 1):
		"""
		return list of optimized slope factor, intercept, and R2 value
		@param label_type: "index" or "label"
		@param label_x
		@param label_y
		@param deg: Degree of the fitting polynomial
		@return self: [polynomial_coefficients, residuals, r2]
		"""
		x = None
		y = None
		if label_type == "index":
			x = self._datas[self._datas.columns[label_x]]
			y = self._datas[self._datas.columns[label_y]]
		elif label_type == "label":
			x = self._datas[label_x]
			y = self._datas[label_y]
		else:
			sys.stderr.write("ERROR: undefined label_type at get_factor() in DataGroup class.\n")
			sys.exit(1)

		coeff = np.polyfit(x, y, deg)
		r1 = np.corrcoef(x, y)[0,1]
		r2 = r1 ** 2
		return coeff.tolist() + [r1, r2]


	def save_csv(self, output_file):
		"""
		save to csv file
		@param output_file: output file path
		@return self
		"""
		self._datas.to_csv(output_file, header = True, index = True)
		return self


# =============== main =============== #
# if __name__ == '__main__':
# 	main()
