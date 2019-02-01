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


	def rename_label(self, direction, old_name, new_name):
		"""
		rename label
		@param direction: "row" or "column"
		@param old_name: oldname
		@param new_name: new name
		@return self
		"""
		if direction == "row":
			self._datas = self._datas.rename(index = {old_name : new_name})
		elif direction == "column":
			self._datas = self._datas.rename(columns = {old_name : new_name})
		else:
			sys.stderr.write("ERROR: undefined direction at rename_label() in DataGroup().\n")
			sys.exit(1)
		return self


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


	def get_diff(self, label_type, label_x, label_y):
		"""
		return difference^2
		@param label_type: "index" or "label"
		@param label_x
		@param label_y
		@return difference^2 and sum of difference^2
		"""
		diff = (self._datas[label_x].values - self._datas[label_y].values) ** 2
		return [diff, np.sum(diff)]




	def save_csv(self, output_file):
		"""
		save to csv file
		@param output_file: output file path
		@return self
		"""
		self._datas.to_csv(output_file, header = True, index = True)
		return self



class DataGroup2:
	""" DataGroup class (Ver. 2) """
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
