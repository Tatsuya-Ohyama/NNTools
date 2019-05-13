#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pickle

from classes.Parameter import Parameter


# =============== class =============== #
class Sequence:
	""" Sequence class """
	def __init__(self, name):
		# member variables
		self._name = ""
		self._sequence = ""
		self._energy_type = ""
		self._is_self_complement = False
		self._base_pair = {}
		self._parameter_list = []


		# initiation
		self.set_name(name)


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
		@return self (for chain method)
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


	def append_parameter(self, parameter):
		"""
		append parameter to parameter list
		@param parameter: parameter name
		@return self
		"""
		self._parameter_list.append(parameter)
		return self


	def set_parameter_type(self, parameter_list):
		"""
		set parameter list
		@param parameter_list: list of parameters
		@return self
		"""
		self._parameter_list = parameter_list
		return self


	def set_energy_type(self, energy_type):
		"""
		set energy_type
		@param energy_type: dH, dS, dG
		@return self
		"""
		self._energy_type = energy_type
		return self


	def set_sequence(self, sequence, base_pairs):
		"""
		set sequence method
		@param sequence: sequence
		@param base_pairs: dict for base pairs
		@return: self
		"""
		if sequence is not None:
			self._sequence = list(sequence)

			misstype = [base for base in self._sequence if base not in base_pairs.keys()]
			if len(misstype) != 0:
				sys.stderr.write("ERROR: misstype of base: {0}.\n".format(misstype))
				sys.exit(1)

			complement = [base_pairs[base] for base in self._sequence]
			self._is_self_complement = self._sequence == list(reversed(complement))
		return self


	def get_name(self):
		"""
		return name
		@return name
		"""
		return self._name


	def get_sequence(self, sequence_type = "list"):
		"""
		return sequence
		@param sequence_type: "list" or "string" (Default: list)
		@return sequence
		"""
		if sequence_type == "string":
			return "".join(self._sequence)
		elif sequence_type == "list":
			return self._sequence
		else:
			sys.stderr.write("ERROR: undefined sequence_type at get_sequence() in Sequence class.\n")
			sys.exit(1)


	def is_complement(self):
		"""
		return self complement or not
		@return True or False
		"""
		return self._is_self_complement


	def get_energy_type(self):
		"""
		return energy_type
		@return energy_type
		"""
		return self._energy_type


	def get_freq(self, parameter_types, base_pairs):
		"""
		return pair frequency
		@param parameter_types: list for parameter types
		@param base_pairs: dict for base pairs
		@return pair frequency list
		"""
		freq = {param: 0 for param in parameter_types}

		init_param = [param for param in parameter_types if "init" in param]
		init_param_data = [param.replace("init", "").replace("_", "") for param in init_param]
		for idx, param in enumerate(init_param_data):
			if param == "":
				# single initiation
				freq[init_param[idx]] += 1
			elif self._sequence[0] in list(param):
				# multiple initiation
				freq[init_param[idx]] += 1

		end_param = [param for param in parameter_types if "end" in param]
		end_param_data = [param.replace("end", "").replace("_", "") for param in end_param]
		for idx, param in enumerate(end_param_data):
			if param == "":
				# single end
				freq[end_param[idx]] += 1
			elif self._sequence[-1] in list(param):
				# multiple end
				freq[end_param[idx]] += 1

		if "symmetry" in freq.keys() and self._is_self_complement:
			freq["symmetry"] += 1

		if "length" in freq.keys():
			freq["length"] = len(self._sequence)

		for base_idx in range(len(self._sequence) - 1):
			pair_forward = [self._sequence[base_idx], self._sequence[base_idx + 1]]
			pair_reverse = [base_pairs[self._sequence[base_idx]], base_pairs[self._sequence[base_idx + 1]]]
			pair_type = "/".join(["".join(pair_forward), "".join(pair_reverse)])
			if pair_type not in self._parameter_list:
				pair_type = "/".join(["".join(reversed(pair_reverse)), "".join(reversed(pair_forward))])
			freq[pair_type] += 1

		return freq


	def get_energy(self, obj_parameter, base_pairs):
		"""
		return energy value
		@param base_pairs: dict for base pairs
		@param obj_parameter: Parameter object
		@return energy_value
		"""
		# calculate energy
		energy = 0.0

		init_param = [param for param in obj_parameter.get_parameter().keys() if "init" in param]
		init_param_data = [param.replace("init", "").replace("_", "") for param in init_param]
		for idx, param in enumerate(init_param_data):
			if param == "":
				# single initiation
				energy += obj_parameter.get_parameter(init_param[idx])[0]
			elif self._sequence[0] in list(param):
				# multiple initiation
				energy += obj_parameter.get_parameter(init_param[idx])[0]

		end_param = [param for param in obj_parameter.get_parameter().keys() if "end" in param]
		end_param_data = [param.replace("end", "").replace("_", "") for param in end_param]
		for idx, param in enumerate(end_param_data):
			if param == "":
				# single end
				energy += obj_parameter.get_parameter(end_param[idx])[0]
			elif self._sequence[-1] in list(param):
				# multiple end
				energy += obj_parameter.get_parameter(end_param[idx])[0]

		for base_idx in range(len(self._sequence) - 1):
			pair_forward = [self._sequence[base_idx], self._sequence[base_idx + 1]]
			pair_reverse = [base_pairs[self._sequence[base_idx]], base_pairs[self._sequence[base_idx + 1]]]
			pair_type = "/".join(["".join(pair_forward), "".join(pair_reverse)])
			if pair_type not in self._parameter_list:
				pair_type = "/".join(["".join(reversed(pair_reverse)), "".join(reversed(pair_forward))])

			energy += obj_parameter.get_parameter(pair_type)[0]

		if self._is_self_complement and "symmetry" in obj_parameter.get_parameter().keys():
			energy += obj_parameter.get_parameter("symmetry")[0]

		if "length" in obj_parameter.get_parameter().keys():
			energy += obj_parameter.get_parameter("length")[0] * len(self._sequence)

		return energy


# =============== main =============== #
# if __name__ == '__main__':
# 	main()
