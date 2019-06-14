#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pickle
import re

from classes.Parameter import Parameter


# =============== class =============== #
class Sequence:
	""" Sequence class """
	def __init__(self, name):
		# member variables
		self._name = ""
		self._sequence = ""
		self._complement = ""
		self._is_self_complement = False

		self._cache_parameter_types = []
		self._cache_base_pairs = {}
		self._cache_freq = {}


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

			self._complement = [base_pairs[base] for base in self._sequence]
			self._is_self_complement = self._sequence == list(reversed(self._complement))
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


	def get_freq(self, parameter_types, base_pairs, flag_cache = True):
		"""
		return pair frequency
		@param parameter_types: list for parameter types or parameter object
		@param base_pairs: dict for base pairs
		@param flag_cache: use cache data (Default: True)
		@return pair frequency list
		"""
		if type(parameter_types) == Parameter:
			parameter_types = parameter_types.get_parameter(data_type = "name")

		if flag_cache and parameter_types == self._cache_parameter_types and base_pairs == self._cache_base_pairs:
			# flag_cache is True and condition is the same => use cache
			return self._cache_freq

		else:
			# build or rebuild parameter list
			self._cache_parameter_types = parameter_types
			self._cache_base_pairs = base_pairs
			self._cache_freq = {param: 0 for param in parameter_types}

			for param in parameter_types:
				# special penalty
				if param.startswith("init"):
					# initiation parameter
					list_init = [self._sequence[0] + self._complement[0], self._complement[0] + self._sequence[0]]
					if param.replace("init_", "") in list_init:
						self._cache_freq[param] += 1

				elif param.startswith("length"):
					# length parameter
					self._cache_freq[param] = len(self._sequence)

				elif param.startswith("symmetry"):
					# symmetry parameter
					if self._is_self_complement:
						self._cache_freq[param] += 1

				elif param.startswith("re:"):
					# regexp parameter
					re_exps = param.replace("re:", "").split("/")
					flag_match = 0
					for regexp, sequence in zip(re_exps, ["".join(self._sequence), "".join(self._complement)]):
						if re.search(regexp, sequence):
							flag_match += 1
					if flag_match == 2:
						# match sequence and complementary sequence
						self._cache_freq[param] += 1

				elif "/" not in param:
					sys.stderr.write("ERROR: undefined parameter at get_freq() in Sequence class ({0}).\n".format(param))
					sys.exit(1)

			for base_idx in range(len(self._sequence) - 1):
				pair_type = "/".join(["".join(self._sequence[base_idx : base_idx + 2]), "".join(self._complement[base_idx : base_idx + 2])])
				if pair_type not in parameter_types:
					# if not in parameter_type, generate reversed parameter_type
					pair_type = "/".join(["".join(reversed(self._complement[base_idx : base_idx + 2])), "".join(reversed(self._sequence[base_idx : base_idx + 2]))])
				self._cache_freq[pair_type] += 1

			return self._cache_freq


	def get_energy(self, obj_parameter, base_pairs):
		"""
		return energy value
		@param base_pairs: dict for base pairs
		@param obj_parameter: Parameter object
		@return energy_value
		"""
		# calculate energy
		energy = 0.0
		freq = self.get_freq(obj_parameter, base_pairs)
		energy = sum([cnt_pair * obj_parameter.get_parameter(parameter_type)[0] for parameter_type, cnt_pair in freq.items()])

		return energy


# =============== main =============== #
# if __name__ == '__main__':
# 	main()
