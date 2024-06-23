#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re

from mods.parameter import Parameter


# =============== class =============== #
class Sequence:
	""" Sequence class """
	def __init__(self, name):
		# member variables
		self._name = ""
		self._sequence = []
		self._complement = []
		self._is_self_complement = False

		self._cache_parameter_types = []
		self._cache_freq = {}

		# initiation
		self.set_name(name)


	@property
	def name(self):
		return self._name

	@property
	def sequence(self):
		return self._sequence

	@property
	def complement(self):
		return self._complement

	@property
	def is_complement(self):
		return self._is_self_complement


	def set_name(self, name):
		"""
		Method to set name

		Args:
			name (str): name

		Returns:
			self
		"""
		self._name = name
		return self


	def set_sequence(self, sequence, complement=None):
		"""
		Method to set sequence

		Args:
			sequence (str): sequence
			complement (str,optional): complement sequence

		Returns:
			self
		"""
		self._sequence = list(sequence)

		if complement is not None:
			self._complement = list(complement)
			self._is_self_complement = self._sequence == list(reversed(self._complement))

		return self


	def generate_complement(self, base_pairs):
		"""
		Method to generate complement sequence from base pair information

		Args:
			base_pairs (dict): base pairs

		Returns:
			self
		"""
		misstype = [base for base in self._sequence if base not in base_pairs.keys()]
		if len(misstype) != 0:
			sys.stderr.write("ERROR: Misstype of base: {0} in `{1}`.\n".format(misstype, self._name))
			sys.exit(1)

		self._complement = [base_pairs[base] for base in self._sequence]
		self._is_self_complement = self._sequence == list(reversed(self._complement))

		return self


	def get_sequence(self, data_type="list"):
		"""
		Method to return sequence

		Args:
			data_type (str, optional): "list" or "string" (Default: "list")

		Returns:
			str: sequence
		"""
		if data_type == "string":
			return "".join(self._sequence)
		elif data_type == "list":
			return self._sequence
		else:
			sys.stderr.write("ERROR: Undefined data_type at get_sequence() in Sequence class.\n")
			sys.exit(1)


	def get_complement(self, data_type="list"):
		"""
		Method to return complement sequence

		Args:
			data_type (str, optional): "list" or "string" (Default: "list")

		Returns:
			list or str: sequence
		"""
		if data_type == "string":
			return "".join(self._complement)
		elif data_type == "list":
			return self._complement
		else:
			sys.stderr.write("ERROR: Undefined sequence_type at get_sequence() in Sequence class.\n")
			sys.exit(1)


	def get_freq(self, obj_parameter, flag_cache=True):
		"""
		Method to return pair frequency

		Args:
			obj_parameter(Parameter object): parameter object
			flag_cache (bool, optional): use cache data (Default: True)

		Returns:
			list: pair frequency
		"""
		parameter_types = obj_parameter.get_parameter(data_type="name")

		if flag_cache and parameter_types == self._cache_parameter_types:
			# flag_cache is True and condition is the same => use cache
			return self._cache_freq

		else:
			# build parameter list
			self._cache_parameter_types = parameter_types
			self._cache_freq = {param: 0 for param in parameter_types}

			sequence = self.get_sequence("string")
			complement = self.get_complement("string")

			# check parameter_types
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
					obj_match1 = re.search(re_exps[0], sequence)
					obj_match2 = re.search(re_exps[1], complement)

					if obj_match1 and obj_match2:
						# match sequence and complementary sequence
						self._cache_freq[param] += 1

				elif param.startswith("reg:"):
					# regexp parameter for repeated
					re_exps = param.replace("reg:", "").split("/")
					flag_match = []
					for regexp, sequence in zip(re_exps, ["".join(self._sequence), "".join(self._complement)]):
						flag_match.append([x.span() for x in re.finditer(regexp, sequence)])

					if flag_match[0] == flag_match[1]:
						# match sequence and complementary sequence
						self._cache_freq[param] = len(flag_match[0])

				elif "/" in param:
					# search in correct order
					query_seq1, query_seq1c = param.split("/", 1)
					length_pattern = len(query_seq1)
					pairs = ["/".join([sequence[i:i+length_pattern], complement[i:i+length_pattern]]) for i in range(len(sequence)-length_pattern+1)]
					self._cache_freq[param] += pairs.count(param)

					# search in reverse order
					if obj_parameter.one_direction:
						continue

					query_seq2, query_seq2c = "".join(list(reversed(query_seq1c))), "".join(list(reversed(query_seq1)))
					if query_seq1 == query_seq2:
						# skip reverse order search at symmetry parameter
						continue

					pairs_reverse = [v[::-1] for v in pairs]
					self._cache_freq[param] += pairs_reverse.count(param)

				else:
					sys.stderr.write("ERROR: Undefined parameter at get_freq() in Sequence class ({0}).\n".format(param))
					sys.exit(1)

			return self._cache_freq


	def get_energy(self, obj_parameter, data_type="raw"):
		"""
		Method to return energy value

		Args:
			obj_parameter (objParameter): Paramete object
			data_type (str, optional): raw ([value, minimum value with error, maximum value with error]) or fix (value, error+/-) (Default: "raw")

		Returns:
			float: energy value
		"""
		# calculate energy
		energy = 0.0
		freq = self.get_freq(obj_parameter)
		energy = sum([cnt_pair * obj_parameter.get_parameter(parameter_type, data_type)[0] for parameter_type, cnt_pair in freq.items()])

		return energy
