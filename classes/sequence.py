#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re

from classes.parameter import Parameter
<<<<<<< HEAD
=======

>>>>>>> NNcalcTool.py


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


	@property
	def name(self):
		return self._name

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


	def set_sequence(self, sequence, base_pairs):
		"""
		Method to set sequence

		Args:
			sequence (str): sequence
			base_pairs (dict): base pairs

		Returns:
			self
		"""
		if sequence is not None:
			self._sequence = list(sequence)

			misstype = [base for base in self._sequence if base not in base_pairs.keys()]
			if len(misstype) != 0:
				sys.stderr.write("ERROR: misstype of base: {0} in `{1}`.\n".format(misstype, self._name))
				sys.exit(1)

			self._complement = [base_pairs[base] for base in self._sequence]
			self._is_self_complement = self._sequence == list(reversed(self._complement))
		return self


	def get_sequence(self, sequence_type="list"):
		"""
		Method to return sequence

		Args:
			sequence_type (str, optional): "list" or "string" (Default: "list")

		Returns:
			str: sequence
		"""
		if sequence_type == "string":
			return "".join(self._sequence)
		elif sequence_type == "list":
			return self._sequence
		else:
			sys.stderr.write("ERROR: undefined sequence_type at get_sequence() in Sequence class.\n")
			sys.exit(1)


	def get_freq(self, parameter_types, base_pairs, flag_cache=True):
		"""
		Method to return pair frequency

		Args:
			parameter_types (list or objParameter): list for parameter types or objParameter
			base_pairs (dict): base pairs
			flag_cache (bool, optional): use cache data (Default: True)

		Returns:
			list: pair frequency
		"""
		if type(parameter_types) == Parameter:
			parameter_types = parameter_types.get_parameter(data_type="name")

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

				elif param.startswith("reg:"):
					# regexp parameter for repeated
					re_exps = param.replace("reg:", "").split("/")
					flag_match = []
					for regexp, sequence in zip(re_exps, ["".join(self._sequence), "".join(self._complement)]):
						flag_match.append([x.span() for x in re.finditer(regexp, sequence)])

					if flag_match[0] == flag_match[1]:
						# match sequence and complementary sequence
						self._cache_freq[param] = len(flag_match[0])

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
		Method to return energy value

		Args:
			obj_parameter (objParameter): Paramete object
			base_pairs (dict): base pairs

		Returns:
			float: energy value
		"""
		# calculate energy
		energy = 0.0
		freq = self.get_freq(obj_parameter, base_pairs)
		energy = sum([cnt_pair * obj_parameter.get_parameter(parameter_type)[0] for parameter_type, cnt_pair in freq.items()])

		return energy
