#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os


# =============== function =============== #
def check_overwrite(file):
	"""
	function to check for overwrite
	@param file(str): target file path
	@return None
	"""
	if os.path.exists(file):
		# Warn if file exists
		sys.stderr.write("WARN: %s exists. Overwrite it? (y/N): " % file)
		sys.stderr.flush()
		user = sys.stdin.readline().strip()

		if user != "y":
			# When there is no permission to overwrite
			sys.exit(0)
		else:
			# If there is permission to overwrite, delete it in advance considering conflicts
			os.remove(file)



def check_exist(path, mode, flag_exit = True):
	"""
	function to check for file existence
	@param path(str): target file path
	@param mode(int): 1(existence) / 2(existence for file) / 3(existence for dir)
	@param flag_exit(bool): Exit if not present (Default: True)
	@param (bool) or exit(None)
	"""
	if mode == 1:
		if not os.path.exists(path):
			sys.stderr.write("ERROR: No such path (%s)\n" % path)
			if flag_exit:
				sys.exit(1)
			else:
				return False
	elif mode == 2:
		if not os.path.isfile(path):
			sys.stderr.write("ERROR: No such file (%s)\n" % path)
			if flag_exit:
				sys.exit(1)
			else:
				return False
	elif mode == 3:
		if not os.path.isdir(path):
			sys.stderr.write("ERROR: No such directory (%s)\n" % path)
			if flag_exit:
				sys.exit(1)
			else:
				return False
	else:
		sys.stderr.write("ERROR: Subroutine error: Not specified mode\n")
		sys.exit(1)
	return True


# =============== main =============== #
# if __name__ == '__main__':
# 	main()
