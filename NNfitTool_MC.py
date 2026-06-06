#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
NNfitTool_MC.py -  Fitting program for parameters of Nearest Neighbor method using Least square solver and Monte Carlo method
"""


# =============== modules =============== #
import sys, signal
sys.dont_write_bytecode = True
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
import re
import os
import numpy as np
import pandas as pd
from scipy.optimize import least_squares



# =============== constant =============== #
PROGRAM_NAME = os.path.basename(sys.argv[0])
EPS = 1e-4
SEED = 1234
VERSION = "1.0"
TARGETS = ["dH", "dS", "dG"]



# =============== function =============== #
def read_exp(exp_file):
	"""
	Function to read experimental data

	Args:
		exp_file (str): .csv file with sequence and experimental values

	Returns:
		pd.DataFrame: experimental dataframe
	"""
	df = pd.read_csv(exp_file, header=0, index_col=False)

	list_col_sequence = [v for v in df.columns if v.lower().startswith("sequence")]
	rename_map = {}
	if len(list_col_sequence) == 2:
		rename_map[list_col_sequence[0]] = "Sequence A"
		rename_map[list_col_sequence[1]] = "Sequence B"

	elif len(list_col_sequence) == 1:
		rename_map[list_col_sequence[0]] = "Sequence A"

	df = df.rename(columns=rename_map)

	if "Comment" not in df.columns:
		df.loc[:,"Comment"] = ""

	if len(list_col_sequence) == 0:
		raise ValueError("Experimental .csv file must contain either 'Sequence' or both 'Sequence A' and 'Sequence B'.")

	elif len(list_col_sequence) == 1:
		df.loc[:,"Sequence B"] = None

	required_numeric = [
		"dH", "dH (error)",
		"dS", "dS (error)",
		"dG", "dG (error)",
	]
	missing = [c for c in required_numeric if c not in df.columns]
	if missing:
		raise ValueError(f"Missing columns in experimental CSV: {missing}")

	df.loc[:,"Comment"] = df.loc[:,"Comment"].astype(str)
	df.loc[:,"Sequence A"] = df.loc[:,"Sequence A"].astype(str)

	if len(list_col_sequence) == 1:
		df.loc[:,"Sequence B"] = df.loc[:,"Sequence B"].where(df["Sequence B"].notna(), None)
		df.loc[:,"Sequence B"] = df.loc[:,"Sequence B"].apply(lambda x: str(x) if x is not None else None)

	for c in required_numeric:
		df.loc[:,c] = pd.to_numeric(df.loc[:,c], errors="coerce")

	return df.loc[:,[
		"Comment", "Sequence A", "Sequence B",
		"dH", "dH (error)",
		"dS", "dS (error)",
		"dG", "dG (error)",
	]]


def read_params(param_file):
	"""
	Function to read parameter data

	Args:
		param_file (str): .csv file with parameter file

	Returns:
		pd.DataFrame: parameter dataframe
	"""
	df = pd.read_csv(param_file, header=0, index_col=0).dropna(axis=1, how="all")
	rename_map = {
		"Direction (dH)": "dH (direction)",
		"Direction (dS)": "dS (direction)",
		"Direction (dG)": "dG (direction)",
	}
	df = df.rename(columns=rename_map)
	required = [
		"dH", "dS", "dG",
		"dH (change)", "dS (change)", "dG (change)",
		"dH (direction)", "dS (direction)", "dG (direction)",
	]
	missing = [c for c in required if c not in df.columns]
	if missing:
		raise ValueError(f"Missing columns in parameter CSV: {missing}")
	for c in ["dH", "dS", "dG"]:
		df.loc[:,c] = pd.to_numeric(df.loc[:,c], errors="coerce")
	for c in ["dH (change)", "dS (change)", "dG (change)"]:
		df.loc[:,c] = df.loc[:,c].astype(str).str.lower().map({"true": True, "false": False})
	for c in ["dH (direction)", "dS (direction)", "dG (direction)"]:
		df.loc[:,c] = pd.to_numeric(df.loc[:,c], errors="coerce").fillna(0).astype(int)
	return df


def get_basepairs(params):
	"""
	Function to get base pairs

	Args:
		params (list): parameter list

	Returns:
		dict: base pair list
	"""
	basepairs = {}
	for param in params:
		if param.startswith(("init", "length", "symmetry", "re:", "reg:")):
			continue
		if "/" not in param:
			continue
		seq1, seq2 = param.split("/", 1)
		if len(seq1) != len(seq2):
			continue
		for b1, b2 in zip(seq1, seq2):
			if b1 in basepairs and basepairs[b1] != b2:
				raise ValueError(f"Conflicting base pair mapping for {b1}: {basepairs[b1]} vs {b2}")
			if b2 in basepairs and basepairs[b2] != b1:
				raise ValueError(f"Conflicting base pair mapping for {b2}: {basepairs[b2]} vs {b1}")
			basepairs[b1] = b2
			basepairs[b2] = b1
	return basepairs


def get_complement(sequence, basepairs):
	"""
	Function to get complement sequence using basepairs

	Args:
		sequence (str): sequence
		basepairs (dict): base pair list

	Raises:

	Returns:
		list: complement sequence list
	"""
	seq = list(sequence)
	try:
		return [basepairs[b] for b in seq]
	except KeyError as e:
		raise ValueError(f"Unknown base in sequence {sequence}: {e.args[0]}")


def count_param(sequence, complement, param):
	"""
	Function to count matching parameter

	Args:
		sequence (str): sequence
		complement (str): complement sequence
		param (str): a Nearest neighbor parameter (pattern)

	Returns:
		int: number of matched pattern
	"""
	if param == "symmetry":
		return int(sequence == "".join(reversed(complement)))

	if param.startswith("re:"):
		patterns = param.replace("re:", "").split("/", 1)
		if len(patterns) != 2:
			return 0
		return int(re.search(patterns[0], sequence) is not None and re.search(patterns[1], complement) is not None)

	if param.startswith("reg:"):
		patterns = param.replace("reg:", "").split("/", 1)
		if len(patterns) != 2:
			return 0
		pos1 = [obj_match.span()[0] for obj_match in re.finditer(patterns[0], sequence)]
		pos2 = [obj_match.span()[0] for obj_match in re.finditer(patterns[1], complement)]
		return len(pos1) if pos1 == pos2 else 0

	if "/" not in param:
		return 0

	query1, query2 = param.split("/", 1)
	L = len(query1)
	pairs = [sequence[i:i+L] + "/" + complement[i:i+L] for i in range(len(sequence) - L + 1)]
	cnt = pairs.count(param)
	rev_param = query2[::-1] + "/" + query1[::-1]
	if rev_param != param:
		cnt += pairs.count(rev_param)
	return cnt


def get_freq_matrix(sequences_A, sequences_B, params):
	"""
	Function to count matched patterns for sequences

	Args:
		sequences_A (list): sequence A list
		sequences_B (list): sequence B list
		params (list): parameter (pattern) list

	Returns:
		pd.DataFrame: row:sequence, col:pattern_name, dtype=float
	"""
	rows = []
	for seq, comp in zip(sequences_A, sequences_B):
		row = {p: count_param(seq, comp, p) for p in params}
		rows.append(row)
	return pd.DataFrame(rows, index=sequences_A, columns=params).fillna(0).astype(float)


def build_bounds(df_params, target):
	"""
	Function to build bounds

	Args:
		df_params (pd.DataFrame): parameter dataframe
		target (str): dH, dS, or dG

	Returns:
		free_idx: index for change flag
		x0(np.array): base value array
		lb(np.array): lower bound array
		ub(np.array): upper bound array
	"""
	vals = df_params.loc[:,target].to_numpy(float)	# parameter values for target
	change = df_params.loc[:,f"{target} (change)"].to_numpy(bool)	# change flag
	direction = df_params.loc[:,f"{target} (direction)"].to_numpy(int)	# direction option
	free_idx = np.where(change)[0]	# index for change flag with True
	x0 = vals[free_idx].copy()		# value with change flag
	lb = np.full_like(x0, -np.inf, dtype=float)	# change range for lower bound
	ub = np.full_like(x0, np.inf, dtype=float)	# change range for upper bound
	for j, idx in enumerate(free_idx):
		if direction[idx] == -1:
			ub[j] = vals[idx]
		elif direction[idx] == 1:
			lb[j] = vals[idx]
	return free_idx, x0, lb, ub


def unpack_theta(x_free, df_params, target, free_idx):
	theta = df_params[target].to_numpy(float).copy()
	theta[free_idx] = x_free
	return theta


def fit_one_target(df_exp, df_params, target, rng=None):
	"""
	Function to fitting for one target (dH, dS, or dG)

	Args:
		df_exp (pd.DataFrame): experiment dataframe
		df_params (pd.DataFrame): parameter dataframe
		target (str): dH, dS, or dG
		rng (generator, optional): random generator (Default: None)

	Returns:
		dict: fitting result for one target
	"""
	if rng is None:
		rng = np.random.default_rng()

	params = df_params.index.tolist()
	X = get_freq_matrix(df_exp["Sequence A"].tolist(), df_exp["Sequence B"].tolist(), params).to_numpy(float)
	y = df_exp[target].to_numpy(float)
	sigma = df_exp[f"{target} (error)"].to_numpy(float)	# sigma = error
	sigma = np.where(np.isfinite(sigma) & (sigma > 0), sigma, EPS)	# sigma 内の NaN, -inf, +inf を EPS で置換

	free_idx, x0, lb, ub = build_bounds(df_params, target)
	fixed_mask = np.ones(len(params), dtype=bool)
	fixed_mask[free_idx] = False
	theta_fixed = df_params[target].to_numpy(float)[fixed_mask]
	X_free = X[:, free_idx]
	X_fixed = X[:, fixed_mask]
	y_eff = y - X_fixed @ theta_fixed

	def residuals(x):
		pred = X_free @ x
		return (pred - y_eff) / sigma

	res = least_squares(residuals, x0=x0, bounds=(lb, ub), method="trf")
	theta = unpack_theta(res.x, df_params, target, free_idx)
	pred = X @ theta
	resid = y - pred
	chi2 = np.sum((resid / sigma) ** 2)
	dof = max(len(y) - len(free_idx), 1)
	reduced_chi2 = chi2 / dof
	return {
		"theta": theta,
		"pred": pred,
		"resid": resid,
		"chi2": chi2,
		"reduced_chi2": reduced_chi2,
		"success": bool(res.success),
		"message": res.message,
		"nfev": int(res.nfev),
		"X": X,
		"sigma": sigma,
	}


def sample_experimental_values(df_exp, targets, rng, mode="normal"):
	"""
	Function to generate one sample value set with error range

	Args:
		df_exp (pd.DataFrame): experiment dataframe
		targets (str): dH, dS, or dG
		rng (generator): random value generator
		mode (str, optional): generate value based on "normal" distribution or "uniform" (Default: "normal")

	Returns:
		pd.DataFrame: dataframe with the same format as experiment dataframe
	"""
	df_sample = df_exp.copy()
	for target in targets:
		y = df_sample.loc[:,target].to_numpy(float)
		e = df_sample.loc[:,f"{target} (error)"].to_numpy(float)
		e = np.where(np.isfinite(e) & (e > 0), e, EPS)
		if mode == "normal":
			y_new = rng.normal(loc=y, scale=e)
		elif mode == "uniform":
			y_new = rng.uniform(low=y - e, high=y + e)
		else:
			raise ValueError("mode must be 'normal' or 'uniform'")
		df_sample[target] = y_new
	return df_sample


def fit_by_monte_carlo(df_exp, df_params, n_mc=100, noise_model="normal", seed=SEED, calc_dS=True, temperature=None):
	"""
	Function to monte carlo fitting

	Args:
		df_exp (pd.Dataframe): experiment dataframe
		df_params (pd.Dataframe): parameter dataframe
		n_mc (int, optional): monte carlo iteration (Default: 100)
		noise_model (str, optional): generate value based on "normal" distribution or "uniform" (Default: "normal")
		seed (int, optional): random seed (Default: SEED)
		calc_dS (bool): calculate dS from dH, dG, and temperature (Default: False)
		temperature (float): temperature

	Returns:
		pd.DataFrame: fitting parameter
		pd.DataFrame: summary log
		pd.DataFrame: iteration log
	"""
	rng = np.random.default_rng(seed)
	records = []
	fit_logs = []

	targets = ["dH", "dG"]
	if calc_dS:
		targets = ["dH", "dS", "dG"]

	for i in range(n_mc):
		df_MC = sample_experimental_values(df_exp, targets, rng, mode=noise_model)

		for target in targets:
			fit = fit_one_target(df_MC, df_params, target)

			fit_logs.append({
				"mc_index": i,
				"target": target,
				"chi2": fit["chi2"],
				"reduced_chi2": fit["reduced_chi2"],
				"success": fit["success"],
				"nfev": fit["nfev"],
				"message": fit["message"],
			})

			for p, val in zip(df_params.index, fit["theta"]):
				records.append({
					"mc_index": i,
					"Parameter": p,
					"target": target,
					"value": val,
				})

	df_samples = pd.DataFrame(records)
	if len(targets) == 2:
		df_tmp = df_samples.copy()
		df_dH = df_tmp.loc[df_tmp.loc[:,"target"] == "dH",["mc_index", "Parameter", "value"]].rename(columns={"value": "dH"})
		df_dG = df_tmp.loc[df_tmp.loc[:,"target"] == "dG",["mc_index", "Parameter", "value"]].rename(columns={"value": "dG"})

		df_dS = df_dH.merge(df_dG, on=["mc_index", "Parameter"], how="inner")
		df_dS.loc[:,"target"] = "dS"
		df_dS.loc[:,"value"] = (df_dS.loc[:,"dH"] - df_dS.loc[:,"dG"]) / temperature * 1000.0

		df_samples = pd.concat([df_tmp, df_dS.loc[:,["mc_index", "Parameter", "target", "value"]]], ignore_index=True)

	summary_rows = []
	for (param, target), g in df_samples.groupby(["Parameter", "target"]):
		vals = g["value"].to_numpy(float)
		summary_rows.append({
			"Parameter": param,
			"target": target,
			"value": float(np.mean(vals)),
			"error": float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0,
			"q2.5": float(np.quantile(vals, 0.025)),
			"q50": float(np.quantile(vals, 0.5)),
			"q97.5": float(np.quantile(vals, 0.975)),
		})

		df_samples_wide = (df_samples.pivot(index=["mc_index", "target"], columns="Parameter", values="value").reset_index())
		df_fitting_log = pd.DataFrame(fit_logs).merge(
			df_samples_wide,
			on=["mc_index", "target"],
			how="left"
		)

	df_fitting_params = pd.DataFrame(summary_rows).set_index("Parameter", drop=True)
	return df_fitting_params, df_fitting_log


def fit_nominal(df_exp, df_params, calc_dS=True, temperature=None):
	"""
	Function to fitting values without error

	Args:
		df_exp (pd.DataFrame): experiment dataframe
		df_params (pd.DataFrame): parameter dataframe
		calc_dS (bool): calculate dS from dH, dG, and temperature (Default: False)
		temperature (float): temperature (Default: 310.15)

	Returns:
		pd.DataFrame: fitting parameter
		pd.DataFrame: log
	"""
	rows = []
	fit_logs = []

	targets = ["dH", "dG"]
	if calc_dS:
		targets = ["dH", "dS", "dG"]

	for target in targets:
		fit = fit_one_target(df_exp, df_params, target)

		for p, val in zip(df_params.index, fit["theta"]):
			rows.append({
				"Parameter": p,
				"target": target,
				"value": val,
				"error": np.nan,
				"q2.5": np.nan,
				"q50": np.nan,
				"q97.5": np.nan,
			})

		fit_logs.append({
			"target": target,
			"chi2": fit["chi2"],
			"reduced_chi2": fit["reduced_chi2"],
			"success": fit["success"],
			"nfev": fit["nfev"],
			"message": fit["message"],
		})

	df_fitting_params = pd.DataFrame(rows).set_index("Parameter", drop=True)
	if len(targets) == 2:
		df_tmp = df_fitting_params.copy()
		df_dH = df_tmp.loc[df_tmp.loc[:,"target"] == "dH","value"].rename("dH")
		df_dG = df_tmp.loc[df_tmp.loc[:,"target"] == "dG","value"].rename("dG")

		df_dS = pd.concat([df_dH, df_dG], axis=1)
		df_dS.loc[:,"target"] = "dS"
		df_dS.loc[:,"value"] = (df_dS.loc[:,"dH"] - df_dS.loc[:,"dG"]) / temperature * 1000.0
		df_dS.loc[:,"error"] = np.nan
		df_dS.loc[:,"q2.5"] = np.nan
		df_dS.loc[:,"q50"] = np.nan
		df_dS.loc[:,"q97.5"] = np.nan

		df_dS = df_dS.loc[:,["target", "value", "error", "q2.5", "q50", "q97.5"]]
		df_fitting_params = pd.concat([df_tmp, df_dS])

	return df_fitting_params, pd.DataFrame(fit_logs)


def predict_from_parameters(df_exp, df_params, targets):
	"""
	Function to fitting values without error

	Args:
		df_exp (pd.DataFrame): experiment dataframe
		df_params (pd.DataFrame): parameter dataframe
		targets (str): dH, dS, or dG

	Returns:
		pd.DataFrame: prediction result
	"""
	rows = []

	for target in targets:
		sub = df_params.loc[df_params.loc[:,"target"] == target,:].copy()
		theta_map = dict(zip(sub.index, sub["value"]))

		params = list(theta_map.keys())
		X = get_freq_matrix(
			df_exp["Sequence A"].tolist(),
			df_exp["Sequence B"].tolist(),
			params,
		)

		theta = np.array([theta_map[p] for p in params], dtype=float)
		pred = X.to_numpy(float) @ theta
		exp = df_exp[target].to_numpy(float)
		resid = exp - pred

		rows.append(pd.DataFrame({
			"Comment": df_exp["Comment"],
			"Sequence A": df_exp["Sequence A"],
			"Sequence B": df_exp["Sequence B"],
			"target": target,
			"exp": exp,
			"pred": pred,
			"residual": resid,
		}))

	return pd.concat(rows, ignore_index=True)


def confirm_overwrite(target_file, force_overwrite):
	"""
	Function to confirm for overwrite

	Args:
		target_file (str): target file path
		force_overwrite (bool): overwrite forcibly
	"""

	if not os.path.exists(target_file):
		return None

	if force_overwrite:
		os.remove(target_file)
		return None

	sys.stderr.write(f"WARN: {target_file} exists. Overwrite it? (y/N): ")
	sys.stderr.flush()
	user = sys.stdin.readline().lower().strip()

	if user != "y":
		sys.exit(0)

	else:
		os.remove(target_file)



# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Fitting program for parameters of Nearest Neighbor method using Least square solver and Monte Carlo method")
	parser.add_argument("-x", dest="EXP_FILE", metavar="EXP.csv", required=True, 	help = """sequence and experimental value file
column: Label, Sequence, dH, dH(error), dS, dS(error), dG, and dG(error)
""")
	parser.add_argument("-r", dest="PARAM_FILE", metavar="REF_PARAM.csv", required=True,
	help = """referenced parameter values
column: `Parameter`, `dH`, `dS`, `dG`, space, `dH (change)`, `dS (change)`, `dG (change)`, space, `dH (Direction)`, `dS (Direction)`, and `dG (Direction)`
`Parameter`: AA/TT, GC/CG, or etc, and re: (regexp) and reg: (counted pattern by regexp)
e.g., "re:^A/^T" and "re:^T/^A" (initial parameter for A/T (both specification require))
e.g., "reg:.*?G.*?/.*?C.*?" (number of G/C pair parameter)
e.g., "reg:./." (length parameter)
""")
	parser.add_argument("-o", dest="OUTPUT_FILE", metavar="OUTPUT.csv", required=True, help="output file")
	parser.add_argument("--mc", dest="MC_ITERATION", type=int, default=100, help="Number of Monte Carlo replicates")
	parser.add_argument("-S", dest="FLAG_SEPARATE", action="store_true", default=False, help="Separately calculate dS (Default: OFF (dS is calculated by Gibbs free energy equation))")
	parser.add_argument("-T", dest="TEMPERATURE", metavar="TEMPERATURE", type=float, default=None, required="-S" not in sys.argv, help="temperature for experimental data (Default: 310.15)")
	parser.add_argument("-e", "--error",  dest="FLAG_ERROR", action="store_true", default=False, help="consider with experimental value with error using least square and MC; return the fitting parameters of least square with not -e option")
	parser.add_argument("--noise", dest="NOISE", choices=["normal", "uniform"], default="normal", help="Noise model for sampling experimental values")
	parser.add_argument("--seed", dest="SEED", type=int, default=SEED, help="Random seed")
	parser.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite forcibly")
	parser.add_argument("--debug_dir", dest="DEBUG_DIR", help="debug directory")
	args = parser.parse_args()

	# read input files
	df_exp = read_exp(args.EXP_FILE)
	df_params = read_params(args.PARAM_FILE)
	query_sequenceB_empty = df_exp.loc[:,"Sequence B"].isna()
	if query_sequenceB_empty.sum() != 0:
		basepairs = get_basepairs(df_params["Parameter"].tolist())
		df_exp.loc[query_sequenceB_empty,"Sequence B"] = df_exp.loc[query_sequenceB_empty,"Sequence A"].apply(lambda x: "".join(get_complement(x, basepairs)))

	if "dS" in TARGETS:
		med = np.nanmedian(np.abs(df_exp["dS"].to_numpy(float)))
		if med > 1.0:
			df_exp = df_exp.copy()
			df_exp["dS"] = df_exp["dS"] / 1000.0
			df_exp["dS (error)"] = df_exp["dS (error)"] / 1000.0

	# fitting
	df_fitting_params = None
	df_fitting_summary = None
	if args.FLAG_ERROR:
		df_fitting_params, df_fitting_summary = fit_by_monte_carlo(df_exp, df_params, n_mc=args.MC_ITERATION, noise_model=args.NOISE, seed=args.SEED, calc_dS=args.FLAG_SEPARATE, temperature=args.TEMPERATURE)

	else:
		df_fitting_params, df_fitting_summary = fit_nominal(df_exp, df_params, calc_dS=args.FLAG_SEPARATE, temperature=args.TEMPERATURE)

	df_pred = predict_from_parameters(df_exp, df_fitting_params, TARGETS)


	# output (debug)
	if args.DEBUG_DIR is not None:
		if not os.path.isdir(args.DEBUG_DIR):
			os.makedirs(args.DEBUG_DIR)

		output_file = os.path.join(args.DEBUG_DIR, "fitting_paramters.csv")
		confirm_overwrite(output_file, args.FLAG_OVERWRITE)
		df_fitting_params.to_csv(output_file, index=True)

		output_file = os.path.join(args.DEBUG_DIR, "fitting_summary.csv")
		confirm_overwrite(output_file, args.FLAG_OVERWRITE)
		df_fitting_summary.to_csv(output_file, index=False)

		output_file = os.path.join(args.DEBUG_DIR, "args.txt")
		confirm_overwrite(output_file, args.FLAG_OVERWRITE)
		with open(output_file, "w", encoding="utf-8") as obj_output:
			for k, v in vars(args).items():
				obj_output.write(f"{k}: {v}\n")


	# output
	df_output_input_section = pd.DataFrame(["<< Input >>"])
	df_output_input_main = pd.DataFrame([
		["Program", f"{PROGRAM_NAME} (Ver. {VERSION})"],
		["Experimental data", args.EXP_FILE],
		["Reference parameter", args.PARAM_FILE],
		["MC iteration", args.MC_ITERATION],
		["Separate calculation", args.FLAG_SEPARATE],
		["Temperature", args.TEMPERATURE],
		["with_error", args.FLAG_ERROR],
		["Noise", args.NOISE],
		["Seed", args.SEED],
	], columns=range(2))
	df_output_input_separator = pd.DataFrame([""])

	df_output_params_main = df_params.copy()
	df_output_params_main.insert(0, "Parameter", df_params.index)
	df_output_params_header = pd.DataFrame([df_output_params_main.columns])
	df_output_params_main.index = range(df_output_params_main.shape[0])
	df_output_params_main.columns = range(df_output_params_main.shape[1])
	df_output_params_separator = pd.DataFrame([[""],[""]])

	df_output_fitting_params_section = pd.DataFrame(["<< Results >>"])
	df_output_fitting_params_main = pd.DataFrame(index=df_params.index)
	df_output_fitting_params_main.loc[:,"Parameter"] = df_params.index
	for target in TARGETS:
		query = df_fitting_params.loc[:,"target"] == target
		for col_name_source, col_name_output in zip(["value", "error"], [target, f"{target} (error)"]):
			df_extract = df_fitting_params.loc[query,col_name_source].copy().fillna(0.0)
			df_extract.name = col_name_output
			df_output_fitting_params_main = pd.concat([df_output_fitting_params_main, df_extract], axis=1)
	df_output_fitting_params_header = pd.DataFrame([df_output_fitting_params_main.columns])
	df_output_fitting_params_main.index = range(df_output_fitting_params_main.shape[0])
	df_output_fitting_params_main.columns = range(df_output_fitting_params_main.shape[1])
	df_output_fitting_params_separator = pd.DataFrame([""])

	df_output_sequence_section = pd.DataFrame(["<< Sequence >>"])
	df_output_sequence_main = df_exp.loc[:,["Comment","Sequence A","Sequence B"]].copy()
	df_output_sequence_main.columns = ["Comment","Sequence A","Sequence B"]
	df_output_sequence_main.index = df_output_sequence_main.loc[:,"Sequence A"]
	df_output_sequence_exp = pd.DataFrame(index=df_output_sequence_main.loc[:,"Sequence A"])
	df_output_sequence_pred = pd.DataFrame(index=df_output_sequence_main.loc[:,"Sequence A"])
	df_output_sequence_diff = pd.DataFrame(index=df_output_sequence_main.loc[:,"Sequence A"])
	for target in TARGETS:
		query = df_pred.loc[:,"target"] == target

		df_extract_exp = df_pred.loc[query,"exp"].copy()
		df_extract_exp.index = df_pred.loc[query,"Sequence A"]
		df_extract_exp.name = f"Exp. ({target})"
		df_output_sequence_exp = pd.concat([df_output_sequence_exp, df_extract_exp], axis=1)

		df_extract_pred = df_pred.loc[query,"pred"].copy()
		df_extract_pred.index = df_pred.loc[query,"Sequence A"]
		df_extract_pred.name = f"Pred. ({target})"
		df_output_sequence_pred = pd.concat([df_output_sequence_pred, df_extract_pred], axis=1)

		df_extract_diff = df_pred.loc[query,"residual"].copy()
		df_extract_diff.index = df_pred.loc[query,"Sequence A"]
		df_extract_diff.name = f"Diff. ({target})"
		df_output_sequence_diff = pd.concat([df_output_sequence_diff, df_extract_diff], axis=1)

	df_output_sequence_main = pd.concat([df_output_sequence_main, df_output_sequence_exp, df_output_sequence_pred, df_output_sequence_diff], axis=1)
	df_freq_matrix = get_freq_matrix(df_exp.loc[:,"Sequence A"].tolist(), df_exp.loc[:,"Sequence B"].tolist(), df_params.index.tolist())
	columns = df_freq_matrix.columns
	df_output_sequence_main = pd.concat([df_output_sequence_main, df_freq_matrix], axis=1)
	df_output_sequence_main = df_output_sequence_main.astype({v: int for v in columns})
	df_output_sequence_header = pd.DataFrame([df_output_sequence_main.columns])
	df_output_sequence_main.columns = range(df_output_sequence_main.shape[1])
	df_output_sequence_main.index = range(df_output_sequence_main.shape[0])

	df_output = pd.concat([
		df_output_input_section,
		df_output_input_main,
		df_output_input_separator,
		df_output_params_header,
		df_output_params_main,
		df_output_params_separator,
		df_output_fitting_params_section,
		df_output_fitting_params_header,
		df_output_fitting_params_main,
		df_output_fitting_params_separator,
		df_output_sequence_section,
		df_output_sequence_header,
		df_output_sequence_main,
	], ignore_index=True)
	confirm_overwrite(args.OUTPUT_FILE, args.FLAG_OVERWRITE)
	df_output.to_csv(args.OUTPUT_FILE, header=False, index=False)

