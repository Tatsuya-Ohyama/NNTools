#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%% ==============================
# modules
import sys, signal
sys.dont_write_bytecode = True
signal.signal(signal.SIGINT, signal.SIG_DFL)

import pandas as pd
import numpy as np
from scipy.optimize import lsq_linear



#%% ==============================
# constant
NN_TYPES = [
	"AA", "AT", "TA", "CA", "GT",
	"CT", "GA", "CG", "GC", "GG"
]



#%% ==============================
# function
def read_exp(exp_file):
	"""
	Function to read experimental file

	Args:
		exp_file (str): experimental file

	Returns:
		pd.DataFrame
	"""
	df = pd.read_csv(exp_file, header=0, index_col=False)
	df = df.astype({
		"Comment": str,
		"Sequence": str,
		"dH": float,
		"dH (error)": float,
		"dS": float,
		"dS (error)": float,
		"dG": float,
		"dG (error)": float,
	})
	return df


def read_params(param_file):
	"""
	Function to read parameter file

	Args:
		param_file (str): parameter file

	Returns:
		pd.DataFrame
	"""
	df = pd.read_csv(param_file, header=0, index_col=False)
	df = df.dropna(axis=1)
	df = df.astype({
		"Parameter": str,
		"dH": float,
		"dS": float,
		"dG": float,
		"dH (change)": bool,
		"dS (change)": bool,
		"dG (change)": bool,
		"Direction (dH)": int,
		"Direction (dS)": int,
		"Direction (dG)": int,
	})
	return df


def build_design_matrix(sequences, nn_types):
	"""
	NN カウント + INIT 列を含むデザイン行列を作成。
	列順: [NN_TYPES..., INIT]
	"""
	n_seq = len(sequences)
	n_param = len(nn_types)
	X = np.zeros((n_seq, n_param), dtype=float)

	for i, seq in enumerate(sequences):
		seq = seq.upper().replace(" ", "")
		for j in range(len(seq) - 1):
			dinuc = seq[j:j+2]
			if dinuc in nn_types:
				k = nn_types.index(dinuc)
				X[i, k] += 1.0

	return X


def fit_nn_wls_fixed_and_bounded(
	sequences,
	y_exp,
	y_err,
	fixed_params=None,   # 例: {"AA": -1.0, "INIT": 0.0}
	lower_bounds=None,   # 例: {"AA": -np.inf, "AT": -np.inf, "INIT": 0.0}
	upper_bounds=None,   # 例: {"AA": 0.0, "AT": 0.0}
	nn_types=NN_TYPES
):
	"""
	固定パラメータ + 境界付き（符号制約など）WLS。
	fixed_params, lower_bounds, upper_bounds を None にすると
	「全部自由・無制約」として動作する。
	"""
	y = np.asarray(y_exp, dtype=float)
	sigma = np.asarray(y_err, dtype=float)
	X_full = build_design_matrix(sequences, nn_types)

	# パラメータ名を X_full の列数から決める
	# 列順: NN_TYPES..., INIT
	n_par = X_full.shape[1]
	names = list(nn_types) + ["INIT"]  # ここで len(names) == n_par のはず

	# None → {} に変換
	fixed_params = {} if fixed_params is None else dict(fixed_params)
	lower_bounds = {} if lower_bounds is None else dict(lower_bounds)
	upper_bounds = {} if upper_bounds is None else dict(upper_bounds)

	# --- 固定パラメータのマスクと値 ---
	fixed_mask = np.zeros(n_par, dtype=bool)
	fixed_values = np.zeros(n_par, dtype=float)

	for name, val in fixed_params.items():
		if name not in names:
			raise ValueError(f"Unknown parameter name in fixed_params: {name}")
		idx = names.index(name)
		fixed_mask[idx] = True
		fixed_values[idx] = float(val)

	free_mask = ~fixed_mask

	# X_full を固定／自由に分解
	# fixed_params が空 {} の場合、fixed_mask はすべて False で、X_fixed は列数 0 の行列になる
	X_fixed = X_full[:, fixed_mask]
	X_free  = X_full[:, free_mask]

	# 固定分を右辺に移す
	if X_fixed.shape[1] > 0:
		y_eff = y - X_fixed @ fixed_values[fixed_mask]
	else:
		y_eff = y.copy()

	# --- WLS → 通常 LS へ ---
	w = 1.0 / (sigma ** 2)
	sqrt_w = np.sqrt(w)

	A = X_free * sqrt_w[:, None]  # (n_samples, n_free_params)
	b = y_eff * sqrt_w

	n_free = A.shape[1]

	# 自由パラメータ用の境界（全自由の場合も含む）
	lb_free = np.full(n_free, -np.inf, dtype=float)
	ub_free = np.full(n_free, +np.inf, dtype=float)

	# free_names を names と free_mask から構成
	free_names = [name for name, is_free in zip(names, free_mask) if is_free]

	for i, name in enumerate(free_names):
		if name in lower_bounds:
			lb_free[i] = lower_bounds[name]
		if name in upper_bounds:
			ub_free[i] = upper_bounds[name]

	# --- 境界付き線形最小二乗 ---
	if n_free > 0:
		res = lsq_linear(A, b, bounds=(lb_free, ub_free), method="trf")
		beta_free = res.x

		# 共分散近似（境界を無視した線形近似）
		JtJ = A.T @ A
		cov_free = np.linalg.pinv(JtJ)
		beta_free_err = np.sqrt(np.diag(cov_free))
	else:
		# 全部固定されている場合
		beta_free = np.zeros(0, dtype=float)
		cov_free = np.zeros((0, 0), dtype=float)
		beta_free_err = np.zeros(0, dtype=float)

	# --- フルパラメータへ復元 ---
	beta_full = np.zeros(n_par, dtype=float)
	beta_err_full = np.zeros(n_par, dtype=float)

	# 固定
	beta_full[fixed_mask] = fixed_values[fixed_mask]
	beta_err_full[fixed_mask] = 0.0

	# 自由
	beta_full[free_mask] = beta_free
	beta_err_full[free_mask] = beta_free_err

	# 予測・残差・χ²
	y_pred = X_full @ beta_full
	residuals = y - y_pred
	chi2 = np.sum(((y - y_pred) / sigma) ** 2)

	# 自由度: データ数 - rank(A)（free が 0 の場合は rank 0）
	if n_free > 0:
		rank_free = np.linalg.matrix_rank(A)
	else:
		rank_free = 0
	dof = len(y) - rank_free
	red_chi2 = chi2 / dof if dof > 0 else np.nan

	return beta_full, beta_err_full, cov_free, residuals, y_pred, chi2, red_chi2, names




#%% ==============================
# main
sequences = [
	"ATGC",
	"AATT",
	"CGCG",
	"GGAT",
	"TATA",
]
y_exp = [-5.2, -4.8, -6.0, -5.5, -3.9]     # ΔG°
y_err = [0.2, 0.1, 0.3, 0.2, 0.15]         # ±σ

# 例1: INIT を 0 に固定し、AA, AT は必ず負（安定化）とする
fixed = {"INIT": 0.0}
lower = {"AA": -np.inf, "AT": -np.inf}
upper = {"AA": 0.0,     "AT": 0.0}

fixed = None
lower = None
upper = None

beta, beta_err, cov_free, res, y_pred, chi2, red_chi2, names = fit_nn_wls_fixed_and_bounded(
		sequences,
		y_exp,
		y_err,
		fixed_params=fixed,
		lower_bounds=lower,
		upper_bounds=upper
	)

print("Fitted parameters (ΔG° ± 1σ):")
for name, val, err in zip(names, beta, beta_err):
	print(f"{name:>4s}: {val:7.3f} ± {err:6.3f} kcal/mol")
print(f"\nchi2 = {chi2:.2f}, reduced chi2 = {red_chi2:.2f}")


