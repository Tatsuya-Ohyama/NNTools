#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Graph class using matplotlib
"""

import sys
import matplotlib.pyplot as plt



# =============== constant (user) =============== #
FIG_SIZE = [600, 500]
DPI = 80
FONT_FAMILY = "Noto Sans CJK JP"
FONT_SIZE = 22
FONT_COLOR = "#000000"
BACKGROUND_COLOR = "#FFFFFF"
BACKGROUND_ALPHA = 1.0
BORDER_COLOR = "#000000"
TIC_X = [0, 50, 10]
TIC_Y = [0, 10, 1]
TICK_DIRECTION = "inout"
TICK_LENGTH = 6.0
TICK_PAD = 5.0
TICK_ROTATION_X = 90
TICK_ROTATION_Y = 0




# =============== function =============== #
# init function
def matplot_init(font_family=FONT_FAMILY, font_size=FONT_SIZE, tick_direction=TICK_DIRECTION, tick_length=TICK_LENGTH, tick_pad=TICK_PAD):
	"""
	matplotlib 起動時に初期化する関数 (フォントキャッシュクリア＆デフォルト値設定)

	Args:
		font_family (str, optional): フォントの種類 (Default: `FONT_FAMILY`)
		font_size (float, optional): フォントのサイズ (Default: `FONT_SIZE`)
		tick_direction (str, optional): 目盛りの向き `in`, `out` or `inout` (Default: `TICK_DIRECTION`)
		tick_length (float, optional): 目盛りの長さ (Default: `TICK_LENGTH`)
		tick_pad (float, optional): 目盛りと目盛りラベル間の距離 (Default: `TICK_PAD`)
	"""
	import matplotlib
	import os
	import shutil
	if os.path.isdir(matplotlib.get_cachedir()):
		for name in os.listdir(matplotlib.get_cachedir()):
			if os.path.isfile(name):
				os.remove(name)

			elif os.path.isdir(name):
				shutil.rmtree(name)

	# グラフのデフォルト値の設定
	plt.rcParams['axes.linewidth'] = 2.0
	plt.rcParams['axes.axisbelow'] = True
	plt.rcParams['font.family'] = font_family
	plt.rcParams['font.size'] = font_size
	plt.rcParams['font.weight'] = "medium"
	plt.rcParams['xtick.direction'] = tick_direction
	plt.rcParams['xtick.major.size'] = tick_length
	plt.rcParams['xtick.major.pad'] = tick_pad
	plt.rcParams['ytick.direction'] = tick_direction
	plt.rcParams['ytick.major.size'] = tick_length
	plt.rcParams['ytick.major.pad'] = tick_pad
	plt.rcParams.update({"mathtext.default": "regular"})



# =============== class =============== #
# =============== Graph2D =============== #
class Graph2D:
	def __init__(self, width=FIG_SIZE[0], height=FIG_SIZE[1], background_color=BACKGROUND_COLOR, background_alpha=BACKGROUND_ALPHA):
		self._obj_fig = None		# fig オブジェクト
		self._obj_ax = None			# カレント ax オブジェクト
		self._list_obj_ax = []		# ax オブジェクトのリスト

		self.__make_obj_fig(width, height, background_color, background_alpha)


	@property
	def fig(self):
		return self._obj_fig

	@property
	def ax(self):
		return self._obj_ax

	@property
	def axs(self):
		return self._list_obj_ax


	def __make_obj_fig(self, width, height, background_color, background_alpha):
		"""
		matplotlib の fig オブジェクトを生成するメソッド

		Args:
			width (int, optional): fig の幅 (px) (Default: FIG_SIZE[0])
			height (int, optional): fig の高さ (px) (Default: FIG_SIZE[1])
			background_color (str): 背景色
			background_alpha (float): 背景色のアルファ値
		"""
		self._obj_fig = plt.figure(1, figsize=(width / DPI, height / DPI), dpi=DPI)
		self._obj_fig.patch.set_facecolor(background_color)		# 図全体の背景色
		self._obj_fig.patch.set_alpha(background_alpha)			# 図全体の背景色のアルファ


	def append_ax(self, row, col, ax_idx, background_color=BACKGROUND_COLOR, background_alpha=BACKGROUND_ALPHA, border_color=BORDER_COLOR):
		"""
		ax オブジェクトを追加する関数

		Args:
			row (int): fig オブジェクト内に配置する行数
			col (int): fig オブジェクト内に配置する列数
			ax_idx (int): fig オブジェクト内の位置インデックス
			background_color (str, optional): 背景色 (Default: BACKGROUND_COLOR)
			background_alpha (float, optional): 背景色のアルファ値 (Default: BACKGROUND_ALPHA)
			border_color (str), optional): 枠の色 (Default: BORDER_COLOR)

		Returns:
			self
		"""
		# ax オブジェクト作成
		ax = self._obj_fig.add_subplot(row, col, ax_idx)
		ax.patch.set_facecolor(background_color)		# 各グラフの背景色
		ax.patch.set_alpha(background_alpha)			# 各グラフの背景色のアルファ
		ax.spines["right"].set_color(border_color)		# 枠の有無および色
		ax.spines["left"].set_color(border_color)		# 枠の有無および色
		ax.spines["top"].set_color(border_color)		# 枠の有無および色
		ax.spines["bottom"].set_color(border_color)		# 枠の有無および色

		self._list_obj_ax.append(ax)
		self._obj_ax = self._list_obj_ax[-1]
		return self


	def change_ax(self, idx):
		"""
		カレントの ax オブジェクトを変更するメソッド

		Args:
			idx (int): axs リストのインデックス (0〜n)

		Returns:
			self
		"""
		try:
			self._obj_ax = self._list_obj_ax[idx]
		except IndexError:
			sys.stderr.write("ERROR: invalid index.\n")
		return self


	def set_title(self, title, font_size=FONT_SIZE):
		"""
		現在の ax オブジェクトのグラフにタイトルを設定するメソッド

		Args:
			title (str): グラフタイトル
			font_size (int, optional): フォントサイズ (Default: `FONT_SIZE`)

		Returns:
			self
		"""
		self._obj_ax.set_title(title, fontsize=font_size)
		return self


	def set_label(self, axis, label, font_size=FONT_SIZE, color=FONT_COLOR):
		"""
		現在の ax オブジェクトのグラフに x 軸のラベルを設定するメソッド

		Args:
			axis (str): `x` or `y`: 設定する軸
			label (str): ラベル
			font_size (float, optional): フォントサイズ (Default: `FONT_SIZE`)
			color (str, optional): ラベルの色 (Default: `FONT_COLOR`)

		Returns:
			self
		"""
		if font_size is None:
			font_size = FONT_SIZE

		if axis.lower() == "x":
			self._obj_ax.set_xlabel(label, fontsize=font_size, weight="medium", color=color)
		elif axis.lower() == "y":
			self._obj_ax.set_ylabel(label, fontsize=font_size, weight="medium", color=color)
		return self


	def set_range(self, axis, start, end):
		"""
		現在の ax オブジェクトのグラフの軸に対する描画範囲を設定するメソッド

		Args:
			axis (str): `x` or `y`: 設定する軸
			start (float or None): 描画開始値
			end (float or None): 描画終了値

		Returns:
			self
		"""
		if axis.lower() == "x":
			if start is None or end is None:
				self._obj_ax.relim()
				self._obj_ax.autoscale_view()
				self._obj_ax.set_xlim(auto=True)

			if start is not None:
				self._obj_ax.set_xlim(xmin=start)

			if end is not None:
				self._obj_ax.set_xlim(xmax=end)

		elif axis.lower() == "y":
			if start is None or end is None:
				self._obj_ax.relim()
				self._obj_ax.autoscale_view()
				self._obj_ax.set_ylim(auto=True)

			if start is not None:
				self._obj_ax.set_ylim(ymin=start)

			if end is not None:
				self._obj_ax.set_ylim(ymax=end)

		else:
			sys.stderr.write("ERROR: invalid axis.\n")
			sys.exit(1)
		return self


	def set_ticks(self, axis, list_ticks, list_labels=None, format=None, font_size=FONT_SIZE, color=FONT_COLOR, tick_direction=TICK_DIRECTION, tick_length=TICK_LENGTH, tick_pad=TICK_PAD, rotation=None):
		"""
		現在の ax オブジェクトのグラフの目盛りを設定をするメソッド

		Args:
			axis (str): `x` or `y`: 設定する軸
			list_ticks (list): [x1, ...]
			list_labels (list): [tic_label1, ...]
			format (str, optional): 軸ラベルのフォーマット
			font_size (float or str): フォントサイズ (Default: `FONT_SIZE`)
			color (str): 目盛りの色 (Default: `FONT_COLOR`)
			tick_direction(str): 目盛りの向き (Default: `TICK_DIRECTION`)
			tick_length (float): 目盛りの長さ (Default: `TICK_LENGTH`)
			tick_pad (float): 目盛りと目盛りラベルの間隔 (Default: `TICK_PAD`)
			rotation (float): 目盛りラベルの角度 (Default: `TICK_ROTATION_X` or `TICK_ROTATION_Y`)

		Returns:
			self
		"""
		if axis.lower() == "x":
			if rotation is None:
				rotation = TICK_ROTATION_X
			self._obj_ax.tick_params(axis="x", which="major", length=tick_length, direction=tick_direction, labelsize=font_size, pad=tick_pad, rotation=rotation)
			self._obj_ax.set_xticks(list_ticks)
			if list_labels is not None:
				self._obj_ax.set_xticklabels(list_labels, fontsize=font_size, color=color)
			if format is not None:
				self._obj_ax.xaxis.set_major_formatter(plt.FormatStrFormatter(format))

		elif axis.lower() == "y":
			if rotation is None:
				rotation = TICK_ROTATION_Y
			self._obj_ax.tick_params(axis="y", which="major", length=tick_length, direction=tick_direction, labelsize=font_size, pad=tick_pad, rotation=rotation)
			self._obj_ax.set_yticks(list_ticks)
			if list_labels is not None:
				self._obj_ax.set_yticklabels(list_labels, fontsize=font_size, color=color)
			if format is not None:
				self._obj_ax.yaxis.set_major_formatter(plt.FormatStrFormatter(format))

		else:
			sys.stderr.write("ERROR: invalid axis.\n")
			sys.exit(1)
		return self


	def add_zeroaxis_bar(self, axis, v=0, color=BORDER_COLOR):
		"""
		現在の ax オブジェクトのグラフの特定の軸の特定の x 値に線を加えるメソッド

		Args:
			axis (str): `x` or `y`: 設定する軸
			v (float): 線を追加する位置 (Default: 0)
			color (str): 線の色 (Default: `BORDER_COLOR`)

		Returns:
			self
		"""
		if axis.lower() == "x":
			self._obj_ax.axvline(x=v, color=color)
		elif axis.lower() == "y":
			self._obj_ax.axhline(y=v, color=color)
		else:
			sys.stderr.write("ERROR: invalid axis.\n")
			sys.exit(1)
		return self


	def add_grid(self, axis, color=BORDER_COLOR):
		"""
		現在の ax オブジェクトのグラフの x 軸にグリッド線を追加するメソッド

		Args:
			axis (str): `x` or `y`: 設定する軸

		Returns:
			self
		"""
		if axis.lower() == "x":
			self._obj_ax.grid(True, axis="x", color=color, linewidth=1, dashes=(5,2.5))
		elif axis.lower() == "y":
			self._obj_ax.grid(True, axis="y", color=color, linewidth=1, dashes=(5,2.5))
		return self


	def set_legend(self, show=True, col=1, loc="best", font_size=FONT_SIZE, color=FONT_COLOR, facecolor=BACKGROUND_COLOR):
		"""
		現在の ax オブジェクトのグラフの凡例を設定するメソッド

		Args:
			show (bool, optional): 凡例を表示するかのフラグ (Default: True)
			col (int, optional): 凡例の列数 (Default: 1)
			loc (str, optional): 位置 (Default: "best")
			font_size (float, optional): フォントサイズ (Default: `FONT_SIZE`)
			color (str, optional): フォントの色 (Default: `FONT_COLOR`)

		Returns:
			self
		"""
		if show:
			self._obj_ax.legend(ncol=col, fontsize=font_size, labelcolor=color, loc=loc, facecolor=facecolor)
		else:
			self._obj_ax.get_legend().remove()
		return self


	def make_graph(self, interface, output_file=None, no_tight_layout=False):
		"""
		グラフオブジェクトからグラフの表示、あるいは保存するメソッド

		Args:
			interface (str): "GUI" or "CUI"
			output_file (str, optional): 保存するグラフの .png ファイルのパス (Default: None)
		"""
		if no_tight_layout == False:
			self._obj_fig.tight_layout()

		if output_file is not None:
			plt.savefig(output_file)

		if interface == "GUI":
			plt.show()

		elif interface != "CUI":
			sys.stderr.write("ERROR: invalid interface.\n")
			sys.exit(1)

		plt.close()



class DataSet:
	"""
	データクラス
	"""
	def __init__(self, name, input_file=None, output_file=None, legend=None):
		self._name = None
		self._input_file = None
		self._output_file = None
		self._legend = None

		self._name = name
		self.set_input_file(input_file)
		self.set_output_file(output_file)
		self.set_legend(legend)


	@property
	def name(self):
		return self._name

	@property
	def input_file(self):
		return self._input_file

	@property
	def output_file(self):
		if self._output_file is None:
			self._output_file = "{0}.png".format(self._name)
		return self._output_file

	@property
	def legend(self):
		if self._legend is None:
			self._legend = self._name
		return self._legend


	def set_input_file(self, input_file):
		"""
		入力ファイルを設定するメソッド

		Args:
			input_file (str): input file

		Returns:
			self
		"""
		self._input_file = input_file
		return self


	def set_output_file(self, output_file):
		"""
		出力ファイルを設定するメソッド

		Args:
			output_file (str): output file

		Returns:
			self
		"""
		self._output_file = output_file
		return self


	def set_legend(self, legend):
		"""
		凡例名を設定するメソッド

		Args:
			legend (str): legend

		Returns:
			self
		"""
		self._legend = legend
		return self
