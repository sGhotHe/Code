####################################################################################
# INTRODUCTION:
# This code is to do plot several Sensitivity Analysis feature
# Created by Hebs at 23/6/7/9:09
# Contact: hebishuo@pku.edu.cn
####################################################################################

default_paras = [] # default parameters list
default_paras.append('n')
default_paras.append('nI')
default_paras.append('nI2')
default_paras.append('nBC')
default_paras.append('kBC')
default_paras.append('PNSD')
default_paras.append('MS')
default_paras.append('VD')
default_paras.append('CT')
default_paras.append('kappa')
default_paras.append('kappaI')
default_paras.append('kappaI2')
default_paras.append('rhoBC')
default_paras.append('BCPNSD')
default_paras.append('BCPMSD')
default_paras.append('BCI')
default_paras.append('BCAE')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import time

import calNRSA
import calANOVA
import calCorr
import calCons

def set_label_color(paras, NRSA, labels):
	'''
	This function is to set label different color by judging whether its value is small.
	Using NRSA result to judge. Set big value black and small gray.
	input:
		paras	: parameters name lise, array, string
		NRSA	: result of Nominal Range Sensitivity Analysis, array, float
		labels	: labels to set colors, array, string
	output:
		colors	: paras label list, array, string
	'''
	paras = np.array(paras)[np.where(np.array(NRSA)>0.1*np.max(NRSA))]
	colors = []
	for i in range(len(labels)):
		if labels[i] in paras:
			colors.append('k')
		else:
			colors.append('lightgrey')
	
	return colors

def sort_arrays(array1, array2):
	'''
	This function is to sort array1 and array2 by array1 size
	designed by ChatGPT
	input:
		array1		: array 1, array, float
		array2		: array 2, array
	output:
		sorted_1	: sorted array 1, sort from largets to smallest, array, float
		sorted_2	: sorted array 2, array
	'''
	pairs = zip(array1, array2)
	sorted_pairs = sorted(pairs, key=lambda x: x[0], reverse=True)
	sorted_1, sorted_2 = zip(*sorted_pairs)
	
	return np.array(sorted_1), sorted_2

def linear_color(color1, color2, n):
	'''
	This function is to use two colors to get n colors, using linear interpolation
	designed by ChatGPT
	input:
		color1	: color 1, string
		color2	: color 2, string
		n	: color numbers, must >= 2
	output:
		colors	: color list, array
	'''
	rgb1 = mcolors.to_rgb(color1)
	rgb2 = mcolors.to_rgb(color2)
	
	if n==1:
		return rgb1
	
	r_step = (rgb2[0]-rgb1[0]) / (n-1)
	g_step = (rgb2[1]-rgb1[1]) / (n-1)
	b_step = (rgb2[2]-rgb1[2]) / (n-1)
	
	colors = []
	
	for i in range(n):
		r = rgb1[0] + i * r_step
		g = rgb1[1] + i * g_step
		b = rgb1[2] + i * b_step
		colors.append((r, g, b))
	
	return colors

def cal_lim(array, **args):
	'''
	This function is to calculate array plot lim
	input:
		array   : origin data
		**upper : upper blank rate left for plot, default 0.1
		**down  : down blank rate left for plot, default 0.1
		**mx    : lim max value, defualt np.max(array) * upper
		**mn    : lim min value, defualt np.min(array) * down
	output:
		lim     : plot plt.ylim()
	'''
	if 'upper' in args:
		upper = args['upper']
	else:
		upper = 0.1
	if 'down' in args:
		down = args['down']
	else:
		down = 0.1
	if 'mx' in args:
		mx = args['mx']
		mx_mark = True
	else:
		mx_mark = False
	if 'mn' in args:
		mn = args['mn']
		mn_mark = True
	else:
		mn_mark = False
	
	ma = np.max(array)
	mi = np.min(array)
	up = (ma-mi) * upper
	dn = (ma-mi) * down
	lim_max = ma + up
	lim_min = mi - dn
	
	if mx_mark:
		lim_max = mx
	if mn_mark:
		lim_min = mn
	
	lim = [lim_min, lim_max]
	return lim

def plot(dtime, run_num, rate, par_range, **args):
	'''
	This function is to plot bar feature
	input:
		dtime		: runMie.py output time, string, example: '230101'
		run_num	: running times number, int
		rate	: linear change rate, float
		par_range	: parameters change range rate, float
		**path		: runMie.py output save path, string, default 'output/Mie/'
		**paras	: parameters list, array, string, default default_paras
		**K		: ANOVA group number, int, default 10
		**save		: save flag, bool, default False
		**save_path	: save path, string, default 'figure/SA/
		**debug	: output flag, bool, default False
	output:
		NRSA, ANOVA, Correlation Coefficient and Constraint Factor of AOD, SSA and g
		rank bar feature
	'''
	if 'path' in args:
		path = args['path']
	else:
		path = 'output/Mie/'
	if 'paras' in args:
		paras = args['paras']
	else:
		paras = default_paras
	if 'K' in args:
		K = args['K']
	else:
		K = 10
	if 'save' in args:
		save = args['save']
	else:
		save = False
	if 'save_path' in args:
		save_path = args['save_path']
	else:
		save_path = 'figure/SA/'
	if 'debug' in args:
		debug = args['debug']
	else:
		debug = False
	
	# calculate 4 SA results
	
	if debug:
		print('calculating...')
	
	NRSA = calNRSA.cal(run_num, rate, paras=paras, debug=debug)
	NRSA = abs(np.mean(NRSA, axis=2))
	if debug:
		print('NRSA done')
	
	ANOVA = abs(calANOVA.cal(dtime, paras=paras, K=K, path=path))
	if debug:
		print('ANOVA done')
	Corr_infos = calCorr.cal(dtime, paras=paras, path=path, debug=debug)
	Corr = abs(Corr_infos['PRCC'])
	if debug:
		print('Correlation Coefficient done')
	Cons_infos = calCons.cal(dtime, par_range, paras=paras, path=path)
	Cons = Cons_infos['Cons']  # in DY*(DX/X)^-1
	run_num = Cons_infos['run_num']
	if debug:
		print('Constraint Factors done')
	
	# plot
	
	if debug:
		print('done\nplotting...')
	
	outputs = ['AOD', 'SSA', 'g']
	out_num = len(outputs)
	
	fs = 13 # fontsize
	fw = 'bold' # fontweight
	width = 0.5 # bar width
	
	plt.figure(figsize=(20,6))
	plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05, wspace=0.08, hspace=0.3)
	colors = linear_color('k', 'lightgray', len(Cons[0]))
	
	for i in range(out_num):
		# do parameter labels change
		for j in range(len(paras)):
			if paras[j]=='rhoBC':
				paras[j] = '\u03C1$_{LAC}$'
			elif paras[j]=='n':
				paras[j] = 'CRI$_{shell}$'
			elif paras[j]=='nBC':
				paras[j] = 'n$_{core}$'
			elif paras[j]=='kBC':
				paras[j] = 'k$_{core}$'
			elif paras[j]=='BCI':
				paras[j] = 'LACI'
			elif paras[j]=='BCAE':
				paras[j] = 'AE$_{LAC}$'
			elif paras[j]=='kappa':
				paras[j] = '\u03BA'
			elif paras[j]=='kappaI':
				paras[j] = '\u03BAI'
			elif paras[j]=='kappaI2':
				paras[j] = '\u03BAI2'
			elif paras[j]=='BCPNSD':
				paras[j] = 'PNSD$_{BC}$'
			elif paras[j]=='BCPMSD':
				paras[j] = 'PMSD$_{BC}$'
		
		sorted_NRSA,	sorted_NRSA_paras	= sort_arrays(NRSA[i,:], paras)
		sorted_ANOVA,	sorted_ANOVA_paras	= sort_arrays(ANOVA[i,:], paras)
		sorted_Corr,	sorted_Corr_paras	= sort_arrays(Corr[i,:], paras)
		sorted_Cons,	sorted_Cons_paras	= sort_arrays(Cons[i,:], paras)
		
		NRSA_label_colors = set_label_color(sorted_NRSA_paras, sorted_NRSA, sorted_NRSA_paras)
		ANOVA_label_colors = set_label_color(sorted_NRSA_paras, sorted_NRSA, sorted_ANOVA_paras)
		Corr_label_colors = set_label_color(sorted_NRSA_paras, sorted_NRSA, sorted_Corr_paras)
		Cons_label_colors = set_label_color(sorted_NRSA_paras, sorted_NRSA, sorted_Cons_paras)
		
		# NRSA
		
		plt.subplot(4,3,1+i)
		plt.title(outputs[i], fontsize=fs+2, fontweight=fw)
		plt.bar(sorted_NRSA_paras, sorted_NRSA, color=colors, width=width)
		plt.ylim(cal_lim(NRSA[i,:], upper=0.35, mn=0))
		plt.yticks([])
		if i==0:
			plt.ylabel('dy/dx', fontsize=fs, fontweight=fw)
		plt.xticks(fontsize=fs, fontweight=fw)
		
		ax = plt.gca()
		[t.set_color(j) for (j,t) in zip(NRSA_label_colors, ax.xaxis.get_ticklabels())]
		if i==0:
			plt.text(0.01, 0.8, '(a) OAT', fontsize=fs, fontweight=fw, transform=ax.transAxes)
		
		# ANOVA
		
		plt.subplot(4,3,4+i)
		plt.bar(sorted_ANOVA_paras, sorted_ANOVA, color=colors, width=width)
		plt.ylim(cal_lim(ANOVA[i,:], upper=0.35, mn=0))
		plt.yticks([])
		if i==0:
			plt.ylabel('F-value', fontsize=fs, fontweight=fw)
		plt.xticks(fontsize=fs, fontweight=fw)
		
		ax = plt.gca()
		[t.set_color(j) for (j,t) in zip(ANOVA_label_colors, ax.xaxis.get_ticklabels())]
		if i==0:
			plt.text(0.01, 0.8, '(b) ANOVA', fontsize=fs, fontweight=fw, transform=ax.transAxes)
		
		# Correlation Coefficient
		
		plt.subplot(4,3,7+i)
		plt.bar(sorted_Corr_paras, sorted_Corr, color=colors, width=width)
		plt.ylim(cal_lim(Corr[i,:], upper=0.35, mn=0))
		plt.yticks([])
		if i==0:
			plt.ylabel('CC', fontsize=fs, fontweight=fw)
		plt.xticks(fontsize=fs, fontweight=fw)
		
		ax = plt.gca()
		[t.set_color(j) for (j,t) in zip(Corr_label_colors, ax.xaxis.get_ticklabels())]
		if i==0:
			plt.text(0.01, 0.8, '(c) PCC', fontsize=fs, fontweight=fw, transform=ax.transAxes)
		
		# Constraint Factor
		
		plt.subplot(4,3,10+i)
		plt.bar(sorted_Cons_paras, sorted_Cons*100, color=colors, width=width)
		plt.ylim(cal_lim(Cons[i,:]*100, upper=0.35, mn=0))
		plt.yticks(fontsize=fs, fontweight=fw)
		if i==0:
			plt.ylabel('dy/(dx/x),10$^{-2}$', fontsize=fs, fontweight=fw)
		plt.xticks(fontsize=fs, fontweight=fw)
		
		ax = plt.gca()
		[t.set_color(j) for (j,t) in zip(Cons_label_colors, ax.xaxis.get_ticklabels())]
		if i==0:
			plt.text(0.01, 0.8, '(d) CP', fontsize=fs, fontweight=fw, transform=ax.transAxes)
			plt.text(0.7, 0.8, 'run num = '+str(run_num), fontsize=fs, fontweight=fw, transform=ax.transAxes)
		plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':', axis='y')
		
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_SA.pdf')
	else:
		plt.show()

if __name__ == '__main__':
	paras = [] # default parameters list
	paras.append('n')
	#paras.append('nI')
	#paras.append('nI2')
	paras.append('nBC')
	paras.append('kBC')
	paras.append('PNSD')
	paras.append('MS')
	#paras.append('VD')
	paras.append('CT')
	paras.append('kappa')
	#paras.append('kappaI')
	#paras.append('kappaI2')
	paras.append('rhoBC')
	#paras.append('BCPNSD')
	#paras.append('BCPMSD')
	#paras.append('BCI')
	paras.append('BCAE')
	
	plot('230604', 2, 0.3, 0.05, paras=paras, debug=True, save=False)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
