####################################################################################
# INTRODUCTION:
# This code is to plot the difference of OAT and CF
# Created by Hebs at 23/9/18/15:38
# Contact: hebishuo@pku.edu.cn
####################################################################################

default_paras = [] # default parameters list
default_paras.append('n')
default_paras.append('nH1')
default_paras.append('nBC')
default_paras.append('kBC')
default_paras.append('PNSD')
default_paras.append('MS')
default_paras.append('MSH1')
default_paras.append('MSH2')
default_paras.append('CT')
default_paras.append('CTH1')
default_paras.append('kappa')
default_paras.append('kappaH1')
default_paras.append('rhoBC')
default_paras.append('BCPNSD')
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
	
	return sorted_1, sorted_2

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
		run_num		: running times number, int
		rate		: linear change rate, float
		par_range	: parameters change range rate, float
		**path		: runMie.py output save path, string, default 'output/Mie/'
		**paras		: parameters list, array, string, default default_paras
		**save		: save flag, bool, default False
		**save_path	: save path, string, default 'figure/SA/
		**debug	: output flag, bool, default False
	output:
		OAT and Constraint Factor of AOD, SSA and g
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
	
	# calculate 2 SA results and 1 vs result
	
	if debug:
		print('calculating...')
	
	NRSA = calNRSA.cal(run_num, rate, paras=paras, debug=debug)
	NRSA = abs(np.mean(NRSA, axis=2))
	if debug:
		print('NRSA done')
	
	Cons_infos = calCons.cal(dtime, par_range, paras=paras, path=path)
	Cons = Cons_infos['Cons']
	run_num = Cons_infos['run_num']
	if debug:
		print('Constraint Factors done')
	
	vs = abs((NRSA - Cons) / Cons)
	if debug:
		print('vs done')
	
	# plot
	
	if debug:
		print('done\nplotting...')
	
	outputs = ['AOD', 'SSA', 'g']
	letters = ['(a) ', '(b) ', '(c) ']
	out_num = len(outputs)
	par_num = len(paras)
	
	fs = 13 # fontsize
	fw = 'bold' # fontweight
	
	plt.figure(figsize=(8,6))
	plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.1, wspace=0.05, hspace=0.1)
	colors = linear_color('k', 'lightgrey', 2)
	fs = 13 # fontsize
	fw = 'bold' # fontweight
	
	# do parameter labels change
	for j in range(par_num):
		if paras[j]=='n':
			paras[j] = 'CRI$_{shell}$'
		if paras[j]=='rhoBC':
			paras[j] = 'rhoLAC'
		elif paras[j]=='nBC':
			paras[j] = 'nLAC'
		elif paras[j]=='kBC':
			paras[j] = 'kLAC'
		elif paras[j]=='BCI':
			paras[j] = 'LACI'
		elif paras[j]=='BCAE':
			paras[j] = 'AELAC'
		elif paras[j]=='kappa':
			paras[j] = chr(954)
	
	plt.title('OAT vs CF', fontsize=fs+7, fontweight=fw)
	for i in range(out_num):
		plt.subplot(311+i)
		for j in range(par_num):
			if j==0 and i==0:
				plt.bar(j+1-0.1, NRSA[i,j], width=0.2, color=colors[0], label='OAT')
				plt.bar(j+1+0.1, Cons[i,j], width=0.2, color=colors[1], label='CP')
				plt.legend(fontsize=fs+2)
			else:
				plt.bar(j+1-0.1, NRSA[i,j], width=0.2, color=colors[0])
				plt.bar(j+1+0.1, Cons[i,j], width=0.2, color=colors[1])
			plt.text(j+1-0.17, max(NRSA[i,j],Cons[i,j])+0.001, str(round(vs[i,j]*100)/100), fontsize=fs+3, fontweight=fw, color='gray')
		if i==(out_num-1):
			plt.xticks(np.arange(5)+1, paras, fontsize=fs+5, fontweight=fw)
			plt.ylabel('sensitivity', fontsize=fs+5, fontweight=fw)
			plt.yticks([0,0.04,0.08],[0,0.04,0.08],fontsize=fs, fontweight=fw)
		else:
			plt.xticks([])
			plt.yticks(fontsize=fs, fontweight=fw)
		plt.ylim(0,max(np.max(NRSA[i]),np.max(Cons[i]))*1.8)
		plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':', axis='y')
		#plt.ylim(cal_lim(NRSA[i,:], upper=0.35, mn=0))
		#plt.yticks([])
		#plt.ylabel('k', fontsize=fs, fontweight=fw)
		#plt.xticks(fontsize=fs, fontweight=fw, rotation=315)
		
		ax = plt.gca()
		plt.text(0.03, 0.8, letters[i]+outputs[i], fontsize=fs+5, fontweight=fw, transform=ax.transAxes)
		
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_OATvsCF.pdf')
	else:
		plt.show()
	plt.close()
	
if __name__ == '__main__':
	paras = [] # default parameters list
	paras.append('n')
	#paras.append('nI')
	#paras.append('nI2')
	#paras.append('nBC')
	#paras.append('kBC')
	paras.append('PNSD')
	paras.append('MS')
	#paras.append('VD')
	paras.append('CT')
	paras.append('kappa')
	#paras.append('kappaI')
	#paras.append('kappaI2')
	#paras.append('rhoBC')
	#paras.append('BCPNSD')
	#paras.append('BCPMSD')
	#paras.append('BCI')
	#paras.append('BCAE')
	
	plot('230604', 2, 0.3, 0.05, paras=paras, debug=True, save=False)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
