####################################################################################
# INTRODUCTION:
# This code is to plot correlation coefficient figure
# Created by Hebs at 22/10/5/9:18
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.stats import rankdata

import calCorr
import calMie

def cal_lim(array, **args):
	'''
	This function is to calculate array plot lim
	input:
		array   : origin data
		**upper : upper blank rate left for plot, default 0.1
		**down  : down blank rate left for plot, default 0.1
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
	mn = np.min(array)
	mx = np.max(array)
	up = (mx-mn) * upper
	dn = (mx-mn) * down
	lim = [mn-dn, mx+up]
	return lim

def normalize(data):
	'''
	This function is to normalize data
	input:
		data	: input data, array, float
	output:
		data	: output normalized data, array, float
	'''
	# find the maximum data 'max',
	# to do: data = data / 'max'
	m = data[0]
	for i in range(len(data)):
		if data[i]>m:
			m = data[i]
	
	data = data / m
	return data

def plot_residual(dtime, **args):
	'''
	This function is to compare difference with CC, RCC, PCC and PRCC parameter and output raletive relationship, plot in figure
	input:
		dtime       : data time, string
		**save      : save flag, bool, default False
		**save_path : figure save path, string, defualt 'figure/Corr/'
	output:
		figure
	'''
	if 'save' in args:
		save = args['save']
	else:
		save = False
	if 'save_path' in args:
		save_path = args['save_path']
	else:
		save_path = 'figure/Corr/'
	
	parameters = []
	parameters.append('n')
	parameters.append('nI')
	parameters.append('nBC')
	parameters.append('kBC')
	parameters.append('PNSD')
	parameters.append('MS')
	parameters.append('VD')
	parameters.append('CT')
	parameters.append('kappa')
	parameters.append('kappaI')
	parameters.append('rhoBC')
	parameters.append('BCPNSD')
	parameters.append('BCPMSD')
	parameters.append('BCI')
	parameters.append('BCAE')
	par_num = len(parameters)
	
	outputs = ['AOD', 'SSA', 'g']
	out_num = len(outputs)
	
	path = 'output/Mie/' + dtime + '/all/'
	data = np.load(path+'infos.npy', allow_pickle=True)
	run_num = len(data)
	
	paras = np.zeros((par_num,run_num)) # all parameters
	outs = np.zeros((out_num,run_num))
	
	for i in range(par_num):
		for j in range(run_num):
			paras[i,j] = data[j][parameters[i]+'_rate']
	
	for i in range(out_num):
		for j in range(run_num):
			outs[i,j] = data[j][outputs[i]]
	
	r_paras = np.zeros(paras.shape)
	r_outs = np.zeros(outs.shape)
	
	res_paras = np.zeros((out_num, par_num, run_num))
	res_outs = np.zeros((out_num, par_num, run_num))
	
	res_r_paras = np.zeros((out_num, par_num, run_num))
	res_r_outs = np.zeros((out_num, par_num, run_num))
	
	for i in range(par_num):
		r_paras[i] = rankdata(paras[i])
	
	for i in range(out_num):
		r_outs[i] = rankdata(outs[i])
		for j in range(par_num):
			res_paras[i,j], res_outs[i,j] = calCorr.cal_residual(paras, outs[i], j)
			res_r_paras[i,j], res_r_outs[i,j] = calCorr.cal_residual(r_paras, r_outs[i], j)
	
	fontsize = 20
	
	for i in range(par_num):
		plt.figure(figsize=(12,9))
		plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.1, hspace=0.1)
		
		# CC
		
		plt.subplot(341)
		plt.scatter(paras[i], outs[0], s=5, c='deepskyblue')
		plt.xticks([])
		plt.yticks([])
		plt.title('CC', fontsize=fontsize, fontweight='bold')
		plt.ylabel('AOD', fontsize=fontsize, fontweight='bold')
		plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
		
		plt.subplot(345)
		plt.scatter(paras[i], outs[1], s=5, c='deepskyblue')
		plt.xticks([])
		plt.yticks([])
		plt.ylabel('SSA', fontsize=fontsize, fontweight='bold')
		plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
		
		plt.subplot(349)
		plt.scatter(paras[i], outs[2], s=5, c='deepskyblue')
		plt.xticks([])
		plt.yticks([])
		plt.xlabel(parameters[i], fontsize=fontsize, fontweight='bold')
		plt.ylabel('g', fontsize=fontsize, fontweight='bold')
		plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
		
		# RCC
		
		plt.subplot(342)
		plt.scatter(r_paras[i], r_outs[0], s=5, c='deepskyblue')
		plt.xticks([])
		plt.yticks([])
		plt.title('RCC', fontsize=fontsize, fontweight='bold')
		plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
		
		plt.subplot(346)
		plt.scatter(r_paras[i], r_outs[1], s=5, c='deepskyblue')
		plt.xticks([])
		plt.yticks([])
		plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
		
		plt.subplot(3,4,10)
		plt.scatter(r_paras[i], r_outs[2], s=5, c='deepskyblue')
		plt.xticks([])
		plt.yticks([])
		plt.xlabel('rank '+parameters[i], fontsize=fontsize, fontweight='bold')
		plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
		
		# PCC
		
		plt.subplot(343)
		plt.scatter(res_paras[0,i], res_outs[0,i], s=5, c='deepskyblue')
		plt.xticks([])
		plt.yticks([])
		plt.title('PCC', fontsize=fontsize, fontweight='bold')
		plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
		
		plt.subplot(347)
		plt.scatter(res_paras[1,i], res_outs[1,i], s=5, c='deepskyblue')
		plt.xticks([])
		plt.yticks([])
		plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
		
		plt.subplot(3,4,11)
		plt.scatter(res_paras[2,i], res_outs[2,i], s=5, c='deepskyblue')
		plt.xticks([])
		plt.yticks([])
		plt.xlabel('res '+parameters[i], fontsize=fontsize, fontweight='bold')
		plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
		
		# PRCC
		
		plt.subplot(344)
		plt.scatter(res_r_paras[0,i], res_r_outs[0,i], s=5, c='deepskyblue')
		plt.xticks([])
		plt.yticks([])
		plt.title('PRCC', fontsize=fontsize, fontweight='bold')
		plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
		
		plt.subplot(348)
		plt.scatter(res_r_paras[1,i], res_r_outs[1,i], s=5, c='deepskyblue')
		plt.xticks([])
		plt.yticks([])
		plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
		
		plt.subplot(3,4,12)
		plt.scatter(res_r_paras[2,i], res_r_outs[2,i], s=5, c='deepskyblue')
		plt.xticks([])
		plt.yticks([])
		plt.xlabel('res rank '+parameters[i], fontsize=fontsize, fontweight='bold')
		plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
		
		if save:
			plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_'+parameters[i]+'.pdf')
		else:
			plt.show()
		plt.close()

def plot_corr(dtime, **args):
	'''
	This function is to plot correlation coefficient ranking bar figure
	input:
		dtime       : data time, string
		**save      : save flag, bool, default False
		**save_path : figure save path, string, defualt 'figure/Corr/'
	output:
	'''
	if 'save' in args:
		save = args['save']
	else:
		save = False
	if 'save_path' in args:
		save_path = args['save_path']
	else:
		save_path = 'figure/Corr/'
	
	parameters = []
	parameters.append('n')
	parameters.append('nI')
	parameters.append('nBC')
	parameters.append('kBC')
	parameters.append('PNSD')
	parameters.append('MS')
	parameters.append('VD')
	parameters.append('CT')
	parameters.append('kappa')
	parameters.append('kappaI')
	parameters.append('rhoBC')
	parameters.append('BCPNSD')
	parameters.append('BCPMSD')
	parameters.append('BCI')
	parameters.append('BCAE')
	
	outputs = ['AOD', 'SSA', 'g']
	
	path1 = 'output/Mie/' + dtime + '/all/'
	CC, CC_P, PCC, PCC_P, RCC, RCC_P, PRCC, PRCC_P = calCorr.cal(path1)
	data1 = np.load(path1+'infos.npy', allow_pickle=True)
	AOD1 = np.zeros(len(data1))
	SSA1 = np.zeros(len(data1))
	g1 = np.zeros(len(data1))
	for i in range(len(data1)):
		AOD1[i] = data1[i]['AOD']
		SSA1[i] = data1[i]['SSA']
		g1[i] = data1[i]['g']
	
	val = np.zeros((len(outputs),len(parameters)))
	
	for i in range(len(parameters)):
		path2 = 'output/Mie/' + dtime + '/' + parameters[i] + '/'
		data2 = np.load(path2+'infos.npy', allow_pickle=True)
		AOD2 = np.zeros(len(data2))
		SSA2 = np.zeros(len(data2))
		g2 = np.zeros(len(data2))
		for j in range(len(data2)):
			AOD2[j] = data2[j]['AOD']
			SSA2[j] = data2[j]['SSA']
			g2[j] = data2[j]['g']
		val[0,i] = (np.nanstd(AOD1)-np.nanstd(AOD2)) / np.nanstd(AOD1)
		val[1,i] = (np.nanstd(SSA1)-np.nanstd(SSA2)) / np.nanstd(SSA1)
		val[2,i] = (np.nanstd(g1)-np.nanstd(g2)) / np.nanstd(g1)
	
	'''
	for i in range(3):
		val[i] = normalize(val[i])
		print(val[i])
	'''
	
	idx = np.ones(len(parameters)) * (len(parameters)+1)
	width = 0.35
	fontsize=14
	
	plt.figure(figsize=(12,10))
	plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.1, hspace=0.4)
	for i in range(len(outputs)):
		plt.subplot(311+i)
		x = idx - rankdata(abs(val[i]), method='ordinal')
		# output ranking
		# to output sensitivity, change code into:
		# y = val[i]
		# need to normalize data
		
		plt.bar(x-width/2, rankdata(abs(val[i]), method='ordinal'), color='deepskyblue', width=width, label='constraint factor')
		plt.bar(x+width/2, rankdata(abs(PRCC[i]), method='ordinal'), color='powderblue', width=width, label='correlation coefficient')
		for j in range(len(val[i])):
			if round(abs(val[i])[j],3) > 0.01:
				plt.text(x[j], max(rankdata(abs(val[i]), method='ordinal')[j],rankdata(abs(PRCC[i]), method='ordinal')[j])+0.5, round(abs(val[i])[j],2), fontsize=fontsize, fontweight='bold', ha='center')
		plt.ylim(0, rankdata(abs(val[i]), method='ordinal')[0]*1.2)
		plt.xticks(x, parameters, rotation=45, fontweight='bold', fontsize=fontsize)
		plt.yticks([])
		plt.ylabel(outputs[i], fontweight='bold', fontsize=fontsize)
		if i == 0:
			plt.legend()
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_all.pdf')
	else:
		plt.show()
	plt.close()

if __name__ == '__main__':
	plot_residual('230604', save=False)
	#plot_corr('230604', save=False)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
