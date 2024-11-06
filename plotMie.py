####################################################################################
# INTRODUCTION:
# This code is to plot calMie results
# Created by Hebs at 22/9/9/14:45
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

default_outputs = [] # default outputs list
default_outputs.append('AOD')
default_outputs.append('SSA')
default_outputs.append('g')

import numpy as np
import matplotlib.pyplot as plt
import sys
import time

import calCons
import calOAT

def read_rate(parameter, path):
	'''
	This function is to read change rate from infos.npy
	input:
		parameter : parameter name, string
		path      : infos.npy store path, string
	output:
		rate      : parameter change rate, in float, list
	'''
	name = parameter + '_rate'
	data = np.load(path+'infos.npy', allow_pickle=True)
	if name not in data[0]:
		print('No such parameter ' + parameter + '. Please check')
		sys.exit()
	rate = np.zeros(len(data))
	for i in range(len(data)):
		rate[i] = data[i][name]
	return rate

def cal_PDF(data, num):
	'''
	This funciton is to calculating data PDF
	input:
		data   : data to statistic, numpy array
		num    : bin num, int
	output:
		x      : bin, numpy array
		PDF    : posibility distribution function, numpy array
	'''
	data_min = np.nanmin(data)
	data_max = np.nanmax(data)
	x = np.arange(data_min, data_max, (data_max-data_min)/num)
	PDF = np.zeros(len(x))
	for y in data:
		i = 0
		while i<num-1:
			if y>=x[i] and y<x[i+1]:
				PDF[i] = PDF[i] + 1
				break
			i = i + 1
		if y==data_max:
			PDF[-1] = PDF[-1] + 1
	return x, PDF

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

def plot(parameter, num, dtime, **args):
	'''
	This function is to plot all Mie parameter change result
	input:
		parameter    : parameter name, string
		num          : bin number, int
		dtime        : data time, string
		**save       : whether to save figure, boolean, default False
		**save_path  : figure save path, string, default 'figure/Mie/'
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
		save_path = 'figure/Mie/'
	
	data1_path = 'output/Mie/' + dtime + '/all/'
	data2_path = 'output/Mie/' + dtime + '/' + parameter + '/'
	
	print('loading...')
	
	data1 = np.load(data1_path+'infos.npy', allow_pickle=True)
	AOD1 = np.zeros(len(data1))
	SSA1 = np.zeros(len(data1))
	g1 = np.zeros(len(data1))
	for i in range(len(data1)):
		AOD1[i] = data1[i]['AOD']
		SSA1[i] = data1[i]['SSA']
		g1[i] = data1[i]['g']
	rate1 = read_rate(parameter, data1_path)
	
	data2 = np.load(data2_path+'infos.npy', allow_pickle=True)
	AOD2 = np.zeros(len(data2))
	SSA2 = np.zeros(len(data2))
	g2 = np.zeros(len(data2))
	for i in range(len(data2)):
		AOD2[i] = data2[i]['AOD']
		SSA2[i] = data2[i]['SSA']
		g2[i] = data2[i]['g']
	rate2 = read_rate(parameter, data2_path)
	
	AOD1 = AOD1 * 1000 # for better show
	AOD2 = AOD2 * 1000
	
	print('done')
	print('calculating...')
	
	# in calculating, do this:
	# count the PDF of nrf and parameter rate
	
	rate1_x, rate1_y = cal_PDF(rate1, num*2)
	width_rate1 = (rate1_x[-1]-rate1_x[0]) / (num*2)
	rate2_x, rate2_y = cal_PDF(rate2, num)
	width_rate2 = (rate2_x[-1]-rate2_x[0]) / num
	AOD1_x, AOD1_y = cal_PDF(AOD1, num+3)
	width_AOD1 = (AOD1_x[-1]-AOD1_x[0]) / (num+3)
	AOD2_x, AOD2_y = cal_PDF(AOD2, num)
	width_AOD2 = (AOD2_x[-1]-AOD2_x[0]) / num
	SSA1_x, SSA1_y = cal_PDF(SSA1, num+3)
	width_SSA1 = (SSA1_x[-1]-SSA1_x[0]) / (num+3)
	SSA2_x, SSA2_y = cal_PDF(SSA2, num)
	width_SSA2 = (SSA2_x[-1]-SSA2_x[0]) / num
	g1_x, g1_y = cal_PDF(g1, num+3)
	width_g1 = (g1_x[-1]-g1_x[0]) / (num+3)
	g2_x, g2_y = cal_PDF(g2, num)
	width_g2 = (g2_x[-1]-g2_x[0]) / num
	
	val_AOD1 = str(round(np.nanmean(AOD1)*10)/10)+' \u00B1 '+str(round(np.nanstd(AOD1)*100)/100)
	val_AOD2 = str(round(np.nanmean(AOD2)*10)/10)+' \u00B1 '+str(round(np.nanstd(AOD2)*100)/100)
	val_SSA1 = str(round(np.nanmean(SSA1)*100)/100)+' \u00B1 '+str(round(np.nanstd(SSA1)*1000)/1000)
	val_SSA2 = str(round(np.nanmean(SSA2)*100)/100)+' \u00B1 '+str(round(np.nanstd(SSA2)*1000)/1000)
	val_g1 = str(round(np.nanmean(g1)*100)/100)+' \u00B1 '+str(round(np.nanstd(g1)*1000)/1000)
	val_g2 = str(round(np.nanmean(g2)*100)/100)+' \u00B1 '+str(round(np.nanstd(g2)*1000)/1000)
	
	print('done')
	print('plotting...')
	
	fontsize = 20
	
	plt.figure(figsize=(12,15))
	plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1, wspace=0.3, hspace=0.3)
	
	# plot three sub figure
	# fig1: rate PDF
	# fig2: nrf PDF
	# fig3: rate and nrf 2 dimontional distribution
	# if possible, combine all fig into one
	
	# AOD
	
	plt.subplot(341)
	plt.bar(rate2_x, rate2_y/max(rate2_y), color='deepskyblue', width=width_rate2, linewidth=0.5, label='5%')
	plt.bar(rate1_x, rate1_y/max(rate1_y), color='', edgecolor='k', width=width_rate1, label='30%')
	plt.xlabel(parameter, fontsize=12, fontweight='bold')
	plt.ylabel('Relative Frequency', fontsize=12, fontweight='bold')
	plt.xticks(fontsize=fontsize, fontweight='bold')
	plt.yticks(fontsize=fontsize, fontweight='bold')
	plt.legend()
	
	plt.subplot(342)
	plt.bar(AOD2_x, AOD2_y/max(AOD2_y), color='deepskyblue', width=width_AOD2)
	plt.bar(AOD1_x, AOD1_y/max(AOD1_y), color='', edgecolor='k', width=width_AOD1)
	plt.xlabel('AOD, $\\times 10^{-3}$', fontsize=12, fontweight='bold')
	plt.xticks(fontsize=fontsize, fontweight='bold')
	plt.yticks(fontsize=fontsize, fontweight='bold')
	
	plt.subplot(322)
	plt.scatter(rate1, AOD1, s=30, c='k', label='30%: '+val_AOD1)
	plt.scatter(rate2, AOD2, s=30, c='', edgecolors='deepskyblue', linewidths=1, label='5%: '+val_AOD2)
	plt.xlabel(parameter+' rate', fontsize=12, fontweight='bold')
	plt.ylabel('AOD, $\\times 10^{-3}$', fontsize=12, fontweight='bold')
	plt.xticks(fontsize=fontsize, fontweight='bold')
	plt.yticks(fontsize=fontsize, fontweight='bold')
	plt.legend()
	plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
	
	# SSA
	
	plt.subplot(345)
	plt.bar(rate2_x, rate2_y/max(rate2_y), color='deepskyblue', width=width_rate2, linewidth=0.5, label='5%')
	plt.bar(rate1_x, rate1_y/max(rate1_y), color='', edgecolor='k', width=width_rate1, label='30%')
	plt.xlabel(parameter, fontsize=fontsize, fontweight='bold')
	plt.ylabel('Relative Frequency', fontsize=12, fontweight='bold')
	plt.xticks(fontsize=fontsize, fontweight='bold')
	plt.yticks(fontsize=fontsize, fontweight='bold')
	plt.legend()
	
	plt.subplot(346)
	plt.bar(SSA2_x, SSA2_y/max(SSA2_y), color='deepskyblue', width=width_SSA2)
	plt.bar(SSA1_x, SSA1_y/max(SSA1_y), color='', edgecolor='k', width=width_SSA1)
	plt.xlabel('SSA', fontsize=fontsize, fontweight='bold')
	plt.xticks(fontsize=fontsize, fontweight='bold')
	plt.yticks(fontsize=fontsize, fontweight='bold')
	
	plt.subplot(324)
	plt.scatter(rate1, SSA1, s=30, c='k', label='30%: '+val_SSA1)
	plt.scatter(rate2, SSA2, s=30, c='', edgecolors='deepskyblue', linewidths=1, label='5%: '+val_SSA2)
	plt.xlabel(parameter+' rate', fontsize=fontsize, fontweight='bold')
	plt.ylabel('SSA', fontsize=fontsize, fontweight='bold')
	plt.xticks(fontsize=fontsize, fontweight='bold')
	plt.yticks(fontsize=fontsize, fontweight='bold')
	plt.legend()
	plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
	
	# g
	
	plt.subplot(349)
	plt.bar(rate2_x, rate2_y/max(rate2_y), color='deepskyblue', width=width_rate2, linewidth=0.5, label='5%')
	plt.bar(rate1_x, rate1_y/max(rate1_y), color='', edgecolor='k', width=width_rate1, label='30%')
	plt.xlabel(parameter, fontsize=fontsize, fontweight='bold')
	plt.ylabel('Relative Frequency', fontsize=fontsize, fontweight='bold')
	plt.xticks(fontsize=fontsize, fontweight='bold')
	plt.yticks(fontsize=fontsize, fontweight='bold')
	plt.legend()
	
	plt.subplot(3,4,10)
	plt.bar(g2_x, g2_y/max(g2_y), color='deepskyblue', width=width_g2)
	plt.bar(g1_x, g1_y/max(g1_y), color='', edgecolor='k', width=width_g1)
	plt.xlabel('g', fontsize=fontsize, fontweight='bold')
	plt.xticks(fontsize=fontsize, fontweight='bold')
	plt.yticks(fontsize=fontsize, fontweight='bold')
	
	plt.subplot(326)
	plt.scatter(rate1, g1, s=30, c='k', label='30%: '+val_g1)
	plt.scatter(rate2, g2, s=30, c='', edgecolors='deepskyblue', linewidths=1, label='5%: '+val_g2)
	plt.xlabel(parameter+' rate', fontsize=fontsize, fontweight='bold')
	plt.ylabel('g', fontsize=fontsize, fontweight='bold')
	plt.xticks(fontsize=fontsize, fontweight='bold')
	plt.yticks(fontsize=fontsize, fontweight='bold')
	plt.legend()
	plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
	
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_'+parameter+'.pdf')
	else:
		plt.show()
	plt.close()

def plot_sub(dtime, **args):
	'''
	This function is to plot all 2-D sub fig
	input:
		dtime       : data time, string
		**save      : save flag, bool, default False
		**save_path : figure save path, string, defualt 'figure/Mie/'
		**paras     : parameters list for plotting, array, string, default default_paras
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
		save_path = 'figure/Mie/'
	if 'paras' in args:
		paras = args['paras']
	else:
		paras = default_paras
	
	path = 'output/Mie/' + dtime + '/'
	
	data1_path = path + 'all/'
	data1 = np.load(data1_path+'infos.npy', allow_pickle=True)
	run_num = len(data1)
	AOD1 = np.zeros(len(data1))
	SSA1 = np.zeros(len(data1))
	g1 = np.zeros(len(data1))
	for i in range(len(data1)):
		AOD1[i] = data1[i]['AOD']
		SSA1[i] = data1[i]['SSA']
		g1[i] = data1[i]['g']
	
	Cons = calCons.cal(dtime, paras=paras)
	
	plt.figure(figsize=(16,8))
	plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.05, wspace=0.05, hspace=0.05)
	
	fontsize = 13
	AOD_text_height = 0.9
	SSA_text_height = 0.8
	g_text_height = 0.9
	AOD_text2_height = 0.05
	SSA_text2_height = -0.25
	g_text2_height = -0.05
	AOD_yticks = [0.000,0.002,0.004,0.006,0.008]
	AOD_yticks2 = [0,2,4,6,8]
	AOD_yticks_empty = ['','','','','']
	SSA_yticks = [0.8,0.85,0.9,0.95,1]
	SSA_yticks_empty = ['','','','','']
	g_yticks = [0.14,0.15,0.16,0.17,0.18]
	g_yticks_empty = ['','','','','']
	AOD_ylim = cal_lim(AOD1, upper=0.2, down=0.1)
	SSA_ylim = cal_lim(SSA1, upper=0.2, down=-0.65)
	g_ylim = cal_lim(g1, upper=0.2, down=-0.1)
	
	# AOD
	
	for i in range(len(paras)):
		data2_path = path + paras[i] + '/'
		data2 = np.load(data2_path+'infos.npy', allow_pickle=True)
		rate1 = read_rate(paras[i], data1_path)
		rate2 = read_rate(paras[i], data2_path)
		AOD2 = np.zeros(len(data2))
		for j in range(len(data2)):
			AOD2[j] = data2[j]['AOD']
		
		val_AOD = round(Cons[0,i]*10000) / 100
		
		plt.subplot(3,6,i+1)
		plt.scatter(rate1, AOD1, s=5, c='k', label='30%')
		plt.scatter(rate2, AOD2, s=10, c='', edgecolors='deepskyblue', linewidths=0.5, label='5%')
		plt.xlim(0.5,1.5)
		#plt.ylim(AOD_ylim)
		if i==0:
			plt.xlabel('Change Rate', fontsize=fontsize, fontweight='bold')
			plt.ylabel('AOD, $\\times10^{-3}$', fontsize=fontsize, fontweight='bold')
			ax = plt.gca()
			ax.xaxis.set_ticks_position('top')
			ax.xaxis.set_label_position('top')
			plt.xticks([0.7,1,1.3], ['-30%','1','30%'], fontsize=fontsize, fontweight='bold')
			#plt.yticks(AOD_yticks, AOD_yticks2, fontsize=fontsize, fontweight='bold')
			plt.legend(loc='upper right', fontsize=fontsize-3)
			plt.text(0.83, AOD_text2_height*AOD_ylim[1]+(1-AOD_text2_height)*AOD_ylim[0], 'run num = '+str(run_num), fontsize=fontsize, weight='bold')
		else:
			plt.xticks([0.7,1,1.3], ['','',''], fontsize=fontsize, fontweight='bold')
			#plt.yticks(AOD_yticks, AOD_yticks_empty, fontsize=fontsize, fontweight='bold')
		plt.text(0.55, AOD_text_height*AOD_ylim[1]+(1-AOD_text_height)*AOD_ylim[0], '('+chr(97+i)+') '+paras[i]+': '+str(val_AOD)+'%', fontsize=fontsize, weight='bold')
		plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
		print(round((i+1)/len(paras)*100), '% done...')
	
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_allAOD.pdf')
	else:
		plt.show()
	plt.close()
	
	# SSA
	
	plt.figure(figsize=(12,10))
	plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.05, wspace=0.05, hspace=0.05)
	
	for i in range(len(paras)):
		data2_path = path + paras[i] + '/'
		data2 = np.load(data2_path+'infos.npy', allow_pickle=True)
		rate1 = read_rate(paras[i], data1_path)
		rate2 = read_rate(paras[i], data2_path)
		SSA2 = np.zeros(len(data2))
		for j in range(len(data2)):
			SSA2[j] = data2[j]['SSA']
		
		val_SSA = round(Cons[1,i]*10000) / 100
		
		plt.subplot(3,6,i+1)
		plt.scatter(rate1, SSA1, s=5, c='k', label='30%')
		plt.scatter(rate2, SSA2, s=10, c='', edgecolors='deepskyblue', linewidths=0.5, label='5%')
		plt.xlim(0.5,1.5)
		#plt.ylim(SSA_ylim)
		if i==0:
			plt.xlabel('Change Rate', fontsize=fontsize, fontweight='bold')
			plt.ylabel('SSA', fontsize=fontsize, fontweight='bold')
			ax = plt.gca()
			ax.xaxis.set_ticks_position('top')
			ax.xaxis.set_label_position('top')
			plt.xticks([0.7,1,1.3], ['-30%','1','30%'], fontsize=fontsize, fontweight='bold')
			#plt.yticks(SSA_yticks, SSA_yticks, fontsize=fontsize, fontweight='bold')
			plt.legend(loc='upper right', fontsize=fontsize-3)
			plt.text(0.83, SSA_text2_height*SSA_ylim[1]+(1-SSA_text2_height)*SSA_ylim[0], 'run num = '+str(run_num), fontsize=fontsize, weight='bold')
		else:
			plt.xticks([0.7,1,1.3], ['','',''], fontsize=fontsize, fontweight='bold')
			#plt.yticks(SSA_yticks, SSA_yticks_empty, fontsize=fontsize, fontweight='bold')
		plt.text(0.55, SSA_text_height*SSA_ylim[1]+(1-SSA_text_height)*SSA_ylim[0], '('+chr(97+i)+') '+paras[i]+': '+str(val_SSA)+'%', fontsize=fontsize, weight='bold')
		plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
		print(round((i+1)/len(paras)*100), '% done...')
	
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_allSSA.pdf')
	else:
		plt.show()
	plt.close()
	
	# g
	
	plt.figure(figsize=(12,10))
	plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.05, wspace=0.05, hspace=0.05)
	
	for i in range(len(paras)):
		data2_path = path + paras[i] + '/'
		data2 = np.load(data2_path+'infos.npy', allow_pickle=True)
		rate1 = read_rate(paras[i], data1_path)
		rate2 = read_rate(paras[i], data2_path)
		g2 = np.zeros(len(data2))
		for j in range(len(data2)):
			g2[j] = data2[j]['g']
		
		val_g = round(Cons[2,i]*10000) / 100
		
		plt.subplot(3,6,i+1)
		plt.scatter(rate1, g1, s=5, c='k', label='30%')
		plt.scatter(rate2, g2, s=10, c='', edgecolors='deepskyblue', linewidths=0.5, label='5%')
		plt.xlim(0.5,1.5)
		#plt.ylim(g_ylim)
		if i==0:
			plt.xlabel('Change Rate', fontsize=fontsize, fontweight='bold')
			plt.ylabel('g', fontsize=fontsize, fontweight='bold')
			ax = plt.gca()
			ax.xaxis.set_ticks_position('top')
			ax.xaxis.set_label_position('top')
			plt.xticks([0.7,1,1.3], ['-30%','1','30%'], fontsize=fontsize, fontweight='bold')
			#plt.yticks(g_yticks, g_yticks, fontsize=fontsize, fontweight='bold')
			plt.legend(loc='upper right', fontsize=fontsize-3)
			plt.text(0.83, g_text2_height*g_ylim[1]+(1-g_text2_height)*g_ylim[0], 'run num = '+str(run_num), fontsize=fontsize, weight='bold')
		else:
			plt.xticks([0.7,1,1.3], ['','',''], fontsize=fontsize, fontweight='bold')
			#plt.yticks(g_yticks, g_yticks_empty, fontsize=fontsize, fontweight='bold')
		plt.text(0.55, g_text_height*g_ylim[1]+(1-g_text_height)*g_ylim[0], '('+chr(97+i)+') '+paras[i]+': '+str(val_g)+'%', fontsize=fontsize, weight='bold')
		plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
		print(round((i+1)/len(paras)*100), '% done...')
	
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_allg.pdf')
	else:
		plt.show()
	plt.close()

def norm(x, loc, scale):
	'''
	This function is normal distribution function
	input:
		x     : x, float
		loc   : normal distribution location parameter, float
		scale : normal distribution scale parameter, float
	output:
		y     : f(x), float
	'''
	y = 1 / np.sqrt(2*np.pi) / scale * np.exp(-(x-loc)**2/2/scale**2)
	return y

def norm_cdf(x, loc, scale, **args):
	'''
	This function is normal distribution cumulative distribution function
	input:
		x     : x, float
		loc   : normal distribution location parameter, float
		scale : normal distribution scale parameter, float
		**dx  : steps, float, default 1e-3
	output:
		y     : f(x), float
	'''
	if 'dx' in args:
		dx = args['dx']
	else:
		dx = 1e-2
	
	x = np.arange(-1e2, x, dx)
	y = 0
	for i in range(len(x)):
		dy = norm(x[i], loc, scale) * dx
		y += dy
	return y

def plot_OATvsCP(dtime, par_range, **args):
	'''
	This function is to plot OAT and CP sensitivity results
	input:
		dtime		: runDARF.py output time, string, example: '230101'
		par_range	: parameter change rate range, float
		**paras		: parameters list, array, string, default default_paras
		**outputs	: outputs list, array, string, default default_outputs
		**save		: save flag, bool, default False
		**save_path	: save path, string, default 'figure/DARF/
	output:
		figure
	'''
	if 'paras' in args:
		paras = args['paras']
	else:
		paras = default_paras
	if 'outputs' in args:
		outputs = args['outputs']
	else:
		outputs = default_outputs
	if 'save' in args:
		save = args['save']
	else:
		save = False
	if 'save_path' in args:
		save_path = args['save_path']
	else:
		save_path = 'figure/DARF/'
	
	# read data
	
	OAT = calOAT.cal(dtime, outputs=outputs, paras=paras)['OAT']
	CP = calCons.cal(dtime, par_range, outputs=outputs, paras=paras)['Cons']
	
	vs = abs((OAT-CP)/CP)
	vs[np.isnan(vs)] = -1 # turn np.nan to 1
	vs[np.isinf(vs)] = -1 # turn np.inf to 1
	out_num = len(outputs)
	par_num = len(paras)
	
	# change parameter names
	'''
	for i in range(par_num):
		if paras[i]=='n':
			paras[i] = 'n_shell'
		if paras[i]=='nI':
			paras[i] = 'n_shell H1'
		if paras[i]=='nI2':
			paras[i] = 'n_shell H2'
		if paras[i]=='nBC':
			paras[i] = 'n_BC'
		if paras[i]=='kBC':
			paras[i] = 'k_BC'
		#if paras[i]=='MS':
		#	paras[i] = 'Mixing state'
		#if paras[i]=='VD':
		#	paras[i] = 'Vertical distribution'
		#if paras[i]=='CT':
		#	paras[i] = 'Coating thickness'
		if paras[i]=='kappaI':
			paras[i] = 'kappa H1'
		if paras[i]=='kappaI2':
			paras[i] = 'kappa H2'
		if paras[i]=='rhoBC':
			paras[i] = '\u03C1_BC'
		#if paras[i]=='BCAE':
		#	paras[i] = 'BC absorbing enhancement'
		if paras[i]=='amb':
		#	paras[i] = 'Atmospheric parameters'
			paras[i] = 'AEP'
	'''
	# set plot parameters
	
	fs = 8 # fontsize
	fw = 'bold' # fontweight
	width = 0.4 # bar width
	lw = 2 # line width
	rotation = 45 # x labels rotaion angle
	
	ylim = np.zeros(out_num)
	for i in range(out_num):
		ylim[i] = max(np.max(OAT[i]),np.max(CP[i]))
	
	plt.figure(figsize=(10,6))
	plt.subplots_adjust(left=0.1, right=0.9, top=0.95, bottom=0.2, wspace=0.1, hspace=0.15)
	
	# start plotting
	
	for i in range(out_num):
		ax = plt.subplot(out_num,1,i+1)
		# main plot function
		for j in range(par_num):
			if i==0 and j==0:
				plt.bar(j+1-width/2, CP[i,j], width=width, label='CP', color='k')
				plt.bar(j+1+width/2, OAT[i,j], width=width, label='OAT', color='gray')
				plt.legend(fontsize=fs, loc='upper right')
			else:
				plt.bar(j+1-width/2, CP[i,j], width=width, color='k')
				plt.bar(j+1+width/2, OAT[i,j], width=width, color='gray')
			if vs[i,j]!=-1: # -1 means no sensitivity
				plt.text(j+1-width*0.5, max(CP[i,j],OAT[i,j])+ylim[i]*0.05, str(round(vs[i,j]*100))+'%', fontsize=fs, fontweight=fw, color='gray')
			
		# y lim
		plt.ylim([0,ylim[i]*1.3])
		# y label
		if i==0:
			plt.ylabel('dAOD/dx%', fontsize=fs, fontweight=fw)
		elif i==1:
			plt.ylabel('dSSA/dx%', fontsize=fs, fontweight=fw)
		else:
			plt.ylabel('dg/dx%', fontsize=fs, fontweight=fw)
		
		# x label
		if i==out_num-1:
			plt.xticks(np.arange(par_num)+1, paras, fontsize=fs, fontweight=fw, rotation=rotation, horizontalalignment='right')
		else:
			plt.xticks([])
		ax.yaxis.tick_right()
	
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_SA.pdf')
	else:
		plt.show()

if __name__ == '__main__':
	paras = []
	paras.append('n')
	paras.append('nH1')
	paras.append('nBC')
	paras.append('kBC')
	paras.append('PNSD')
	paras.append('MS')
	paras.append('MSH1')
	paras.append('MSH2')
	paras.append('CT')
	paras.append('CTH1')
	paras.append('kappa')
	paras.append('kappaH1')
	paras.append('rhoBC')
	paras.append('BCPNSD')
	paras.append('BCAE')
	
	#plot(paras[0], 5, '230604')
	#plot_sub('230604', paras=['n', 'nI2', 'kappa'], save=False)
	plot_OATvsCP('240603', 0.05, save=False)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
