####################################################################################
# INTRODUCTION:
# This code is to plot all output result to compare parameter's sensitivity
# Created by Hebs at 22/6/21/9:33
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import matplotlib.pyplot as plt
import sys
import time

import read11

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

def plot_old(parameter, b, s, **args):
	'''
	This function is plot main function
	input:
		parameter    : parameter name, string
		b            : big range, percent
		s            : small range, percent
		**data_path  : output data storage path, string, default 'output/all/'
		**infos_path : infos.npy storage path, string, default './'
		**save       : whether to save figure, boolean, default False
		**save_path  : figure save path, string, default 'figure/all/'
	output:
		figure
	'''
	if 'data_path' in args:
		data_path = args['data_path']
	else:
		data_path = 'output/all/'
	if 'infos_path' in args:
		infos_path = args['infos_path']
	else:
		infos_path = ''
	if 'save' in args:
		save = args['save']
	else:
		save = False
	if 'save_path' in args:
		save_path = args['save_path']
	else:
		save_path = 'figure/all/'
	
	print('loading...')
	data = np.load(infos_path+'infos.npy', allow_pickle=True)
	num = len(data)
	nrf = np.zeros(num)
	'''
	for i in range(num):
		fn = data_path + str(i) + '.txt'
		try:
			iout11 = read11.read11(filename=fn)
			nrf[i] = iout11['fxdn'][0]-iout11['fxup'][0]
		except:
			nrf[i] = np.nan
	'''
	rate = read_rate(parameter, infos_path)
	
	print('done')
	print('calculating...')
	
	rate_slct = []
	nrf_slct = []
	
	for i in range(num):
		if rate[i] >= (1-s/100) and rate[i] <= (1+s/100):
			rate_slct.append(rate[i])
			#nrf_slct.append(nrf[i])
	
	num_slct = len(rate_slct)
	'''
	val1 = str(round(np.nanmean(nrf))) + ' \u00B1 ' + str(round(np.nanstd(nrf)*10)/10)
	val2 = str(round(np.nanmean(nrf_slct))) + ' \u00B1 ' + str(round(np.nanstd(nrf_slct)*10)/10)
	'''
	print('done')
	print('plotting...')
	
	plt.figure(figsize=(15,5))
	plt.subplots_adjust(left=0.1, right=0.87, top=0.99, bottom=0.1, wspace=0.25)
	'''
	plt.subplot(122)
	plt.scatter(rate, nrf, s=30, c='k', label='30%: '+val1)
	plt.scatter(rate_slct, nrf_slct, s=70, c='r', label='10%: '+val2)
	plt.xlabel(parameter+' rate', fontsize=12, fontweight='bold')
	plt.ylabel('Net Radiative Flux, $W/m^2$', fontsize=12, fontweight='bold')
	plt.xticks(fontsize=12, fontweight='bold')
	plt.yticks(fontsize=12, fontweight='bold')
	plt.legend(bbox_to_anchor=(1.37,1))
	plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
	'''
	x1 = (np.arange(21)/10-1) * b / 100 + 1
	x2 = (np.arange(21)/10-1) * s / 100 + 1
	y1 = np.zeros(21)
	y2 = np.zeros(21)
	xticks = [1-b/100,1-b/200,1,1+b/200,1+b/100]
	for i in range(num):
		ii = round((rate[i]-min(rate)) / (max(rate)-min(rate)) * 20)
		y1[ii] = y1[ii] + 1
	for i in range(num_slct):
		ii = round((rate_slct[i]-min(rate_slct)) / (max(rate_slct)-min(rate_slct)) * 20)
		y2[ii] = y2[ii] + 1
	
	plt.subplot(141)
	plt.bar(x1, y1/max(y1), color='w', edgecolor='k', width=0.029)
	#plt.bar(x2, y2/max(y2), color='r', edgecolor='k', width=0.01, linewidth=0.5)
	plt.xlabel(parameter, fontsize=12, fontweight='bold')
	plt.ylabel('Relative Frequency', fontsize=12, fontweight='bold')
	plt.xticks(xticks, xticks, fontsize=12, fontweight='bold')
	plt.yticks(fontsize=12, fontweight='bold')
	'''
	x1 = np.arange(np.nanmin(nrf), np.nanmax(nrf), (np.nanmax(nrf)-np.nanmin(nrf))/21)
	if x1.shape[0] == 22:
		x1 = x1[:-1]
	x2 = np.arange(np.nanmin(nrf_slct), np.nanmax(nrf_slct), (np.nanmax(nrf_slct)-np.nanmin(nrf_slct))/21)
	if x2.shape[0] == 22:
		x2 = x2[:-1]
	y1 = np.zeros(21)
	y2 = np.zeros(21)
	for i in range(num):
		try:
			ii = round((nrf[i]-np.nanmin(nrf)) / (np.nanmax(nrf)-np.nanmin(nrf)) * 20)
			y1[ii] = y1[ii] + 1
		except:
			ii = 0
	for i in range(num_slct):
		try:
			ii = round((nrf_slct[i]-np.nanmin(nrf_slct)) / (np.nanmax(nrf_slct)-np.nanmin(nrf_slct)) * 20)
			y2[ii] = y2[ii] + 1
		except:
			ii = 0
	plt.subplot(142)
	plt.bar(x1, y1/max(y1), color='w', edgecolor='k', width=5)
	plt.bar(x2, y2/max(y2), color='r', edgecolor='k', width=5)
	plt.xlabel('Net Radiative Flux', fontsize=12, fontweight='bold')
	plt.xticks([])
	plt.yticks(fontsize=12, fontweight='bold')
	'''
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_'+parameter+'.pdf')
	plt.show()
	plt.close()

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

def read_nrf(data_path, num):
	'''
	This function is to read net radiative flux
	input:
		data_path : data path, string
		num       : bin number, int
	output:
		nrf       : net radiative flux, array, W/m2
	'''
	data = np.load(data_path+'infos.npy', allow_pickle=True)
	num = len(data)
	nrf = np.zeros(num)
	
	fn0 = 'output/all/20220825_zero/zero.txt'
	iout11 = read11.read11(filename=fn0)
	nrf0 = iout11['fxdn'][0]-iout11['fxup'][0]
	
	for i in range(num):
		fn = data_path + str(i) + '.txt'
		try:
			iout11 = read11.read11(filename=fn)
			nrf[i] = nrf0 - (iout11['fxdn'][0]-iout11['fxup'][0])
		except:
			nrf[i] = np.nan
	
	return nrf

def plot(parameter, num, **args):
	'''
	This function is to plot all change result
	input:
		parameter    : parameter name, string
		num          : bin number, int
		**data1_path : output data storage path, string, default 'output/all/'
		**data2_path : output data storage path, string, default 'output/all/'
		**save       : whether to save figure, boolean, default False
		**save_path  : figure save path, string, default 'figure/all/'
	output:
		figure
	'''
	if 'data1_path' in args:
		data1_path = args['data1_path']
	else:
		data1_path = 'output/all/'
	if 'data2_path' in args:
		data2_path = args['data2_path']
	else:
		data2_path = 'output/all/'
	if 'save' in args:
		save = args['save']
	else:
		save = False
	if 'save_path' in args:
		save_path = args['save_path']
	else:
		save_path = 'figure/all/'
	
	print('loading...')
	
	data1 = np.load(data1_path+'infos.npy', allow_pickle=True)
	rate1 = read_rate(parameter, data1_path)
	num1 = len(data1)
	nrf1 = np.zeros(num1)
	
	data2 = np.load(data2_path+'infos.npy', allow_pickle=True)
	rate2 = read_rate(parameter, data2_path)
	num2 = len(data2)
	nrf2 = np.zeros(num2)
	
	fn0 = 'output/all/20220825_zero/zero.txt'
	iout11 = read11.read11(filename=fn0)
	nrf0 = iout11['fxdn'][0]-iout11['fxup'][0]
	
	for i in range(num1):
		fn1 = data1_path + str(i) + '.txt'
		try:
			iout11 = read11.read11(filename=fn1)
			nrf1[i] = nrf0 - (iout11['fxdn'][0]-iout11['fxup'][0])
		except:
			nrf1[i] = np.nan
	
	for i in range(num2):
		fn2 = data2_path + str(i) + '.txt'
		try:
			iout11 = read11.read11(filename=fn2)
			nrf2[i] = nrf0 - (iout11['fxdn'][0]-iout11['fxup'][0])
		except:
			nrf2[i] = np.nan
	
	print('done')
	print('calculating...')
	
	# in calculating, do this:
	# count the PDF of nrf and parameter rate
	
	rate1_x, rate1_y = cal_PDF(rate1, num*2)
	rate2_x, rate2_y = cal_PDF(rate2, num)
	nrf1_x, nrf1_y = cal_PDF(nrf1, num+3)
	nrf2_x, nrf2_y = cal_PDF(nrf2, num)
	
	val1 = str(round(np.nanmean(nrf1)*10)/10) + ' \u00B1 ' + str(round(np.nanstd(nrf1)*10)/10)
	val2 = str(round(np.nanmean(nrf2)*10)/10) + ' \u00B1 ' + str(round(np.nanstd(nrf2)*10)/10)
	
	print('done')
	print('plotting...')
	
	plt.figure(figsize=(15,5))
	plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.25)
	
	# plot three sub figure
	# fig1: rate PDF
	# fig2: nrf PDF
	# fig3: rate and nrf 2 dimontional distribution
	# if possible, combine all fig into one
	
	plt.subplot(141)
	plt.bar(rate2_x, rate2_y/max(rate2_y), color='r', width=0.024, linewidth=0.5)
	plt.bar(rate1_x, rate1_y/max(rate1_y), color='', edgecolor='k', width=0.06)
	plt.xlabel(parameter, fontsize=12, fontweight='bold')
	plt.ylabel('Relative Frequency', fontsize=12, fontweight='bold')
	plt.xticks(fontsize=12, fontweight='bold')
	plt.yticks(fontsize=12, fontweight='bold')
	
	plt.subplot(142)
	plt.bar(nrf2_x, nrf2_y/max(nrf2_y), color='r', width=10.5)
	plt.bar(nrf1_x, nrf1_y/max(nrf1_y), color='', edgecolor='k', width=13)
	plt.xlabel(parameter, fontsize=12, fontweight='bold')
	plt.xlabel('Net Radiative Flux', fontsize=12, fontweight='bold')
	plt.xticks(fontsize=12, fontweight='bold')
	plt.yticks(fontsize=12, fontweight='bold')
	
	plt.subplot(122)
	plt.scatter(rate1, nrf1, s=30, c='k', label='50%: '+val1)
	plt.scatter(rate2, nrf2, s=30, c='', edgecolors='r', linewidths=1, label='10%: '+val2)
	plt.xlabel(parameter+' rate', fontsize=12, fontweight='bold')
	plt.ylabel('Net Radiative Flux, $W/m^2$', fontsize=12, fontweight='bold')
	plt.xticks(fontsize=12, fontweight='bold')
	plt.yticks(fontsize=12, fontweight='bold')
	plt.legend()
	plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
	
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_'+parameter+'.pdf')
	plt.show()
	plt.close()

def cal(parameter, s, **args):
	'''
	This function is to calculate parameter's sensitivity
	input:
		parameter    : parameter name, string
		s            : small range, percent
		**data_path  : output data storage path, string, default 'output/all/'
		**infos_path : infos.npy storage path, string, default './'
	output:
		DARF change, percent
	'''
	if 'data_path' in args:
		data_path = args['data_path']
	else:
		data_path = 'output/all/'
	if 'infos_path' in args:
		infos_path = args['infos_path']
	else:
		infos_path = 'output/all/'
	
	data = np.load(infos_path+'infos.npy', allow_pickle=True)
	num = len(data)
	nrf = np.zeros(num)
	for i in range(num):
		fn = data_path + str(i) + '.txt'
		try:
			iout11 = read11.read11(filename=fn)
			nrf[i] = iout11['fxdn'][0]-iout11['fxup'][0]
		except:
			nrf[i] = np.nan
	rate = read_rate(parameter, 'output/all/')
	
	rate_slct = []
	nrf_slct = []
	
	for i in range(num):
		if rate[i] >= (1-s/100) and rate[i] <= (1+s/100):
			rate_slct.append(rate[i])
			nrf_slct.append(nrf[i])
	
	change = (np.nanstd(nrf)-np.nanstd(nrf_slct)) / np.nanstd(nrf) * 100
	change = round(change*100) / 100
	return change

def plot_sub(num, **args):
	'''
	This function is to plot all 2-D sub fig
	input:
		num          : bin number, int
	'''
	if 'save' in args:
		save = args['save']
	else:
		save = False
	if 'save_path' in args:
		save_path = args['save_path']
	else:
		save_path = 'figure/all/'
	
	parameters = []
	parameters.append('AOD')
	parameters.append('SSA')
	parameters.append('PF')
	parameters.append('albedo')
	parameters.append('n')
	parameters.append('k')
	parameters.append('nBC')
	parameters.append('kBC')
	parameters.append('PNSD')
	parameters.append('CT')
	parameters.append('MS')
	parameters.append('VD')
	parameters.append('kappa')
	parameters.append('BCPNSD')
	parameters.append('BCAE')
	
	kappa_path = 'output/all/20220731_10%kappa/'
	CT_path = 'output/all/20220731_10%CT/' # ugly
	MS_path = 'output/all/20220731_10%MS/' # ugly
	VD_path = 'output/all/20220731_10%VD/'
	BCAE_path = 'output/all/20220731_10%BCAE/'
	PNSD_path = 'output/all/20220731_10%PNSD/'
	BCPNSD_path = 'output/all/20220731_10%BCPNSD/'
	n_path = 'output/all/20220731_10%n/'
	k_path = 'output/all/20220731_10%k/' # ugly
	kBC_path = 'output/all/20220731_10%kBC/'
	nBC_path = 'output/all/20220731_10%nBC/'
	PF_path = 'output/all/20220731_10%PF/' # ugly
	albedo_path = 'output/all/20220731_10%albedo/'
	AOD_path = 'output/all/20220815_10%AOD/' # ugly
	SSA_path = 'output/all/20220815_10%SSA/'
	all_path = 'output/all/20220818_all/'
	
	plt.figure(figsize=(12,10))
	plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.05, wspace=0.1, hspace=0.1)
	
	for i in range(len(parameters)):
		data1_path = all_path
		data2_path = eval(parameters[i]+'_path')
		rate1 = read_rate(parameters[i], data1_path)
		rate2 = read_rate(parameters[i], data2_path)
		nrf1 = read_nrf(data1_path, num)
		nrf2 = read_nrf(data2_path, num)
		val = round((np.nanstd(nrf1)-np.nanstd(nrf2))/np.nanstd(nrf1)*1000)/10
		plt.subplot(4,4,i+1)
		plt.scatter(rate1, nrf1, s=5, c='k', label='50%')
		plt.scatter(rate2, nrf2, s=10, c='', edgecolors='r', linewidths=0.5, label='10%')
		plt.xlim(0.25,1.75)
		if i==0:
			plt.xlabel('Change Rate', fontsize=12, fontweight='bold')
			plt.ylabel('Net Radiative Flux, $W/m^2$', fontsize=12, fontweight='bold')
			ax = plt.gca()
			ax.xaxis.set_ticks_position('top')
			ax.xaxis.set_label_position('top')
			plt.xticks([0.5,1,1.5], ['-50%',0,'+50%'], fontsize=12, fontweight='bold')
			plt.yticks([-50,0,50,100], [-50,0,50,100], fontsize=12, fontweight='bold')
			plt.legend(loc='upper right')
		else:
			plt.xticks([0.5,1,1.5], ['','',''], fontsize=12, fontweight='bold')
			plt.yticks([-50,0,50,100], ['','','',''], fontsize=12, fontweight='bold')
		plt.text(0.3, 100, '('+chr(97+i)+') '+parameters[i]+': '+str(val)+'%', fontsize=12, weight='bold')
		plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
		print(round((i+1)/len(parameters)*100), '% done...')
	
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_all.pdf')
	plt.show()
	plt.close()

if __name__ == '__main__':
	'''
	data = np.load('infos.npy', allow_pickle=True)
	for key in data[0].keys():
		if data[0][key] and key != 'Mie_infos':
			key = key[:-5]
			#print(key, ':', cal(key, 10), '%')
			plot(key, 30, 10, data_path='output/all/', infos_path='')
	'''
	kappa_path = 'output/all/20220731_10%kappa/'
	CT_path = 'output/all/20220731_10%CT/' # ugly
	MS_path = 'output/all/20220731_10%MS/' # ugly
	VD_path = 'output/all/20220731_10%VD/'
	BCAE_path = 'output/all/20220731_10%BCAE/'
	PNSD_path = 'output/all/20220731_10%PNSD/'
	BCPNSD_path = 'output/all/20220731_10%BCPNSD/'
	n_path = 'output/all/20220731_10%n/'
	k_path = 'output/all/20220731_10%k/' # ugly
	kBC_path = 'output/all/20220731_10%kBC/'
	nBC_path = 'output/all/20220731_10%nBC/'
	PF_path = 'output/all/20220731_10%PF/' # ugly
	albedo_path = 'output/all/20220731_10%albedo/'
	AOD_path = 'output/all/20220815_10%AOD/' # ugly
	SSA_path = 'output/all/20220815_10%SSA/'
	'''
	parameter = 'n'
	path = eval(parameter+'_path')
	plot(parameter, 10, data1_path='output/all/20220818_all/', data2_path=path, save=True)
	'''
	plot_sub(10, save=True)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
