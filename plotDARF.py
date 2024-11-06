####################################################################################
# INTRODUCTION:
# This code is to plot runDARF outputs results
# Created by Hebs at 24/3/17/10:27
# Contact: hebishuo@pku.edu.cn
####################################################################################
'''
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
default_paras.append('BCAE')
default_paras.append('amb')
default_paras.append('albedo')
'''
default_paras = [] # default parameters list
'''
default_paras.append('n')
default_paras.append('nH1')
default_paras.append('nBC')
default_paras.append('kBC')
default_paras.append('PNSD')
default_paras.append('MS')
default_paras.append('MSH1')
default_paras.append('MSH2')
default_paras.append('VD')
default_paras.append('CT')
default_paras.append('CTH1')
default_paras.append('kappa')
default_paras.append('kappaH1')
default_paras.append('rhoBC')
default_paras.append('BCPNSD')
default_paras.append('BCAE')
default_paras.append('amb')
default_paras.append('albedo')
'''
default_paras.append('n')
default_paras.append('PNSD')
default_paras.append('kappa')

default_paras.append('nBC')
default_paras.append('kBC')
default_paras.append('rhoBC')
default_paras.append('BCPNSD')

default_paras.append('MS')
default_paras.append('CT')
default_paras.append('BCAE')

default_paras.append('VD')
default_paras.append('amb')
default_paras.append('albedo')

default_paras.append('nH1')
default_paras.append('kappaH1')
default_paras.append('CTH1')
default_paras.append('MSH1')
default_paras.append('MSH2')

default_outputs = [] # default outputs list
default_outputs.append('RF_top')
default_outputs.append('RF_bot')
default_outputs.append('dtau_bot')
default_outputs.append('waer_bot')
default_outputs.append('g_bot')

import numpy as np
import matplotlib.pyplot as plt
import time
import calSA
import calCorr
import calANOVA

def plot(dtime, par_range, **args):
	'''
	This function is to plot sensitivity results
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
	
	OAT = calSA.read_OAT(dtime, outputs=outputs, paras=paras)['OAT']
	CP = calSA.cal_CP(dtime, par_range, outputs=outputs, paras=paras)['CP']
	
	vs = abs((OAT-CP)/CP)
	vs[np.isnan(vs)] = -1 # turn np.nan to 1
	vs[np.isinf(vs)] = -1 # turn np.inf to 1
	out_num = len(outputs)
	par_num = len(paras)
	
	# change parameter names
	
	for i in range(par_num):
		if paras[i]=='n':
			paras[i] = 'n shell'
		if paras[i]=='nH1':
			paras[i] = 'n H1'
		if paras[i]=='nBC':
			paras[i] = 'n BC'
		if paras[i]=='kBC':
			paras[i] = 'k BC'
		#if paras[i]=='MS':
		#	paras[i] = 'Mixing state'
		if paras[i]=='MSH1':
			paras[i] = 'MS H1'
		#if paras[i]=='VD':
		#	paras[i] = 'Vertical distribution'
		#if paras[i]=='CT':
		#	paras[i] = 'Coating thickness'
		if paras[i]=='CTH1':
			paras[i] = 'CT H1'
		if paras[i]=='kappaH1':
			paras[i] = 'kappa H1'
		if paras[i]=='rhoBC':
			paras[i] = '\u03C1 BC'
		#if paras[i]=='BCAE':
		#	paras[i] = 'BC absorbing enhancement'
		if paras[i]=='amb':
		#	paras[i] = 'Atmospheric parameters'
			paras[i] = 'AEP'
	
	# set plot parameters
	
	fs = 8 # fontsize
	fw = 'bold' # fontweight
	width = 0.4 # bar width
	lw = 2 # line width
	rotation = 45 # x labels rotaion angle
	colors = ['k','k','k','r','r','r','r','b','b','b','g','g','g','g','g']
	colors = ['k','k','k','r','r','r','r','b','b','b','y','y','y','g','g','g','g','g']
	
	ylim = np.zeros(out_num)
	for i in range(out_num):
		ylim[i] = max(np.max(OAT[i]),np.max(CP[i]))
	
	plt.figure(figsize=(8,out_num+1))
	plt.subplots_adjust(left=0.1, right=0.9, top=0.95, bottom=0.2, wspace=0.1, hspace=0.15)
	
	# start plotting
	
	for i in range(out_num):
		ax = plt.subplot(out_num,1,i+1)
		# main plot function
		for j in range(par_num):
			if i==0 and j==0:
				plt.bar(j+1-width/2, CP[i,j], width=width, label='CP', color=colors[j])
				plt.bar(j+1+width/2, OAT[i,j], width=width, label='OAT', edgecolor=colors[j], color='w')
				plt.legend(fontsize=fs, loc='upper right')
			else:
				plt.bar(j+1-width/2, CP[i,j], width=width, color=colors[j])
				plt.bar(j+1+width/2, OAT[i,j], width=width, edgecolor=colors[j], color='w')
			if vs[i,j]!=-1: # -1 means no sensitivity
				plt.text(j+1-width*0.7, max(CP[i,j],OAT[i,j])+ylim[i]*0.05, str(round(vs[i,j]*100))+'%', fontsize=fs*0.8, fontweight=fw, color='gray')
			
		# y lim
		plt.ylim([0,ylim[i]*1.3])
		# y tick
		ax.yaxis.tick_right()
		plt.yticks(fontsize=fs, fontweight=fw)
		# y label
		if outputs[i]=='RF_top':
			plt.ylabel('top RF, dRF/dx%', fontsize=fs, fontweight=fw)
		elif outputs[i]=='RF_bot':
			plt.ylabel('bottom RF, dRF/dx%', fontsize=fs, fontweight=fw)
		elif outputs[i]=='dtau_bot':
			plt.ylabel('dAOD/dx%', fontsize=fs, fontweight=fw)
		elif outputs[i]=='waer_bot':
			plt.ylabel('dSSA/dx%', fontsize=fs, fontweight=fw)
		else:
			plt.ylabel('dg/dx%', fontsize=fs, fontweight=fw)
		
		# x label
		if i==out_num-1:
			ax=plt.gca()
			plt.xticks(np.arange(par_num)+1, paras, fontsize=fs, fontweight=fw, rotation=rotation, horizontalalignment='right')
			for xtick, color in zip(ax.get_xticklabels(), colors):
				xtick.set_color(color)
			# by chatGPT
		else:
			plt.xticks([])
	
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_SA.pdf')
	else:
		plt.show()

def plot2(dtime, par_range, **args):
	'''
	This function is to plot aerosol optical parameters' sensitivity changed with height
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
	
	CP_infos = calSA.cal_CP3(dtime, par_range, outputs=outputs, paras=paras)
	CP = CP_infos['CP']
	nn = CP_infos['nn']
	out_num = len(outputs)
	par_num = len(paras)
	print(CP)
	print(paras)
	
	# change parameter names
	
	for i in range(par_num):
		if paras[i]=='n':
			paras[i] = 'n_shell'
		if paras[i]=='nI':
			paras[i] = 'n_shell I'
		if paras[i]=='nI2':
			paras[i] = 'n_shell I2'
		if paras[i]=='nBC':
			paras[i] = 'n_BC'
		if paras[i]=='kBC':
			paras[i] = 'k_BC'
		if paras[i]=='MS':
			paras[i] = 'Mixing state'
		if paras[i]=='VD':
			paras[i] = 'Vertical distribution'
		if paras[i]=='CT':
			paras[i] = 'Coating thickness'
		if paras[i]=='kappaI':
			paras[i] = 'kappa I'
		if paras[i]=='kappaI2':
			paras[i] = 'kappa I2'
		if paras[i]=='rhoBC':
			paras[i] = '\u03C1_BC'
		if paras[i]=='BCAE':
			paras[i] = 'BC absorbing enhancement'
		if paras[i]=='amb':
			paras[i] = 'Atmospheric parameters'
			
	
	# set plot parameters
	
	fs = 6 # fontsize
	fw = 'bold' # fontweight
	width = 0.4 # bar width
	lw = 2 # line width
	rotation = 45 # x labels rotaion angle
	
	ylim = np.zeros(out_num)
	for i in range(out_num):
		ylim[i] = np.max(CP[i])
	
	plt.figure(figsize=(12,6))
	plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.2, wspace=0.1, hspace=0.15)
	
	# start plotting
	
	for i in range(out_num):
		# main plot function
		for j in range(nn):
			plt.subplot(nn,out_num,i+1+j*out_num)
			for k in range(par_num-1): # delete albedo
				plt.bar(k+1, CP[i,j,k], width=width, color='k')
				#plt.text(k+1-width*0.7, CP[i,j,k]+ylim[i]*0.25, str(round(CP[i,j,k]*100)/100), fontsize=fs, fontweight=fw, color='k')
			# title
			if i==0 and j==0:
				plt.title('AOD', fontsize=fs+2, fontweight=fw)
			if i==1 and j==0:
				plt.title('SSA', fontsize=fs+2, fontweight=fw)
			if i==2 and j==0:
				plt.title('g', fontsize=fs+2, fontweight=fw)
			# y lim
			plt.ylim([0,ylim[i]*1.2])
			# y label
			if i==0 and j==0:
				plt.ylabel('dy/dx', fontsize=fs+1, fontweight=fw)
			# y ticks
			plt.yticks([])
			# x label
			if j==nn-1:
				plt.xticks(np.arange(par_num-1)+1, paras[:-1], fontsize=fs, fontweight=fw, rotation=rotation, horizontalalignment='right')
			else:
				plt.xticks([])
	
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_SA.pdf')
	else:
		plt.show()

def plot3(dtime, par_range, **args):
	'''
	This function is to plot aerosol optical parameters' sensitivity changed with height
	figure x-axis and y-axis are value and height
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
	
	CP_infos = calSA.cal_CP3(dtime, par_range, outputs=outputs)
	CP = CP_infos['CP'] # top to bottom
	nn = CP_infos['nn']
	h = [0.0000E+00, 2.0748E-01, 4.7365E-01, 
	8.8059E-01, 1.5248E+00, 2.5132E+00, 
	3.9613E+00, 5.9916E+00, 8.7329E+00, 
	1.2320E+01, 1.6892E+01, 2.2594E+01, 
	2.9574E+01, 3.7986E+01, 4.7986E+01] # from bottom to top
	h = np.array(h)
	h = h[::-1] # from top to bottom
	out_num = len(outputs)
	par_num = len(paras)
	
	# CP need to normalize
	
	for i in range(out_num):
		for j in range(par_num):
			CP[i,:,j] = CP[i,:,j] / sum(CP[i,:,j])
	
	# change parameter names
	
	for i in range(par_num):
		if paras[i]=='n':
			paras[i] = 'n_shell'
		if paras[i]=='nI':
			paras[i] = 'n_shell I'
		if paras[i]=='nI2':
			paras[i] = 'n_shell I2'
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
			paras[i] = 'kappa I'
		if paras[i]=='kappaI2':
			paras[i] = 'kappa I2'
		if paras[i]=='rhoBC':
			paras[i] = '\u03C1_BC'
		#if paras[i]=='BCAE':
		#	paras[i] = 'BC absorbing enhancement'
		if paras[i]=='amb':
		#	paras[i] = 'Atmospheric parameters'
			paras[i] = 'AEP'
	
	# change output names
	
	for i in range(out_num):
		if outputs[i]=='dtau':
			outputs[i] = 'AOD'
		if outputs[i]=='waer':
			outputs[i] = 'SSA'
	
	# set plot parameters
	
	fs = 10 # fontsize
	fw = 'bold' # fontweight
	width = 0.4 # bar width
	lw = 2 # line width
	rotation = 45 # x labels rotaion angle
	linestyle = ['-','--','-.',':']
	
	ylim = np.zeros(out_num)
	for i in range(out_num):
		ylim[i] = np.max(CP[i])
	
	plt.figure(figsize=(10,3))
	plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.2, wspace=0.15, hspace=0.1)
	
	# start plotting
	
	for i in range(out_num):
		plt.subplot(1, out_num, i+1)
		for j in range(par_num):
			# main function
			plt.plot(CP[i,:,j], h, label=paras[j], lw=lw)
		# title
		plt.title(outputs[i], fontweight=fw, fontsize=fs)
		# x label
		if i==0:
			plt.xlabel('Sensitivity', fontweight=fw, fontsize=fs)
		plt.xticks([0,0.1,0.2,0.3], fontweight=fw, fontsize=fs)
		plt.xlim(0,0.3)
		# y label
		if i==0:
			plt.ylabel('Height, km', fontweight=fw, fontsize=fs)
			plt.yticks(fontweight=fw, fontsize=fs)
			plt.legend(fontsize=fs*0.9)
		else:
			plt.yticks([])
		plt.ylim(0)
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_SA.pdf')
	else:
		plt.show()

def plot_all(dtime, par_range, **args):
	'''
	This function is to plot aerosol optical parameters' sensitivity 
	with all 4 SA method, including OAT, ANOVA, PRCC and CP
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

	OAT = calSA.read_OAT(dtime, outputs=outputs, paras=paras)['OAT']
	ANOVA = np.abs(calANOVA.cal(dtime, path='output/DARF/', outputs=outputs, paras=paras))
	PRCC = np.abs(calCorr.cal(dtime, path='output/DARF/', outputs=outputs, paras=paras)['PRCC'])
	CP = calSA.cal_CP(dtime, par_range, outputs=outputs, paras=paras)['CP']
	# need to be nomalized
	for i in range(len(outputs)):
		OAT[i]		 = OAT[i]	 / sum(OAT[i])
		ANOVA[i]	 = ANOVA[i]	 / sum(ANOVA[i])
		PRCC[i]		 = PRCC[i]	 / sum(PRCC[i])
		CP[i]		 = CP[i]	 / sum(CP[i])
	
	out_num = len(outputs)
	par_num = len(paras)
	
	# change parameter names
	
	for i in range(par_num):
		if paras[i]=='n':
			paras[i] = 'n_shell'
		if paras[i]=='nI':
			paras[i] = 'n_shell I'
		if paras[i]=='nI2':
			paras[i] = 'n_shell I2'
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
			paras[i] = 'kappa I'
		if paras[i]=='kappaI2':
			paras[i] = 'kappa I2'
		if paras[i]=='rhoBC':
			paras[i] = '\u03C1_BC'
		#if paras[i]=='BCAE':
		#	paras[i] = 'BC absorbing enhancement'
		if paras[i]=='amb':
		#	paras[i] = 'Atmospheric parameters'
			paras[i] = 'AEP'
	
	# set plot parameters
	
	fs = 8 # fontsize
	fw = 'bold' # fontweight
	width = 0.2 # bar width
	lw = 2 # line width
	rotation = 45 # x labels rotaion angle
	
	ylim = np.zeros(out_num)
	for i in range(out_num):
		ylim[i] = max(np.max(OAT[i]),np.max(ANOVA[i]),np.max(PRCC[i]),np.max(CP[i]))
	
	plt.figure(figsize=(10,4))
	plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.4, wspace=0.1, hspace=0.15)
	
	# start plotting
	
	for i in range(out_num):
		plt.subplot(out_num,1,i+1)
		# main plot function
		for j in range(par_num):
			if i==0 and j==0:
				plt.bar(j+1-width*1.5, OAT[i,j], width=width, label='OAT', color='k')
				plt.bar(j+1-width*0.5, ANOVA[i,j], width=width, label='ANOVA', color='dimgray')
				plt.bar(j+1+width*0.5, PRCC[i,j], width=width, label='PRCC', color='gray')
				plt.bar(j+1+width*1.5, CP[i,j], width=width, label='CP', color='lightgray')
				plt.legend(fontsize=fs-1)
			else:
				plt.bar(j+1-width*1.5, OAT[i,j], width=width, color='k')
				plt.bar(j+1-width*0.5, ANOVA[i,j], width=width, color='dimgray')
				plt.bar(j+1+width*0.5, PRCC[i,j], width=width, color='gray')
				plt.bar(j+1+width*1.5, CP[i,j], width=width, color='lightgray')
		
		# y lim
		plt.ylim([0,ylim[i]*1.1])
		# y label
		if i==0:
			plt.ylabel('top RF, dRF/dx', fontsize=fs, fontweight=fw)
		elif i==1:
			plt.ylabel('bottom RF, dRF/dx', fontsize=fs, fontweight=fw)
		'''
		if i==0:
			plt.ylabel('dAOD/dx, $10^{-1}$', fontsize=fs-2, fontweight=fw)
		elif i==1:
			plt.ylabel('dSSA/dx, $10^{-1}$', fontsize=fs-2, fontweight=fw)
		else:
			plt.ylabel('dg/dx, $10^{-1}$', fontsize=fs-2, fontweight=fw)
		'''
		# x label
		if i==out_num-1:
			plt.xticks(np.arange(par_num)+1, paras, fontsize=fs, fontweight=fw, rotation=rotation, horizontalalignment='right')
		else:
			plt.xticks([])
	
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_SA.pdf')
	else:
		plt.show()

def plot4(dtime, par_range, **args):
	'''
	This function is to plot sensitivity rank results
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
	
	CP = calSA.cal_CP(dtime, par_range, outputs=outputs, paras=paras)['CP']
	
	out_num = len(outputs)
	par_num = len(paras)
	
	# sort data
	
	weight = CP[0] / sum(CP[0]) + CP[1] / sum(CP[1])
	combined = zip(weight, CP[0], CP[1], paras)
	sorted_combined = sorted(combined, key=lambda x: x[0], reverse=True)
	weight, CP[0], CP[1], paras = zip(*sorted_combined)
	paras = np.array(paras) # from chatGPT
	
	# change parameter names
	
	for i in range(par_num):
		if paras[i]=='n':
			paras[i] = 'n_shell'
		if paras[i]=='nI':
			paras[i] = 'n_shell I'
		if paras[i]=='nI2':
			paras[i] = 'n_shell I2'
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
			paras[i] = 'kappa I'
		if paras[i]=='kappaI2':
			paras[i] = 'kappa I2'
		if paras[i]=='rhoBC':
			paras[i] = '\u03C1_BC'
		#if paras[i]=='BCAE':
		#	paras[i] = 'BC absorbing enhancement'
		if paras[i]=='amb':
		#	paras[i] = 'Atmospheric parameters'
			paras[i] = 'AEP'
	
	# set plot parameters
	
	fs = 10 # fontsize
	fw = 'bold' # fontweight
	width = 0.4 # bar width
	lw = 2 # line width
	rotation = 45 # x labels rotaion angle
	
	plt.figure(figsize=(8,3))
	plt.subplots_adjust(left=0.1, right=0.9, top=0.95, bottom=0.25, wspace=0.1, hspace=0.1)
	
	# start plotting
	
	# main plot function
	ax1 = plt.subplot(111)
	ax2 = ax1.twinx()
	for i in range(par_num):
		if i==0 :
			a = ax1.bar(i+1-width/2, CP[0,i], width=width, label='top RF', color='k')
			b = ax2.bar(i+1+width/2, CP[1,i], width=width, label='bottom RF', color='gray')
			plt.legend([a,b], [c.get_label() for c in [a,b]], fontsize=fs) # from chatGPT
		else:
			ax1.bar(i+1-width/2, CP[0,i], width=width, color='k')
			ax2.bar(i+1+width/2, CP[1,i], width=width, color='gray')
	# y lim
	ax1.set_ylim(0, np.max(CP[0])*1.15)
	ax2.set_ylim(0, np.max(CP[1])*1.15)
	# y label
	ax1.tick_params(axis='y', labelsize=fs)
	ax1.set_ylabel('Top RF, dRF/dx%', fontsize=fs, fontweight=fw)
	ax2.tick_params(axis='y', labelsize=fs)
	ax2.set_ylabel('Bottom RF, dRF/dx%', fontsize=fs, fontweight=fw)
	'''
	plt.ylabel('Sensitivity, $10^{-1}$', fontsize=fs-2, fontweight=fw)
	'''
	# x label
	ax1.set_xticks(np.arange(par_num)+1, paras, fontsize=fs, fontweight=fw, rotation=rotation, horizontalalignment='right')
	
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_SA.pdf')
	else:
		plt.show()

def plot5(dtime, par_range, **args):
	'''
	This function is to plot CP sensitivity results only
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
	
	CP = calSA.cal_CP(dtime, par_range, outputs=outputs, paras=paras)['CP']
	
	out_num = len(outputs)
	par_num = len(paras)
	
	# sort data
	'''
	weight = CP[0] / sum(CP[0]) + CP[1] / sum(CP[1]) + CP[2] / sum(CP[2])
	combined = zip(weight, CP[0], CP[1], CP[2], paras)
	sorted_combined = sorted(combined, key=lambda x: x[0], reverse=True)
	weight, CP[0], CP[1], CP[2], paras = zip(*sorted_combined)
	paras = np.array(paras) # from chatGPT
	'''
	# change parameter names
	
	for i in range(par_num):
		if paras[i]=='n':
			paras[i] = 'n shell'
		if paras[i]=='nH1':
			paras[i] = 'n H1'
		if paras[i]=='nBC':
			paras[i] = 'n BC'
		if paras[i]=='kBC':
			paras[i] = 'k BC'
		#if paras[i]=='MS':
		#	paras[i] = 'Mixing state'
		if paras[i]=='MSH1':
			paras[i] = 'MS H1'
		#if paras[i]=='VD':
		#	paras[i] = 'Vertical distribution'
		#if paras[i]=='CT':
		#	paras[i] = 'Coating thickness'
		if paras[i]=='CTH1':
			paras[i] = 'CT H1'
		if paras[i]=='kappaH1':
			paras[i] = 'kappa H1'
		if paras[i]=='rhoBC':
			paras[i] = '\u03C1 BC'
		#if paras[i]=='BCAE':
		#	paras[i] = 'BC absorbing enhancement'
		if paras[i]=='amb':
		#	paras[i] = 'Atmospheric parameters'
			paras[i] = 'AEP'
	
	# set plot parameters
	
	fs = 9 # fontsize
	fw = 'bold' # fontweight
	width = 0.5 # bar width
	lw = 2 # line width
	rotation = 45 # x labels rotaion angle
	#colors = ['k','k','k','r','r','r','r','b','b','b','g','g','g','g','g']
	colors = ['k','k','k','r','r','r','r','b','b','b','y','y','y','g','g','g','g','g']
	
	plt.figure(figsize=(7,out_num*2))
	plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.17, wspace=0.1, hspace=0.3)
	
	# start plotting
	
	for i in range(out_num):
		
		# plot n shell for its high sensitivity
		ax = plt.subplot2grid((out_num*3,1),(i*3,0),colspan=1,rowspan=1)
		for j in range(par_num):
			plt.bar(j+1, CP[i,j], width=width, color=colors[j])
			#plt.text(j+1-width*0.7, CP[i,j], str(round(CP[i,j]*100)/100), fontsize=fs*0.8, fontweight=fw, color='gray')
		kwargs = dict(marker=[(-1,-0.85),(1,0.85)],markersize=5,linestyle='none',color='k',mec='k',mew=1,clip_on=False)
		ax.plot([0,1],[0,0],transform=ax.transAxes,**kwargs)
		ax.spines['bottom'].set_visible(False)
		
		plt.ylim(round(np.max(CP[i])*1.1-np.max(CP[i,1:])*1.1*0.5), round(np.max(CP[i])*1.1))
		plt.yticks([round(np.max(CP[i])*1.1-np.max(CP[i,1:])*1.1*0.5), round(np.max(CP[i])*1.1)],fontsize=fs, fontweight=fw)
		plt.xticks([])
		
		# plot other parameters
		ax = plt.subplot2grid((out_num*3,1),(i*3+1,0),colspan=1,rowspan=2)
		
		# main plot function
		for j in range(par_num):
			plt.bar(j+1, CP[i,j], width=width, color=colors[j])
			#plt.text(j+1-width*0.7, CP[i,j], str(round(CP[i,j]*100)/100), fontsize=fs*0.8, fontweight=fw, color='gray')
		ax.plot([0,1],[1,1],transform=ax.transAxes,**kwargs)
		ax.spines['top'].set_visible(False)
		
		# y lim
		plt.ylim(0,round(np.max(CP[i,1:])*1.1,3))
		
		# y label
		if outputs[i]=='RF_top':
			plt.ylabel('           top RF, dRF/dx%', labelpad=10+7, fontsize=fs, fontweight=fw)
		elif outputs[i]=='RF_bot':
			plt.ylabel('           bottom RF, dRF/dx%', labelpad=10, fontsize=fs, fontweight=fw)
		elif outputs[i]=='dtau_bot':
			plt.ylabel('           dAOD/dx%', labelpad=10, fontsize=fs, fontweight=fw)
		elif outputs[i]=='waer_bot':
			plt.ylabel('           dSSA/dx%', labelpad=10, fontsize=fs, fontweight=fw)
		else:
			plt.ylabel('           dg/dx%', labelpad=10, fontsize=fs, fontweight=fw)
		
		# y ticks
		plt.yticks([0,round(np.max(CP[i,1:])*1.1)/3,round(np.max(CP[i,1:])*1.1)/3*2,round(np.max(CP[i,1:])*1.1)], ['0',str(round(np.max(CP[i,1:])*1.1/3)),str(round(np.max(CP[i,1:])*1.1/3*2)),str(round(np.max(CP[i,1:])*1.1))], fontsize=fs, fontweight=fw)
		
		# x label
		if i==out_num-1:
			ax=plt.gca()
			plt.xticks(np.arange(par_num)+1, paras, fontsize=fs, fontweight=fw, rotation=rotation, horizontalalignment='right')
			for xtick, color in zip(ax.get_xticklabels(), colors):
				xtick.set_color(color)
		else:
			plt.xticks([])
	
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_SA.pdf')
	else:
		plt.show()

def plot6(dtime, par_range, **args):
	'''
	This function is to plot CP sensitivity results only
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
	
	CP = calSA.cal_CP(dtime, par_range, outputs=outputs, paras=paras)['CP']
	
	out_num = len(outputs)
	par_num = len(paras)
	
	# sort data
	'''
	weight = CP[0] / sum(CP[0]) + CP[1] / sum(CP[1]) + CP[2] / sum(CP[2])
	combined = zip(weight, CP[0], CP[1], CP[2], paras)
	sorted_combined = sorted(combined, key=lambda x: x[0], reverse=True)
	weight, CP[0], CP[1], CP[2], paras = zip(*sorted_combined)
	paras = np.array(paras) # from chatGPT
	'''
	# change parameter names
	
	for i in range(par_num):
		if paras[i]=='n':
			paras[i] = 'n shell'
		if paras[i]=='nH1':
			paras[i] = 'n H1'
		if paras[i]=='nBC':
			paras[i] = 'n BC'
		if paras[i]=='kBC':
			paras[i] = 'k BC'
		#if paras[i]=='MS':
		#	paras[i] = 'Mixing state'
		if paras[i]=='MSH1':
			paras[i] = 'MS H1'
		#if paras[i]=='VD':
		#	paras[i] = 'Vertical distribution'
		#if paras[i]=='CT':
		#	paras[i] = 'Coating thickness'
		if paras[i]=='CTH1':
			paras[i] = 'CT H1'
		if paras[i]=='kappaH1':
			paras[i] = 'kappa H1'
		if paras[i]=='rhoBC':
			paras[i] = '\u03C1 BC'
		#if paras[i]=='BCAE':
		#	paras[i] = 'BC absorbing enhancement'
		if paras[i]=='amb':
		#	paras[i] = 'Atmospheric parameters'
			paras[i] = 'AEP'
	
	# set plot parameters
	
	fs = 9 # fontsize
	fw = 'bold' # fontweight
	width = 0.5 # bar width
	lw = 2 # line width
	rotation = 45 # x labels rotaion angle
	colors = ['k','k','k','r','r','r','r','b','b','b','g','g','g','g','g']
	
	plt.figure(figsize=(7,out_num*2))
	plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.12, wspace=0.1, hspace=0.3)
	
	# start plotting
	
	for i in range(out_num):
		
		# plot n shell for its high sensitivity
		ax = plt.subplot2grid((out_num*3,1),(i*3,0),colspan=1,rowspan=1)
		for j in range(par_num):
			plt.bar(j+1, CP[i,j], width=width, color=colors[j])
		kwargs = dict(marker=[(-1,-0.85),(1,0.85)],markersize=5,linestyle='none',color='k',mec='k',mew=1,clip_on=False)
		ax.plot([0,1],[0,0],transform=ax.transAxes,**kwargs)
		ax.spines['bottom'].set_visible(False)
		
		plt.ylim(round(np.max(CP[i])*1.1-np.max(CP[i,1:])*1.1*0.5,2), round(np.max(CP[i])*1.1,2))
		plt.yticks([round(np.max(CP[i])*1.1-np.max(CP[i,1:])*1.1*0.5,2), round(np.max(CP[i])*1.1,2)],fontsize=fs, fontweight=fw)
		plt.xticks([])
		
		# plot other parameters
		ax = plt.subplot2grid((out_num*3,1),(i*3+1,0),colspan=1,rowspan=2)
		
		# main plot function
		for j in range(par_num):
			plt.bar(j+1, CP[i,j], width=width, color=colors[j])
		ax.plot([0,1],[1,1],transform=ax.transAxes,**kwargs)
		ax.spines['top'].set_visible(False)
		
		# y lim
		plt.ylim(0,round(np.max(CP[i,1:])*1.1,3))
		
		# y label
		if outputs[i]=='RF_top':
			plt.ylabel('           top RF, dRF/dx%', labelpad=10+7, fontsize=fs, fontweight=fw)
		elif outputs[i]=='RF_bot':
			plt.ylabel('           bottom RF, dRF/dx%', labelpad=10, fontsize=fs, fontweight=fw)
		elif outputs[i]=='dtau_bot':
			plt.ylabel('           dAOD/dx%', labelpad=10, fontsize=fs, fontweight=fw)
		elif outputs[i]=='waer_bot':
			plt.ylabel('           dSSA/dx%', labelpad=10, fontsize=fs, fontweight=fw)
		else:
			plt.ylabel('           dg/dx%', labelpad=10, fontsize=fs, fontweight=fw)
		
		# y ticks
		plt.yticks([0,round(np.max(CP[i,1:])*1.1,3)/3,round(np.max(CP[i,1:])*1.1,3)/3*2,round(np.max(CP[i,1:])*1.1,3)], ['0',str(round(np.max(CP[i,1:])*1.1/3,2)),str(round(np.max(CP[i,1:])*1.1/3*2,2)),str(round(np.max(CP[i,1:])*1.1,2))], fontsize=fs, fontweight=fw)
		
		# x label
		if i==out_num-1:
			ax=plt.gca()
			plt.xticks(np.arange(par_num)+1, paras, fontsize=fs, fontweight=fw, rotation=rotation, horizontalalignment='right')
			for xtick, color in zip(ax.get_xticklabels(), colors):
				xtick.set_color(color)
		else:
			plt.xticks([])
	
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_SA.pdf')
	else:
		plt.show()

def plot7(dtime, par_range, **args):
	'''
	This function is to plot CP sensitivity and non-linear proportion
	input:
		dtime		: runDARF.py output time, string, example: '230101'
		par_range	: parameter change rate range, float
		**paras		: parameters list, array, string, default default_paras
		**outputs	: outputs list, array, string, default default_outputs
		**save		: save flag, bool, default False
		**save_path	: save path, string, default 'figure/DARF/
	output:
		x axis for sensitivity
		y axis for factor list, different color for different factor group
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
	
	CP = calSA.cal_CP(dtime, par_range, outputs=outputs, paras=paras)['CP']
	OAT = calSA.read_OAT(dtime, outputs=outputs, paras=paras)['OAT']
	vs = abs((OAT-CP)/CP)
	#vs[np.isnan(vs)] = 0
	#vs[np.isinf(vs)] = 0
	#print(np.sum(vs*CP/np.sum(CP)))
	vs[np.isnan(vs)] = -1 # turn np.nan to -1
	vs[np.isinf(vs)] = -1 # turn np.inf to -1
	
	out_num = len(outputs)
	par_num = len(paras)
	
	# sort data
	'''
	weight = CP[0] / sum(CP[0]) + CP[1] / sum(CP[1]) + CP[2] / sum(CP[2])
	combined = zip(weight, CP[0], CP[1], CP[2], paras)
	sorted_combined = sorted(combined, key=lambda x: x[0], reverse=True)
	weight, CP[0], CP[1], CP[2], paras = zip(*sorted_combined)
	paras = np.array(paras) # from chatGPT
	'''
	# change parameter names
	
	for i in range(par_num):
		if paras[i]=='n':
			paras[i] = 'CRI$_{shell,dry}$'
		if paras[i]=='nH1':
			paras[i] = 'CRI$_{shell,dry}$ H1'
		if paras[i]=='PNSD':
			paras[i] = 'PNSD$_{dry}$'
		if paras[i]=='nBC':
			paras[i] = 'n$_{LAC}$'
		if paras[i]=='kBC':
			paras[i] = 'k$_{LAC}$'
		if paras[i]=='BCPNSD':
			paras[i] = 'LACPNSD'
		if paras[i]=='MS':
			paras[i] = 'Mixing state'
		if paras[i]=='MSH1':
			paras[i] = 'Mixing state H1'
		if paras[i]=='MSH2':
			paras[i] = 'Mixing state H2'
		if paras[i]=='VD':
			paras[i] = 'Vertical distribution'
		if paras[i]=='CT':
			paras[i] = 'Coating thickness'
		if paras[i]=='CTH1':
			paras[i] = 'Coating thickness H1'
		if paras[i]=='kappaH1':
			paras[i] = 'kappa H1'
		if paras[i]=='rhoBC':
			paras[i] = '\u03C1$_{LAC}$'
		if paras[i]=='BCAE':
			paras[i] = 'LAC absorbing enhancement'
		if paras[i]=='amb':
			paras[i] = 'RH'
	
	# set plot parameters
	
	figsize=(out_num*4, par_num*0.35)
	fs = 9 # fontsize
	fw = 'bold' # fontweight
	height = 0.5 # bar height
	lw = 2 # line width
	rotation = 0 # x labels rotaion angle
	edgecolor = 'k' # bar edge color
	linewidth = 1.2 # bar edge line width
	kwargs = dict(marker=[(-1,-1),(1,1)],markersize=8,linestyle='none',color='k',mec='k',mew=1,clip_on=False)
	bbox_props = dict(boxstyle='square', fc='w', lw=0.5)
	#colors = ['k','k','k','r','r','r','r','b','b','b','y','y','y','g','g','g','g','g']
	colors = ['k','k','k','r','r','r','r','b','b','b','gold','gold','gold']
	
	# start plotting
	
	plt.figure(figsize=figsize)
	plt.subplots_adjust(left=0.3, right=0.95, top=0.95, bottom=0.15, wspace=0.3, hspace=0.1)
	
	# start plotting
	
	for i in range(out_num):
		
		# plot n shell for its high sensitivity
		ax = plt.subplot2grid((1,out_num*3),(0,i*3+2),colspan=1,rowspan=1)
		for j in range(par_num):
			plt.barh(j+1, CP[i,-j-1], height=height, color=colors[-j-1], edgecolor=edgecolor, linewidth=linewidth)
			plt.barh(j+1, 0, left=CP[i,-j-1], height=height*1.5, edgecolor=edgecolor, linewidth=linewidth*1.5)
			plt.barh(j+1, OAT[i,-j-1]-CP[i,-j-1], left=CP[i,-j-1], height=height, color='w', edgecolor=edgecolor, linewidth=linewidth)
			#plt.text(j+1-height*0.7, CP[i,j], str(round(CP[i,j]*100)/100), fontsize=fs*0.8, fontweight=fw, color='gray')
		ax.plot([0,0],[1,0],transform=ax.transAxes,**kwargs)
		ax.spines['left'].set_visible(False)
		
		max_num_1 = max(np.max(CP[i]),np.max(OAT[i]))
		max_num_2 = max(np.max(CP[i,1:]),np.max(OAT[i,1:]))
		
		xlim = (round(max_num_1*1.1-max_num_2*1.1*0.5,2), round(max_num_1*1.1,2))
		plt.xlim(xlim)
		plt.xticks([xlim[0],xlim[1]], [round(xlim[0]),round(xlim[1])], fontsize=fs, fontweight=fw)
		plt.yticks([])
		
		# plot other parameters
		ax = plt.subplot2grid((1,out_num*3),(0,i*3),colspan=2,rowspan=1)
		
		# main plot function
		for j in range(par_num):
			plt.barh(j+1, CP[i,-j-1], height=height, color=colors[-j-1], edgecolor=edgecolor, linewidth=linewidth)
			plt.barh(j+1, 0, left=CP[i,-j-1], height=height*1.5, edgecolor=edgecolor, linewidth=linewidth*1.5)
			plt.barh(j+1, OAT[i,-j-1]-CP[i,-j-1], left=CP[i,-j-1], height=height, color='w', edgecolor=edgecolor, linewidth=linewidth)
			#plt.text(max(CP[i,-j-1],OAT[i,-j-1])+height*0.7, j+1-height*0.3, str(round(CP[i,j]*100)/100), fontsize=fs*0.8, fontweight=fw, color='gray')
			if i==out_num-1:
				plt.text(0.7, 0.05, 'color bar: CP\nwhite bar: OAT-CP', fontsize=fs, fontweight=fw, bbox=bbox_props, transform=ax.transAxes,)
		ax.plot([1,1],[1,0],transform=ax.transAxes,**kwargs)
		ax.spines['right'].set_visible(False)
		
		# x label
		if outputs[i]=='RF_top':
			plt.xlabel('                    RF$_{ari,top}$ sensitivity, dRF/dx%', labelpad=10, fontsize=fs*1.2, fontweight=fw)
		elif outputs[i]=='RF_bot':
			plt.xlabel('                    RF$_{ari,bottom}$ sensitivity, dRF/dx%', labelpad=10, fontsize=fs*1.2, fontweight=fw)
		elif outputs[i]=='dtau_bot':
			plt.xlabel('                    AOD sensitivity, dAOD/dx%', labelpad=10, fontsize=fs*1.2, fontweight=fw)
		elif outputs[i]=='waer_bot':
			plt.xlabel('                    SSA sensitivity, dSSA/dx%', labelpad=10, fontsize=fs*1.2, fontweight=fw)
		else:
			plt.xlabel('                    g sensitivity, dg/dx%', labelpad=10, fontsize=fs*1.2, fontweight=fw)
		
		# x lim
		plt.xlim(0,round(max_num_2*1.1,3))
		
		# x ticks
		
		plt.xticks([0,max_num_2*1.1/3,max_num_2*1.1/3*2,max_num_2*1.1], 
		['0',str(round(max_num_2*1.1/3)),str(round(max_num_2*1.1/3*2)),str(round(max_num_2*1.1))], 
		fontsize=fs, fontweight=fw)
		
		# y label
		if i==0:
			ax=plt.gca()
			plt.yticks(np.arange(par_num)+1, paras[::-1], fontsize=fs, fontweight=fw, rotation=rotation, horizontalalignment='right')
			for ytick, color in zip(ax.get_yticklabels(), colors[::-1]):
				ytick.set_color(color)
		else:
			plt.yticks([])
	
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_SA.pdf')
	else:
		plt.show()

def plot8(dtime, par_range, **args):
	'''
	This function is to plot CP sensitivity and non-linear proportion
	input:
		dtime		: runDARF.py output time, string, example: '230101'
		par_range	: parameter change rate range, float
		**paras		: parameters list, array, string, default default_paras
		**outputs	: outputs list, array, string, default default_outputs
		**save		: save flag, bool, default False
		**save_path	: save path, string, default 'figure/DARF/
	output:
		x axis for sensitivity
		y axis for factor list, different color for different factor group
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
	
	CP = calSA.cal_CP(dtime, par_range, outputs=outputs, paras=paras)['CP']
	OAT = calSA.read_OAT(dtime, outputs=outputs, paras=paras)['OAT']
	vs = abs((OAT-CP)/CP)
	vs[np.isnan(vs)] = -1 # turn np.nan to -1
	vs[np.isinf(vs)] = -1 # turn np.inf to -1
	
	out_num = len(outputs)
	par_num = len(paras)
	
	# sort data
	'''
	weight = CP[0] / sum(CP[0]) + CP[1] / sum(CP[1]) + CP[2] / sum(CP[2])
	combined = zip(weight, CP[0], CP[1], CP[2], paras)
	sorted_combined = sorted(combined, key=lambda x: x[0], reverse=True)
	weight, CP[0], CP[1], CP[2], paras = zip(*sorted_combined)
	paras = np.array(paras) # from chatGPT
	'''
	# change parameter names
	
	for i in range(par_num):
		if paras[i]=='n':
			paras[i] = 'CRI$_{shell,dry}$'
		if paras[i]=='nH1':
			paras[i] = 'CRI$_{shell,dry}$ H1'
		if paras[i]=='nBC':
			paras[i] = 'n$_{LAC}$'
		if paras[i]=='kBC':
			paras[i] = 'k$_{LAC}$'
		if paras[i]=='BCPNSD':
			paras[i] = 'LACPNSD'
		if paras[i]=='MS':
			paras[i] = 'Mixing state'
		if paras[i]=='MSH1':
			paras[i] = 'Mixing state H1'
		if paras[i]=='MSH2':
			paras[i] = 'Mixing state H2'
		if paras[i]=='VD':
			paras[i] = 'Vertical profile'
		if paras[i]=='CT':
			paras[i] = 'Coating thickness'
		if paras[i]=='CTH1':
			paras[i] = 'Coating thickness H1'
		if paras[i]=='kappaH1':
			paras[i] = 'kappa H1'
		if paras[i]=='rhoBC':
			paras[i] = '\u03C1$_{LAC}$'
		if paras[i]=='BCAE':
			paras[i] = 'LAC absorbing enhancement'
		if paras[i]=='amb':
			paras[i] = 'RH'
	
	# set plot parameters
	
	figsize=(out_num*4,par_num*0.35)
	fs = 9 # fontsize
	fw = 'bold' # fontweight
	height = 0.5 # bar height
	lw = 2 # line width
	rotation = 0 # x labels rotaion angle
	edgecolor = 'k' # bar edge color
	linewidth = 1.2 # bar edge line width
	kwargs = dict(marker=[(-1,-1),(1,1)],markersize=8,linestyle='none',color='k',mec='k',mew=1,clip_on=False)
	bbox_props = dict(boxstyle='square', fc='w', lw=0.5)
	#colors = ['k','k','k','r','r','r','r','b','b','b','g','g','g','g','g']
	colors = ['k','k','k','r','r','r','r','b','b','b']
	
	# start plotting
	
	plt.figure(figsize=figsize)
	plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.2, wspace=0.3, hspace=0.1)
	
	# start plotting
	
	for i in range(out_num):
		
		# plot n shell for its high sensitivity
		ax = plt.subplot2grid((1,out_num*3),(0,i*3+2),colspan=1,rowspan=1)
		for j in range(par_num):
			plt.barh(j+1, CP[i,-j-1], height=height, color=colors[-j-1], edgecolor=edgecolor, linewidth=linewidth)
			plt.barh(j+1, 0, left=CP[i,-j-1], height=height*1.5, edgecolor=edgecolor, linewidth=linewidth*1.5)
			plt.barh(j+1, OAT[i,-j-1]-CP[i,-j-1], left=CP[i,-j-1], height=height, color='w', edgecolor=edgecolor, linewidth=linewidth)
			#plt.text(j+1-height*0.7, CP[i,j], str(round(CP[i,j]*100)/100), fontsize=fs*0.8, fontweight=fw, color='gray')
		ax.plot([0,0],[1,0],transform=ax.transAxes,**kwargs)
		ax.spines['left'].set_visible(False)
		
		max_num_1 = max(np.max(CP[i]),np.max(OAT[i]))
		max_num_2 = max(np.max(CP[i,1:]),np.max(OAT[i,1:]))
		
		xlim = (round(max_num_1*1.1-max_num_2*1.1*0.5,2), round(max_num_1*1.1,2))
		plt.xlim(xlim)
		plt.xticks([xlim[0],xlim[1]], [round(xlim[0]*100),round(xlim[1]*100)], fontsize=fs, fontweight=fw)
		plt.yticks([])
		
		# plot other parameters
		ax = plt.subplot2grid((1,out_num*3),(0,i*3),colspan=2,rowspan=1)
		
		# main plot function
		for j in range(par_num):
			plt.barh(j+1, CP[i,-j-1], height=height, color=colors[-j-1], edgecolor=edgecolor, linewidth=linewidth)
			plt.barh(j+1, 0, left=CP[i,-j-1], height=height*1.5, edgecolor=edgecolor, linewidth=linewidth*1.5)
			plt.barh(j+1, OAT[i,-j-1]-CP[i,-j-1], left=CP[i,-j-1], height=height, color='w', edgecolor=edgecolor, linewidth=linewidth)
			#plt.text(j+1-height*0.7, CP[i,j], str(round(CP[i,j]*100)/100), fontsize=fs*0.8, fontweight=fw, color='gray')
			if i==0:
				plt.text(0.8, 0.07, 'color bar: CP\nwhite bar: OAT-CP', fontsize=fs, fontweight=fw, bbox=bbox_props, transform=ax.transAxes,)
		ax.plot([1,1],[1,0],transform=ax.transAxes,**kwargs)
		ax.spines['right'].set_visible(False)
		
		# x label
		if outputs[i]=='RF_top':
			plt.xlabel('                    RF$_{ari,top}$ sensitivity, dRF/dx%', labelpad=10, fontsize=fs*1.2, fontweight=fw)
		elif outputs[i]=='RF_bot':
			plt.xlabel('                    RF$_{ari,bottom}$ sensitivity, dRF/dx%', labelpad=10, fontsize=fs*1.2, fontweight=fw)
		elif outputs[i]=='dtau_bot':
			plt.xlabel('                    AOD sensitivity, $10^{-2}$dAOD/dx%', labelpad=10, fontsize=fs*1.2, fontweight=fw)
		elif outputs[i]=='waer_bot':
			plt.xlabel('                    SSA sensitivity, $10^{-2}$dSSA/dx%', labelpad=10, fontsize=fs*1.2, fontweight=fw)
		else:
			plt.xlabel('                    g sensitivity, $10^{-2}$dg/dx%', labelpad=10, fontsize=fs*1.2, fontweight=fw)
		
		# x lim
		plt.xlim(0,round(max_num_2*1.1,3))
		
		# x ticks
		
		plt.xticks([0,max_num_2*1.1/3,max_num_2*1.1/3*2,max_num_2*1.1], 
		['0',round(max_num_2*1.1/3*100),round(max_num_2*1.1/3*2*100),round(max_num_2*1.1*100)], 
		fontsize=fs, fontweight=fw)
		
		# y label
		if i==0:
			ax=plt.gca()
			plt.yticks(np.arange(par_num)+1, paras[::-1], fontsize=fs, fontweight=fw, rotation=rotation, horizontalalignment='right')
			for ytick, color in zip(ax.get_yticklabels(), colors[::-1]):
				ytick.set_color(color)
		else:
			plt.yticks([])
	
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_SA.pdf')
	else:
		plt.show()

if __name__ == '__main__':
	paras = [] # default parameters list
	paras.append('n')
	paras.append('PNSD')
	paras.append('kappa')
	
	paras.append('nBC')
	paras.append('kBC')
	paras.append('rhoBC')
	paras.append('BCPNSD')
	
	paras.append('MS')
	paras.append('CT')
	paras.append('BCAE')
	
	paras.append('VD')
	paras.append('amb')
	paras.append('albedo')
	'''
	paras.append('nH1')
	paras.append('kappaH1')
	paras.append('CTH1')
	paras.append('MSH1')
	paras.append('MSH2')
	'''
	AOP_paras = []
	AOP_paras.append('n')
	AOP_paras.append('PNSD')
	AOP_paras.append('kappa')
	
	AOP_paras.append('nBC')
	AOP_paras.append('kBC')
	AOP_paras.append('rhoBC')
	AOP_paras.append('BCPNSD')
	
	AOP_paras.append('MS')
	AOP_paras.append('CT')
	AOP_paras.append('BCAE')
	'''
	AOP_paras.append('nH1')
	AOP_paras.append('kappaH1')
	AOP_paras.append('CTH1')
	AOP_paras.append('MSH1')
	AOP_paras.append('MSH2')
	'''
	
	#plot('240610', 0.05, save=False, outputs=['RF_top','RF_bot'], paras=['n','nH1','kappa','kappaH1','CT','CTH1','MS','MSH1','MSH2'])
	#plot('240610', 0.05, save=False, outputs=['RF_top','RF_bot'], paras=paras)
	#plot('240610', 0.05, save=False, outputs=['dtau_bot','waer_bot','g_bot'], paras=paras)
	#plot2('240610', 0.05, save=False, outputs=['dtau','waer','g'])
	#plot_all('240610', 0.05, save=False, outputs=['RF_top','RF_bot'], paras=paras)
	#plot3('240610', 0.05, save=False, outputs=['dtau','waer','g'], paras=['n','PNSD','kappa','MS','CT'])
	#plot4('240610', 0.05, save=False, outputs=['RF_top','RF_bot'], paras=paras)
	#plot5('240610', 0.05, save=False, outputs=['RF_top','RF_bot'], paras=paras)
	#plot6('240610', 0.05, save=False, outputs=['dtau_bot','waer_bot','g_bot'], paras=paras)
	plot7('241010', 0.05, save=False, outputs=['RF_top','RF_bot'], paras=paras)
	plot8('241010', 0.05, save=False, outputs=['dtau_bot','waer_bot','g_bot'], paras=AOP_paras)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
