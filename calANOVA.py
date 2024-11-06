####################################################################################
# INTRODUCTION:
# This code is to do ANOVA and cal F-test result
# Created by Hebs at 23/5/8/9:08
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
default_paras.append('BCPMSD')
default_paras.append('BCI')
default_paras.append('BCAE')
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

import numpy as np

def cal_ANOVA(paras, y, **args):
	'''
	This function is to calculate ANOVA and F-test
	input:
		paras	: parameters matrix, array, in shape (par_num, run_num)
		y	: output list, array
		**K	: group number, int, default 10
	output:
		ANOVA	: result of analysis of variance's F-test, array, in shape (par_num)
	'''
	if 'K' in args:
		K = args['K']
	else:
		K = 10
	
	ANOVA = np.zeros(len(paras))
	
	for i in range(len(paras)):
		if np.mean(paras[i])==0:
			ANOVA[i] = 0
			continue
		# to do ANOVA, 1st: divide data in groups; 2rd: do F-test
		boundaries = np.percentile(paras[i], np.linspace(0,100,K+1))
		groups = np.digitize(paras[i], boundaries)
		grouped_data = [y[groups == i] for i in range(1,11)] # by ChatGPT
		########################################
		#print(grouped_data)
		#Y = np.array(grouped_data)
		# for different python version
		########################################
		Y = grouped_data
		
		# now calculate F value
		F1 = 0
		F2 = 0
		for j in range(K):
			F1 = F1 + len(Y[j]) * (np.mean(Y[j])-np.mean(y))**2 / (K-1)
			for k in range(len(Y[j])):
				F2 = F2 + (Y[j][k]-np.mean(Y[j]))**2 / (len(y)-K)
		
		ANOVA[i] = F1 / F2
	
	return ANOVA

def cal(dtime, **args):
	'''
	input:
		dtime		: runMie.py output time, string, example: '230101'
		**path		: data path, string, defualt 'output/Mie/'
		**K		: group number, int, default 10
		**paras	: parameters list, array, string, default default_paras
	output:
		ANOVA		: result of analysis of variance's F-test, array
	'''
	if 'path' in args:
		path = args['path']
	else:
		path = 'output/Mie/'
	if 'K' in args:
		K = args['K']
	else:
		K = 10
	if 'paras' in args:
		paras = args['paras']
	else:
		paras = default_paras
	if 'outputs' in args:
		outputs = args['outputs']
	else:
		outputs = ['AOD', 'SSA', 'g']
	
	path = path + dtime + '/all/'
	par_num = len(paras)
	out_num = len(outputs)
	
	data = np.load(path+'infos.npy', allow_pickle=True)
	run_num = len(data)
	
	all_paras = np.zeros((par_num,run_num)) # all parameters
	outs = np.zeros((out_num,run_num))
	
	for i in range(par_num):
		for j in range(run_num):
			all_paras[i,j] = data[j][paras[i]+'_rate']
	
	for i in range(out_num):
		for j in range(run_num):
			outs[i,j] = data[j][outputs[i]]
	
	ANOVA = np.zeros((out_num,par_num))
	
	for i in range(out_num):
		ANOVA[i] = cal_ANOVA(all_paras, outs[i], K=K)
	
	return ANOVA

if __name__ == '__main__':
	#ANOVA = cal('230530', K=10, paras=['n','nI','rhoBC'])
	ANOVA = cal('240310', path='output/DARF/', K=10, outputs=['RF_top','RF_bot'])
	print(ANOVA)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
