####################################################################################
# INTRODUCTION:
# This code is to calculationg correlation coefficient
# Created by Hebs at 22/10/4/16:21
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
from scipy.stats import pearsonr, rankdata
from numpy.linalg import lstsq

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

def pearson_corr(paras, y):
	'''
	pearson correlation coefficient calculating function
	input:
		paras : parameters matrix, array, in shape (par_num, run_num)
		y     : output list, array
	output:
		CC    : Spearman correlation coefficient, array
		CC_P  : Spearman correlation coefficient p-value, array
	'''
	CC = np.zeros(len(paras))
	CC_P = np.zeros(len(paras))
	
	for i in range(len(CC)):
		if np.mean(paras[i])==0:
			CC[i] = 0
			CC_P[i] = 99999 # mask value
		else:
			CC[i] = pearsonr(paras[i],y)[0]
			CC_P[i] = pearsonr(paras[i],y)[1]
	
	return CC, CC_P

def partial_corr(paras, y):
	'''
	partial correlation coefficient calculating function
	input:
		paras : parameters matrix, array, in shape (par_num, run_num)
		y     : output list, array
	output:
		PCC   : Spearman correlation coefficient, array
		PCC_P : Spearman correlation coefficient p-value, array
	'''
	C = np.vstack((paras,y)).T
	PCC = np.zeros(len(paras))
	PCC_P = np.zeros(len(paras))
	
	for i in range(len(paras)):
		if np.mean(paras[i])==0:
			PCC[i] = 0
			PCC_P[i] = 99999 # mask value
		else:
			idx = np.ones(len(paras)+1, dtype=bool)
			idx[i] = False
			idx[-1] = False
			beta_i = lstsq(C[:,idx], C[:,i])[0]
			beta_y = lstsq(C[:,idx], C[:,-1])[0]
			res_i = C[:,i] - C[:,idx].dot(beta_i)
			res_y = C[:,-1] - C[:,idx].dot(beta_y)
			PCC[i] = pearsonr(res_i, res_y)[0]
			PCC_P[i] = pearsonr(res_i, res_y)[1]
	
	return PCC, PCC_P

def cal_residual(paras, y, i):
	'''
	This function is to calculate residual between parameter and output
	input:
		paras  : parameters matrix, array, in shape (par_num, run_num)
		y      : output list, array
		i      : ith parameter, int
	output:
		res_i  : ith parameter residual
		res_y  : output residual
	'''
	C = paras.T
	PCC = np.zeros(len(paras))
	PCC_P = np.zeros(len(paras))
	
	idx = np.ones(len(paras), dtype=np.bool)
	idx[i] = False
	beta_i = lstsq(C[:,idx], C[:,i])[0]
	beta_y = lstsq(C[:,idx], y)[0]
	res_i = C[:,i] - C[:,idx].dot(beta_i)
	res_y = y - C[:,idx].dot(beta_y)
	
	return res_i, res_y

def cal_correlation(paras, y, **args):
	'''
	This function is to calculate Pearson correlation coefficient between two parameters
	input:
		paras  : parameters matrix, array, in shape (par_num, run_num)
		y      : output list, array
		**debug: output flag, bool, default False
	output:
		CC     : Pearson correlation coefficient, float, array
		RCC    : rank correlation coefficient, float, array
		PCC    : partial correlation coefficient, float, array
		PRCC   : partial rank correlation coefficient, float, array
		those coefficient all have P value, in name XX_P
	'''
	if 'debug' in args:
		debug = args['debug']
	else:
		debug = False
	
	if debug:
		print('calculating Pearson correlation coefficient...')
	CC, CC_P = pearson_corr(paras, y)
	
	if debug:
		print('calculating partial correlation coefficient...')
	PCC, PCC_P = partial_corr(paras, y)
	
	rank_paras = np.zeros(paras.shape)
	for i in range(len(paras)):
		if np.mean(paras[i])==0:
			rank_paras[i] = np.zeros(len(paras[i]))
		else:
			rank_paras[i] = rankdata(paras[i])
	rank_y = rankdata(y)
	
	if debug:
		print('calculating ranking correlation coefficient...')
	RCC, RCC_P = pearson_corr(rank_paras, rank_y)
	
	if debug:
		print('calculating partial ranking correlation coefficient...')
	PRCC, PRCC_P = partial_corr(rank_paras, rank_y)
	
	return CC, CC_P, PCC, PCC_P, RCC, RCC_P, PRCC, PRCC_P

def cal(dtime, **args):
	'''
	input:
		dtime       : runMie.py output time, string, example: '230101'
		**save      : save flag, bool, default False
		**save_path : save path, string, default 'output/CC/'
		**debug     : output flag, bool, default False
		**paras     : parameters list, array, string, default default_paras
	output:
		in dictionary
		CC          : Pearson correlation coefficient, float, array
		RCC         : rank correlation coefficient, float, array
		PCC         : partial correlation coefficient, float, array
		PRCC        : partial rank correlation coefficient, float, array
		those coefficient all have P value
	'''
	if 'path' in args:
		path = args['path']
	else:
		path = 'output/Mie/'
	if 'save' in args:
		save = args['save']
	else:
		save = False
	if 'save_path' in args:
		save_path = args['save_path']
	else:
		save_path = 'output/CC/'
	if 'debug' in args:
		debug = args['debug']
	else:
		debug = False
	if 'paras' in args:
		paras = args['paras']
	else:
		paras = default_paras
	
	par_num = len(paras)
	outputs = ['AOD', 'SSA', 'g']
	out_num = len(outputs)
	path = path + dtime + '/all/'
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
	
	CC = np.zeros((out_num,par_num))
	CC_P = np.zeros((out_num,par_num))
	RCC = np.zeros((out_num,par_num))
	RCC_P = np.zeros((out_num,par_num))
	PCC = np.zeros((out_num,par_num))
	PCC_P = np.zeros((out_num,par_num))
	PRCC = np.zeros((out_num,par_num))
	PRCC_P = np.zeros((out_num,par_num))
	
	for i in range(out_num):
		CC[i], CC_P[i], PCC[i], PCC_P[i], RCC[i], RCC_P[i], PRCC[i], PRCC_P[i] = cal_correlation(all_paras, outs[i], debug=debug)
	
	if save:
		info = dict(CC=CC, CC_P=CC_P, PCC=PCC, PCC_P=PCC_P, RCC=RCC, RCC_P=RCC_P, PRCC=PRCC, PRCC_P=PRCC_P, paras=paras, outputs=outputs)
		np.save(save_path + 'infos.npy', info)
	
	Corr = dict(CC=CC, CC_P=CC_P, PCC=PCC, PCC_P=PCC_P, RCC=RCC, RCC_P=RCC_P, PRCC=PRCC, PRCC_P=PRCC_P)
	return Corr

if __name__ == '__main__':
	Corr = cal('230604', debug=False, paras=['n','VD'])
	print(Corr['CC'][0])
	print(Corr['CC_P'][0])
	print(Corr['PCC'][0])
	print(Corr['PCC_P'][0])
	print(Corr['RCC'][0])
	print(Corr['RCC_P'][0])
	print(Corr['PRCC'][0])
	print(Corr['PRCC_P'][0])
	
	'''
	idx = np.ones(len(default_paras)) * (len(default_paras)+1)
	for i in range(len(PCC[2])):
		print(round(PCC[2,i],3),end=' ')
	print('')
	print(idx-rankdata(abs(PCC[2])))
	print(default_paras)
	
	for i in range(len(PRCC[2])):
		print(round(PRCC[2,i],3),end=' ')
	print('')
	print(idx-rankdata(abs(PRCC[2])))
	print(default_paras)
	'''
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
