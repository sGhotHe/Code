####################################################################################
# INTRODUCTION:
# This code is to calculate co-parameters' Constraint Factor 
# Sensitivity Analysis's result
# Created by Hebs at 23/9/16/10:39
# Contact: hebishuo@pku.edu.cn
####################################################################################
# warning
# this function may have serious math problem
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
default_paras.append('co')

import numpy as np

def cal(dtime, **args):
	'''
	This function is to calculate Constraint Factor Sensitivity Analysis use runMie.py's result
	input:
		dtime		: runMie.py output time, string, example: '230101'
		**path		: runMie.py output save path, string, default 'output/Mie/'
		**paras		: parameters list, array, string, default default_paras
	output:
		in dictionary
		paras		: parameters list, array, string
		run_num	: running number, int
		Cons		: Constraint Factor Sensitivity Analysis result, array, float,
					in shape (3, paras_num), 3 for AOD, SSA and g respectively
		Stds		: Constraint Factor standard deviation, array, float, 
					in shape (3, paras_num+1), 3 for AOD, SSA and g respectively, 
					1 for all rate
	'''
	if 'path' in args:
		path = args['path']
	else:
		path = 'output/Mie/'
	if 'paras' in args:
		paras = args['paras']
	else:
		paras = default_paras
	
	path = path + dtime + '/'
	
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
	
	par_num = len(default_paras)
	Stds = np.zeros((3, par_num+1))
	Stds[0,0] = np.nanstd(AOD1)
	Stds[1,0] = np.nanstd(SSA1)
	Stds[2,0] = np.nanstd(g1)
	
	Cons = np.zeros((3, par_num))
	
	for i in range(par_num):
		if default_paras[i] in paras:
			data2_path = path + default_paras[i] + '/'
			data2 = np.load(data2_path+'infos.npy', allow_pickle=True)
			AOD2 = np.zeros(len(data2))
			SSA2 = np.zeros(len(data2))
			g2 = np.zeros(len(data2))
			for j in range(len(data2)):
				AOD2[j] = data2[j]['AOD']
				SSA2[j] = data2[j]['SSA']
				g2[j] = data2[j]['g']
			Stds[0,i+1] = np.nanstd(AOD2)
			Stds[1,i+1] = np.nanstd(SSA2)
			Stds[2,i+1] = np.nanstd(g2)
			Cons[0,i] = (np.nanstd(AOD1)-np.nanstd(AOD2)) / np.nanstd(AOD1)
			Cons[1,i] = (np.nanstd(SSA1)-np.nanstd(SSA2)) / np.nanstd(SSA1)
			Cons[2,i] = (np.nanstd(g1)-np.nanstd(g2)) / np.nanstd(g1)
		else:
			Stds[:,i+1] = np.ones(3) * 99999
			Cons[:,i] = np.ones(3) * 99999 # mask value
	
	Stds = Stds[np.where(Stds!=99999)].reshape(3,-1)
	Cons = Cons[np.where(Cons!=99999)].reshape(3,-1)
	infos = dict(run_num=run_num, paras=paras, Stds=Stds, Cons=Cons)
	return infos
	
if __name__ == '__main__':
	Cons = cal('230916_n_kappa', paras=['n','kappa','co'])
	print(Cons['Cons'])
	print(Cons['Stds'])
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
