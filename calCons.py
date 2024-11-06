####################################################################################
# INTRODUCTION:
# This code is to calculate Constraint Factor Sensitivity Analysis's result
# Created by Hebs at 23/6/7/9:21
# Contact: hebishuo@pku.edu.cn
####################################################################################

default_range = 0.3
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

def cal(dtime, par_range, **args):
	'''
	This function is to calculate Constraint Factor Sensitivity Analysis use runMie.py's result
	input:
		dtime		: runMie.py output time, string, example: '230101'
		par_range	: parameter change rate range, float
		**path		: runMie.py output save path, string, default 'output/Mie/'
		**paras		: parameters list, array, string, default default_paras
	output:
		in dictionary
		paras		: parameters list, array, string
		run_num	: running number, int
		Cons		: Constraint Factor Sensitivity Analysis result, array, float,
				  in shape (3, paras_num), 3 for AOD, SSA and g respectively
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
	
	# read aerosol optical parameters
	data1_path = path + 'all/'
	data1 = np.load(data1_path+'infos.npy', allow_pickle=True)
	run_num = len(data1)
	par_num = len(default_paras)
	
	AOD1 = np.zeros(len(data1))
	SSA1 = np.zeros(len(data1))
	g1 = np.zeros(len(data1))
	
	for i in range(len(data1)):
		AOD1[i] = data1[i]['AOD']
		SSA1[i] = data1[i]['SSA']
		g1[i] = data1[i]['g']
	
	# read factor values
	x = np.load(data1_path+'paras.npy', allow_pickle=True).item()
	
	Cons = np.zeros((3, par_num))
	# according to error propagation, sensitivity would be:
	# S = d(sigma_y^2) / d(sigma_x^2)
	# calulate sigma_x for every factor first
	
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
			
			Cons[0,i] = np.sqrt(abs((np.nanvar(AOD1)-np.nanvar(AOD2))/0.33**2/(default_range**2-par_range**2)))
			Cons[1,i] = np.sqrt(abs((np.nanvar(SSA1)-np.nanvar(SSA2))/0.33**2/(default_range**2-par_range**2)))
			Cons[2,i] = np.sqrt(abs((np.nanvar(g1)-np.nanvar(g2))/0.33**2/(default_range**2-par_range**2)))
		else:
			Cons[:,i] = np.ones(3) * 99999 # mask value
	
	Cons = Cons[np.where(Cons!=99999)].reshape(3,-1)
	infos = dict(run_num=run_num, paras=paras, Cons=Cons)
	return infos

if __name__ == '__main__':
	Cons = cal('240604', 0.05)
	print(Cons['Cons'])
	print(default_paras)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
