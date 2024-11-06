####################################################################################
# INTRODUCTION:
# This code is to calculate runMie outputs' OAT sensitivity
# Created by Hebs at 24/6/3/9:14
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

default_outputs = [] # default outputs list
default_outputs.append('AOD')
default_outputs.append('SSA')
default_outputs.append('g')

import numpy as np
import os
import runMie
import readTaizhou
import calPNSD

def run(dtime, par_range, **args):
	'''
	This function is to run OAT of runMie.py
	input:
		dtime		: runMie.py output time, string, example: '230101'
		par_range	: parameter change rate range, float
		**path		: runMie.py output save path, string, default 'output/Mie/'
		**debug		: debug flat, bool, default False
	output:
		runMie.py output
	'''
	if 'path' in args:
		path = args['path']
	else:
		path = 'output/Mie/'
	if 'debug' in args:
		debug = args['debug']
	else:
		debug = False
	
	path = path + dtime + '/OAT/'
	
	# read data
	
	sp2 = readTaizhou.read_Taizhou('data/sp2/Taizhou.npy')
	Dps = sp2['Dps']
	DBC = sp2['DBC']
	DBCps = sp2['DBCps']
	PNSD = sp2['PNSD']
	DMASP2 = sp2['DMASP2']
	
	# turn time sequent data to mean data
	PNSD = np.nanmean(PNSD, axis=0)
	DMASP2 = np.nanmean(DMASP2, axis=0)
	# turn DMASP2 to BCPNSD,
	# due to DBCps has different dlogDBCps by bin, have to do special treatment
	BCPNSD = np.zeros(DMASP2.shape) # turn dn/dlogDBC to dn/dlogDBCps/dlogDBC
	for i in range(len(DBCps)):
		dlogDBC = calPNSD.cal_dlnDp(DBC)
		if i<len(DBCps)-1:
			dlogDBCps_i = np.log10(DBCps[i+1]/DBCps[i])
		else:
			dlogDBCps_i = np.log10(DBCps[i]/DBCps[i-1])
		BCPNSD[i] = DMASP2[i] / dlogDBCps_i
	
	clean_PNSD = PNSD / 50 # to fit Beijing aerosol number distribution level
	dirty_PNSD = PNSD / 10
	clean_BCPNSD = BCPNSD / 10
	dirty_BCPNSD = BCPNSD / 2
	
	# set parameters running list
	
	parameters = default_paras
	
	par_num = len(default_paras)
	paras_rate = np.zeros((par_num, par_num, 2))
	
	for i in range(par_num):
		for j in range(par_num):
			if i==j:
				paras_rate[i,j] = np.array([1-par_range, 1+par_range])
			else:
				paras_rate[i,j] = np.ones(2)
	
	# make directory for data saving
	
	os.system('mkdir '+path)
	os.system('mkdir '+path+'all/')
	for i in range(par_num):
		fn = path + parameters[i]
		os.system('mkdir '+fn)
	
	# running
	
	if debug:
		print('done\ncalculating...')
	
	for i in range(par_num):
		infos = []
		output_path = path + parameters[i] + '/'
		for j in range(2):
			info, paras = runMie.run(1000, [525], 70, MSH2x=5, 
			n_rate			= paras_rate[i,0,j], 
			nH1_rate		= paras_rate[i,1,j], 
			nBC_rate		= paras_rate[i,2,j], 
			kBC_rate		= paras_rate[i,3,j], 
			PNSD_rate		= paras_rate[i,4,j], 
			MS_rate			= paras_rate[i,5,j], 
			MSH1_rate		= paras_rate[i,6,j], 
			MSH2_rate		= paras_rate[i,7,j], 
			CT_rate			= paras_rate[i,8,j], 
			CTH1_rate		= paras_rate[i,9,j], 
			kappa_rate		= paras_rate[i,10,j], 
			kappaH1_rate	= paras_rate[i,11,j], 
			rhoBC_rate		= paras_rate[i,12,j], 
			BCPNSD_rate		= paras_rate[i,13,j], 
			BCAE_rate		= paras_rate[i,14,j], 
			Dps=Dps, DBC=DBC, DBCps=DBCps, clean_PNSD=clean_PNSD, dirty_PNSD=dirty_PNSD, clean_BCPNSD=clean_BCPNSD, dirty_BCPNSD=dirty_BCPNSD, debug=debug)
			infos.append(info)
			np.save(path+parameters[i]+'/infos.npy', infos)
			if i==0 and j==0:
				np.save(path+'all/paras.npy', dict(paras=paras,paras_rate=paras_rate))
			if debug:
				print(round((i*2+j+1)/par_num/2*100), '% done...')
	
	if debug:
		print('done')

def cal(dtime, **args):
	'''
	This function is to read OAT running results
	input:
		dtime		: runMie.py output time, string, example: '230101'
		**path		: runMie.py output save path, string, default 'output/Mie/'
		**paras		: parameters list, array, string, default default_paras
		**outputs	: outputs list, array, string, default default_outputs
	output:
		in dictionary
		paras	: parameters list, array, string
		outputs	: outputs list, array, string
		OAT		: OAT reasults, array, in shape(len(outputs),len(paras))
	'''
	if 'path' in args:
		path = args['path']
	else:
		path = 'output/Mie/'
	if 'paras' in args:
		paras = args['paras']
	else:
		paras = default_paras
	if 'outputs' in args:
		outputs = args['outputs']
	else:
		outputs = default_outputs
	
	path = path + dtime + '/OAT/'
	
	# read data
	
	paras_data = np.load(path+'all/paras.npy', allow_pickle=True).item()
	paras_rate = paras_data['paras_rate']
	par_num = len(default_paras)
	out_num = len(outputs)
	
	OAT = np.zeros((out_num,par_num))
	
	for i in range(par_num):
		if default_paras[i] in paras:
			data_path = path + default_paras[i] + '/'
			data = np.load(data_path+'infos.npy', allow_pickle=True)
			
			for j in range(out_num):
				OAT[j,i] = abs((data[1][outputs[j]]-data[0][outputs[j]])/(paras_rate[i,i,1]-paras_rate[i,i,0])) + 1e-9
	
	OAT = OAT[np.where(OAT!=0)].reshape(out_num,-1)
	infos = dict(paras=paras, outputs=outputs, OAT=OAT)
	return infos

if __name__ == '__main__':
	run('240603', 0.1, debug=True)
	infos = cal('240603')
	print(infos['OAT'])
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
