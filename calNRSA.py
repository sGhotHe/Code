####################################################################################
# INTRODUCTION:
# This code is to calculate Nominal Range Sensitivity Analysis
# Created by Hebs at 23/6/5/15:43
# Contact: hebishuo@pku.edu.cn
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

import numpy as np

import runMie
import readTaizhou

def cal(run_num, par_range, **args):
	'''
	This function is to calculate Nominal Range Sensitivity Analysis
	input:
		run_num	: running times number, int
		par_range	: parameters change range, float
		**paras	: parameters list, array, string, default default_paras
		**debug	: debug flat, bool, default False
	output:
		NRSA		: result of Nominal Range Sensitivity Analysis, array, in shape 
				  (3, par_num, run_num-1)
	'''
	if 'paras' in args:
		paras = args['paras']
	else:
		paras = default_paras
	if 'debug' in args:
		debug = args['debug']
	else:
		debug = False
	
	#read data
	
	if debug:
		print('reading...')
	
	sp2 = readTaizhou.read_Taizhou('data/sp2/Taizhou.npy')
	Dps = sp2['Dps']
	DBC = sp2['DBC']
	DBCps = sp2['DBCps']
	PNSD = sp2['PNSD']
	BCPNSD = sp2['DMASP2']
	N = np.nansum(PNSD, axis=1)
	NBC = np.nansum(np.nansum(BCPNSD, axis=2), axis=1)
	
	N_max = 0
	NBC_max = 0
	i_max = 0
	iBC_max = 0
	
	for i in range(len(N)):
		if N[i] > N_max:
			N_max = N[i]
			i_max = i
		if NBC[i] > NBC_max:
			NBC_max = NBC[i]
			iBC_max = i
	
	clean_PNSD = np.nanmean(PNSD, axis=0) / 5
	dirty_PNSD = np.nanmean(PNSD, axis=0)
	clean_BCPNSD = np.nanmean(BCPNSD, axis=0) / 5
	dirty_BCPNSD = np.nanmean(BCPNSD, axis=0) 
	
	data = np.load('data/BCMD/BCPMSD.npy', allow_pickle=True).item()
	DBC2 = data['DBC']
	BCPMSD = data['BCPMSD']
	
	clean_BCPMSD = BCPMSD / 5
	dirty_BCPMSD = BCPMSD
	
	# set parameters running list
	
	par_num = len(default_paras)
	paras_rate = np.zeros((par_num, par_num, run_num))
	
	for i in range(par_num):
		for j in range(par_num):
			if i==j and default_paras[i] in paras:
				paras_rate[i,j] = np.arange(1-par_range, 1+par_range, 2*par_range/(run_num-1))
			else:
				paras_rate[i,j] = np.ones(run_num)
	
	# running
	
	if debug:
		print('done\ncalculating...')
	
	NRSA = np.zeros((3, par_num, run_num-1)) # 3 for AOD, SSA and g
	
	for i in range(par_num):
		if default_paras[i] in paras:
			AOD = np.zeros(run_num)
			SSA = np.zeros(run_num)
			g = np.zeros(run_num)
			para_rate = np.zeros((par_num, run_num))
			for j in range(run_num):
				info, para = runMie.run(2000, 100, [525], nI2x=2, kappaI2x=2,
				n_rate=	paras_rate[i,0,j], 
				nI_rate=	paras_rate[i,1,j], 
				nI2_rate=	paras_rate[i,2,j],
				nBC_rate=	paras_rate[i,3,j], 
				kBC_rate=	paras_rate[i,4,j], 
				PNSD_rate=	paras_rate[i,5,j], 
				MS_rate=	paras_rate[i,6,j], 
				VD_rate=	paras_rate[i,7,j], 
				CT_rate=	paras_rate[i,8,j], 
				kappa_rate=	paras_rate[i,9,j], 
				kappaI_rate=	paras_rate[i,10,j], 
				kappaI2_rate=	paras_rate[i,11,j], 
				rhoBC_rate=	paras_rate[i,12,j], 
				BCPNSD_rate=	paras_rate[i,13,j], 
				BCPMSD_rate=	paras_rate[i,14,j], 
				BCI_rate=	paras_rate[i,15,j], 
				BCAE_rate=	paras_rate[i,16,j], 
				Dps=Dps, DBC=DBC, DBCps=DBCps, clean_PNSD=clean_PNSD, dirty_PNSD=dirty_PNSD, clean_BCPNSD=clean_BCPNSD, dirty_BCPNSD=dirty_BCPNSD, DBC2=DBC2, clean_BCPMSD=clean_BCPMSD, dirty_BCPMSD=dirty_BCPMSD, debug=False)
				AOD[j] = info['AOD']
				SSA[j] = info['SSA']
				g[j] = info['g']
				for k in range(par_num):
					para_rate[k,j] = info[default_paras[k]+'_rate']
				if debug:
					print(round((i*run_num+j+1)/par_num/run_num*100), '% done...')
			
			NRSA[0,i] = (AOD[1:]-AOD[:-1]) / (para_rate[i,1:]-para_rate[i,:-1])
			NRSA[1,i] = (SSA[1:]-SSA[:-1]) / (para_rate[i,1:]-para_rate[i,:-1])
			NRSA[2,i] = (g[1:]-g[:-1]) / (para_rate[i,1:]-para_rate[i,:-1])
		else:
			NRSA[:,i] = np.ones((3,run_num-1)) * 99999 # mask value
			if debug:
				print(round((i+1)/par_num*100), '% done...')
	
	if debug:
		print('done')
	
	NRSA = NRSA[np.where(NRSA!=99999)].reshape(3,-1,run_num-1)
	return NRSA

if __name__ == '__main__':
	NRSA = cal(2, 0.1, debug=True)
	print(NRSA)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
