####################################################################################
# INTRODUCTION:
# This code is to read DMA-SP2 data
# Created by Hebs at 21/11/2/10:03
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np

def read_Taizhou(fn):
	'''
	This function is to raed Taizhou data
	input:
		fn       : SP2 data file path, string
	output:
		Taizhou  : Taizhou data
			Dps     : partical diameter distribution, array, nm
			PNSD    : partical number concentration size distribution, array, cm^-3
			DBC     : BC core diameter distribution, array, nm
			DBCps   : BC core have partical diameter distribution, array, nm
			DMASP2  : BC size-resolved particle number size distribution, array, dn/dlnDBC
			BCPNSD  : BC partical number size distribution, array, dn/dlnDBC
			BCPVSD  : BC partical volume size distribution, array, d(um^3/cm^3)/dlnDBC
			ksca    : scattering coefficient, Mm^-1
			wl_sca  : ksca wave length, nm
			kabs    : absorbing coefficient, Mm^-1
			wl_abs  : kabs wave length, nm
	'''
	data = np.load(fn, allow_pickle=True, encoding='latin1')
	Dps = data.item()['Dp']
	PNSD = data.item()['PNSDs']
	DBCps = data.item()['DpSp2']
	DBC = data.item()['DcSp2']
	DMASP2 = data.item()['dN']
	BCPNSD = data.item()['BCpnsd']
	BCPVSD = data.item()['BCpvsd']
	ksca = data.item()['SCA']
	wl_sca = data.item()['WvNehp']
	kabs = data.item()['ABSpass']
	wl_abs = data.item()['WvPass']
	jul = data.item()['jul']
	sp2 = dict(Dps=Dps, PNSD=PNSD, DBCps=DBCps, DBC=DBC, DMASP2=DMASP2, BCPNSD=BCPNSD, BCPVSD=BCPVSD, ksca=ksca, wl_sca=wl_sca, kabs=kabs, wl_abs=wl_abs, jul=jul)
	return sp2

if __name__ == '__main__':
	data = np.load('data/sp2/Taizhou.npy', allow_pickle=True, encoding='latin1')
	keys = []
	for key in data.item().keys():
		keys.append(key)
	keys.sort()
	for key in keys:
		print(key)
	DBC = data.item()['DcSp2']
	BCpnsd = data.item()['BCpnsd']
	print(DBC.shape)
	'''
	Dps = data.item()['Dps']
	print(Dps.shape)
	PNSD = data.item()['PNSDs']
	print(PNSD.shape)
	BC_Wv = data.item()['BC_Wv']
	print(BC_Wv)
	BC = data.item()['BC']
	print(BC)
	time = data.item()['time']
	print(time)
	
	item = data.item()['WvPass']
	print(item.shape)
	print(item)
	
	sp2 = read_SP2('data/sp2/allTaizhouData_2.npy')
	print(sp2['BCPNSD'].shape)
	
	data = np.load('data/sp2/Taizhou.npy', allow_pickle=True, encoding='latin1').item()
	keys = []
	for key in data.keys():
		keys.append(key)
	keys.sort()
	for key in keys:
		print(key, data[key].shape)
	
	
	x = np.arange(len(PNSD))
	import matplotlib.pyplot as plt
	plt.figure(figsize=(12,6))
	plt.plot(x, v_BCPNSD*1e-9, label='BCpnsd', lw=10, alpha=0.5)
	plt.plot(x, v_DMAsp2*2e-9, label='doubled dN')
	plt.plot(x, v_DMAsp2corrected*2e-9, label='doubled dNcorrected')
	plt.xticks([], fontsize=15, fontweight='bold')
	plt.xlabel('time', fontsize=15, fontweight='bold')
	plt.ylim(0, 240)
	plt.ylabel('BC total volume, $um^3/cm^3$', fontsize=15, fontweight='bold')
	plt.yticks(fontsize=15, fontweight='bold')
	plt.legend(fontsize=15)
	plt.show()
	plt.close()
	plt.figure(figsize=(8,6))
	plt.scatter(v_BCPNSD*1e-9, v_DMAsp2*1e-9, label='dN')
	plt.scatter(v_BCPNSD*1e-9, v_DMAsp2corrected*1e-9, label='dNcorrected')
	plt.plot(np.arange(200), np.arange(200)/2, lw=2, ls='--', label='1:2 line', c='green')
	plt.xlim(0,200)
	plt.ylim(0,100)
	plt.xlabel('BCpnsd total volume, $um^3/cm^3$', fontsize=15, fontweight='bold')
	plt.xticks([0, 40, 80, 120, 160, 200], [0, 40, 80, 120, 160, 200], fontsize=15, fontweight='bold')
	plt.ylabel('BC total volume, $um^3/cm^3$', fontsize=15, fontweight='bold')
	plt.yticks(fontsize=15, fontweight='bold')
	plt.legend(fontsize=15)
	plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
	plt.show()
		
	sp2 = read_Taizhou('data/sp2/Taizhou.npy')
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
	clean_PNSD = np.nanmean(PNSD, axis=0)
	dirty_PNSD = PNSD[i_max]
	clean_BCPNSD = np.nanmean(BCPNSD, axis=0)
	dirty_BCPNSD = BCPNSD[iBC_max]
	
	import matplotlib.pyplot as plt
	plt.plot(np.arange(len(DBCps)),np.nansum(clean_BCPNSD,axis=1))
	plt.plot(np.arange(len(DBCps)),np.nansum(dirty_BCPNSD,axis=1))
	plt.show()
	plt.close()
	'''
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
