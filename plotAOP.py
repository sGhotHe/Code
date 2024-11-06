####################################################################################
# INTRODUCTION:
# This code is to use pymiescatt to plot AOP factors' sensitivity
# Created by Hebs at 24/6/20/9:08
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import matplotlib.pyplot as plt
import PyMieScatt as ps

import readTaizhou
import calPNSD
import calNI
import time

def plot(RH, wl, **args):
	'''
	This function is to plot PyMieScatt output factor sensitivity,
	factors including: n, kappa, PNSD
	input:
		RH		: relative humidity, percent, float
		wl		: wavelength, nm, float
		**rate	: factor change rate, float, default 1.1
		**save	: save flag, bool, default False
	output:
		figure
	'''
	if 'rate' in args:
		rate = args['rate']
	else:
		rate = 1.01
	if 'save' in args:
		save = args['save']
	else:
		save = False
	if 'save_path' in args:
		save_path = args['save_path']
	else:
		save_path = 'figure/DARF/'
	
	# read data
	
	sp2 = readTaizhou.read_Taizhou('data/sp2/Taizhou.npy')
	Dps = sp2['Dps']
	DBC = sp2['DBC']
	DBCps = sp2['DBCps']
	PNSD = sp2['PNSD']
	DMASP2 = sp2['DMASP2'] # dn/dlogDBC
	
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
	
	clean_PNSD = PNSD / 1000 # to fit Beijing aerosol number distribution level
	dirty_PNSD = PNSD / 20
	clean_BCPNSD = BCPNSD / 500
	dirty_BCPNSD = BCPNSD / 100
	
	def cal_PNSD(PNSD_rate):
		return clean_PNSD * PNSD_rate + dirty_PNSD * (1-PNSD_rate)
	
	PNSD = cal_PNSD(0.5)
	PNSD_new = cal_PNSD(0.5*rate)
	BCPNSD = clean_BCPNSD * 0.5 + dirty_BCPNSD * 0.5
	
	nShell = 1.58 # from Beijing University observation
	nBC = 1.67 # from G.Zhao
	kBC = 0.67 # from G.Zhao
	kappa = 0.21537370768611558 # from Beijing University observation
	MS = 0.7
	
	# calculate sensitivity
	
	S_n_aod = np.zeros(len(Dps))
	S_n_ssa = np.zeros(len(Dps))
	S_n_g = np.zeros(len(Dps))
	S_PNSD_aod = np.zeros(len(Dps))
	S_PNSD_ssa = np.zeros(len(Dps))
	S_PNSD_g = np.zeros(len(Dps))
	S_kappa_aod = np.zeros(len(Dps))
	S_kappa_ssa = np.zeros(len(Dps))
	S_kappa_g = np.zeros(len(Dps))
	'''
	for i in range(len(Dps)):
		DBC_i = Dps[i] / 5
		
		D_i, m_i = kappaKohler.RH2D(Dps[i], DBC_i, mBC, n+1e-9j, kappa, RH, wl)
		MieQCoreShell = ps.MieQCoreShell(mBC, m_i, wl, DBC_i, D_i)
		D_i, m_i = kappaKohler.RH2D(Dps[i], 0, mBC, n+1e-9j, kappa, RH, wl)
		MieQ = ps.MieQ(m_i, wl, D_i)
		aod = MieQCoreShell[1] * MS + MieQ[1] * (1-MS)
		ssa = aod / (MieQCoreShell[0]*MS+MieQ[0]*(1-MS))
		g = MieQCoreShell[3] * MS + MieQ[3] * (1-MS)
		
		D_i, m_i = kappaKohler.RH2D(Dps[i], DBC_i, mBC, n*rate+1e-9j, kappa, RH, wl)
		MieQCoreShell_n = ps.MieQCoreShell(mBC, m_i, wl, DBC_i, D_i)
		D_i, m_i = kappaKohler.RH2D(Dps[i], 0, mBC, n*rate+1e-9j, kappa, RH, wl)
		MieQ_n = ps.MieQ(m_i, wl, D_i)
		aod_n = MieQCoreShell_n[1] * MS + MieQ_n[1] * (1-MS)
		ssa_n = aod_n / (MieQCoreShell_n[0]*MS+MieQ_n[0]*(1-MS))
		g_n = MieQCoreShell_n[3] * MS + MieQ_n[3] * (1-MS)
		
		D_i, m_i = kappaKohler.RH2D(Dps[i], DBC_i, mBC, n+1e-9j, kappa*rate, RH, wl)
		MieQCoreShell_kappa = ps.MieQCoreShell(mBC, m_i, wl, DBC_i, D_i)
		D_i, m_i = kappaKohler.RH2D(Dps[i], 0, mBC, n+1e-9j, kappa*rate, RH, wl)
		MieQ_kappa = ps.MieQ(m_i, wl, D_i)
		aod_kappa = MieQCoreShell_kappa[1] * MS + MieQ_kappa[1] * (1-MS)
		ssa_kappa = aod_kappa / (MieQCoreShell_kappa[0]*MS+MieQ_kappa[0]*(1-MS))
		g_kappa = MieQCoreShell_kappa[3] * MS + MieQ_kappa[3] * (1-MS)
		
		S_n_aod[i] = (aod_n-aod) / ((rate-1)*100) * PNSD[i]
		S_n_ssa[i] = (ssa_n-ssa) / ((rate-1)*100) * PNSD[i]
		S_n_g[i] = (g_n-g) / ((rate-1)*100) * PNSD[i]
		
		S_PNSD_aod[i] = aod / ((rate-1)*100) * (PNSD[i]-PNSD_new[i])
		S_PNSD_ssa[i] = ssa / ((rate-1)*100) * (PNSD[i]-PNSD_new[i])
		S_PNSD_g[i] = g / ((rate-1)*100) * (PNSD[i]-PNSD_new[i])
		
		S_kappa_aod[i] = (aod_kappa-aod) / ((rate-1)*100) * PNSD[i]
		S_kappa_ssa[i] = (ssa_kappa-ssa) / ((rate-1)*100) * PNSD[i]
		S_kappa_g[i] = (g_kappa-g) / ((rate-1)*100) * PNSD[i]
	'''
	'''
	kext, waer, g = calNI.cal_Mie5(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, 0, kappa, 0, 0, 0, 1, 0, RH, wl, 6, 1, 0, 1, angularResolution=30)
	kext_n, waer_n, g_n = calNI.cal_Mie5(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell*rate, 0, kappa, 0, 0, 0, 1, 0, RH, wl, 6, 1, 0, 1, angularResolution=30)
	kext_PNSD, waer_PNSD, g_PNSD = calNI.cal_Mie5(Dps, PNSD_new, DBCps, DBC, BCPNSD, kBC, nBC, nShell, 0, kappa, 0, 0, 0, 1, 0, RH, wl, 6, 1, 0, 1, angularResolution=30)
	kext_kappa, waer_kappa, g_kappa = calNI.cal_Mie5(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, 0, kappa*rate, 0, 0, 0, 1, 0, RH, wl, 6, 1, 0, 1, angularResolution=30)
	infos = dict(kext=kext, waer=waer, g=g, kext_n=kext_n, waer_n=waer_n, g_n=g_n, kext_PNSD=kext_PNSD, waer_PNSD=waer_PNSD, g_PNSD=g_PNSD, kext_kappa=kext_kappa, waer_kappa=waer_kappa, g_kappa=g_kappa)
	np.save('Mie_infos.npy', infos)
	'''
	infos = np.load('Mie_infos.npy', allow_pickle=True).item()
	kext = infos['kext']
	waer = infos['waer']
	g = infos['g']
	kext_n = infos['kext_n']
	waer_n = infos['waer_n']
	g_n = infos['g_n']
	kext_PNSD = infos['kext_PNSD']
	waer_PNSD = infos['waer_PNSD']
	g_PNSD = infos['g_PNSD']
	kext_kappa = infos['kext_kappa']
	waer_kappa = infos['waer_kappa']
	g_kappa = infos['g_kappa']
	
	S_n_aod = (kext_n-kext) / (rate-1) / 100 # in dy/dx%
	S_PNSD_aod = -(kext_PNSD-kext) / (rate-1) / 100
	S_kappa_aod = (kext_kappa-kext) / (rate-1) / 100
	S_n_ssa = (waer_n-waer) / (rate-1) / 100
	S_PNSD_ssa = -(waer_PNSD-waer) / (rate-1) / 100
	S_kappa_ssa = (waer_kappa-waer) / (rate-1) / 100
	S_n_g = (g_n-g) / (rate-1) / 100
	S_PNSD_g = (g_PNSD-g) / (rate-1) / 100
	S_kappa_g = (g_kappa-g) / (rate-1) / 100
	
	# start plotting
	
	plt.figure(figsize=(10,3))
	plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.2, wspace=0.3, hspace=0.3)
	fs = 10
	fw = 'bold'
	lw = 2.5
	
	ax = plt.subplot(131)
	plt.plot(Dps, S_n_aod, label='$S_{n\ shell}$', linewidth=lw)
	plt.plot(Dps, S_PNSD_aod, label='$S_{PNSD\ dry}$', linewidth=lw)
	plt.plot(Dps, S_kappa_aod, label='$S_{kappa}$', linewidth=lw)
	plt.text(0.05, 0.85, '\u03bb = 525nm', fontsize=fs, fontweight=fw, color='k', transform=ax.transAxes)
	plt.title('kext', fontsize=fs, fontweight=fw)
	plt.xscale('log')
	plt.xlabel('Dps, nm', fontsize=fs, fontweight=fw)
	plt.xticks(fontsize=fs, fontweight=fw)
	plt.ylabel('Sensitivity, dy/dx%', fontsize=fs, fontweight=fw)
	plt.yticks(fontsize=fs, fontweight=fw)
	plt.legend()
	
	plt.subplot(132)
	plt.plot(Dps, S_n_ssa, linewidth=lw)
	plt.plot(Dps, S_PNSD_ssa, linewidth=lw)
	plt.plot(Dps, S_kappa_ssa, linewidth=lw)
	plt.title('ksca / kext', fontsize=fs, fontweight=fw)
	plt.xscale('log')
	plt.xticks(fontsize=fs, fontweight=fw)
	plt.yticks(fontsize=fs, fontweight=fw)
	
	plt.subplot(133)
	plt.plot(Dps, S_n_g, linewidth=lw)
	plt.plot(Dps, S_PNSD_g, linewidth=lw)
	plt.plot(Dps, S_kappa_g, linewidth=lw)
	plt.title('g', fontsize=fs, fontweight=fw)
	plt.xscale('log')
	plt.xticks(fontsize=fs, fontweight=fw)
	plt.yticks(fontsize=fs, fontweight=fw)
	
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_Mie.pdf')
	else:
		plt.show()

if __name__ == '__main__':
	plot(70, 525, save=False)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
