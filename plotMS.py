####################################################################################
# INTRODUCTION:
# This code is to plot how BC and kappa mixing state affact particles' optical properties figure
# Created by Hebs at 23/2/23/9:57
# Contact: hebishuo@pku.edu.cn
####################################################################################

import sys
import os
import time
import numpy as np
import matplotlib.pyplot as plt

import calMie
import calBC
import calRH
import readTaizhou
import readAtms

def cal_nI(min_n, max_n, nI, **args):
	'''
	This function is to calculate optical properties and 30% n change result stantard deviation optical properties data
	input:
		min_n :	minimum complex refractive index real part, float
		max_n :	maximum complex refractive index real part, float
		nI :		maximum complex refractive index real part change slope, float
		**dx :		x axis bin size, float, default 0.1
		**dy :		y axis bin size, float, default 0.1
		**debug :	debug flag, bool, default False
	output:
		n :		complex refractive index real part data, array, float
		nI :		complex refractive index real part mixing state data, array, float
		AOD :		aerosol optical depth, two dimentional array, float
		SSA :		single scattering albedo, two dimentional array, float
		g :		asymmetry factor, two dimentional array, float
		std_AOD :	standard deviation aerosol optical depth, two dimentional array, 					float
		std_SSA :	standard deviation single scattering albedo, two dimentional 				array, float
		std_g :	standard deviation asymmetry factor, two dimentional array, float
	'''
	if 'dx' in args:
		dx = args['dx']
	else:
		dx = 0.01
	if 'dy' in args:
		dy = args['dy']
	else:
		dy = 0.01
	if 'debug' in args:
		debug = args['debug']
	else:
		debug = False
	
	# read in defualt data
	
	fn = 'data/sp2/Taizhou.npy'
	if os.path.exists(fn):
		if debug:
			print('reading BCPNSD data...')
		BCPNSD_exist = True
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
		
		clean_PNSD = np.nanmean(PNSD, axis=0) / 10 # cm-3/dlnDps
		dirty_PNSD = np.nanmean(PNSD, axis=0)
		clean_BCPNSD = np.nanmean(BCPNSD, axis=0) / 10 # cm-3/dlogDBC
		dirty_BCPNSD = np.nanmean(BCPNSD, axis=0)
	else:
		print('no default PNSD data at', fn, ', please check')
		sys.exit()
	
	PNSD = clean_PNSD * 0.5 + dirty_PNSD * 0.5
	BCPNSD = clean_BCPNSD * 0.5 + dirty_BCPNSD * 0.5
	wl = [525]
	Z = 2000
	dZ = 100
	MS = 0.7
	BCPNSD = clean_BCPNSD * 0.5 + dirty_BCPNSD * 0.5
	BCI = calBC.cal_BCI(DBC, DBCps, BCPNSD)
	VD = 28 / (28+39)
	rhoBC = 0.95
	nBC, kBC = calBC.rhoBC2mBC(rhoBC)
	kappa = 0.21537370768611558
	kappaI = 0
	#nShell = 1.58
	#nI = 0
	atms = readAtms.read()
	z = atms['z'] * 1e3
	dz = z[:-1] - z[1:]
	wh = atms['wh']
	t = atms['t']
	RH = np.zeros(len(wh))
	for i in range(len(RH)):
		RH[i] = calRH.wh2RH(t[i], wh[i])
		if RH[i]<1:
			RH[i] = 1
	for i in range(len(wl)-1):
		if (wl[i+1]-wl[i])<0:
			print('wl(k) < wl(k+1) doesn\'t meet. Please check.')
			sys.exit()
	if Z > z[0]:
		print('Too high height. Please check.')
		sys.exit()
	i = 0
	while i < len(z)-1:
		if z[i]>=Z and z[i+1]<=Z:
			break
		i = i + 1
	RH = RH[i]
	BCAE = 1
	CT = 1
	PF = 1
	
	# start running
	
	n = np.arange(min_n, max_n, dx)
	nI = np.arange(-1*nI, nI, dy)
	AOD = np.zeros((len(n),len(nI)))
	SSA = np.zeros((len(n),len(nI)))
	g = np.zeros((len(n),len(nI)))
	std_AOD = np.zeros((len(n),len(nI)))
	std_SSA = np.zeros((len(n),len(nI)))
	std_g = np.zeros((len(n),len(nI)))
	
	for i in range(len(n)):
		for j in range(len(nI)):
			AOD[i,j], SSA[i,j], g[i,j] = calMie.Mie(Dps, PNSD, DBCps, DBC, BCPNSD, nBC, kBC, n[i], nI[j], kappa, kappaI, RH, wl[0], Z, dZ, BCAE, CT, VD, PF, debug=debug)
			AOD1, SSA1, g1 = calMie.Mie(Dps, PNSD, DBCps, DBC, BCPNSD, nBC, kBC, n[i]*0.7, nI[j], kappa, kappaI, RH, wl[0], Z, dZ, BCAE, CT, VD, PF, debug=debug)
			AOD2, SSA2, g2 = calMie.Mie(Dps, PNSD, DBCps, DBC, BCPNSD, nBC, kBC, n[i]*1.3, nI[j], kappa, kappaI, RH, wl[0], Z, dZ, BCAE, CT, VD, PF, debug=debug)
			std_AOD[i,j] = abs(AOD1-AOD2) / 6 / AOD[i,j]
			std_SSA[i,j] = abs(SSA1-SSA2) / 6 / SSA[i,j]
			std_g[i,j] = abs(g1-g2) / 6 / g[i,j]
			# for normal distribution, 6 * sigma = 99.73% x axis interval
	
	return n, nI, AOD, SSA, g, std_AOD, std_SSA, std_g

def cal_kappaI(min_kappa, max_kappa, kappaI, **args):
	'''
	This function is to calculate optical properties and 30% kappa change result stantard deviation optical properties data
	input:
		min_kappa :	minimum kappa, float
		max_kappa :	maximum kappa, float
		kappaI :	maximum kappa change slope, float
		**dx :		x axis bin size, float, default 0.1
		**dy :		y axis bin size, float, default 0.1
		**debug :	debug flag, bool, default False
	output:
		kappa :	kappa data, array, float
		kappaI :	kappa mixing inhomogeneity data, array, float
		AOD :		aerosol optical depth, two dimentional array, float
		SSA :		single scattering albedo, two dimentional array, float
		g :		asymmetry factor, two dimentional array, float
		std_AOD :	standard deviation aerosol optical depth, two dimentional array, 					float
		std_SSA :	standard deviation single scattering albedo, two dimentional 				array, float
		std_g :	standard deviation asymmetry factor, two dimentional array, float
	'''
	if 'dx' in args:
		dx = args['dx']
	else:
		dx = 0.01
	if 'dy' in args:
		dy = args['dy']
	else:
		dy = 0.01
	if 'debug' in args:
		debug = args['debug']
	else:
		debug = False
	
	# read in defualt data
	
	fn = 'data/sp2/Taizhou.npy'
	if os.path.exists(fn):
		if debug:
			print('reading BCPNSD data...')
		BCPNSD_exist = True
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
		
		clean_PNSD = np.nanmean(PNSD, axis=0) / 10 # cm-3/dlnDps
		dirty_PNSD = np.nanmean(PNSD, axis=0)
		clean_BCPNSD = np.nanmean(BCPNSD, axis=0) / 10 # cm-3/dlogDBC
		dirty_BCPNSD = np.nanmean(BCPNSD, axis=0)
	else:
		print('no default PNSD data at', fn, ', please check')
		sys.exit()
	
	PNSD = clean_PNSD * 0.5 + dirty_PNSD * 0.5
	BCPNSD = clean_BCPNSD * 0.5 + dirty_BCPNSD * 0.5
	wl = [525]
	Z = 2000
	dZ = 100
	MS = 0.7
	BCPNSD = clean_BCPNSD * 0.5 + dirty_BCPNSD * 0.5
	BCI = calBC.cal_BCI(DBC, DBCps, BCPNSD)
	VD = 28 / (28+39)
	rhoBC = 0.95
	nBC, kBC = calBC.rhoBC2mBC(rhoBC)
	#kappa = 0.21537370768611558
	#kappaI = 0
	nShell = 1.58
	nI = 0
	atms = readAtms.read()
	z = atms['z'] * 1e3
	dz = z[:-1] - z[1:]
	wh = atms['wh']
	t = atms['t']
	RH = np.zeros(len(wh))
	for i in range(len(RH)):
		RH[i] = calRH.wh2RH(t[i], wh[i])
		if RH[i]<1:
			RH[i] = 1
	for i in range(len(wl)-1):
		if (wl[i+1]-wl[i])<0:
			print('wl(k) < wl(k+1) doesn\'t meet. Please check.')
			sys.exit()
	if Z > z[0]:
		print('Too high height. Please check.')
		sys.exit()
	i = 0
	while i < len(z)-1:
		if z[i]>=Z and z[i+1]<=Z:
			break
		i = i + 1
	RH = RH[i]
	BCAE = 1
	CT = 1
	PF = 1
	
	# start running
	
	kappa = np.arange(min_kappa, max_kappa, dx)
	kappaI = np.arange(-1*kappaI, kappaI, dy)
	AOD = np.zeros((len(kappa),len(kappaI)))
	SSA = np.zeros((len(kappa),len(kappaI)))
	g = np.zeros((len(kappa),len(kappaI)))
	std_AOD = np.zeros((len(kappa),len(kappaI)))
	std_SSA = np.zeros((len(kappa),len(kappaI)))
	std_g = np.zeros((len(kappa),len(kappaI)))
	
	for i in range(len(kappa)):
		for j in range(len(kappaI)):
			AOD[i,j], SSA[i,j], g[i,j] = calMie.Mie(Dps, PNSD, DBCps, DBC, BCPNSD, nBC, kBC, nShell, nI, kappa[i], kappaI[j], RH, wl[0], Z, dZ, BCAE, CT, VD, PF, debug=debug)
			AOD1, SSA1, g1 = calMie.Mie(Dps, PNSD, DBCps, DBC, BCPNSD, nBC, kBC, nShell, nI, kappa[i]*0.7, kappaI[j], RH, wl[0], Z, dZ, BCAE, CT, VD, PF, debug=debug)
			AOD2, SSA2, g2 = calMie.Mie(Dps, PNSD, DBCps, DBC, BCPNSD, nBC, kBC, nShell, nI, kappa[i]*1.3, kappaI[j], RH, wl[0], Z, dZ, BCAE, CT, VD, PF, debug=debug)
			std_AOD[i,j] = abs(AOD1-AOD2) / 6 / AOD[i,j]
			std_SSA[i,j] = abs(SSA1-SSA2) / 6 / SSA[i,j]
			std_g[i,j] = abs(g1-g2) / 6 / g[i,j]
			# for normal distribution, 6 * sigma = 99.73% x axis interval
	
	return kappa, kappaI, AOD, SSA, g, std_AOD, std_SSA, std_g

def plot_kappaI(kappa, kappaI, AOD, SSA, g, std_AOD, std_SSA, std_g, **args):
	'''
	This function is to plot kappa mixing state effect,
	it's a two dimentional figure, x axis is kappa and y axis is kappaI, 
	whith aerosol optical depth AOD, single scattering albedo SSA or 
	asymmetry factor g color map, for that is the most influenced properties
	input:
		kappa :	kappa data, array, float
		kappaI :	kappa mixing state slope data, array, float
		AOD :		aerosol optical depth data, 2-d array, float
		SSA :		single scattering albedo data, 2-d array, float
		g :		asymmetry factor data, 2-d array, float
		std_AOD :	standard aerosol optical depth data for 50% kappa change, 
				2-d array, float
		std_SSA :	standard single scattering albedo data for 50% kappa change, 
				2-d array, float
		std_g :	standard asymmetry factor data for 50% kappa change, 2-d array, 					float
		**debug :	debug flag, bool, default False
		**save :	whether to save figure, bool, default False
		**save_path :	figure save path, string, default 'figure/MS/'
	output:
		figure 1 :	x kappa, y kappaI and colored AOD, SSA and g
		figure 2 :	x kappa, y kappaI and colored standard deviation AOD, SSA and g
	'''
	if 'debug' in args:
		debug = args['debug']
	else:
		debug = False
	if 'save' in args:
		save = args['save']
	else:
		save = False
	if 'save_path' in args:
		save_path = args['save_path']
	else:
		save_path = 'figure/MS/'
	
	left = 0.2
	bottom = 0.1
	right = 0.9
	top = 0.99
	wspace = 0.145
	hspace = 0.11
	dw = 0.146
	dh = 0.03
	dy = 0.01
	lw = 3
	fs = 12
	fw = 'bold'
	font = dict(weight=fw, size=fs)
	
	# plot AOD
	
	fig = plt.figure(figsize=(8,6))
	plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
	
	axs = plt.axes([left, bottom+(top-bottom)/2+dy, right-left, (top-bottom)/2-dy])
	ax = plt.pcolormesh(kappa, kappaI, AOD.T*1000, cmap='coolwarm')
	cb = plt.colorbar(ax, extend='both', pad=0.05)
	cb.ax.tick_params(labelsize=fs)
	cb.set_label('AOD, $\\times10^{-3}$', fontdict=font)
	plt.xticks([], [])
	plt.yticks(fontsize=fs, fontweight=fw)
	plt.ylabel('Kappa inhomogeneity', fontsize=fs, fontweight=fw)
	plt.text(0.02, 0.9, '(a)', fontsize=fs, fontweight=fw, fontstyle='italic', transform=axs.transAxes)
	
	axs = plt.axes([left, bottom, right-left, (top-bottom)/2-dy])
	ax = plt.pcolormesh(kappa, kappaI, std_AOD.T*100, cmap='coolwarm', vmax = 1.5)
	cb = plt.colorbar(ax, extend='both', pad=0.05)
	cb.ax.tick_params(labelsize=fs)
	#cb.set_ticks([0.6,0.7,0.8,0.9])
	#cb.set_ticklabels(['0.6','0.7','0.8','0.9'])
	cb.set_label('$Std_{AOD}$, %', fontdict=font)
	plt.xticks(fontsize=fs, fontweight=fw)
	plt.xlabel('Kappa', fontsize=fs, fontweight=fw)
	plt.yticks(fontsize=fs, fontweight=fw)
	plt.ylabel('Kappa inhomogeneity', fontsize=fs, fontweight=fw)
	plt.text(0.02, 0.9, '(b)', fontsize=fs, fontweight=fw, fontstyle='italic', transform=axs.transAxes)
	
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_kappaI_AOD.pdf')
	else:
		plt.show()
	plt.close()
	
	# plot SSA
	
	fig = plt.figure(figsize=(8,6))
	plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
	
	axs = plt.axes([left, bottom+(top-bottom)/2+dy, right-left, (top-bottom)/2-dy])
	ax = plt.pcolormesh(kappa, kappaI, SSA.T, cmap='coolwarm')
	cb = plt.colorbar(ax, extend='both', pad=0.05)
	cb.ax.tick_params(labelsize=fs)
	cb.set_label('SSA', fontdict=font)
	plt.xticks([], [])
	plt.yticks(fontsize=fs, fontweight=fw)
	plt.ylabel('Kappa inhomogeneity', fontsize=fs, fontweight=fw)
	plt.text(0.02, 0.9, '(a)', fontsize=fs, fontweight=fw, fontstyle='italic', transform=axs.transAxes)
	
	axs = plt.axes([left, bottom, right-left, (top-bottom)/2-dy])
	ax = plt.pcolormesh(kappa, kappaI, std_SSA.T*100, cmap='coolwarm', vmax=0.11)
	cb = plt.colorbar(ax, extend='both', pad=0.05)
	cb.ax.tick_params(labelsize=fs)
	#cb.set_ticks([0.6,0.7,0.8,0.9])
	#cb.set_ticklabels(['0.6','0.7','0.8','0.9'])
	cb.set_label('$Std_{SSA}$, %', fontdict=font)
	plt.xticks(fontsize=fs, fontweight=fw)
	plt.xlabel('Kappa', fontsize=fs, fontweight=fw)
	plt.yticks(fontsize=fs, fontweight=fw)
	plt.ylabel('Kappa inhomogeneity', fontsize=fs, fontweight=fw)
	plt.text(0.02, 0.9, '(b)', fontsize=fs, fontweight=fw, fontstyle='italic', transform=axs.transAxes)
	
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_kappaI_SSA.pdf')
	else:
		plt.show()
	plt.close()
	
	# plot g
	
	fig = plt.figure(figsize=(8,6))
	plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
	
	axs = plt.axes([left, bottom+(top-bottom)/2+dy, right-left, (top-bottom)/2-dy])
	ax = plt.pcolormesh(kappa, kappaI, g.T, cmap='coolwarm')
	cb = plt.colorbar(ax, extend='both', pad=0.05)
	cb.ax.tick_params(labelsize=fs)
	cb.set_label('G', fontdict=font)
	plt.xticks([], [])
	plt.yticks(fontsize=fs, fontweight=fw)
	plt.ylabel('Kappa inhomogeneity', fontsize=fs, fontweight=fw)
	plt.text(0.02, 0.9, '(a)', fontsize=fs, fontweight=fw, fontstyle='italic', transform=axs.transAxes)
	
	axs = plt.axes([left, bottom, right-left, (top-bottom)/2-dy])
	ax = plt.pcolormesh(kappa, kappaI, std_g.T*100, cmap='coolwarm', vmin=0.55, vmax=0.95)
	cb = plt.colorbar(ax, extend='both', pad=0.05)
	cb.ax.tick_params(labelsize=fs)
	cb.set_ticks([0.6,0.7,0.8,0.9])
	cb.set_ticklabels(['0.6','0.7','0.8','0.9'])
	cb.set_label('$Std_g$, %', fontdict=font)
	plt.xticks(fontsize=fs, fontweight=fw)
	plt.xlabel('Kappa', fontsize=fs, fontweight=fw)
	plt.yticks(fontsize=fs, fontweight=fw)
	plt.ylabel('Kappa inhomogeneity', fontsize=fs, fontweight=fw)
	plt.text(0.02, 0.9, '(b)', fontsize=fs, fontweight=fw, fontstyle='italic', transform=axs.transAxes)
	
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_kappaI_g.pdf')
	else:
		plt.show()
	plt.close()

def plot_nI(n, nI, AOD, SSA, g, std_AOD, std_SSA, std_g, **args):
	'''
	This function is to plot n mixing state effect,
	it's a two dimentional figure, x axis is n and y axis is nI, 
	whith aerosol optical depth AOD, single scattering albedo SSA or 
	asymmetry factor g color map, for that is the most influenced properties
	input:
		n :		n data, array, float
		nI :		n mixing state slope data, array, float
		AOD :		aerosol optical depth data, 2-d array, float
		SSA :		single scattering albedo data, 2-d array, float
		g :		asymmetry factor data, 2-d array, float
		std_AOD :	standard aerosol optical depth data for 30% n change, 
				2-d array, float
		std_SSA :	standard single scattering albedo data for 30% n change, 
				2-d array, float
		std_g :	standard asymmetry factor data for 30% n change, 2-d array, 					float
		**debug :	debug flag, bool, default False
		**save :	whether to save figure, bool, default False
		**save_path :	figure save path, string, default 'figure/MS/'
	output:
		figure 1 :	x n, y nI and colored AOD, SSA and g
		figure 2 :	x n, y nI and colored standard deviation AOD, SSA and g
	'''
	if 'debug' in args:
		debug = args['debug']
	else:
		debug = False
	if 'save' in args:
		save = args['save']
	else:
		save = False
	if 'save_path' in args:
		save_path = args['save_path']
	else:
		save_path = 'figure/MS/'
	
	left = 0.2
	bottom = 0.1
	right = 0.9
	top = 0.99
	wspace = 0.145
	hspace = 0.11
	dw = 0.146
	dh = 0.03
	dy = 0.01
	lw = 3
	fs = 12
	fw = 'bold'
	font = dict(weight=fw, size=fs)
	
	# plot AOD
	
	fig = plt.figure(figsize=(8,6))
	plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
	
	axs = plt.axes([left, bottom+(top-bottom)/2+dy, right-left, (top-bottom)/2-dy])
	ax = plt.pcolormesh(n, nI, AOD.T*1000, cmap='coolwarm')
	cb = plt.colorbar(ax, extend='both', pad=0.05)
	cb.ax.tick_params(labelsize=fs)
	cb.set_label('AOD, $\\times10^{-3}$', fontdict=font)
	plt.xticks([], [])
	plt.yticks([-0.02,-0.01,0,0.01,0.02], [-0.02,-0.01,0,0.01,0.02], fontsize=fs, fontweight=fw)
	plt.ylabel('N inhomogeneity', fontsize=fs, fontweight=fw)
	plt.text(0.02, 0.9, '(a)', fontsize=fs, fontweight=fw, fontstyle='italic', transform=axs.transAxes)
	
	axs = plt.axes([left, bottom, right-left, (top-bottom)/2-dy])
	ax = plt.pcolormesh(n, nI, std_AOD.T*100, cmap='coolwarm')
	cb = plt.colorbar(ax, extend='both', pad=0.05)
	cb.ax.tick_params(labelsize=fs)
	#cb.set_ticks([0.6,0.7,0.8,0.9])
	#cb.set_ticklabels(['0.6','0.7','0.8','0.9'])
	cb.set_label('$Std_{AOD}$, %', fontdict=font)
	plt.xticks([0.14,0.15,0.16,0.17, 0.1795], [0.14,0.15,0.16,0.17,0.18], fontsize=fs, fontweight=fw)
	plt.xlabel('N', fontsize=fs, fontweight=fw)
	plt.yticks([-0.02,-0.01,0,0.01,0.02], [-0.02,-0.01,0,0.01,0.02], fontsize=fs, fontweight=fw)
	plt.ylabel('N inhomogeneity', fontsize=fs, fontweight=fw)
	plt.text(0.02, 0.9, '(b)', fontsize=fs, fontweight=fw, fontstyle='italic', transform=axs.transAxes)
	
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_nI_AOD.pdf')
	else:
		plt.show()
	plt.close()
	
	# plot SSA
	
	fig = plt.figure(figsize=(8,6))
	plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
	
	axs = plt.axes([left, bottom+(top-bottom)/2+dy, right-left, (top-bottom)/2-dy])
	ax = plt.pcolormesh(n, nI, SSA.T, cmap='coolwarm')
	cb = plt.colorbar(ax, extend='both', pad=0.05)
	cb.ax.tick_params(labelsize=fs)
	cb.set_label('SSA', fontdict=font)
	plt.xticks([], [])
	plt.yticks([-0.02,-0.01,0,0.01,0.02], [-0.02,-0.01,0,0.01,0.02], fontsize=fs, fontweight=fw)
	plt.ylabel('N inhomogeneity', fontsize=fs, fontweight=fw)
	plt.text(0.02, 0.9, '(a)', fontsize=fs, fontweight=fw, fontstyle='italic', transform=axs.transAxes)
	
	axs = plt.axes([left, bottom, right-left, (top-bottom)/2-dy])
	ax = plt.pcolormesh(n, nI, std_SSA.T*100, cmap='coolwarm')
	cb = plt.colorbar(ax, extend='both', pad=0.05)
	cb.ax.tick_params(labelsize=fs)
	#cb.set_ticks([0.6,0.7,0.8,0.9])
	#cb.set_ticklabels(['0.6','0.7','0.8','0.9'])
	cb.set_label('$Std_{SSA}$, %', fontdict=font)
	plt.xticks([0.14,0.15,0.16,0.17, 0.1795], [0.14,0.15,0.16,0.17,0.18], fontsize=fs, fontweight=fw)
	plt.xlabel('N', fontsize=fs, fontweight=fw)
	plt.yticks([-0.02,-0.01,0,0.01,0.02], [-0.02,-0.01,0,0.01,0.02], fontsize=fs, fontweight=fw)
	plt.ylabel('N inhomogeneity', fontsize=fs, fontweight=fw)
	plt.text(0.02, 0.9, '(b)', fontsize=fs, fontweight=fw, fontstyle='italic', transform=axs.transAxes)
	
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_nI_SSA.pdf')
	else:
		plt.show()
	plt.close()
	
	# plot g
	
	fig = plt.figure(figsize=(8,6))
	plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
	
	axs = plt.axes([left, bottom+(top-bottom)/2+dy, right-left, (top-bottom)/2-dy])
	ax = plt.pcolormesh(n, nI, g.T, cmap='coolwarm')
	cb = plt.colorbar(ax, extend='both', pad=0.05)
	cb.ax.tick_params(labelsize=fs)
	cb.set_label('G', fontdict=font)
	plt.xticks([], [])
	plt.yticks([-0.02,-0.01,0,0.01,0.02], [-0.02,-0.01,0,0.01,0.02], fontsize=fs, fontweight=fw)
	plt.ylabel('N inhomogeneity', fontsize=fs, fontweight=fw)
	plt.text(0.02, 0.9, '(a)', fontsize=fs, fontweight=fw, fontstyle='italic', transform=axs.transAxes)
	
	axs = plt.axes([left, bottom, right-left, (top-bottom)/2-dy])
	ax = plt.pcolormesh(n, nI, std_g.T*100, cmap='coolwarm')
	cb = plt.colorbar(ax, extend='both', pad=0.05)
	cb.ax.tick_params(labelsize=fs)
	#cb.set_ticks([0.6,0.7,0.8,0.9])
	#cb.set_ticklabels(['0.6','0.7','0.8','0.9'])
	cb.set_label('$Std_g$, %', fontdict=font)
	plt.xticks([0.14,0.15,0.16,0.17, 0.1795], [0.14,0.15,0.16,0.17,0.18], fontsize=fs, fontweight=fw)
	plt.xlabel('N', fontsize=fs, fontweight=fw)
	plt.yticks([-0.02,-0.01,0,0.01,0.02], [-0.02,-0.01,0,0.01,0.02], fontsize=fs, fontweight=fw)
	plt.ylabel('N inhomogeneity', fontsize=fs, fontweight=fw)
	plt.text(0.02, 0.9, '(b)', fontsize=fs, fontweight=fw, fontstyle='italic', transform=axs.transAxes)
	
	if save:
		plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_nI_g.pdf')
	else:
		plt.show()
	plt.close()

if __name__ == '__main__':
	'''
	kappa, kappaI, AOD, SSA, g, std_AOD, std_SSA, std_g = cal_kappaI(0.15, 0.25, 0.13, debug=True, dx=0.1/100, dy=0.2/50)
	infos = dict(kappa=kappa, kappaI=kappaI, AOD=AOD, SSA=SSA, g=g, std_AOD=std_AOD, std_SSA=std_SSA, std_g=std_g)
	np.save('output/MS/infos_kappa.npy', infos)
	
	n, nI, AOD, SSA, g, std_AOD, std_SSA, std_g = cal_nI(0.14, 0.18, 0.03, debug=True, dx=0.04/100, dy=0.03/50)
	infos = dict(n=n, nI=nI, AOD=AOD, SSA=SSA, g=g, std_AOD=std_AOD, std_SSA=std_SSA, std_g=std_g)
	np.save('output/MS/infos_n.npy', infos)
	'''
	data = np.load('output/MS/infos_n.npy', allow_pickle=True).item()
	n = data['n']
	nI = data['nI']
	AOD = data['AOD']
	SSA = data['SSA']
	g = data['g']
	std_AOD = data['std_AOD']
	std_SSA = data['std_SSA']
	std_g = data['std_g']
	
	plot_nI(n, nI, AOD, SSA, g, std_AOD, std_SSA, std_g, save=False)
	'''
	data = np.load('output/MS/infos.npy', allow_pickle=True).item()
	kappa = data['kappa']
	kappaI = data['kappaI']
	AOD = data['AOD']
	SSA = data['SSA']
	g = data['g']
	std_AOD = data['std_AOD']
	std_SSA = data['std_SSA']
	std_g = data['std_g']
	
	plot_kappaI(kappa, kappaI, AOD, SSA, g, std_AOD, std_SSA, std_g, save=False)
	'''
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
