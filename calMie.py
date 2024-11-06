####################################################################################
# INTRODUCTION:
# This code is to do some Mie theory calculating
# Created by Hebs at 22/9/7/12:37
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import sys

import calNI
import calPNSD
import phaseFunc

def read_rate(parameter, path):
	'''
	This function is to read change rate from infos.npy
	input:
		parameter : parameter name, string
		path      : infos.npy store path, string
	output:
		rate      : parameter change rate, in float, list
	'''
	name = parameter + '_rate'
	data = np.load(path+'infos.npy', allow_pickle=True)
	if name not in data[0]:
		print('No such parameter ' + parameter + '. Please check')
		sys.exit()
	rate = np.zeros(len(data))
	for i in range(len(data)):
		rate[i] = data[i][name]
	return rate

def cal_PDF(data, num):
	'''
	This funciton is to calculating data PDF
	input:
		data   : data to statistic, numpy array
		num    : bin num, int
	output:
		x      : bin, numpy array
		PDF    : posibility distribution function, numpy array
	'''
	data_min = np.nanmin(data)
	data_max = np.nanmax(data)
	x = np.arange(data_min, data_max, (data_max-data_min)/num)
	PDF = np.zeros(len(x))
	for y in data:
		i = 0
		while i<num-1:
			if y>=x[i] and y<x[i+1]:
				PDF[i] = PDF[i] + 1
				break
			i = i + 1
		if y==data_max:
			PDF[-1] = PDF[-1] + 1
	return x, PDF

def cal_lim(array, **args):
	'''
	This function is to calculate array plot lim
	input:
		array   : origin data
		**upper : upper blank rate left for plot, default 0.1
		**down  : down blank rate left for plot, default 0.1
	output:
		lim     : plot plt.ylim()
	'''
	if 'upper' in args:
		upper = args['upper']
	else:
		upper = 0.1
	if 'down' in args:
		down = args['down']
	else:
		down = 0.1
	mn = np.min(array)
	mx = np.max(array)
	up = (mx-mn) * upper
	dn = (mx-mn) * down
	lim = [mn-dn, mx+up]
	return lim

def Mie2(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, nShellH1, kappa, kappaH1, MS, MSH1, MSH2x, MSH2y,dz, RH, wl, CT, CTH1, BCAE, **args):
	'''
	This function is to use Mie modal to calculate aerosol parameters AOD, SSA and g
	input:
		Dps      : diameter distribution, array, nm
		PNSD     : number distribution, array, dn/dlogDps, cm^-3
		DBCps    : BC particle diameter size, array, nm
		DBC      : BC core diameter size, array, nm
		BCPNSD   : BC particle number concentration, array, dn/dlogDBCps/dlogDBC, cm^-3
		kBC      : BC core imagine part of complex refractive index, float
		nBC      : BC core real part of complex refractive index, float
		n        : shell complex refractive index, float
		nH1      : n mixing state change rate by diameter, float
		kappa    : hygroscopicity parameter, float
		kappaH1  : hygroscopicity parameter mixing state change rate by diameter, float
		MS       : mixing state, float
		MSH1     : mixing state change rate by diameter, float
		MSH2x    : mixing state bin peak number, int
		MSH2y    : mixing state bin peak saperate rate, float
		dz       : thickness, meter, float
		RH       : relative humidity, float, percent
		wl       : wave length, nm, float
		CT       : BC core coating thickness adjust parameter, float
		CTH1     : BC core coating thickness change rate by diameter, float
		BCAE     : BC core absorbing enhancement adjust parameter, float
		**args:
			debug				: debug flag, bool, default False
	output:
		aerosol optical parameters
		AOD		: aerosol optical depth, np.array, in shape (len(wl))
		SSA		: single scattering albedo, np.array, in shape (len(wl))
		g       : asymmetry factor, np.array, in shape (len(wl))
	'''
	if 'debug' in args:
		debug = args['debug']
	else:
		debug = False
	
	if debug:
		print('calculating...')
	
	AOD = np.zeros(len(wl))
	SSA = np.zeros(len(wl))
	g = np.zeros(len(wl))
	
	for i in range(len(wl)):
		kext_i, ksca_i, g[i] = calNI.cal_Mie3(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, nShellH1, kappa, kappaH1, MS, MSH1, MSH2x, MSH2y, RH, wl[i], CT, CTH1, BCAE)
		AOD[i] = kext_i * dz * 1e-6
		SSA[i] = ksca_i / kext_i
		if debug:
			print('AOD:', AOD)
			print('SSA:', SSA)
			print('g:', g)
	
	if debug:
		print('done')
	
	return AOD, SSA, g

def Mie(Dps, PNSD, DBCps, DBC, BCPNSD, nBC, kBC, nShell, nI1, nI2x, nI2y, kappa, kappaI1, kappaI2x, kappaI2y, RH, wl, z, dz, BCAE, CT, VD, PF, **args):
	'''
	This function is to use Mie modal to calculate Mie parameters,
	parameters including: 
		AOD
		SSA
		g
	input:
		Dps          : diameter distribution, array, nm
		PNSD         : number distribution, array, dn/dlogDp
		DBCps        : BC particle diameter size, array, nm
		DBC          : BC core diameter size, array, nm
		BCPNSD       : BC particle number size distribution, array, dn/dlogDBC
		nBC          : BC core real part of complex refractive index
		kBC          : BC core imagine part of complex refractive index
		nShell       : shell complex refractive index
		nI1          : n mixing state change rate by diameter
		nI2x         : n mixing state size bin peak number
		nI2y         : n mixing state size bin peak separate rate
		kappa        : hygroscopicity parameter, float
		kappaI1      : hygroscopicity parameter mixing state change rate by diameter
		kappaI2x     : hygroscopicity parameter mixing state size bin peak number
		kappaI2y     : hygroscopicity parameter mixing state size bin peak separate rate
		RH           : raletive humidity, float, percent
		wl           : wave length, float, nm
		z            : height, float, m
		dz           : thickness, float, m
		BCAE         : BC absorbing enhancement adjustment, float
		CT           : BC coating thickness adjustment, float
		VD           : vertical distribution mixing, 0 for type A and 1 for type B
		PF           : phase function adjustment, float
		**debug      : debug flag, default False, bool
		**ddz        : AOD calculation step, default dz/100
	output:
		AOD    	is the optical depth increment within level i at wavelength k,
				information is specified in top-down order. 

		SSA    	is the single scattering albedo

		g      	are legendre moments of the phase function.
				Note that zeroeth moment is not read, it is assumed to be 1.
	'''
	if 'debug' in args:
		debug = args['debug']
	else:
		debug = False
	if 'ddz' in args:
		ddz = args['ddz']
	else:
		ddz = dz / 10
	
	if debug:
		print('calculating...')
	
	# Before calculating, some factors change the input parameters, 
	# including: coating thickness, mixing state, BC mixing state, BC absorbing enhancement
	
	kext, ksca, g = calNI.cal_Mie(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, nI1, nI2x, nI2y, kappa, kappaI1, kappaI2x, kappaI2y, RH, wl, CT, BCAE)
	g = g * PF
	n = round(dz/ddz)
	H_integral = 0
	for i in range(n):
		H_integral += calPNSD.cal_VD(Dps, PNSD, z+n*ddz, VD) * ddz
	AOD = kext * 1e-6 * H_integral
	SSA = ksca / kext
	if debug:
		print('AOD:', AOD)
		print('SSA:', SSA)
		print('g:', g)
		print('done')
	
	return AOD, SSA, g

def cal_AOD(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, kappa, RH, wl, z, dz, CT, BCAE, VD, **args):
	'''
	This function is to calculate AOD for a small layer of aerosols
	input:
		Dps     : diameter distribution, array, nm
		PNSD    : number distribution, array, dn/dlogDp
		DBCps   : BC particle diameter size, array, nm
		DBC     : BC core diameter size, array, nm
		BCPNSD  : BC particle number concentration, array, cm^-3
		kBC     : BC core imagine part of complex refractive index
		nBC     : BC core real part of complex refractive index
		nShell  : shell complex refractive index
		kappa   : hygroscopicity parameter, float
		RH      : relative humidity, float, percent
		wl      : wave length, nm
		z       : height, m
		dz      : delta height, m
		CT      : BC core coating thickness adjust parameter
		BCAE    : BC core absorbing enhancement adjust parameter
		VD      : vertical distribution mixing, 0 for type A and 1 for type B
		**ddz   : step, default dz/10
	output:
		dtau    : AOD, aerosol optical depth
	'''
	if 'ddz' in args:
		ddz = args['ddz']
	else:
		ddz = 10
	
	kext = calPNSD.cal_kext(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, kappa, RH, wl, CT, BCAE)
	H_integral = 0
	for i in range(ddz):
		ratio = calPNSD.cal_VD(Dps, PNSD, z+i*dz/ddz, VD)
		H_integral += dz / ddz * ratio
	
	dtau = kext * 1e-6 * H_integral
	
	return dtau

def cal_SSA(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, kappa, RH, wl, CT, BCAE):
	'''
	This function is to calculate SSA for a small layer of aerosols
	input:
		Dps     : diameter distribution, array, nm
		PNSD    : number distribution, array, dn/dlogDp
		DBCps   : BC particle diameter size, array, nm
		DBC     : BC core diameter size, array, nm
		BCPNSD  : BC particle number concentration, array, cm^-3
		kBC     : BC core imagine part of complex refractive index
		nBC     : BC core real part of complex refractive index
		nShell  : shell complex refractive index
		kappa   : hygroscopicity parameter, float
		RH      : relative humidity, float, percent
		wl      : wave length, nm
		CT      : BC core coating thickness adjust parameter
		BCAE    : BC core absorbing enhancement adjust parameter
	output:
		waer    : SSA, single scattering albedo
	'''
	kext = calPNSD.cal_kext(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, kappa, RH, wl, CT, BCAE)
	ksca = calPNSD.cal_ksca(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, kappa, RH, wl, CT)
	waer = ksca / kext
	
	return waer

def cal_g(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, kappa, RH, wl, CT, PF):
	'''
	This function is to calculate asymmetry factor in certain height and raletive humidity
	input:
		Dps     : diameter distribution, array, nm
		PNSD    : number distribution, array, dn/dlogDp
		DBCps   : BC particle diameter size, array, nm
		DBC     : BC core diameter size, array, nm
		BCPNSD  : BC particle number concentration, array, cm^-3
		kBC     : BC core imagine part of complex refractive index
		nBC     : BC core real part of complex refractive index
		nShell  : shell complex refractive index
		kappa   : hygroscopicity parameter, float
		RH      : relative humidity, float, percent
		wl      : wave length, nm
		PF      : phase function change rate
	output:
		g       : asymmetry factor
	'''
	g = calPNSD.cal_g(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, kappa, RH, wl, CT)
	g = g * PF
	
	return g
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
