####################################################################################
# INTRODUCTION:
# This code is to test all factor turburlence sensitivity of Mie parameters: AOD, SSA and g
# Created by Hebs at 22/8/30/12:56
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import sys
import os
import random
from pyDOE import lhs
from scipy.stats import norm

import write_atms
import readAtms
import readTaizhou
import calRH
import calBC
import calMie
import calPNSD

###############################################
# some constants
default_mag = 0.3 # Dps resolution adjust magnification
default_range = 0.3 # 30% turbulence
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

main_paras = [] # main factors
main_paras.append('n')
main_paras.append('MS')
main_paras.append('CT')
main_paras.append('kappa')
###############################################

def LHS_norm(num, loc, scale):
	'''
	This function is Latin hypercubic sampling method code for normal distribution
	input:
		num     : sample number, int
		loc     : normal distribution location parameter, float
		scale   : normal distribution scale parameter, float
	output:
		lhd     : parameter list, array
	'''
	lhd = lhs(1, samples=num)
	lhd = norm(loc=loc, scale=scale).ppf(lhd)
	return lhd

def EFAST_norm(num, loc, scale, freq):
	'''
	This function is extended Fourier Amplitude Sensitivity Test method code for normal distribution
	input:
		num     : sample number, int
		loc     : normal distribution location parameter, float
		scale   : normal distribution scale parameter, float
		freq    : sin signal frequency, float
	output:
		x       : parameter list, array
	'''
	x = np.arange(num*1.0)
	
	phi = random.uniform(-np.pi, np.pi)
	for i in range(len(x)):
		omega = 2 * np.pi * freq / num
		x[i] = 0.5 + 1 / np.pi * np.arcsin(np.sin(omega*x[i]+phi))
		x[i] = norm(loc=loc, scale=scale).ppf(x[i])
	return x

def run(dz, wl, RH, **args):
	'''
	This function is the main function of this code, using Monte Carlo method to test all factors' sensitivity of Mie parameters
	Factors including:
		n for real part of complex refractive index
		nH1 for real part of complex refractive index mixing 1 level heterogeneity
		nBC for BC real part of complex refractive index
		kBC for BC imagine part of complex refractive index
		PNSD for number distribution dirty part proportion
		MS for mixing state external part proportion
		MSH1 for mixing state external part proportion 1 level heterogeneity
		MSH2 for mixing state external part proportion 2 level heterogeneity
		VD for vertical distribution type B part proportion
		CT for coating thickness
		CTH1 for coating thickness 1 level heterogeneity
		kappa for hydroscopicity parameter
		kappaH1 for kappa mixing 1 level heterogeneity
		rhoBC for BC density
		BCPNSD for BC number distribution dirty part proportion, in cm-3
		BCAE for BC absorbing enhancement
	input:
		dz    : thickness, in m
		wl    : the wavelength [ wl(k) < wl(k+1) ], array, nm
		RH    : relative humidity, percent, float
		**args including:
			debug, bool
			data_path, string
			Dps, np.array
			clean_PNSD, np.array
			dirty_PNSD, np.array
			clean_BCPNSD, np.array
			dirty_BCPNSD, np.array
			DBC, np.array
			DBCps, np.array
	output:
		infos:
			aerosol optical parameter information, dict
			including:
				AOD   : aerosol optical depth
				SSA   : single scattering albedo
				g     : asymmetry factor
			in form of: Mie_infos=[dtau[:,0],waer[:,0],pmom[:,0,:]],
			and parameter turburlence information, including all parameters' change rate
		paras:
			origin parameter value, dict
	'''
	if 'n_rate' in args:
		n_rate_on = True
		n_rate = args['n_rate']
	else:
		n_rate_on = False
	if 'n_on' in args:
		n_on = args['n_on']
		if 'n_range' in args:
			n_range = args['n_range']
		else:
			n_range = default_range
	else:
		n_on = False
	if 'nH1_rate' in args:
		nH1_rate_on = True
		nH1_rate = args['nH1_rate']
	else:
		nH1_rate_on = False
	if 'nH1_on' in args:
		nH1_on = args['nH1_on']
		if 'nH1_range' in args:
			nH1_range = args['nH1_range']
		else:
			nH1_range = default_range
	else:
		nH1_on = False
	
	if 'nBC_rate' in args:
		nBC_rate_on = True
		nBC_rate = args['nBC_rate']
	else:
		nBC_rate_on = False
	if 'nBC_on' in args:
		nBC_on = args['nBC_on']
		if 'nBC_range' in args:
			nBC_range = args['nBC_range']
		else:
			nBC_range = default_range
	else:
		nBC_on = False
	if 'kBC_rate' in args:
		kBC_rate_on = True
		kBC_rate = args['kBC_rate']
	else:
		kBC_rate_on = False
	if 'kBC_on' in args:
		kBC_on = args['kBC_on']
		if 'kBC_range' in args:
			kBC_range = args['kBC_range']
		else:
			kBC_range = default_range
	else:
		kBC_on = False
	if 'PNSD_rate' in args:
		PNSD_rate_on = True
		PNSD_rate = args['PNSD_rate']
	else:
		PNSD_rate_on = False
	if 'PNSD_on' in args:
		PNSD_on = args['PNSD_on']
		if 'PNSD_range' in args:
			PNSD_range = args['PNSD_range']
		else:
			PNSD_range = default_range
	else:
		PNSD_on = False
	if 'MS_rate' in args:
		MS_rate_on = True
		MS_rate = args['MS_rate']
	else:
		MS_rate_on = False
	if 'MS_on' in args:
		MS_on = args['MS_on']
		if 'MS_range' in args:
			MS_range = args['MS_range']
		else:
			MS_range = default_range
	else:
		MS_on = False
	
	if 'MSH1_rate' in args:
		MSH1_rate_on = True
		MSH1_rate = args['MSH1_rate']
	else:
		MSH1_rate_on = False
	if 'MSH1_on' in args:
		MSH1_on = args['MSH1_on']
		if 'MSH1_range' in args:
			MSH1_range = args['MSH1_range']
		else:
			MSH1_range = default_range
	else:
		MSH1_on = False
	
	if 'MSH2x' in args:
		MSH2x = args['MSH2x']
	else:
		MSH2x = 1
	
	if 'MSH2_rate' in args:
		MSH2_rate_on = True
		MSH2_rate = args['MSH2_rate']
	else:
		MSH2_rate_on = False
	if 'MSH2_on' in args:
		MSH2_on = args['MSH2_on']
		if 'MSH2_range' in args:
			MSH2_range = args['MSH2_range']
		else:
			MSH2_range = default_range
	else:
		MSH2_on = False
	
	if 'CT_rate' in args:
		CT_rate_on = True
		CT_rate = args['CT_rate']
	else:
		CT_rate_on = False
	if 'CT_on' in args:
		CT_on = args['CT_on']
		if 'CT_range' in args:
			CT_range = args['CT_range']
		else:
			CT_range = default_range
	else:
		CT_on = False
	
	if 'CTH1_rate' in args:
		CTH1_rate_on = True
		CTH1_rate = args['CTH1_rate']
	else:
		CTH1_rate_on = False
	if 'CTH1_on' in args:
		CTH1_on = args['CTH1_on']
		if 'CTH1_range' in args:
			CTH1_range = args['CTH1_range']
		else:
			CTH1_range = default_range
	else:
		CTH1_on = False
	
	if 'kappa_rate' in args:
		kappa_rate_on = True
		kappa_rate = args['kappa_rate']
	else:
		kappa_rate_on = False
	if 'kappa_on' in args:
		kappa_on = args['kappa_on']
		if 'kappa_range' in args:
			kappa_range = args['kappa_range']
		else:
			kappa_range = default_range
	else:
		kappa_on = False
	
	if 'kappaH1_rate' in args:
		kappaH1_rate_on = True
		kappaH1_rate = args['kappaH1_rate']
	else:
		kappaH1_rate_on = False
	if 'kappaH1_on' in args:
		kappaH1_on = args['kappaH1_on']
		if 'kappaH1_range' in args:
			kappaH1_range = args['kappaH1_range']
		else:
			kappaH1_range = default_range
	else:
		kappaH1_on = False
	
	if 'rhoBC_rate' in args:
		rhoBC_rate_on = True
		rhoBC_rate = args['rhoBC_rate']
	else:
		rhoBC_rate_on = False
	if 'rhoBC_on' in args:
		rhoBC_on = args['rhoBC_on']
		if 'rhoBC_range' in args:
			rhoBC_range = args['rhoBC_range']
		else:
			rhoBC_range = default_range
	else:
		rhoBC_on = False
	
	if 'BCPNSD_rate' in args:
		BCPNSD_rate_on = True
		BCPNSD_rate = args['BCPNSD_rate']
	else:
		BCPNSD_rate_on = False
	if 'BCPNSD_on' in args:
		BCPNSD_on = args['BCPNSD_on']
		if 'BCPNSD_range' in args:
			BCPNSD_range = args['BCPNSD_range']
		else:
			BCPNSD_range = default_range
	else:
		BCPNSD_on = False
	
	if 'BCAE_rate' in args:
		BCAE_rate_on = True
		BCAE_rate = args['BCAE_rate']
	else:
		BCAE_rate_on = False
	if 'BCAE_on' in args:
		BCAE_on = args['BCAE_on']
		if 'BCAE_range' in args:
			BCAE_range = args['BCAE_range']
		else:
			BCAE_range = default_range
	else:
		BCAE_on = False
	
	if 'debug' in args:
		debug = args['debug']
	else:
		debug = False
	if 'data_path' in args:
		data_path = args['data_path']
	else:
		data_path = 'data/'
	
	if 'Dps' in args:
		Dps_exist = True
		Dps = args['Dps']
	else:
		Dps_exist = False
	if 'clean_PNSD' in args:
		clean_PNSD_exist = True
		clean_PNSD = args['clean_PNSD']
	else:
		clean_PNSD_exist = False
	if 'dirty_PNSD' in args:
		dirty_PNSD_exist = True
		dirty_PNSD = args['dirty_PNSD']
	else:
		dirty_PNSD_exist = False
	if 'DBC' in args:
		DBC_exist = True
		DBC = args['DBC']
	else:
		DBC_exist = False
	if 'DBCps' in args:
		DBCps_exist = True
		DBCps = args['DBCps']
	else:
		DBCps_exist = False
	if 'clean_BCPNSD' in args:
		clean_BCPNSD_exist = True
		clean_BCPNSD = args['clean_BCPNSD']
	else:
		clean_BCPNSD_exist = False
	if 'dirty_BCPNSD' in args:
		dirty_BCPNSD_exist = True
		dirty_BCPNSD = args['dirty_BCPNSD']
	else:
		dirty_BCPNSD_exist = False
	if 'DBC2' in args:
		DBC2_exist = True
		DBC2 = args['DBC2']
	else:
		DBC2_exist = False
	if 'clean_BCPMSD' in args:
		clean_BCPMSD_exist = True
		clean_BCPMSD = args['clean_BCPMSD']
	else:
		clean_BCPMSD_exist = False
	if 'dirty_BCPMSD' in args:
		dirty_BCPMSD_exist = True
		dirty_BCPMSD = args['dirty_BCPMSD']
	else:
		dirty_BCPMSD_exist = False
	
	# read data, including Dps, PNSD, DBC, DBCps and BCPNSD data
	##################################################################################
	
	if not Dps_exist or not clean_PNSD_exist or not dirty_PNSD_exist or not DBCps_exist or not clean_BCPNSD_exist or not dirty_BCPNSD_exist:
		BCPNSD_exist = False
		fn = data_path + 'sp2/Taizhou.npy'
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
			
			#######################################################
			# need more discuss
			clean_PNSD = np.nanmean(PNSD, axis=0) / 10 # cm-3/dlnDps
			dirty_PNSD = np.nanmean(PNSD, axis=0)
			clean_BCPNSD = np.nanmean(BCPNSD, axis=0) / 10 # cm-3/dlogDBC
			dirty_BCPNSD = np.nanmean(BCPNSD, axis=0)
			#######################################################
	
	else:
		BCPNSD_exist = True
	
	if debug:
		print('particle information:')
		#print('Dps:', Dps)
		#print('clean PNSD:', clean_PNSD)
		#print('dirty PNSD:', dirty_PNSD)
		#print('DBCps:', DBCps)
		#print('DBC:', DBC)
		#print('clean BCPNSD:', clean_BCPNSD)
		#print('dirty BCPNSD:', dirty_BCPNSD)
		print('done')
	
	# after this, DBC, DBCps and BCPNSD data is certained
	
	# read mixing state. If not exist, using default data
	
	MS = 0.7 # from Size distribution and mixing state of blank carbon particles during a heavy air pollution episode in Shanghai, X.Gong et al., 2016
	if debug:
		print('mixing state rate :', MS)
	MSH1 = 0.5
	MSH2 = 0.5
	
	BCPNSD = clean_BCPNSD * 0.5 + dirty_BCPNSD * 0.5
	
	# read BC density. If not exist, using default data
	
	rhoBC = 0.95 # g/m3, from Determination of the refractive index of ambient aerosols, G.Zhao et al., 2020
	if debug:
		print('rho BC :', rhoBC, 'g/cm3')
	
	# BC density decides BC complex refractive index
	# mEC = 2.26 + 1.26i, mAir = 1 + 0i, rhoEC = 1.8, rhoAir = 1e-3 -> 0
	nBC, kBC = calBC.rhoBC2mBC(rhoBC)
	if debug:
		print('mBC :', nBC, '+', kBC, 'i')
	
	# read shell/no BC complex refractive index data. If not exist, using default data
	
	nShell = 1.58 # from Beijing University observation
	if debug:
		print('nShell :', nShell)
	nH1 = 0.5
	
	# read kappa data. If not exist, using default data
	
	CTH1 = 0.5
	
	kappa = 0.21537370768611558 # from Beijing University observation
	# according to Petters and Kreidenweis 2007, kappa can change from 0.002 to 0.67
	# set kappaI change range from 0.5 * kappa to 1.5 * kappa, equal to 0.1
	if debug:
		print('kappa :', kappa)
	kappaH1 = 0.5
	
	# save origin factor value
	
	paras = dict(n=nShell, nH1=nH1, nBC=nBC, kBC=kBC, PNSD=0.5, MS=MS, MSH1=MSH1, MSH2=MSH2, CT=1, CTH1=CTH1, kappa=kappa, kappaH1=kappaH1, rhoBC=rhoBC, BCPNSD=0.5, BCAE=1)
	
	# processing data change
	###################################################################################
	
	# distribution change
	
	if PNSD_rate_on:
		ratio = 0.5 * PNSD_rate
		PNSD = clean_PNSD * ratio + dirty_PNSD * (1-ratio)
	elif PNSD_on:
		# mixing two types of the typical factor
		# using two PNSD modal: one for clean, another for polluted
		# when turbulence on, control the mixing ratio of two modal
		# turbulence off, using the clean modal
		PNSD_rate = 1 + random.normalvariate(mu=0,sigma=1/3) * PNSD_range
		# rate from 0 to 2 in 99% posibility
		ratio = 0.5 * PNSD_rate
		PNSD = dirty_PNSD * ratio + clean_PNSD * (1-ratio)
	else:
		ratio = 0.5
		PNSD = dirty_PNSD * ratio + clean_PNSD * (1-ratio)
		PNSD_rate = 0
	if debug:
		print('PNSD : ', PNSD_rate)
	
	######################################################################
	# adjust Dps resolution to simplify calculate pressure
	#print(PNSD, Dps)
	PNSD, Dps = calPNSD.PNSD_Dps_adjust(PNSD, Dps, default_mag)
	#print(PNSD, Dps)
	######################################################################
	
	if BCPNSD_rate_on:
		ratio = 0.5 * BCPNSD_rate
		BCPNSD = dirty_BCPNSD * ratio + clean_BCPNSD * (1-ratio)
	elif BCPNSD_on:
		# one for clean, another for polluted
		BCPNSD_rate = 1 + random.normalvariate(mu=0,sigma=1/3) * BCPNSD_range
		ratio = 0.5 * BCPNSD_rate
		# setting BCPNSD's mixiture
		BCPNSD = dirty_BCPNSD * ratio + clean_BCPNSD * (1-ratio)
	else:
		ratio = 0.5
		BCPNSD = dirty_BCPNSD * ratio + clean_BCPNSD * (1-ratio)
		BCPNSD_rate = 0
	if debug:
		print('BCPNSD : ', BCPNSD_rate)
	
	# MS, rhoBC and BCPMSD to make BCPNSD_ext
	
	if MS_rate_on:
		MS = MS * MS_rate
		if MS > 1:
			MS = 1
	elif MS_on:
		MS_rate = 1 + random.normalvariate(mu=0,sigma=1/3) * MS_range
		MS = MS * MS_rate
		if MS > 1:
			MS = 1
	else:
		MS_rate = 0
	
	if MSH1_rate_on:
		MSH1 = MSH1 * MSH1_rate
	elif MSH1_on:
		MSH1_rate = 1 + random.normalvariate(mu=0,sigma=1/3) * MSH1_range
		MSH1 = MSH1 * MSH1_rate
	else:
		MSH1_rate = 0
		MSH1 = 0
	if debug:
		print('MS heterogeneity 1 : ', MSH1)
	
	if MSH2_rate_on:
		MSH2y = MSH2 * MSH2_rate
	elif MSH2_on:
		MSH2_rate = 1 + random.normalvariate(mu=0,sigma=1/3) * MSH2_range
		MSH2y = MSH2 * MSH2_rate
	else:
		MSH2_rate = 0
		MSH2x = 1
		MSH2y = 0
	if debug:
		print('MS heterogeneity 2 : ', MSH2y)
	
	if rhoBC_rate_on:
		rhoBC = rhoBC * rhoBC_rate
	elif rhoBC_on:
		rhoBC_rate = 1 + random.normalvariate(mu=0,sigma=1/3) * rhoBC_range
		rhoBC = rhoBC * rhoBC_rate
	else:
		rhoBC_rate = 0
	
	# complex refractive index change
	
	if n_rate_on:
		nShell = nShell * n_rate
	elif n_on:
		n_rate = 1 + random.normalvariate(mu=0,sigma=1/3) * n_range
		nShell = nShell * n_rate
	else:
		n_rate = 0
	if debug:
		print('n:', nShell)
	
	if nH1_rate_on:
		nH1 = nH1 * (nH1_rate-1)
		# for input rate, the calculate formula is 1 + (-1~1) * range
		# to get -0.03 ~ 0.03, use the calculation above
	elif nH1_on:
		nH1_rate = 1 + random.normalvariate(mu=0,sigma=1/3) * nH1_range
		# according to Zhao 2020, the maximum n change range is 0.02
		# let 0.1 multiply by nI default range 0.3 get 0.03
		nH1 = nH1 * (nH1_rate-1)
	else:
		nH1_rate = 0
		nH1 = 0
	if debug:
		print('n heterogeneity 1 : ', nH1)
	
	# rhoBC to calculate new nBC and kBC
	
	nBC, kBC = calBC.rhoBC2mBC(rhoBC)
	
	if nBC_rate_on:
		nBC = nBC * nBC_rate
	elif nBC_on:
		nBC_rate = 1 + random.normalvariate(mu=0,sigma=1/3) * nBC_range
		nBC = nBC * nBC_rate
	else:
		nBC_rate = 0
	if debug:
		print('nBC:', nBC)
	
	if kBC_rate_on:
		kBC = kBC * kBC_rate
	elif kBC_on:
		kBC_rate = 1 + random.normalvariate(mu=0,sigma=1/3) * kBC_range
		kBC = kBC * kBC_rate
	else:
		kBC_rate = 0
	if debug:
		print('kBC:', kBC)
	
	# other microphysical properties change
	
	if CT_rate_on:
		CT = CT_rate
	elif CT_on:
		CT_rate = 1 + random.normalvariate(mu=0,sigma=1/3) * CT_range
		CT = CT_rate
	else:
		CT = 1
		CT_rate = 0
	if debug:
		print('Coating thickness:', CT)
	
	if CTH1_rate_on:
		CTH1 = CTH1 * CTH1_rate
	elif CTH1_on:
		CTH1_rate = 1 + random.normalvariate(mu=0,sigma=1/3) * CTH1_range
		CTH1 = CTH1 * CTH1_rate
	else:
		CTH1_rate = 0
		CTH1 = 0
	if debug:
		print('CT heterogeneity 1 : ', CTH1)
	
	if kappa_rate_on:
		kappa = kappa * kappa_rate
	elif kappa_on:
		kappa_rate = 1 + random.normalvariate(mu=0,sigma=1/3) * kappa_range
		kappa = kappa * kappa_rate
	else:
		kappa_rate = 0
	if debug:
		print('Kappa:', kappa)
	
	# read kappa mixing state data. If not exist, using default data
	
	if kappaH1_rate_on:
		kappaH1 = kappaH1 * (kappaH1_rate-1)
	elif kappaH1_on:
		kappaH1_rate = 1 + random.normalvariate(mu=0,sigma=1/3) * kappaH1_range
		# 0.33 multiple by 0.3, max range equal to 0.1
		kappaH1 = kappaH1 * (kappaI_rate-1)
	else:
		kappaH1_rate = 0
		kappaH1 = 0
	if debug:
		print('kappa heterogeneity 1 : ', kappaH1)
	
	if BCAE_rate_on:
		BCAE = BCAE_rate
	elif BCAE_on:
		BCAE_rate = 1 + random.normalvariate(mu=0,sigma=1/3) * BCAE_range
		BCAE = BCAE_rate
	else:
		BCAE = 1
		BCAE_rate = 0
	if debug:
		print('Black carbon absorbing enhancement:', BCAE)
	
	# start running
	##################################################################################
	
	AOD = np.zeros(len(wl))
	SSA = np.zeros(len(wl))
	g = np.zeros(len(wl))
	
	for i in range(len(wl)):
		AOD[i], SSA[i], g[i] = calMie.Mie2(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, nH1, kappa, kappaH1, MS, MSH1, MSH2x, MSH2y, dz, RH, wl, CT, CTH1, BCAE, debug=debug)
		if debug:
			print('AOD :', AOD[i], ',\tSSA :', SSA[i], ',\tg :', g[i])
	
	# write data
	##################################################################################
	'''
	n for real part of complex refractive index
	k for imagine part of complex refractive index
	nBC for BC real part of complex refractive index
	kBC for BC imagine part of complex refractive index
	PNSD for number distribution
	MS for mixing state
	VD for vertical distribution
	CT for coating thickness
	kappa for hydroscopicity parameter
	kappaI for kappa mixing states -> need to discuss
	rhoBC for BC density
	BCPNSD for BC number distribution
	BCPMSD for BC mass distribution -> need to discuss
	BCI for BC mixing inhomogeneity -> need to discuss
	BCAE for BC absorbing enhancement
	'''
	
	infos = dict(AOD=AOD, SSA=SSA, g=g, n_rate=n_rate, nH1_rate=nH1_rate, nBC_rate=nBC_rate, kBC_rate=kBC_rate, PNSD_rate=PNSD_rate, MS_rate=MS_rate, MSH1_rate=MSH1_rate, MSH2_rate=MSH2_rate, CT_rate=CT_rate, CTH1_rate=CTH1_rate, kappa_rate=kappa_rate, kappaH1_rate=kappaH1_rate, rhoBC_rate=rhoBC_rate, BCPNSD_rate=BCPNSD_rate, BCAE_rate=BCAE_rate)
	print(infos)
	return infos, paras

def make_dir(time, **args):
	'''
	This function is to make directory for runMie code output
	input:
		time         : year month and day, string, for example 220921
		**paras      : parameters list, array, in string, default all parameters
	output:
		parameter directory with time
	'''
	if 'paras' in args:
		paras = args['paras']
	else:
		paras = default_paras
	
	import os
	
	path = 'output/Mie/' + time
	os.system('mkdir '+path)
	os.system('mkdir '+path+'/all')
	for i in range(len(paras)):
		fn = path + '/' + paras[i]
		os.system('mkdir '+fn)

def LHS_run(time, run_num, par_range, **args):
	'''
	This function is to circulation run Mie in Latin hypercubic sampling method
	input:
		time      : running time to record
		run_num   : running times number
		par_range : parameters change range
		**debug   : debug flag, bool, default False
	output:
		aerosol optical parameter information, including:
			AOD   : aerosol optical depth
			SSA   : single scattering albedo
			g     : asymmetry factor
	'''
	if 'debug' in args:
		debug = args['debug']
	else:
		debug = False
	
	#read data
	
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
	
	# use LHS produce all parameters for whole run
	
	parameters = default_paras
	
	par_num = len(parameters)
	all_rate = np.zeros((par_num, run_num))
	paras_rate = np.zeros((par_num, par_num, run_num))
	
	# for all run and parameter run
	
	for i in range(par_num):
		all_rate[i] = 1 + LHS_norm(run_num, 0, 1/3).flatten() * default_range
	
	for i in range(par_num):
		for j in range(par_num):
			if i == j:
				paras_rate[i,j] = 1 + LHS_norm(run_num, 0, 1/3).flatten() * par_range
			else:
				paras_rate[i,j] = all_rate[j]
	
	make_dir(time)
	
	#start running
	
	path = 'output/Mie/' + time + '/'
	
	infos = []
	for i in range(run_num):
		info, paras = run(1000, [525], 70, MSH2x=5, 
		n_rate			= all_rate[0,i], 
		nH1_rate		= all_rate[1,i], 
		nBC_rate		= all_rate[2,i], 
		kBC_rate		= all_rate[3,i], 
		PNSD_rate		= all_rate[4,i], 
		MS_rate			= all_rate[5,i], 
		MSH1_rate		= all_rate[6,i], 
		MSH2_rate		= all_rate[7,i], 
		CT_rate			= all_rate[8,i], 
		CTH1_rate		= all_rate[9,i], 
		kappa_rate		= all_rate[10,i], 
		kappaH1_rate	= all_rate[11,i], 
		rhoBC_rate		= all_rate[12,i], 
		BCPNSD_rate		= all_rate[13,i], 
		BCAE_rate		= all_rate[14,i], 
		Dps=Dps, DBC=DBC, DBCps=DBCps, clean_PNSD=clean_PNSD, dirty_PNSD=dirty_PNSD, clean_BCPNSD=clean_BCPNSD, dirty_BCPNSD=dirty_BCPNSD, debug=debug)
		infos.append(info)
		np.save(path+'all/infos.npy', infos)
		if i==0:
			np.save(path+'all/paras.npy', paras)
	
	for i in range(par_num):
		infos = []
		for j in range(run_num):
			info, paras = run(1000, [525], 70, MSH2x=5, 
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

"""
def co_run(time, run_num, par_range, paras, co_paras, **args):
####################################################################################
# warning:
# this function may have problems in math
####################################################################################
	'''
	This function is to use LHS_run to test synergy of several parameters
	input:
		time      : running time to record
		run_num   : running times number
		par_range : parameters change range
		paras     : parameters list, array, string
		co_paras  : parameters list to calculate synergy, array, string, 
		            in shape (co_num, co_paras_num)
		**debug   : debug flag, bool, default False
	output:
		aerosol optical parameter information, including:
			AOD   : aerosol optical depth
			SSA   : single scattering albedo
			g     : asymmetry factor
	'''
	if 'debug' in args:
		debug = args['debug']
	else:
		debug = False
	
	#co_paras must in paras
	for i in range(len(co_paras)):
		for j in range(len(co_paras[i])):
			if co_paras[i][j] not in paras:
				print('Co-paras not in paras. Please check.')
				sys.exit()
	
	#read data
	
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
	
	# use LHS produce all parameters for whole run
	
	par_num = len(paras)
	all_num = len(default_paras)
	co_num = len(co_paras)
	
	all_rate = np.zeros((all_num, run_num))
	paras_rate = np.zeros((par_num, all_num, run_num))
	co_rate = np.zeros((co_num, all_num, run_num))
	
	# for all run, parameter run and co run
	
	for i in range(all_num):
		all_rate[i] = 1 + LHS_norm(run_num, 0, 1/3).flatten() * default_range
	
	for i in range(par_num):
		for j in range(all_num):
			if paras[i] == default_paras[j]:
				paras_rate[i,j] = 1 + LHS_norm(run_num, 0, 1/3).flatten() * par_range
			else:
				paras_rate[i,j] = all_rate[j]
	
	for i in range(co_num):
		for j in range(all_num):
			if default_paras[j] in co_paras[i]: # paras rate would be same as single para rate
				co_rate[i,j] = paras_rate[paras.index(default_paras[j]),j]
			else:
				co_rate[i,j] = all_rate[j]
	
	# sum co_paras' name
	for i in range(co_num):
		sep = '+'
		paras.append(sep.join(co_paras[i]))
	
	make_dir(time, paras=paras)
	
	#start running
	
	path = 'output/Mie/' + time + '/'
	
	infos = []
	for i in range(run_num):
		info = run(2000, 100, [525], nI2x=2, kappaI2x=2,
		n_rate=		all_rate[0,i], 
		nI_rate=	all_rate[1,i], 
		nI2_rate=	all_rate[2,i],
		nBC_rate=	all_rate[3,i], 
		kBC_rate=	all_rate[4,i], 
		PNSD_rate=	all_rate[5,i], 
		MS_rate=	all_rate[6,i], 
		#VD_rate=	all_rate[7,i], 
		CT_rate=	all_rate[8,i], 
		kappa_rate=	all_rate[9,i], 
		kappaI_rate=	all_rate[10,i], 
		kappaI2_rate=	all_rate[11,i], 
		rhoBC_rate=	all_rate[12,i], 
		#BCPNSD_rate=	all_rate[13,i], 
		#BCPMSD_rate=	all_rate[14,i], 
		BCI_rate=	all_rate[15,i], 
		BCAE_rate=	all_rate[16,i], 
		Dps=Dps, DBC=DBC, DBCps=DBCps, clean_PNSD=clean_PNSD, dirty_PNSD=dirty_PNSD, clean_BCPNSD=clean_BCPNSD, dirty_BCPNSD=dirty_BCPNSD, DBC2=DBC2, clean_BCPMSD=clean_BCPMSD, dirty_BCPMSD=dirty_BCPMSD, debug=debug)
		infos.append(info)
		np.save(path+'all/infos.npy', infos) # save info every single calculation
	
	for i in range(par_num):
		infos = []
		for j in range(run_num):
			info = run(2000, 100, [525], nI2x=2, kappaI2x=2,
			n_rate=		paras_rate[i,0,j], 
			nI_rate=	paras_rate[i,1,j], 
			nI2_rate=	paras_rate[i,2,j],
			nBC_rate=	paras_rate[i,3,j], 
			kBC_rate=	paras_rate[i,4,j], 
			PNSD_rate=	paras_rate[i,5,j], 
			MS_rate=	paras_rate[i,6,j], 
			#VD_rate=	paras_rate[i,7,j], 
			CT_rate=	paras_rate[i,8,j], 
			kappa_rate=	paras_rate[i,9,j], 
			kappaI_rate=	paras_rate[i,10,j], 
			kappaI2_rate=	paras_rate[i,11,j], 
			rhoBC_rate=	paras_rate[i,12,j], 
			#BCPNSD_rate=	paras_rate[i,13,j], 
			#BCPMSD_rate=	paras_rate[i,14,j], 
			BCI_rate=	paras_rate[i,15,j], 
			BCAE_rate=	paras_rate[i,16,j], 
			Dps=Dps, DBC=DBC, DBCps=DBCps, clean_PNSD=clean_PNSD, dirty_PNSD=dirty_PNSD, clean_BCPNSD=clean_BCPNSD, dirty_BCPNSD=dirty_BCPNSD, DBC2=DBC2, clean_BCPMSD=clean_BCPMSD, dirty_BCPMSD=dirty_BCPMSD, debug=debug)
			infos.append(info)
			np.save(path+paras[i]+'/infos.npy', infos)
	
	for i in range(co_num):
		infos = []
		for j in range(run_num):
			info = run(2000, 100, [525], nI2x=2, kappaI2x=2,
			n_rate=		co_rate[i,0,j], 
			nI_rate=	co_rate[i,1,j], 
			nI2_rate=	co_rate[i,2,j],
			nBC_rate=	co_rate[i,3,j], 
			kBC_rate=	co_rate[i,4,j], 
			PNSD_rate=	co_rate[i,5,j], 
			MS_rate=	co_rate[i,6,j], 
			#VD_rate=	co_rate[i,7,j], 
			CT_rate=	co_rate[i,8,j], 
			kappa_rate=	co_rate[i,9,j], 
			kappaI_rate=	co_rate[i,10,j], 
			kappaI2_rate=	co_rate[i,11,j], 
			rhoBC_rate=	co_rate[i,12,j], 
			#BCPNSD_rate=	co_rate[i,13,j], 
			#BCPMSD_rate=	co_rate[i,14,j], 
			BCI_rate=	co_rate[i,15,j], 
			BCAE_rate=	co_rate[i,16,j], 
			Dps=Dps, DBC=DBC, DBCps=DBCps, clean_PNSD=clean_PNSD, dirty_PNSD=dirty_PNSD, clean_BCPNSD=clean_BCPNSD, dirty_BCPNSD=dirty_BCPNSD, DBC2=DBC2, clean_BCPMSD=clean_BCPMSD, dirty_BCPMSD=dirty_BCPMSD, debug=debug)
			infos.append(info)
			np.save(path+paras[i-co_num]+'/infos.npy', infos)
"""

if __name__ == '__main__':
	LHS_run('240604', 2000, 0.05, debug=True)
	#co_run('230916_n_kappa', 10000, 0.05, ['n','kappa'], debug=True)
	#co_run('230916_MS_CT', 10000, 0.05, ['MS','CT'], debug=True)
	#co_paras = [['n','kappa'], ['MS','CT']]
	#co_run('230917', 5000, 0.05, main_paras, co_paras, debug=True)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
