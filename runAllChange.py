####################################################################################
# INTRODUCTION:
# This code is to test all factor turburlence sensitivity of DARF in SBDART
# Created by Hebs at 22/5/10/12:15
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import random
import os
import sys

import calRH
import calPNSD
import kappaKohler
import phaseFunc
import readTaizhou
import readAtms
import write_aerosol
import write_albedo
import runGChange

###############################################
# some constants
default_range = 0.5 # 50% turbulence
###############################################

def run(nn, moma, wl, AR, **args):
	'''
	This function is to use Monte Carlo method to test all DARF factors' sensitivity
	Factors are including: 
		aerosol factors:
			AOD
			SSA
			g
			mixing state
			complex refractive index
			number distribution
			vertical distribution
			phase function
			BC complex refractive index
			BC mass distribution
			BC mixing state
			BC absorbing enhancement
			chemical composition
			BC coating thickness
			BC density
			BC form
			BC distribution
			hygroscopicitic factor
		gass factors:
			Rayleigh scattering
			ozone abosorbing
			water abosorbing
		cloud factors:
			AOD
			SSA
		environment factors:
			pressure
			temperature
			raletive humidity
			solar radiative
			albedo
	Input should including:
		factor turburlence range, factors including:
			albedo
			AOD
			SSA
			g
			MS for mixing state
			k for imagine part of complex refractive index
			n for real part of complex refractive index
			PNSD for number distribution
			PF for phase function
			VD for vertical distribution
			rhoBC for BC density
			kBC for BC imagine part of complex refractive index
			nBC for BC real part of complex refractive index
			BCMD for BC mass distribution
			BCMS for BC mixing state
			BCAE for BC absorbing enhancement
			CT for coating thickness
			kappa
		factor turburlence switch
		angular resolution
		debug flag
		data file path
		intput file path
		output file path
	input:
		nn    : number of atmospheric levels for which aerosol information is specified
		moma  : number of phase function moments
		wl    : the wavelength [ wl(k) < wl(k+1) ], array, nm
		AR    : angular resolution for phase function, degree
	'''
	###############################################
	# args reading
	if 'albedo_on' in args:
		albedo_on = args['albedo_on']
		if 'albedo_range' in args:
			albedo_range = args['albedo_range']
		else:
			albedo_range = default_range
	else:
		albedo_on = False
	if 'AOD_on' in args:
		AOD_on = args['AOD_on']
		if 'AOD_range' in args:
			AOD_range = args['AOD_range']
		else:
			AOD_range = default_range
	else:
		AOD_on = False
	if 'SSA_on' in args:
		SSA_on = args['SSA_on']
		if 'SSA_range' in args:
			SSA_range = args['SSA_range']
		else:
			SSA_range = default_range
	else:
		SSA_on = False
	if 'g_on' in args:
		g_on = args['g_on']
		if 'g_range' in args:
			g_range = args['g_range']
		else:
			g_range = default_range
	else:
		g_on = False
	if 'MS_on' in args:
		MS_on = args['MS_on']
		if 'MS_range' in args:
			MS_range = args['MS_range']
		else:
			MS_range = default_range
	else:
		MS_on = False
	if 'k_on' in args:
		k_on = args['k_on']
		if 'k_range' in args:
			k_range = args['k_range']
		else:
			k_range = default_range
	else:
		k_on = False
	if 'n_on' in args:
		n_on = args['n_on']
		if 'n_range' in args:
			n_range = args['n_range']
		else:
			n_range = default_range
	else:
		n_on = False
	if 'PNSD_on' in args:
		PNSD_on = args['PNSD_on']
		if 'PNSD_range' in args:
			PNSD_range = args['PNSD_range']
		else:
			PNSD_range = default_range
	else:
		PNSD_on = False
	if 'VD_on' in args:
		VD_on = args['VD_on']
		if 'VD_range' in args:
			VD_range = args['VD_range']
		else:
			VD_range = default_range
	else:
		VD_on = False
	if 'rhoBC_on' in args:
		rhoBC_on = args['rhoBC_on']
		if 'rhoBC_range' in args:
			rhoBC_range = args['rhoBC_range']
		else:
			rhoBC_range = default_range
	else:
		rhoBC_on = False
	if 'PF_on' in args:
		PF_on = args['PF_on']
		if 'PF_range' in args:
			PF_range = args['PF_range']
		else:
			PF_range = default_range
	else:
		PF_on = False
	if 'kBC_on' in args:
		kBC_on = args['kBC_on']
		if 'kBC_range' in args:
			kBC_range = args['kBC_range']
		else:
			kBC_range = default_range
	else:
		kBC_on = False
	if 'nBC_on' in args:
		nBC_on = args['nBC_on']
		if 'nBC_range' in args:
			nBC_range = args['nBC_range']
		else:
			nBC_range = default_range
	else:
		nBC_on = False
	if 'BCPNSD_on' in args:
		BCPNSD_on = args['BCPNSD_on']
		if 'BCPNSD_range' in args:
			BCPNSD_range = args['BCPNSD_range']
		else:
			BCPNSD_range = default_range
	else:
		BCPNSD_on = False
	if 'BCMD_on' in args:
		BCMD_on = args['BCMD_on']
		if 'BCMD_range' in args:
			BCMD_range = args['BCMD_range']
		else:
			BCMD_range = default_range
	else:
		BCMD_on = False
	if 'BCMS_on' in args:
		BCMS_on = args['BCMS_on']
		if 'BCMS_range' in args:
			BCMS_range = args['BCMS_range']
		else:
			BCMS_range = default_range
	else:
		BCMS_on = False
	if 'BCAE_on' in args:
		BCAE_on = args['BCAE_on']
		if 'BCAE_range' in args:
			BCAE_range = args['BCAE_range']
		else:
			BCAE_range = default_range
	else:
		BCAE_on = False
	if 'CT_on' in args:
		CT_on = args['CT_on']
		if 'CT_range' in args:
			CT_range = args['CT_range']
		else:
			CT_range = default_range
	else:
		CT_on = False
	if 'kappa_on' in args:
		kappa_on = args['kappa_on']
		if 'kappa_range' in args:
			kappa_range = args['kappa_range']
		else:
			kappa_range = default_range
	else:
		kappa_on = False
	if 'debug' in args:
		debug = args['debug']
	else:
		debug = False
	if 'data_path' in args:
		data_path = args['data_path']
	else:
		data_path = 'data/'
	if 'input_path' in args:
		input_path = args['input_path']
	else:
		input_path = 'input/all'
	if 'output_path' in args:
		output_path = args['output_path']
	else:
		output_path = 'output/all/'
	if 'output_name' in args:
		output_name = args['output_name']
	else:
		output_name = 'output'
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
	#################################################
	# initiating
	#################################################
	# calculating
	'''
	for every subroutine, do as below:
		1.judge the factor whether to turbulence
		2.verify the turbulence range
		3.using a time seed to produce a random value among the range
		4.using the random value to revise input file
		5.back up the origin input file
		6.replace the input file with the revised one
		7.do next factor, until all factors are done
		8.run SBDART, output outcomes to output path
		9.replace the origin input file back
	for aerosol.dat, to create a new file, need these Mie parameters:
		PNSD
			mixing state
			BCPNSD or BC mass distribution + BC density
			BC mixing state
		complex radiative index
			kappa
			RH
			BC m
			BC n
			shell/no BC n
			BCPNSD or BC mass distribution + BC density
			BC mixing state
			BC absorption enhancement
			coating thickness
			BC form?
			BC distribution?
		phase function
	'''
	# revise albedo.dat
	if albedo_on:
		if debug:
			print('albedo on')
		rate = 1 + random.normalvariate(mu=0,sigma=1/3) * albedo_range
		if debug:
			print('albedo rate:', rate)
		write_albedo.change(rate)
		albedo_rate = rate
	else:
		albedo_rate = 0
	
	# revise aerosol.dat
	'''
	The basic data input need to be:
		PNSD
		BCPNSD or BC mass distribution + BC density
		BC m
		BC n
		shell/no BC n
		kappa
		RH
	then using Mie model to calculate Mie parameters:
		complex radiative index
		phase function
	then calculating SBDART parameters:
		AOD
		SSA
		g
	finally write a new aerosol.data file
	'''
	# read data
	# read PNSD and BCPNSD data
	if not Dps_exist or not clean_PNSD_exist or not dirty_PNSD_exist or not DBCps_exist or not clean_BCPNSD_exist or not dirty_BCPNSD_exist:
		BCPNSD_exist = False
		fn = data_path + 'sp2/Taizhou.npy'
		if os.path.exists(fn):
			if debug:
				print('reading BCPNSD data...')
			BCPNSD_exist = True
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
			'''
			clean_PNSD = np.nanmean(PNSD, axis=0)
			dirty_PNSD = PNSD[i_max]
			clean_BCPNSD = np.nanmean(BCPNSD, axis=0)
			dirty_BCPNSD = BCPNSD[iBC_max]
			'''
			clean_PNSD = np.nanmean(PNSD, axis=0) / 10
			dirty_PNSD = np.nanmean(PNSD, axis=0) / 2
			clean_BCPNSD = np.nanmean(BCPNSD, axis=0) / 10
			dirty_BCPNSD = np.nanmean(BCPNSD, axis=0) / 2
	
	else:
		BCPNSD_exist = True
			
	if debug:
		print('Dps:', Dps)
		print('clean PNSD:', clean_PNSD)
		print('dirty PNSD:', dirty_PNSD)
		print('DBCps:', DBCps)
		print('DBC:', DBC)
		print('clean BCPNSD:', clean_BCPNSD)
		print('dirty BCPNSD:', dirty_BCPNSD)
		print('done')
	
	# read BCMD data
	fn = data_path + 'BCMD.npy'
	BCMD_exist = False
	if os.path.exists(fn):
		if debug:
			print('reading BCMD data...')
		# read in BCMD data here, clean and dirty saperately
		
		
		
		
		# working zone
		
		
		
		
		BCMD_exist = True
		if debug:
			print('done')
	# if neither BCPNSD nor BCMD data exist, using default data
	if not BCPNSD_exist and not BCMD_exist:
		if debug:
			print('No BCPNSD or BCMD data. Using default data.')
		print('reading default BCPNSD data...')
		if debug:
			print('done')
	# read BCMS data. If not exist, using default data
	fn = data_path + 'BCMS.npy'
	if os.path.exists(fn):
		if debug:
			print('reading BCMS data...')
		if debug:
			print('done')
	else:
		if debug:
			print('No BCMS data. Using default data.')
		print('reading default BCMS data...')
		if debug:
			print('done')
	# read mixing state. If not exist, using default data
	MS = 0.5
	# read vertical distribution type. If not exist, using default data
	VD = 0.5
	# read BC density. If not exist, using default data
	rhoBC = 1.6 # g/m3
	x = 1.25 * (rhoBC -1)
	# turn BCMD and BC density to BCPNSD
	# read BC complex refractive index data. If not exist, using default data
	kBC = 1 + 1.26 * x
	nBC = 1.26 * x + 1
	# raed BC mixing state. If not exist, using default data
	BCMS = 0.5
	# read shell/no BC complex refractive index data. If not exist, using default data
	nShell = 1.58
	# read kappa data. If not exist, using default data
	kappa = 0.21537370768611558
	# read RH data, from top to bottom
	atms = readAtms.read()
	z = atms['z'] * 1e3 # in m
	dz = z[:-1] - z[1:] # dz in each two layers, in m
	wh = atms['wh']
	t = atms['t']
	nz = len(z) # number of atmosphere layers
	# turn water vapor density in g/m3 to RH, in %
	RH = np.zeros(len(wh))
	for i in range(len(RH)):
		RH[i] = calRH.wh2RH(t[i], wh[i])
		if RH[i]<1:
			RH[i] = 1 # minimum RH set to 1
	
	# check the input data validity
	# if nn bigger than atmosphere layer number, exit
	if nn > nz - 1:
		print('Too many aerosol layers. Please check.')
		sys.exit()
	# if wl(k) < wl(k+1) doesn\'t meet, exit
	for i in range(len(wl)-1):
		if (wl[i+1]-wl[i])<0:
			print('wl(k) < wl(k+1) doesn\'t meet. Please check.')
			sys.exit()
	
	# produce random value for each turbulence factor
	'''
	AOD
	SSA
	g
	MS for mixing state
	k for imagine part of complex refractive index
	n for real part of complex refractive index
	PNSD for number distribution
	VD for vertical distribution
	PF for phase function
	kBC for BC imagine part of complex refractive index
	nBC for BC real part of complex refractive index
	BCPNSD
	BCMD for BC mass distribution?
	BCMS for BC mixing state
	BCAE for BC absorbing enhancement
	CT for coating thickness
	kappa
	'''
	# directly change the value of the factor
	if k_on:
		rate = 1 + random.normalvariate(mu=0,sigma=1/3) * k_range
		kBC = kBC * rate
		k_rate = rate
	else:
		k_rate = 0
	if n_on:
		rate = 1 + random.normalvariate(mu=0,sigma=1/3) * n_range
		nBC = nBC * rate
		nShell = nShell * rate
		n_rate = rate
	else:
		n_rate = 0
	if rhoBC_on:
		rate = 1 + random.normalvariate(mu=0,sigma=1/3) * rhoBC_range
		rhoBC = rhoBC * rate
		rhoBC_rate = rate
	else:
		rhoBC_rate = 0
	x = 1.25 * (rhoBC -1)
	# turn BCMD and BC density to BCPNSD
	# read BC complex refractive index data. If not exist, using default data
	kBC = 1 + 1.26 * x
	nBC = 1.26 * x + 1
	if kBC_on:
		rate = 1 + random.normalvariate(mu=0,sigma=1/3) * kBC_range
		kBC = kBC * rate
		kBC_rate = rate
	else:
		kBC_rate = 0
	if nBC_on:
		rate = 1 + random.normalvariate(mu=0,sigma=1/3) * nBC_range
		nBC = nBC * rate
		nBC_rate = rate
	else:
		nBC_rate = 0
	if MS_on:
		rate = 1 + random.normalvariate(mu=0,sigma=1/3) * MS_range
		MS = MS * rate
		if MS > 1:
			MS = 1
		MS_rate = rate
	else:
		MS_rate = 0
	if VD_on:
		rate = 1 + random.normalvariate(mu=0,sigma=1/3) * VD_range
		VD = VD * rate
		if VD > 1:
			VD = 1
		VD_rate = rate
	else:
		VD_rate = 0
	if kappa_on:
		rate = 1 + random.normalvariate(mu=0,sigma=1/3) * kappa_range
		kappa = kappa * rate
		kappa_rate = rate
	else:
		kappa_rate = 0
	
	# mixing two types of the typical factor
	if PNSD_on:
		# using two PNSD modal: one for clean, another for polluted
		# when turbulence on, control the mixing ratio of two modal
		# turbulence off, using the clean modal
		rate = 1 + random.normalvariate(mu=0,sigma=1/3) * PNSD_range
		# rate from 0 to 2 in 99% posibility
		ratio = 0.5 * rate
		PNSD = clean_PNSD * ratio + dirty_PNSD * (1-ratio)
		PNSD_rate = rate
	else:
		PNSD_rate = 0
	if BCPNSD_on:
		# one for clean, another for polluted
		# BCMD and BC density has turned to BCPNSD
		rate = 1 + random.normalvariate(mu=0,sigma=1/3) * BCPNSD_range
		ratio = 0.5 * rate
		# setting BCPNSD's mixiture
		BCPNSD = clean_BCPNSD * ratio + dirty_BCPNSD * (1-ratio)
		BCPNSD_rate = rate
	else:
		BCPNSD_rate = 0
	
	# after PNSD adjust, then do BC micro physical properties adjust
	if CT_on:
		rate = 1 + random.normalvariate(mu=0,sigma=1/3) * CT_range
		CT = rate
		CT_rate = rate
	else:
		CT = 1
		CT_rate = 0
	
	# calculating the parameter that Mie theory needed
	# while calculating, do factors' adjust
	if BCAE_on:
		rate = 1 + random.normalvariate(mu=0,sigma=1/3) * BCAE_range
		BCAE = rate
		BCAE_rate = rate
	else:
		BCAE = 1
		BCAE_rate = 0
	if BCMS_on:
		rate = 1 + random.normalvariate(mu=0,sigma=1/3) * BCMS_range
		BCMS = BCMS * rate
		if BCMS > 1:
			BCMS = 1
		BCMS_rate = rate
	else:
		BCMS_rate = 0
	
	# after Mie calculating, adjust left factors
	if PF_on:
		# according to g, change the phase function's trend: forward or backward
		rate = 1 + random.normalvariate(mu=0,sigma=1/3) * PF_range
		PF = rate
		PF_rate = rate
	else:
		PF = 1
		PF_rate = 0
	
	# do calculating for each wavelength and each layer
	dtau = np.zeros((len(wl), nn)) # for AOD
	waer = np.zeros((len(wl), nn)) # for SSA
	pmom = np.zeros((len(wl), nn, moma)) # for legendre moments of phase function
	for i in range(len(wl)):
		for j in range(nn):
			dtau[i,j], waer[i,j], pmom[i,j] = Mie(Dps, PNSD, DBCps, DBC, BCPNSD, nBC, kBC, nShell, kappa, RH[j], wl[i], z[nz-nn+j], dz[nz-nn+j-1], moma, MS, BCAE, BCMS, CT, VD, AR, PF, debug=debug)
			# need more parameters in function Mie(), including: mixing state, BC absorbing enhancement, coating thickness, etc.
			if debug:
				print(round((i*nn+j+1)/len(wl)/nn*10000)/100, '% done...')
	
	# write a new aerosol.dat file
	write(wl, dtau, waer, pmom)
	
	# after aerosol.dat, do AOD, SSA, g turbulence
	if AOD_on:
		rate = 1 + random.normalvariate(mu=0,sigma=1/3) * AOD_range
		write_aerosol.change_AOD(rate)
		if debug:
			print('AOD rate:', rate)
		AOD_rate = rate
	else:
		AOD_rate = 0
	if SSA_on:
		rate = 1 + random.normalvariate(mu=0,sigma=1/3) * SSA_range
		write_aerosol.change_SSA(rate)
		if debug:
			print('SSA rate:', rate)
		SSA_rate = rate
	else:
		SSA_rate = 0
	if g_on:
		rate = 1 + random.normalvariate(mu=0,sigma=1/3) * g_range
		# do g change, adjust legendre moments of phase function
		# mabye it's same as phase function turbulence?
		g_rate = rate
	else:
		g_rate = 0
	
	#########################################################
	# run SBDART
	os.system('sbdart >'+output_path+output_name+'.txt')
	
	#########################################################
	# change albedo.dat and aerosol.dat back
	write_albedo.change_back()
	write_aerosol.change_back()
	
	#########################################################
	# write factors information into a dict
	factors = dict(albedo_rate=albedo_rate, AOD_rate=AOD_rate, SSA_rate=SSA_rate, g_rate=g_rate, k_rate=k_rate, n_rate=n_rate, kBC_rate=kBC_rate, nBC_rate=nBC_rate, MS_rate=MS_rate, VD_rate=VD_rate, rhoBC_rate=rhoBC_rate, kappa_rate=kappa_rate, PNSD_rate=PNSD_rate, BCPNSD_rate=BCPNSD_rate, CT_rate=CT_rate, BCAE_rate=BCAE_rate, BCMS_rate=BCMS_rate, PF_rate=PF_rate, Mie_infos=[dtau[:,0],waer[:,0],pmom[:,0,:]])
	return factors

def Mie(Dps, PNSD, DBCps, DBC, BCPNSD, nBC, kBC, nShell, kappa, RH, wl, z, dz, moma, MS, BCAE, BCMS, CT, VD, AR, PF, **args):
	'''
	This function is to use Mie modal to calculate Mie parameters,
	parameters including: 
		dtau for AOD
		waer for SSA
		pmom for legendre moments of phase function
	input:
		Dps                 : diameter distribution, array, nm
		PNSD                : number distribution, array, dn/dlogDp
		DBCps               : BC particle diameter size, array, nm
		DBC                 : BC core diameter size, array, nm
		BCPNSD              : BC particle number size distribution, array, dn/dlogDBC
		nBC                 : BC core real part of complex refractive index
		kBC                 : BC core imagine part of complex refractive index
		nShell              : shell complex refractive index
		kappa               : hygroscopicity parameter, float
		RH                  : raletive humidity, float, percent
		wl                  : wave length, float, nm
		z                   : height, float, m
		dz                  : thickness, float, m
		moma                : legendre moments' number of phase function, int
		MS                  : mixing state, 0 for pure internal mix, 
		                      1 for pure external mix, float
		BCAE                : BC absorbing enhancement adjustment, float
		BCMS                : BC mixing state adjustment, float
		CT                  : BC coating thickness adjustment, float
		VD                  : vertical distribution mixing, 0 for type A and 1 for type B
		AR                  : angular resolution, degree
		**debug             : debug flag, default False, bool
	output:
		dtau(i,j)	is the optical depth increment within level i at wavelength k,
				information is specified in top-down order. 

		waer(i,j)	is the single scattering albedo

		pmom(i,j,k)	are legendre moments of the phase function.
				Note that zeroeth moment is not read, it is assumed to be 1.
	'''
	if 'debug' in args:
		debug = args['debug']
	else:
		debug = False
	
	if debug:
		print('calculating...')
	
	# Before calculating, some factors change the input parameters, 
	# including: coating thickness, mixing state, BC mixing state, BC absorbing enhancement
	
	dtau = cal_AOD(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, kappa, RH, wl, z, dz, MS, CT, BCAE, VD)
	if debug:
		print('dtau:', dtau)
	waer = cal_SSA(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, kappa, RH, wl, MS, CT, BCAE)
	if debug:
		print('waer:', waer)
	pmom = np.zeros(moma)
	for i in range(moma):
		pmom[i] = cal_beta(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, kappa, RH, wl, z, i+1, VD, AR, PF)
	if debug:
		print('pmom:', pmom)
	
	# after calculating, some factors would change the value of the outcomes, 
	# including: BC absorbing enhancement
	
	if debug:
		print('done')
	
	return dtau, waer, pmom

def cal_AOD(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, kappa, RH, wl, z, dz, MS, CT, BCAE, VD, **args):
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
		MS      : mixing state, 0 for inertial mix and 1 for external mix
		CT      : BC core coating thickness adjust parameter
		BCAE    : BC core absorbing enhancement adjust parameter
		VD      : vertical distribution mixing, 0 for type A and 1 for type B
		**ddz   : step, default dz/100
	output:
		dtau    : AOD, aerosol optical depth
	'''
	if 'ddz' in args:
		ddz = args['ddz']
	else:
		ddz = dz / 10
	
	n = round(dz/ddz)
	
	PNSD_ext = PNSD * MS # this part have no BC core
	BCPNSD_ext = BCPNSD * MS # this part have no coating
	PNSD_int = PNSD * (1-MS)
	BCPNSD_int = BCPNSD * (1-MS)
	
	kext_int = calPNSD.cal_kext(Dps, PNSD_int, DBCps, DBC, BCPNSD_int, kBC, nBC, nShell, kappa, RH, wl, CT, BCAE)
	kext_ext = calPNSD.Dps2kext(Dps, PNSD_ext, wl, nShell) + calPNSD.Dps2kext(DBC, np.sum(BCPNSD_ext,axis=0), wl, nBC+kBC*1j)
	kext = kext_int + kext_ext
	H_integral = 0
	for i in range(n):
		H_integral += calPNSD.cal_VD(Dps, PNSD, z+n*ddz, VD) * ddz
	
	dtau = kext * 1e-6 * H_integral
	
	return dtau

def cal_SSA(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, kappa, RH, wl, MS, CT, BCAE):
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
		MS      : mixing state, 0 for inertial mix and 1 for external mix
		CT      : BC core coating thickness adjust parameter
		BCAE    : BC core absorbing enhancement adjust parameter
	output:
		waer    : SSA, single scattering albedo
	'''
	PNSD_ext = PNSD * MS # this part have no BC core
	BCPNSD_ext = BCPNSD * MS # this part have no coating
	PNSD_int = PNSD * (1-MS)
	BCPNSD_int = BCPNSD * (1-MS)
	
	kext_int = calPNSD.cal_kext(Dps, PNSD_int, DBCps, DBC, BCPNSD_int, kBC, nBC, nShell, kappa, RH, wl, CT, BCAE)
	kext_ext = calPNSD.Dps2kext(Dps, PNSD_ext, wl, nShell) + calPNSD.Dps2kext(DBC, np.sum(BCPNSD_ext,axis=0), wl, nBC+kBC*1j)
	kext = kext_int + kext_ext
	
	ksca_int = calPNSD.cal_ksca(Dps, PNSD_int, DBCps, DBC, BCPNSD_int, kBC, nBC, nShell, kappa, RH, wl, CT)
	ksca_ext = calPNSD.Dps2k(Dps, PNSD_ext, wl, nShell) + calPNSD.Dps2k(DBC, np.sum(BCPNSD_ext,axis=0), wl, nBC+kBC*1j)
	ksca = ksca_int + ksca_ext
	
	waer = ksca / kext
	
	return waer

def cal_beta(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, kappa, RH, wl, z, n, VD, AR, PF):
	'''
	This function is to calculate Legendre moments in certain height and raletive humidity
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
		n       : moment
		VD      : vertical distribution mixing, 0 for type A and 1 for type B
		AR      : angular resolution, degree
		PF      : phase function change rate
	output:
		beta    : Legendre moments
	'''
	m_BC = nBC + kBC * 1j
	m_shell = nShell
	n_BC = BCPNSD
	
	ratio = calPNSD.cal_VD(Dps, PNSD, z, VD)
	P_bulk, theta = phaseFunc.cal_bulk_phase_function(Dps, PNSD*ratio, DBCps, DBC, n_BC*ratio, m_BC, m_shell, kappa, RH, wl, angularResolution=AR, RHs=np.arange(1/10,21)*5)
	
	af = runGChange.cal_g_change_af(P_bulk, theta, PF)
	P_bulk_new = runGChange.g_change(P_bulk, theta, PF, af=af)
	b = phaseFunc.beta(P_bulk_new, theta, n)
	
	return b

def write(wl, dtau, waer, pmom, **args):
	'''
	This function is to write aerosol.dat
	input:
		wl            : the wavelength [ wl(k) < wl(k+1) ], array, nm
		dtau(i,j)	is the optical depth increment within level i at wavelength k,
				information is specified in top-down order. 

		waer(i,j)	is the single scattering albedo

		pmom(i,j,k)	are legendre moments of the phase function.
				Note that zeroeth moment is not read, it is assumed to be 1.
		**output_path : output path, string, default './'
	output:
		aerosol.dat
	'''
	if 'output_path' in args:
		output_path = args['output_path']
	else:
		output_path = './'
	
	fn = output_path + 'aerosol.dat'
	
	if os.path.exists(fn) and not os.path.exists(fn+'_origin'):
		os.system('mv '+fn+' '+fn+'_origin') # back up
	nn = pmom.shape[1]
	moma = pmom.shape[2]
	
	print('writing...')
	
	with open(fn, 'w') as f:
		f.write(str(nn)+'\t'+str(moma)+'\n')
		for i in range(len(wl)):
			f.write(str(wl[i])+'\n')
			for j in range(nn):
				f.write(str(dtau[i,j])+'\t'+str(waer[i,j])+'\t')
				for k in range(6):
					f.write(str(pmom[i,j,k])+'\t')
				f.write('\n')
	
	print('done')

if __name__ == '__main__':
	#read data once for all run
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
	'''
	clean_PNSD = np.nanmean(PNSD, axis=0)
	dirty_PNSD = PNSD[i_max]
	clean_BCPNSD = np.nanmean(BCPNSD, axis=0)
	dirty_BCPNSD = BCPNSD[iBC_max]
	'''
	clean_PNSD = np.nanmean(PNSD, axis=0) / 10
	dirty_PNSD = np.nanmean(PNSD, axis=0)
	clean_BCPNSD = np.nanmean(BCPNSD, axis=0) / 10
	dirty_BCPNSD = np.nanmean(BCPNSD, axis=0)
	
	#start running
	
	infos = np.load('output/all/20220808_all/infos.npy', allow_pickle=True)
	infos = infos.tolist()
	for i in range(10000):
		factors = run(20, 6, [525], 30, k_on=True, n_on=True, kBC_on=True, nBC_on=True, kappa_on=True, PF_on=True, MS_on=True, VD_on=True, CT_on=True, BCAE_on=True, albedo_on=True, AOD_on=True, SSA_on=True, PNSD_on=True, BCPNSD_on=True, AOD_range=0.1, output_path='output/all/20220818_all/', output_name=str(i+10000), Dps=Dps, DBC=DBC, DBCps=DBCps, clean_PNSD=clean_PNSD, dirty_PNSD=dirty_PNSD, clean_BCPNSD=clean_BCPNSD, dirty_BCPNSD=dirty_BCPNSD, debug=True)
		infos.append(factors)
		np.save('output/all/infos.npy', infos)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
