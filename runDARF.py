####################################################################################
# INTRODUCTION:
# This code is to test all factor turburlence sensitivity of DARF
# Created by Hebs at 24/2/20/14:22
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import sys
import os
import random
import time as tm
from pyDOE import lhs
from scipy.stats import norm
from scipy.stats import uniform
from multiprocessing import Pool

import readTaizhou
import readAtms
import read11
import readAlbedo
import write_atms
import write_aerosol
import write_INPUT
import write_albedo
import calRH
import calBC
import calPNSD
import calDARF

###############################################
# some constants
default_mag = 0.1 # Dps resolution adjust magnification

default_range = 0.3 # 30% turbulence

default_paras = [] # default parameters list
default_paras.append('n')
default_paras.append('nH1')
default_paras.append('nBC')
default_paras.append('kBC')
default_paras.append('PNSD')
default_paras.append('MS')
default_paras.append('MSH1')
default_paras.append('VD')
default_paras.append('CT')
default_paras.append('CTH1')
default_paras.append('kappa')
default_paras.append('kappaH1')
default_paras.append('rhoBC')
default_paras.append('BCPNSD')
default_paras.append('BCAE')
default_paras.append('amb')
default_paras.append('albedo')
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

def LHS_uniform(num, loc, scale):
	'''
	This function is Latin hypercubic sampling method code for uniform distribution
	input:
		num     : sample number, int
		loc     : unifrom distribution location parameter, float
		scale   : uniform distribution scale parameter, float
	output:
		lhd		: parameter list, array
	'''
	lhd = lhs(1, samples=num)
	lhd = uniform(loc=loc, scale=scale).ppf(lhd)
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

def make_dir(time, **args):
	'''
	This function is to make directory for this code output
	input:
		time         : year month and day, string, for example 220921
		**paras      : parameters list, array, in string, default all parameters
	output:
		parameter directory with time, in path output/DARF/
	'''
	if 'paras' in args:
		paras = args['paras']
	else:
		paras = default_paras
	
	path = 'output/DARF/' + time
	os.system('mkdir '+path)
	os.system('mkdir '+path+'/all')
	for i in range(len(paras)):
		fn = path + '/' + paras[i]
		os.system('mkdir '+fn)

def run(wl, nn, moma, **args):
	'''
	This function is to use Monte Carlo method to test all factors' sensitivity of DARF parameters
	input:
		wl		: the wavelength [ wl(k) < wl(k+1) ], np.array, nm
		nn		: number of atmospheric levels for which aerosol information is specified, int
		moma	: number of the phase function's legendre moments, int
		args including:
			factors control switch, XX_rate, XX_rate_on, XX_on and XX_range for each factor XX, 
			XX_rate for change rate set before running, 
			XX_rate use XX_rate_on (bool) to judge existence, default False,
			when XX_rate is set, XX_on is True,
			XX_on for change rate set when running, bool, default False,
			XX_range is needed when XX_on is True, float, default 0.3.
			all factors are:
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
				default_paras.append('BCAE')
				default_paras.append('amb')
				default_paras.append('albedo')
			function control factors, including:
				angularResolution	: phase function angular resolution, degree, default 0.5
				debug				: debug flag, bool, default False
				data_path			: input factor data path, string, default 'data/'
				output_path			: output file storage path, string, default 'output/DARF/'
				output_name			: output file name, array, default ['0.txt', '1.txt', 'atms.dat', 
									  'aerosol.dat', 'albedo.dat']
			Dps					: particle diameter size list, np.array, 
								  use Dps_exist (bool) to judge existence, default False
			clean_PNSD			: clean particle number size distribution, np.array, 
								  use clean_PNSD_exist (bool) to judge existence, default False
			dirty_PNSD			: dirty particle number size distribution, np.array, 
								  use dirty_PNSD_exist (bool) to judge existence, default False
			DBC					: BC-contain particle BC core diameter size distribution, 2-d np.array, 
								  use DBC_exist (bool) to judge existence, default False
			DBCps				: BC-contain particle diameter size list, np.array, 
								  use DBC_exist (bool) to judge existence, default False
			clean_BCPNSD		: clean BC-contain particle number size distribution, np.array, 
								  use clean_BCPNSD_exist (bool) to judge existence, default False
			dirty_BCPNSD		: dirty BC-contain particle number size distribution, np.array, 
								  use dirty_BCPNSD_exist (bool) to judge existence, default False
	output:
		infos	: DARF information, dict
		paras	: origin parameter value, dict
	'''
	
	# how to write this code:
	# use calNI.cal_Mie2(), we can get kext, waer and pmom for every level aerosol for SBDART input,
	# different level have different aerosol number distribution, RH, height and thickness,
	# makes every level has different dtau, waer and pmom.
	# in the running, an atms.dat and an aerosol.dat will form, to run SBDART.
	
	if 'rand' in args:
		rand = args['rand']
	else:
		np.random.seed(int(tm.time()))
		rand = np.random.rand()
	
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
	
	if 'VD_rate' in args:
		VD_rate_on = True
		VD_rate = args['VD_rate']
	else:
		VD_rate_on = False
	if 'VD_on' in args:
		VD_on = args['VD_on']
		if 'VD_range' in args:
			VD_range = args['VD_range']
		else:
			VD_range = default_range
	else:
		VD_on = False
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
	if 'amb_rate' in args:
		amb_rate_on = True
		amb_rate = args['amb_rate']
	else:
		amb_rate_on = False
	if 'amb_on' in args:
		amb_on = args['amb_on']
		if 'amb_range' in args:
			amb_range = args['amb_range']
		else:
			amb_range = default_range
	else:
		amb_on = False
	if 'albedo_rate' in args:
		albedo_rate_on = True
		albedo_rate = args['albedo_rate']
	else:
		albedo_rate_on = False
	if 'albedo_on' in args:
		albedo_on = args['albedo_on']
		if 'albedo_range' in args:
			albedo_range = args['albedo_range']
		else:
			albedo_range = default_range
	else:
		albedo_on = False
	
	if 'debug' in args:
		debug = args['debug']
	else:
		debug = False
	if 'angularResolution' in args:
		angularResolution = args['angularResolution']
	else:
		angularResolution = 0.5
	if 'data_path' in args:
		data_path = args['data_path']
	else:
		data_path = 'data/'
	if 'output_path' in args:
		output_path = args['output_path']
	else:
		output_path = 'output/DARF/'
	if 'output_names' in args:
		output_names = args['output_names']
	else:
		output_names = ['0.txt', '1.txt', 'atms.dat', 'aerosol.dat', 'albedo.dat']
	
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
			'''
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
			'''
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
	
	MSH1 = 1
	MSH2 = 1
	
	# read vertical distribution type. If not exist, using default data
	
	VD = 28 / (28+39) # from Aircraft study of aerosol vertical distributions over Beijing and their optical properties, P.Liu et al., 2009
	if debug:
		print('vertical distribution rate :', VD)
	
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
	
	nH1 = 1
	
	CTH1 = 1
	
	# read kappa data. If not exist, using default data
	
	kappa = 0.21537370768611558 # from Beijing University observation
	# according to Petters and Kreidenweis 2007, kappa can change from 0.002 to 0.67
	# set kappaI change range from 0.5 * kappa to 1.5 * kappa, equal to 0.1
	if debug:
		print('kappa :', kappa)
	
	kappaH1 = 1
	
	# read RH data, from top to bottom
	
	######################################################################################
	# new change 24/2/20:
	# use parameter amb to control ambiant parameters change,
	# the ambiant parameters would be:
	# par_amb = par_summer * amb + par_winter * (1-amb)
	# use new ambiant parameters to write file atms.dat
	######################################################################################
	
	# random number for unique running
	run_file_path = 'run_files/' + str(rand)
	os.system('mkdir '+run_file_path)
	
	# atms.dat change
	
	if amb_rate_on:
		amb = 0.5 * amb_rate
	elif amb_on:
		amb_rate = 1 + random.normalvariate(mu=0,sigma=1/3) * amb_range
		amb = 0.5 * amb_rate
	else:
		amb = 0.5
		amb_rate = 0
	
	write_atms.amb_write('data/era5/data.nc', 40, 116, amb, output=run_file_path+'/atms.dat') # Beijing
	
	# albedo.dat change
	
	write_albedo.write(116, 40, 2019, month=6, day=1, output=run_file_path+'/albedo.dat')
	wl_albedo, albedo = readAlbedo.read(fn=run_file_path+'/albedo.dat')
	albedo = np.array(albedo)
	
	if albedo_rate_on:
		rate = albedo_rate
	elif albedo_on:
		albedo_rate = 1 + random.normalvariate(mu=0,sigma=1/3) * BCPNSD_range
		rate = albedo_rate
	else:
		rate = 1
		albedo_rate = 0
	
	write_albedo.change(rate, fn=run_file_path+'/albedo.dat')
	
	# check the input data validity
	# if wl(k) < wl(k+1) doesn\'t meet, exit
	for i in range(len(wl)-1):
		if (wl[i+1]-wl[i])<0:
			print('wl(k) < wl(k+1) doesn\'t meet. Please check.')
			sys.exit()
	
	# save origin factor value
	
	paras = dict(n=nShell, nH1=nH1, nBC=nBC, kBC=kBC, PNSD=0.5, MS=MS, MSH1=MSH1, MSH2=MSH2, VD=VD, CT=1, CTH1=CTH1, kappa=kappa, kappaH1=kappaH1, rhoBC=rhoBC, BCPNSD=0.5, BCAE=1, amb=0.5, albedo=albedo)
	
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
	
	######################################################################
	# adjust Dps resolution to simplify calculate pressure
	PNSD, Dps = calPNSD.PNSD_Dps_adjust(PNSD, Dps, default_mag)
	new_BCPNSD = np.zeros((len(DBCps),round(len(DBC)*default_mag*2)))
	for i in range(len(DBCps)):
		new_BCPNSD[i], new_DBC = calPNSD.PNSD_Dps_adjust(BCPNSD[i], DBC, default_mag*2)
	BCPNSD = new_BCPNSD
	DBC = new_DBC
	######################################################################
	
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
		MSH1 = MSH1 * (MSH1_rate-1)
	elif MSH1_on:
		MSH1_rate = 1 + random.normalvariate(mu=0,sigma=1/3) * MSH1_range
		MSH1 = MSH1 * (MSH1_rate-1)
	else:
		MSH1_rate = 0
		MSH1 = 0
	if debug:
		print('MS heterogeneity 1 : ', MSH1)
	
	if MSH2_rate_on:
		MSH2y = MSH2 * (MSH2_rate-1)
	elif MSH2_on:
		MSH2_rate = 1 + random.normalvariate(mu=0,sigma=1/3) * MSH2_range
		MSH2y = MSH2 * (MSH2_rate-1)
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
	
	if VD_rate_on:
		VD = VD * VD_rate
		if VD > 1:
			VD = 1
	elif VD_on:
		VD_rate = 1 + random.normalvariate(mu=0,sigma=1/3) * VD_range
		VD = VD * VD_rate
		if VD > 1:
			VD = 1
	else:
		VD_rate = 0
	if debug:
		print('Vertical distribution:', VD)
	
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
		CTH1 = CTH1 * (CTH1_rate-1)
	elif CTH1_on:
		CTH1_rate = 1 + random.normalvariate(mu=0,sigma=1/3) * CTH1_range
		CTH1 = CTH1 * (CTH1_rate-1)
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
	
	dtau, waer, pmom, g = calDARF.DARF(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, nH1, kappa, kappaH1, MS, MSH1, MSH2x, MSH2y, VD, CT, CTH1, BCAE, wl, nn, moma, angularResolution=angularResolution, debug=debug, atms_fn=run_file_path+'/atms.dat')
	'''
	import calMie
	AOD, SSA, g = calMie.Mie(Dps, PNSD, DBCps, DBC, BCPNSD, nBC, kBC, nShell, nI, nI2x, nI2y, kappa, kappaI, kappaI2x, kappaI2y, 60, 525, 2000, 100, BCAE, CT, VD, 1, debug=debug)
	print('AOD ', AOD,' SSA ', SSA,' g ', g)
	'''
	if debug:
		print('bottom dtau :', dtau[0,-1], ',\tbottom waer :', waer[0,-1], ',\tbottom pmom :', pmom[0,-1], ',\tbottom g :', g[0,-1])
	
	# write data into aerosol.dat
	
	write_aerosol.write2(wl, nn, moma, dtau, waer, pmom, output=run_file_path+'/aerosol.dat')
	
	# run SBDART
	
	write_INPUT.write(iaer=0, fn=run_file_path+'/INPUT')
	os.system('cd '+run_file_path+'\nsbdart > 0.txt')
	write_INPUT.write(fn=run_file_path+'/INPUT')
	os.system('cd '+run_file_path+'\nsbdart > 1.txt')
	infos0 = read11.read11(filename=run_file_path+'/0.txt')
	infos1 = read11.read11(filename=run_file_path+'/1.txt')
	RF = (infos1['fxdn']-infos1['fxup']) - (infos0['fxdn']-infos0['fxup'])
	heat1 = infos1['heat']
	heat0 = infos0['heat']
	
	if debug:
		print('RF top = ', RF[0])
		print('RF bottom = ', RF[-1])
	
	# back up the file
	
	os.system('cp '+run_file_path+'/0.txt '+output_path+output_names[0])
	os.system('cp '+run_file_path+'/1.txt '+output_path+output_names[1])
	os.system('cp '+run_file_path+'/atms.dat '+output_path+output_names[2])
	os.system('cp '+run_file_path+'/aerosol.dat '+output_path+output_names[3])
	os.system('cp '+run_file_path+'/albedo.dat '+output_path+output_names[4])
	
	os.system('rm -r '+run_file_path)
	
	infos = dict(RF_top=RF[0], RF_bot=RF[-1], RF=RF, heat1=heat1, heat0=heat0, dtau_bot=dtau[0,-1], dtau=dtau, waer_bot=waer[0,-1], waer=waer, pmom_bot=pmom[0,-1], pmom=pmom, g_bot=g[0,-1], g=g, n_rate=n_rate, nH1_rate=nH1_rate, nBC_rate=nBC_rate, kBC_rate=kBC_rate, PNSD_rate=PNSD_rate, MS_rate=MS_rate, MSH1_rate=MSH1_rate, MSH2_rate=MSH2_rate, VD_rate=VD_rate, CT_rate=CT_rate, CTH1_rate=CTH1_rate, kappa_rate=kappa_rate, kappaH1_rate=kappaH1_rate, rhoBC_rate=rhoBC_rate, BCPNSD_rate=BCPNSD_rate, BCAE_rate=BCAE_rate, amb_rate=amb_rate, albedo_rate=albedo_rate)
	return infos, paras

def LHS_run(time, run_num, par_range, **args):
	'''
	This function is to circulation run Mie in Latin hypercubic sampling method
	input:
		time      : running time to record
		run_num   : running times number
		par_range : parameters change range
		**debug   : debug flag, bool, default False
	output:
		aerosol direct radiative forcing information files, in w/m^2, float
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
	clean_BCPNSD = BCPNSD / 50
	dirty_BCPNSD = BCPNSD / 10
	
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
	
	path = 'output/DARF/' + time + '/'
	
	infos = []
	output_path = path + 'all/'
	for i in range(run_num):
		output_names = ['0_'+str(i)+'.txt', '1_'+str(i)+'.txt', 'atms_'+str(i)+'.dat', 'aerosol_'+str(i)+'.dat', 'albedo_'+str(i)+'.dat']
		info, paras = run([525], 15, 6,
		n_rate			= all_rate[0,i], 
		nH1_rate		= all_rate[1,i], 
		nBC_rate		= all_rate[2,i], 
		kBC_rate		= all_rate[3,i], 
		PNSD_rate		= all_rate[4,i], 
		MS_rate			= all_rate[5,i], 
		MSH1_rate		= all_rate[6,i], 
		VD_rate			= all_rate[7,i], 
		CT_rate			= all_rate[8,i], 
		CTH1_rate		= all_rate[9,i], 
		kappa_rate		= all_rate[10,i], 
		kappaH1_rate	= all_rate[11,i], 
		rhoBC_rate		= all_rate[12,i], 
		BCPNSD_rate		= all_rate[13,i], 
		BCAE_rate		= all_rate[14,i], 
		amb_rate		= all_rate[15,i], 
		albedo_rate		= all_rate[16,i], 
		Dps=Dps, DBC=DBC, DBCps=DBCps, clean_PNSD=clean_PNSD, dirty_PNSD=dirty_PNSD, clean_BCPNSD=clean_BCPNSD, dirty_BCPNSD=dirty_BCPNSD, angularResolution=30, debug=debug, output_path=output_path, output_names=output_names)
		infos.append(info)
		np.save(path+'all/infos.npy', infos)
		if i==0:
			np.save(path+'all/paras.npy', paras)
	
	for i in range(par_num):
		infos = []
		output_path = path + parameters[i] + '/'
		for j in range(run_num):
			output_names = ['0_'+str(j)+'.txt', '1_'+str(j)+'.txt', 'atms_'+str(j)+'.dat', 'aerosol_'+str(j)+'.dat', 'albedo_'+str(j)+'.dat']
			info, paras = run([525], 15, 6,
			n_rate			= paras_rate[i,0,j], 
			nH1_rate		= paras_rate[i,1,j], 
			nBC_rate		= paras_rate[i,2,j], 
			kBC_rate		= paras_rate[i,3,j], 
			PNSD_rate		= paras_rate[i,4,j], 
			MS_rate			= paras_rate[i,5,j], 
			MSH1_rate		= paras_rate[i,6,j], 
			VD_rate			= paras_rate[i,7,j], 
			CT_rate			= paras_rate[i,8,j], 
			CTH1_rate		= paras_rate[i,9,j], 
			kappa_rate		= paras_rate[i,10,j], 
			kappaH1_rate	= paras_rate[i,11,j], 
			rhoBC_rate		= paras_rate[i,12,j], 
			BCPNSD_rate		= paras_rate[i,13,j], 
			BCAE_rate		= paras_rate[i,14,j], 
			amb_rate		= paras_rate[i,15,j], 
			albedo_rate		= paras_rate[i,16,j], 
			Dps=Dps, DBC=DBC, DBCps=DBCps, clean_PNSD=clean_PNSD, dirty_PNSD=dirty_PNSD, clean_BCPNSD=clean_BCPNSD, dirty_BCPNSD=dirty_BCPNSD, angularResolution=30, debug=debug, output_path=output_path, output_names=output_names)
			infos.append(info)
			np.save(path+parameters[i]+'/infos.npy', infos)

def mp_run(parameters):
	'''
	Use Pool.map() to finish multi proccessing
	'''
	# read parameter from parameters
	wl = parameters['wl']
	nn = parameters['nn']
	moma = parameters['moma']
	
	n_rate			= parameters['rates'][0]
	nH1_rate		= parameters['rates'][1]
	nBC_rate		= parameters['rates'][2]
	kBC_rate		= parameters['rates'][3]
	PNSD_rate		= parameters['rates'][4]
	MS_rate			= parameters['rates'][5]
	MSH1_rate		= parameters['rates'][6]
	VD_rate			= parameters['rates'][7]
	CT_rate			= parameters['rates'][8]
	CTH1_rate		= parameters['rates'][9]
	kappa_rate		= parameters['rates'][10]
	kappaH1_rate	= parameters['rates'][11]
	rhoBC_rate		= parameters['rates'][12]
	BCPNSD_rate		= parameters['rates'][13]
	BCAE_rate		= parameters['rates'][14]
	amb_rate		= parameters['rates'][15]
	albedo_rate		= parameters['rates'][16]
	
	Dps=parameters['Dps']
	DBC=parameters['DBC']
	DBCps=parameters['DBCps']
	clean_PNSD=parameters['clean_PNSD']
	dirty_PNSD=parameters['dirty_PNSD']
	clean_BCPNSD=parameters['clean_BCPNSD']
	dirty_BCPNSD=parameters['dirty_BCPNSD']
	
	angularResolution=parameters['angularResolution']
	debug=parameters['debug']
	output_path=parameters['output_path']
	output_names=parameters['output_names']
	rand=parameters['rand']
	
	infos, paras = run(wl, nn, moma, 
	n_rate			= n_rate, 
	nH1_rate		= nH1_rate, 
	nBC_rate		= nBC_rate, 
	kBC_rate		= kBC_rate, 
	PNSD_rate		= PNSD_rate, 
	MS_rate			= MS_rate, 
	MSH1_rate		= MSH1_rate, 
	VD_rate			= VD_rate, 
	CT_rate			= CT_rate, 
	CTH1_rate		= CTH1_rate, 
	kappa_rate		= kappa_rate, 
	kappaH1_rate	= kappaH1_rate, 
	rhoBC_rate		= rhoBC_rate, 
	BCPNSD_rate		= BCPNSD_rate, 
	BCAE_rate		= BCAE_rate, 
	amb_rate		= amb_rate, 
	albedo_rate		= albedo_rate, 
	Dps=Dps, DBC=DBC, DBCps=DBCps, clean_PNSD=clean_PNSD, dirty_PNSD=dirty_PNSD, clean_BCPNSD=clean_BCPNSD, dirty_BCPNSD=dirty_BCPNSD, 
	angularResolution=angularResolution, debug=debug, output_path=output_path, output_names=output_names, rand=rand)
	
	return infos, paras

if __name__ == '__main__':
	#LHS_run('240612', 2, 0.05, debug=True)
	
	# multi processing running
	############################################################
	process_num = 40
	time = '241012'
	run_num = 2000
	############################################################
	par_range = 0.1
	debug = True
	wl = [525]
	nn = 15
	moma = 6
	angularResolution = 30
	
	#read data
	
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
	clean_BCPNSD = BCPNSD / 50
	dirty_BCPNSD = BCPNSD / 10
	
	# use LHS produce all parameters for whole run
	
	parameters = default_paras
	
	par_num = len(parameters)
	all_rate = np.zeros((par_num, run_num))
	paras_rate = np.zeros((par_num, par_num, run_num))
	
	# for all run and parameter run
	
	for i in range(par_num):
		# sigma_uniform = (b-a)/2/sqrt(3) = 1/3, b-a = scale = 2/sqrt(3)
		#all_rate[i] = 1 + LHS_uniform(run_num, -0.577, 1.155).flatten() * default_range
		all_rate[i] = 1 + LHS_norm(run_num, 0, 1/3).flatten() * default_range
	
	for i in range(par_num):
		for j in range(par_num):
			if i == j:
				# sigma_uniform = (b-a)/2/sqrt(3) = 1/3, b-a = scale = 2/sqrt(3)
				paras_rate[i,j] = 1 + LHS_norm(run_num, 0, 1/3).flatten() * par_range
			else:
				paras_rate[i,j] = all_rate[j]
	
	make_dir(time)
	
	#start running
	
	path = 'output/DARF/' + time + '/'
	output_path = path + 'all/'
	np.random.seed(int(tm.time()))
	rand = np.random.rand(run_num)
	
	# create parameters
	
	ps = []
	
	for i in range(run_num):
		output_names = ['0_'+str(i)+'.txt', '1_'+str(i)+'.txt', 'atms_'+str(i)+'.dat', 'aerosol_'+str(i)+'.dat', 'albedo_'+str(i)+'.dat']
		ps.append(dict(wl=wl, nn=nn, moma=moma, rates=all_rate[:,i], Dps=Dps, DBC=DBC, DBCps=DBCps, clean_PNSD=clean_PNSD, dirty_PNSD=dirty_PNSD, clean_BCPNSD=clean_BCPNSD, dirty_BCPNSD=dirty_BCPNSD, angularResolution=angularResolution, debug=debug, output_path=output_path, output_names=output_names, rand=rand[i]))
	
	# run multi processs
	
	with Pool(processes=process_num) as pool:
		results = pool.map(mp_run, ps)
	infos = np.array(results)[:,0]
	paras = np.array(results)[0,1]
	np.save(path+'all/infos.npy', infos)
	np.save(path+'all/paras.npy', paras)
	
	for i in range(par_num):
		output_path = path + parameters[i] + '/'
		ps_i = []
		np.random.seed(int(tm.time()))
		rand = np.random.rand(run_num)
		for j in range(run_num):
			output_names = ['0_'+str(j)+'.txt', '1_'+str(j)+'.txt', 'atms_'+str(j)+'.dat', 'aerosol_'+str(j)+'.dat', 'albedo_'+str(j)+'.dat']
			ps_i.append(dict(wl=wl, nn=nn, moma=moma, rates=paras_rate[i,:,j], Dps=Dps, DBC=DBC, DBCps=DBCps, clean_PNSD=clean_PNSD, dirty_PNSD=dirty_PNSD, clean_BCPNSD=clean_BCPNSD, dirty_BCPNSD=dirty_BCPNSD, angularResolution=angularResolution, debug=debug, output_path=output_path, output_names=output_names, rand=rand[j]))
			
		with Pool(processes=process_num) as pool:
			results = pool.map(mp_run, ps_i)
		infos = np.array(results)[:,0]
		np.save(path+parameters[i]+'/infos.npy', infos)
	print(tm.asctime())
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
