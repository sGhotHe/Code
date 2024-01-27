####################################################################################
# INTRODUCTION:
# This code is to calculate bulk phase function, and to use phase function to calculate Legendre moments
# Created by Hebs at 21/5/28/12:37
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import PyMieScatt as ps
import readTaizhou
import calPNSD
import calBC
import kappaKohler
import os

def cal_phaseFunc(Dp, DBC, m_BC, m_shell, kappa, RH, wl, **args):
	'''
	This function is to calculate phaseFunc, defined as phaseFunc_i = Q_sca_i * 1/4 * pi * Dps_i**2 * P_i(theta)
	input:
		Dp                  : particle diameter, float, nm
		DBC                 : BC core diameter, float, nm
		m_BC                : BC core complex refractive index
		m_shell             : shell complex refractive index
		kappa               : hygroscopicity parameter, float
		RH                  : relative humidity, float, percent
		wl                  : wave length, float, nm
		**angularResolution : angular resolution, degree, default 0.5
	output:
		phaseFunc, array, defined as phaseFunc_i = Q_sca_i * 1/4 * pi * Dps_i**2 * P_i(theta)
	'''
	if 'angularResolution' in args:
		angularResolution = args['angularResolution']
	else:
		angularResolution = 0.5
	
	if Dp>DBC:
		D_RH, m_shell_RH = kappaKohler.RH2D(Dp, DBC, m_BC, m_shell, kappa, RH, wl) # hygroscopic growth
	else:
		D_RH = Dp
		m_shell_RH = m_shell
	if DBC>1:
		sigma_sca = ps.MieQCoreShell(m_BC, m_shell_RH, wl, DBC, D_RH, asCrossSection=True)[1]
	else:
		sigma_sca = ps.MieQ(m_shell_RH, wl, D_RH, asCrossSection=True)[1]
	P = ps.CoreShellScatteringFunction(m_BC, m_shell_RH, wl, DBC, D_RH, angularResolution=angularResolution, normed=True)[3]
	phaseFunc = sigma_sca * P
	return phaseFunc

def write_phaseFuncData(Dps, PNSD, DBCps, DBC, DMASP2, m_BC, m_shell, kappa, wls, **args):
	'''
	This function is to write phaseFuncData.dat for bulk phase function calculation, in different wave length and RH
	input:
		Dps                 : diameter distribution, array, nm
		PNSD                : number distribution, array, dn/dlogDp
		DBCps               : BC particle diameter size, array, nm
		DBC                 : BC core diameter size, array, nm
		DMASP2              : BC particle number size distribution, array, dn/dlogDBC
		m_BC                : BC core complex refractive index
		m_shell             : shell complex refractive index
		kappa               : hygroscopicity parameter, float
		wls                 : wave length list, array, nm
		**angularResolution : angular resolution, degree, default 0.5
		**fn                : data store path, string, default 'data/phaseFuncData/'
		**RHs               : raletive humidity list, array, percent, default 1 to 99
	output:
		phase function data
	'''
	if 'angularResolution' in args:
		angularResolution = args['angularResolution']
	else:
		angularResolution = 0.5
	if 'fn' in args:
		fn = args['fn']
	else:
		fn = 'data/phaseFuncData/'
	if 'RHs' in args:
		RHs = args['RHs']
	else:
		RHs = np.arange(1,100)
	
	# sum_i(Q_sca_i * 1/4 * pi * Dps_i**2 * n_i * P_i(theta)) = P_bulk(theta)
	# sum_theta(P_i(theta)) = 1 for each i
	# except for n_i, other item in equations are only related to Dps
	# set phaseFunc_i = Q_sca_i * 1/4 * pi * Dps_i**2 * P_i(theta)
	# P_bulk(theta) = sum_i(phaseFunk_i * n_i)
	
	# data title format:
	#	phaseFuncData_[wl]
	# data storage format: 
	#	Dps[0]\tDps[1]\t...Dps[-1]
	#	RH[0]\tRH[1]\t...RH[-1]
	#	wl\tangularNumber\tkappa
	#	phaseFunc[0,0,0]\tphaseFunc[0,0,1]\t...phaseFunc[0,0,angularNumber-1]
	#	phaseFunc[0,1,0]\tphaseFunc[0,1,1]\t...phaseFunc[0,1,angularNumber-1]
	#	...
	#	phaseFunc[0,len(Dps)-1,0]\tphaseFunc[0,len(Dps)-1,1]\t... phaseFunc[0,len(Dps)-1,angularNumber-1]
	#	...
	#	...
	#	...
	#	phaseFunc[len(RH)-1,0,0]\t phaseFunc[len(RH)-1,0,1]\t ... phaseFunc[len(RH)-1,0,angularNumber-1]
	#	phaseFunc[len(RH)-1,1,0]\t phaseFunc[len(RH)-1,1,1]\t ... phaseFunc[len(RH)-1,1,angularNumber-1]
	#	...
	#	phaseFunc[len(RH)-1,len(Dps)-1,0]\t phaseFunc[len(RH)-1,len(Dps)-1,1]\t ... phaseFunc[len(RH)-1,len(Dps)-1,angularNumber-1]
	
	# to exceed write and read, use .npy to store data
	# data storage format:
	# 	Dps, RHs, wl, angularNumber, kappa, phaseFunc
	
	# revise at 22/9/19/14:48:
	# METHOD MENTHONED ABOVE HAVE SERIOUS PROBLEM!
	# phase function related to Dps and DBCps distribution both!
	
	angularNumber = round(180/angularResolution)
	phaseFunc = np.zeros((len(wls), len(RHs), len(Dps), angularNumber))
	DMASP2PNSD = calPNSD.DMASP22DMASP2PNSD(DBCps, DMASP2)
	n = calPNSD.PNSD2n(Dps, PNSD, 'd')
	theta = np.arange(angularNumber) * angularResolution * np.pi / 180 # in rad
	
	print('calculating...')
	
	# for core-shell structure, different calculate
	for ii in range(len(wls)):
		for jj in range(len(RHs)):
			for i in range(len(Dps)):
				if n[i]: # PNSD not zero, otherwise no need to calculate
					DMASP2PNSD_i = calPNSD.DMASP2PNSD_Dp(DBCps, DMASP2PNSD, Dps[i])
					DMASP2_i = calPNSD.PNSD2n(Dps, DMASP2PNSD_i, 'e')
					n_BC_i = calPNSD.PNSD2n(DBC, DMASP2_i, 'e')
					for j in range(len(DBC)):
						if n_BC_i[j]: # n_BC not zero
							phaseFunc_BC_ij = cal_phaseFunc(Dps[i], DBC[j], m_BC, m_shell, kappa, RHs[jj], wls[ii], angularResolution=angularResolution)
							# calculate phaseFunc for bulk
							phaseFunc_BC_ij = phaseFunc_BC_ij * n_BC_i[j] / n[i]
							phaseFunc[ii,jj,i] += phaseFunc_BC_ij
					n_noBC_i = n[i] - sum(n_BC_i)
					phaseFunc_noBC_i = cal_phaseFunc(Dps[i], Dps[i], m_shell, m_shell, kappa, RHs[jj], wls[ii], angularResolution=angularResolution)
					phaseFunc_noBC_i = phaseFunc_noBC_i * n_noBC_i / n[i]
					phaseFunc[ii,jj,i] += phaseFunc_noBC_i
			
			print(round((jj+1+ii*len(RHs))/(len(wls)*len(RHs))*1000)/10,'% done...')
		
		print('done\nwriting...')
		data = dict(Dps=Dps, RHs=RHs, angularNumber=angularNumber, kappa=kappa, phaseFunc=phaseFunc)
		np.save(fn+str(wls[ii])+'_'+str(round(kappa*1e3))+'.npy', data)
		print(wls[ii], 'done')
	
	print('all done')

def cal_bulk_phase_function(Dps, PNSD, DBCps, DBC, DMASP2, m_BC, m_shell, kappa, RH, wl, **args):
	'''
	This function is to use size distribution to calculate bulk phase funcion in certain RH
	input:
		Dps                 : diameter distribution, array, nm
		PNSD                : number distribution, array, dn/dlogDp
		DBCps               : BC particle diameter size, array, nm
		DBC                 : BC core diameter size, array, nm
		DMASP2              : BC particle number size distribution, array, dn/dlogDBC
		m_BC                : BC core complex refractive index
		m_shell             : shell complex refractive index
		kappa               : hygroscopicity parameter, float
		RH                  : raletive humidity, float, percent
		wl                  : wave length, float, nm
		**angularResolution : angular resolution, default 0.5, degree
		**fn                : data store file path, default data/phaseFuncData/
		**RHs               : raletive humidity list, array, percent, default np.arange(100), resolution must smaller than 1%
	output:
		P_bulk              : normalized bulk phase function, array
		theta               : angle list, array, rad
	'''
	if 'fn' in args:
		fn = args['fn']
	else:
		fn = 'data/phaseFuncData/'
	if 'angularResolution' in args:
		angularResolution = args['angularResolution']
	else:
		angularResolution = 0.5
	if 'RHs' in args:
		RHs = args['RHs']
	else:
		RHs = np.arange(1,100)
	
	# sum_i(Q_sca_i * 1/4 * pi * Dps_i**2 * n_i * P_i(theta)) = P_bulk(theta)
	# sum_theta(P_i(theta)) = 1 for each i
	# except for n_i, other item in equations are only related to Dps
	# set phaseFunc_i = Q_sca_i * 1/4 * pi * Dps_i**2 * P_i(theta)
	# P_bulk(theta) = sum_i(phaseFunk_i * n_i)
	
	# to exceed the calculation, created a file to store phaseFunc data
	# the following calculate can read from file directly
	# data storage format: 
	#	Dps[0]\tDps[1]\t...Dps[-1]
	#	angularNumber
	#	wl
	#	phaseFunc[0,0]\tphaseFunc[0,1]\t...phaseFunc[0,angularNumber-1]
	#	phaseFunc[1,0]\tphaseFunc[1,1]\t...phaseFunc[1,angularNumber-1]
	#	...
	#	phaseFunc[len(Dps)-1,0]\tphaseFunc[len(Dps)-1,1]\t... phaseFunc[len(Dps)-1,angularNumber-1]
	
	# to further exceed read and write, data stored as .npy
	
	angularNumber = round(180/angularResolution)
	theta = np.arange(angularNumber) * angularResolution * np.pi / 180 # in rad
	n = calPNSD.PNSD2n(Dps, PNSD, index='e')
	
	exist = False # to judge whether file exist
	P_bulk = 0
	
	if os.path.exists(fn+str(wl)+'_'+str(round(kappa*1e3))+'.npy'):
		data = np.load(fn+str(wl)+'_'+str(round(kappa*1e3))+'.npy', allow_pickle=True).item()
		if Dps[0]==data['Dps'][0] and Dps[-1]==data['Dps'][-1] and len(Dps)==len(data['Dps']) and RH>=min(data['RHs']) and RH<=max(data['RHs']) and angularNumber==data['angularNumber']: # file exists
			exist = True
			phaseFunc = data['phaseFunc']
	if exist==False:
		print('No match data, calculating...')
		write_phaseFuncData(Dps, PNSD, DBCps, DBC, DMASP2, m_BC, m_shell, kappa, [wl], angularResolution=angularResolution, fn=fn, RHs=RHs)
		data = np.load(fn+str(wl)+'_'+str(round(kappa*1e3))+'.npy', allow_pickle=True).item()
		phaseFunc = data['phaseFunc']
	
	# ii for RH index
	for i in range(len(RHs)-1):
		if RHs[i]<=RH and RHs[i+1]>RH:
			ii = i
			break
	for i in range(len(Dps)):
		P_bulk += phaseFunc[0,ii,i] * n[i]
	
	P_bulk = P_bulk / sum(P_bulk) # normalized
	return P_bulk, theta

def angle2theta(angles):
	thetas = angles / 180 * np.pi
	return thetas

def P(n, x): # Legendre polinomial, Pn(x)
	if n==0:
		return 1
	elif n==1:
		return x
	p0 = 1
	p1 = x
	for i in range(n-1):
		tmp = ((2*i+3)*x*p1-(i+1)*p0) / (i+2)
		p0 = p1
		p1 = tmp
	return p1

def beta(f, theta, n): # phase function to Legendre moments
	'''
	This function is to use phase function to calculate Legendre moments
	input:
		f     : phase function, angles from small to big
		theta : angle list, angles from small to big, rad
		n     : moment
	output:
		beta  : Legendre moments
	'''
	mus = np.cos(theta)
	a = 0
	b = 0
	for i in range(len(f)-1):
		a += f[i] * P(n, mus[i]) * (mus[i]-mus[i+1])
		b += f[i] * (mus[i]-mus[i+1])
	beta = a / b
	return beta

def cal_beta_old(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, wl, H, n, **args):
	'''
	This function is to use size distribution data to calculate Legendre moments
	input:
		Dps     : diameter distribution, array, nm
		PNSD    : number distribution, array, dn/dlogDp
		DBCps   : BC particle diameter size, array, nm
		DBC     : BC core diameter size, array, nm
		n_BC    : BC particle number concentration, array, cm^-3
		m_BC    : BC core complex refractive index
		m_shell : shell complex refractive index
		wl      : wave length, nm
		H       : height, m
		n       : moment
		**func  : vertical distribution function type, A or B, default B
	output:
		beta    : Legendre moments
	'''
	if 'func' in args:
		func = args['func']
	else:
		func = 'B'
	
	ratio = calBC.cal_vertical_distribution_ratio(n_BC, H, func=func)
	P_bulk, theta = cal_bulk_phase_function(Dps, PNSD*ratio, DBCps, DBC, n_BC*ratio, m_BC, m_shell, wl, angularResolution=10)
	b = beta(P_bulk, theta, n)
	
	return b

def cal_beta(Dps, PNSD, DBCps, DBC, DMASP2, m_BC, m_shell, kappa, RH, wl, z, VD, n, **args):
	'''
	This function is to use size distribution data to calculate Legendre moments
	input:
		Dps                 : diameter distribution, array, nm
		PNSD                : number distribution, array, dn/dlogDp
		DBCps               : BC particle diameter size, array, nm
		DBC                 : BC core diameter size, array, nm
		n_BC                : BC particle number concentration, array, cm^-3
		m_BC                : BC core complex refractive index
		m_shell             : shell complex refractive index
		kappa               : hydroscopic parameter, float
		RH                  : relative humidity, %, float
		wl                  : wave length, nm
		z                   : height, m
		VD                  : vertical distribution mixing, 0 for type A and 1 for type B
		n                   : moment, int
		**angularResolution : angular resolution, default 0.5, degree
		**fn                : data store file path, default data/phaseFuncData/
		**RHs               : raletive humidity list, array, percent, default
	output:
		beta                : Legendre moments
	'''
	if 'fn' in args:
		fn = args['fn']
	else:
		fn = 'data/phaseFuncData/'
	if 'angularResolution' in args:
		angularResolution = args['angularResolution']
	else:
		angularResolution = 0.5
	if 'RHs' in args:
		RHs = args['RHs']
	else:
		RHs = np.arange(1,100)
	
	ratio = calPNSD.cal_VD(Dps, PNSD, z, VD)
	P_bulk, theta = cal_bulk_phase_function(Dps, PNSD*ratio, DBCps, DBC, DMASP2*ratio, m_BC, m_shell, kappa, RH, wl, fn=fn, angularResolution=angularResolution, RHs=RHs)
	b = beta(P_bulk, theta, n)
	
	return b

def cal_g(P, theta):
	'''
	This function is to calculate asymmetry factor g
	input:
		P     : phase function, angles from small to big
		theta : angle list, angles from small to big, rad
	output:
		g     : asymmetry factor
	'''
	mus = np.cos(theta)
	a = 0
	b = 0
	for i in range(len(theta)-1):
		a += mus[i] * P[i] * (theta[i+1]-theta[i])
		b += P[i] * (theta[i+1]-theta[i])
	g = a / b
	return g

def g_change(P, theta, rate, **args):
	'''
	This function is to change asymmetry factor g in certain rate
	input:
		P     : phase function, array
		theta : angles, array, from small to big, rad
		rate  : change delta rate, float
		        ATTENTION: this change rate NOT equal to g change rate, have a little difference
		**af  : adjust factor, to make rate close to g change rate, default 1
	output:
		P_new : changed phase function, array
	'''
	if 'af' in args:
		af = args['af']
	else:
		af = 1
	
	# sum(P) must not change after change
	# use linear change: from forward to backward, change rate from [rate] to -[rate]
	# change rate related to cos(theta)
	
	sum_P = sum(P) # to keep normalized
	P_new = np.zeros(len(P))
	
	for i in range(len(theta)):
		P_new[i] = P[i] * ((rate-1)*np.cos(theta[i])*af+1)
	
	P_new = P_new * sum_P / sum(P)
	
	return P_new

def cal_g_change_af(P, theta, rate, **args):
	'''
	This function is to calculate g_change() parameter af
	input:
		P              : phase function, array
		theta          : angles, array, from small to big, rad
		rate           : change delta rate, float
		**af_min       : minimum of af, default 1
		**af_max       : maximum of af, default 4
		**af_bin_rough : rough bin of af, default 0.2
		**af_bin_fine  : fine bin of af, default 0.01
	output:
		af             : adjust factor, to make rate close to g change rate
	'''
	if 'af_min' in args:
		af_min = args['af_min']
	else:
		af_min = 0.2
	if 'af_max' in args:
		af_max = args['af_max']
	else:
		af_max = 5
	if 'af_bin_rough' in args:
		af_bin_rough = args['af_bin_rough']
	else:
		af_bin_rough = 0.01
	if 'af_bin_fine' in args:
		af_bin_fine = args['af_bin_fine']
	else:
		af_bin_fine = 0.0001
	
	g_old = cal_g(P, theta) * rate # target
	
	af = np.arange(af_min, af_max+af_bin_rough, af_bin_rough)
	ii = 0 # to mark the closest value subscript
	dif = 999999
	
	for i in range(len(af)):
		P_new = g_change(P, theta, rate, af=af[i])
		g_new = cal_g(P_new, theta)
		dg = abs(g_new-g_old)
		if dg<dif:
			dif = dg
			ii = i
	
	af_min_new = af_min + af_bin_rough * (ii-1)
	af_max_new = af_min + af_bin_rough * (ii+1) # select two bin to do fine calculate
	af = np.arange(af_min_new, af_max_new+af_bin_fine, af_bin_fine)
	ii = 0
	
	for i in range(len(af)):
		P_new = g_change(P, theta, rate, af=af[i])
		g_new = cal_g(P_new, theta)
		dg = abs(g_new-g_old)
		if dg<dif:
			dif = dg
			ii = i
	
	return af[ii]

def cal_beta_change_RH(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, kappa, RH, wl, H, n, rate, **args):
	'''
	This function is to use size distribution data and change rate to calculate Legendre moments change in certain RH
	input:
		Dps     : diameter distribution, array, nm
		PNSD    : number distribution, array, dn/dlogDp
		DBCps   : BC particle diameter size, array, nm
		DBC     : BC core diameter size, array, nm
		n_BC    : BC particle number concentration, array, cm^-3
		m_BC    : BC core complex refractive index
		m_shell : shell complex refractive index
		kappa   : hygroscopicity parameter, float
		RH      : relative humidity, float, percent
		wl      : wave length, nm
		H       : height, m
		n       : moment
		rate    : change rate, float
		**func  : vertical distribution function type, A or B, default B
	output:
		beta    : Legendre moments
	'''
	if 'func' in args:
		func = args['func']
	else:
		func = 'B'
	
	ratio = calBC.cal_vertical_distribution_ratio(n_BC, H, func=func)
	P_bulk, theta = cal_bulk_phase_function_RH(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, kappa, RH, wl, angularResolution=10)
	g = cal_g(P_bulk,theta)
	af = cal_g_change_af(P_bulk, theta, rate)
	P_bulk_new = g_change(P_bulk, theta, rate, af=af)
	b = beta(P_bulk_new, theta, n)
	
	return b

if __name__ == '__main__':
	'''
	print('loading...')
	data = readTaizhou.read_Taizhou('data/sp2/Taizhou.npy')
	Dps = data['Dps']
	PNSD = data['PNSD'][1001]
	DBCps = data['DBCps']
	DBC = data['DBC']
	DMASP2 = data['DMASP2'][1001]
	BCPNSD = data['BCPNSD'][1001]
	kabs = data['kabs'][1001]
	ksca = data['ksca'][1001]
	wl_abs = data['wl_abs']
	wl_sca = data['wl_sca']
	
	print('done')
	m_BC = 1.532+0.521j
	m_shell = 1.45+1e-7j
	kappa = 0.2
	wls = [440, 500, 870, 1640]
	write_phaseFuncData(Dps, PNSD, DBCps, DBC, DMASP2, m_BC, m_shell, kappa, wls, angularResolution=10)
	print(cal_bulk_phase_function(Dps, PNSD, DBCps, DBC, DMASP2, m_BC, m_shell, kappa, 30, 440), angularResolution=10)
	'''
	data = np.load('data/phaseFuncData/525_62.npy', allow_pickle=True).item()
	print(data['angularNumber'])
	phaseFunc = data['phaseFunc']
	print(phaseFunc.shape)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
