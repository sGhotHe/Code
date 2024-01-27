####################################################################################
# INTRODUCTION:
# This code is to do some PNSD calculating
# Created by Hebs at 21/10/25/13:14
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import PyMieScatt as ps
import sys
import kappaKohler

def PNSD_Dps_adjust(PNSD, Dps, mag):
	'''
	This function is to adjust PNSD Dps magnification (or resolution)
	input:
		PNSD	: particle number size distribution, array, dn/dlogDp
		Dps	: diameter size distribution, array, nm
		mag	: magnification adjust for Dps, float
	output:
		PNSD_new	: adjusted particle number size distribution, array, dn/dlogDp
		Dps_new	: adjusted diameter size distribution, array, nm
	'''
	logDps = np.log(Dps)
	Dps_num = len(Dps) * mag
	logDps_new = np.arange(logDps[0], logDps[-1], (logDps[-1]-logDps[0])/Dps_num)
	Dps_new = np.exp(logDps_new)
	PNSD_new = np.zeros(len(Dps_new))
	for i in range(len(PNSD_new)):
		PNSD_new[i] = PNSD_Dp(Dps, PNSD, Dps_new[i]) / mag
	return PNSD_new, Dps_new

def cal_dlnDp(Dps):
	'''
	This function is to calculate dlnDp in PNSD
	input:
		Dps    : diameter distribution, array, nm
	output:
		dlnDp  : dlnDp in PNSD
	'''
	dlnDp = (np.log10(Dps[-1])-np.log10(Dps[0]))/(len(Dps)-1)
	return dlnDp

def cal_dlogDp(Dps):
	'''
	This function is to calculate dlogDp in PNSD
	input:
		Dps    : diameter distribution, array, nm
	output:
		dlogDp  : dlnDp in PNSD
	'''
	dlogDp = (np.log(Dps[-1])-np.log(Dps[0]))/(len(Dps)-1)
	return dlogDp

def PNSD2n(Dps, PNSD, index):
	'''
	This function is to turn PNSD (dN/dlnDp) to number concentration (N/V)
	
	|---Dp0---|---Dp1---|---Dp2---|---Dp3---|---Dp4---|........|---Dpx---|
	b0        b1        b2        b3        b4        b5       bx        bx+1
	
	input:
		Dps        : diameter distribution, array, nm
		PNSD       : partical concentration, array, dn/dlogDp
		index      : dlogDp or dlnDp, string, 'e' or 'd'
	output:
		n          : number concentration, array, cm-3
	'''
	if index=='e':
		n = PNSD * cal_dlogDp(Dps)
	elif index=='d':
		n = PNSD * cal_dlnDp(Dps)
	else:
		print('Invalid input for index :', index, '. Please check.')
		sys.eixt()
	return n

def PNSD2m(Dps, PNSD, rho, index):
	'''
	This function is to turn PNSD (dN/dlnDp) to mass concentration (g/cm3)
	input:
		Dps        : diameter distribution, array, nm
		PNSD       : partical concentration, array, dn/dlogDp
		rho        : partical density, float, g/cm3
		index      : dlogDp or dlnDp, string, 'e' or 'd'
	output:
		m          : mass concentration, array, ug/m3
	'''
	n = PNSD2n(Dps, PNSD, index)
	m = np.zeros(len(Dps))
	for i in range(len(Dps)):
		m[i] = np.pi / 6 * rho * (Dps[i]*1e-7)**3 * n[i] * 1e12
	return m

def n2PNSD(Dps, n, index):
	'''
	This function is to turn number concentration (N/V) to PNSD (dN/dlnDp)
	
	|---Dp0---|---Dp1---|---Dp2---|---Dp3---|---Dp4---|........|---Dpx---|
	b0        b1        b2        b3        b4        b5       bx        bx+1
	
	input:
		Dps        : diameter distribution, array, nm
		n          : number concentration, array, cm-3
		index      : dlogDp or dlnDp, string, 'e' or 'd'
	output:
		PNSD       : partical concentration, array, dn/dlnDp
	'''
	if index=='e':
		PNSD = n / cal_dlogDp(Dps)
	elif index=='d':
		PNSD = n / cal_dlnDp(Dps)
	else:
		print('Invalid input for index :', index, '. Please check.')
		sys.eixt()
	return PNSD

def PNSD_Dp(Dps, PNSD, Dp):
	'''
	This function is to use interpolate to calculate PNSD in certain Dp
	input:
		Dps        : diameter distribution, array, nm
		PNSD       : partical concentration, array, dn/dlnDp
		Dp         : partical diameter, float, nm
	output:
		PNSD_Dp    : partical concentration in certain Dp, float, dn/dlnDp
	'''
	Dps_min = np.exp(np.log(Dps[0])-cal_dlogDp(Dps))
	Dps_max = np.exp(np.log(Dps[-1])+cal_dlogDp(Dps))
	if Dp<=Dps_min or Dp>=Dps_max:
		return 0
	Dps_new = np.append(np.append(Dps_min, Dps), Dps_max)
	PNSD_new = np.append(np.append(0, PNSD), 0)
	i = 0
	while i<len(Dps_new)-1:
		if Dp>Dps_new[i] and Dp<=Dps_new[i+1]:
			break
		i += 1
	PNSD_Dp = PNSD_new[i] + (PNSD_new[i+1]-PNSD_new[i]) * np.log(Dp/Dps_new[i]) / np.log(Dps_new[i+1]/Dps_new[i])
	return PNSD_Dp

def DMASP22DMASP2PNSD(DBCps, DMASP2):
	'''
	This function is to turn DMA-SP2 data from number concentration (N/V) to PNSD(dN/dlogDp)
	input:
		DBCps    : diameter distribution, array, nm
		DMASP2   : DMA-SP2 data, array in shape(len(DBCps), len(DBC)), dn/dlnDBC
	output:
		DMASP2PNSD : DMA-SP2 data partical number size distribution, array in shape(len(DBCps), len(DBC)), dn/dlnDBCps/dlnDBC
	'''
	DMASP2PNSD = np.zeros(DMASP2.shape)
	for i in range(len(DBCps)-1):
		DMASP2PNSD[i] = DMASP2[i] / np.log(DBCps[i+1]/DBCps[i])
	DMASP2PNSD[-1] = DMASP2[-1] / np.log(DBCps[-1]/DBCps[-2])
	return DMASP2PNSD

def DMASP2PNSD_Dp(DBCps, DMASP2PNSD, Dp):
	'''
	This function is to use interpolate to calculate DMASP2PNSD in certain Dp
	input:
		DBCps      : diameter distribution, array, nm
		DMASP2PNSD : DMA-SP2 data partical number size distribution, array in shape(len(DBCps), len(DBC)), dn/dlogDBCps/dlogDBC
		Dp         : partical diameter, float, nm
	output:
		PNSD_Dp    : partical concentration in certain Dp, float, dn/dlogDp
	'''
	DBCps_min = np.exp(np.log(DBCps[0])-np.log(DBCps[1]/DBCps[0]))
	DBCps_max = np.exp(np.log(DBCps[-1])+np.log(DBCps[-1]/DBCps[-2]))
	if Dp<=DBCps_min or Dp>=DBCps_max:
		return np.zeros(len(DMASP2PNSD[0]))
	DBCps_new = np.append(np.append(DBCps_min, DBCps), DBCps_max)
	DMASP2PNSD_new = np.vstack((np.vstack((np.zeros(len(DMASP2PNSD[0])), DMASP2PNSD)), np.zeros(len(DMASP2PNSD[0]))))
	i = 0
	while i<len(DBCps_new)-1:
		if Dp>DBCps_new[i] and Dp<=DBCps_new[i+1]:
			break
		i += 1
	DMASP2PNSD_Dp = DMASP2PNSD_new[i] + (DMASP2PNSD_new[i+1]-DMASP2PNSD_new[i]) * np.log(Dp/DBCps_new[i]) / np.log(DBCps_new[i+1]/DBCps_new[i])
	return DMASP2PNSD_Dp

def cal_vertical_distribution(Dps, PNSD, H, **args):
	'''
	This function is to calculate aerosol vertical distribution
	input:
		Dps    : diameter distribution, array, nm
		PNSD   : number distribution, array, dn/dlogDp
		H      : height, m
		**func : vertical distribution function type, A or B, default B
			 type A: Na = N0, if H<H_PBL
			         Na = N0 - k*(H-H_PBL), if H_PBL <= H < H_PBL+H_LD
			         Na = (N0-k*H_LD)*exp(-(H-H_LD-H_PBL)/H_p), if H >= H_PBL + H_LD
			         H_p=1800m, H_PBL=1000m, H_LD=200m, k=28cm^-3m^-1
			 type B: Na = N0exp(-H/H_p)
			         H_p=758m
	output:
		PNSD : number distribution in z, array, dn/dlogDp
	'''
	if 'func' in args:
		func = args['func']
	else:
		func = 'B'
	
	n = sum(PNSD2n(Dps, PNSD, index='d'))
	
	if func=='B':
		H_p = 758 # m
		PNSD_H = PNSD * np.exp(-H/H_p)
	else:
		H_p = 1800 # m
		H_PBL = 1000 # m
		H_LD = 200 # m
		k = 28 # cm^-3m^-1, for N0=8000
		k = k * PNSD2n(Dps,PNSD, index='d') / 8000 # adapt coefficient k for PNSD
		
		if H<H_PBL:
			n_H = n
		elif H<(H_PBL+H_LD):
			n_H = n - k * (H-H_PBL)
		else:
			n_H = (n-k*H_LD) * np.exp(-(H-H_LD-H_PBL)/H_p)
		PNSD_H = PNSD * n_H / n
	return PNSD_H

def Dps2SSA(Dps, PNSD, wl, m):
	'''
	This function is to turn Dps to single scattering albedo SSA
	input:
		Dps       : diameter distribution, 1-d array, nm
		PNSD      : number distribution, 1-d array, dN/dlogDp, cm-3
		wl        : wave length, nm
		m         : complex refractive index
	output:
		SSA       : single scattering albedo
	'''
	bext = 0
	bsca = 0
	
	for i in range(len(Dps)):
		Qext = ps.MieQ(m, wl, Dps[i])[0]
		Qsca = ps.MieQ(m, wl, Dps[i])[1]
		bext += 1/4 * np.pi * Dps[i]**2 * Qext * PNSD[i]
		bsca += 1/4 * np.pi * Dps[i]**2 * Qsca * PNSD[i]
	
	SSA = bsca / bext
	return SSA

def cal_SSA(Dps, PNSD, wl, H, m, **args):
	'''
	This function is to calculate single scattering albedo in certern height
	input:
		Dps    : diameter distribution, array, nm
		PNSD   : number distribution, array, dn/dlogDp
		wl     : wave length, nm
		H      : height, m
		m      : complex rafractive index
		**func : vertical distribution function type, A or B, default B
	output:
		SSA    : single scattering albedo
	'''
	if 'func' in args:
		func = args['func']
	else:
		func = 'B'
	
	PNSD_H = cal_vertical_distribution(Dps, PNSD, H, func=func)
	SSA = Dps2SSA(Dps, PNSD_H, wl, m)
	return SSA

def Dps2k(Dps, PNSD, wl, m):
	'''
	This function is to turn Dps to scattering coefficient k
	input:
		Dps       : diameter distribution, 1-d array, nm
		PNSD      : number distribution, 1-d array, dN/dlogDp, cm-3
		wl        : wave length, nm
		m         : complex refractive index, default value: 1+0j
	output:
		k         : scattering coefficient, Mm^-1
	'''
	n = PNSD2n(Dps, PNSD, 'd')
	k = 0
	
	for i in range(len(Dps)):
		Qsca = ps.MieQ(m, wl, Dps[i])[1]
		Dps_i = Dps[i] * 1e-9 # turn nm to m
		sigma_sca = Qsca * 1/4 * np.pi * Dps_i**2
		k += sigma_sca * n[i] * 1e6 # turn n[i] from cm^-3 to m^-3
	k = k * 1e6 # turn m^-1 to Mm^-1
	return k

def Dps2kext(Dps, PNSD, wl, m):
	'''
	This function is to turn Dps to extinct coefficient k
	input:
		Dps       : diameter distribution, 1-d array, nm
		PNSD      : number distribution, 1-d array, dN/dlogDp, cm-3
		wl        : wave length, nm
		m         : complex refractive index, default value: 1+0j
	output:
		kext      : extince coefficient, Mm^-1
	'''
	n = PNSD2n(Dps, PNSD, 'd')
	kext = 0
	
	for i in range(len(Dps)):
		Qsca = ps.MieQ(m, wl, Dps[i])[0]
		Dps_i = Dps[i] * 1e-9 # turn nm to m
		sigma_sca = Qsca * 1/4 * np.pi * Dps_i**2
		kext += sigma_sca * n[i] * 1e6 # turn n[i] from cm^-3 to m^-3
	kext = kext * 1e6 # turn m^-1 to Mm^-1
	return kext

def Dps2g(Dps, PNSD, wl, m):
	'''
	This function is to turn Dps to asymmetry factor g
	input:
		Dps       : diameter distribution, 1-d array, nm
		PNSD      : number distribution, 1-d array, dN/dlogDp, cm-3
		wl        : wave length, nm
		m         : complex refractive index
	output:
		g         : asymmetry factor
	'''
	n = PNSD2n(Dps, PNSD, 'd')
	g = 0
	
	for i in range(len(Dps)):
		g_i = ps.MieQ(m, wl, Dps[i])[3]
		g += g_i * n[i] / np.nansum(n)
	return g

def cal_ksca_old(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, kappa, RH, wl, CT):
	'''
	This function is to use Dps, PNSD, RH and CT to calculate scattering coefficient
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
	output:
		ksca    : scattering coefficient, Mm-1
	'''
	DBCps_min = DBCps[0]
	DBCps_max = DBCps[-1]
	Dps_min = Dps[0]
	Dps_max = Dps[-1]
	n = PNSD2n(Dps, PNSD, 'e')
	m_BC = nBC + kBC * 1j
	m_shell = nShell
	n_BC = BCPNSD
	
	PNSD_fit = np.zeros(len(DBCps))
	n_fit = np.zeros(len(DBCps))
	n_noBC = np.zeros(len(DBCps))
	ksca = 0
	
	for i in range(len(Dps)):
		if Dps[i]<DBCps_min or Dps[i]>DBCps_max:
			# use PNSD to calculate ksca
			D_i, m_shell_RH_i = kappaKohler.RH2D(Dps[i], 0, m_BC, m_shell, kappa, RH, wl)
			Qsca_noBCi = ps.MieQ(m_shell_RH_i, wl, D_i)[1]
			ksca += Qsca_noBCi * 1/4 * np.pi * (D_i*1e-9)**2 * n[i]*1e6
	
	# use n_BC to calculate ksca
	for i in range(len(DBCps)):
		# first need to fit PNSD to DBCps
		j = 0 # to mark the Dps subscript
		while j<len(Dps)-1:
			if Dps[j]<=DBCps[i] and Dps[j+1]>DBCps[i]:
				break
			j += 1
		if j<len(Dps)-1:
			# DBCps[i] is between Dps[j] and Dps[j+1]
			PNSD_fit[i] = PNSD[j] + (PNSD[j+1]-PNSD[j]) * (np.log(DBCps[i])-np.log(Dps[j])) / (np.log(Dps[j+1])-np.log(Dps[j])) # log linear interpolation
			n_fit[i] = PNSD_fit[i] * cal_dlogDp(DBCps)
			# DBC must smaller than DBCps, otherwise n_BC[i,j] must be zero
			for j in range(len(DBC)):
				if DBC[j]>DBCps[i]:
					n_BC[i,j] = 0
			# if n_BC[i] is bigger than n_fit, n_BC[i] must be adjusted
			# n_BC[i] total number must not bigger than n_fit
			if sum(n_BC[i])>n_fit[i]:
				# and total n_BC[i] must not bigger than n_fit
				n_BC[i] = n_BC[i] * n_fit[i] / sum(n_BC[i])
				# n_noBC[i] = 0
			else:
				n_noBC[i] = n_fit[i] - sum(n_BC[i])
			# calculate no BC core partical ksca
			D_i, m_shell_RH_i = kappaKohler.RH2D(DBCps[i], 0, m_BC, m_shell, kappa, RH, wl)
			Qsca_noBCi = ps.MieQ(m_shell_RH_i, wl, D_i)[1]
			ksca += Qsca_noBCi * 1/4 * np.pi * (D_i*1e-9)**2 * n_noBC[i]*1e6
		# else:
		# DBCps[i] is either bigger than Dps_max or smaller than Dps_min
		# in this situation, n_fit[i] and n_noBC[i] is zero
		# just jump over
		
		for j in range(len(DBC)):
			if DBC[j]<=DBCps[i]:
				D_ij, m_shell_RH_ij = kappaKohler.RH2D(DBCps[i], DBC[j], m_BC, m_shell, kappa, RH, wl)
				Qsca_BCij = ps.MieQCoreShell(m_BC, m_shell_RH_ij, wl, DBC[j], D_ij*CT)[1]
				ksca += Qsca_BCij * 1/4 * np.pi * (D_ij*1e-9)**2 * n_BC[i,j]*1e6
	
	ksca = ksca * 1e6 # Mm^-1
	return ksca

def cal_kext_old(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, kappa, RH, wl, CT, BCAE):
	'''
	This function is to use Dps, PNSD, RH, CT and BCAE to calculate extinct coefficient
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
		kext    : extinct coefficient, Mm-1
	'''
	DBCps_min = DBCps[0]
	DBCps_max = DBCps[-1]
	Dps_min = Dps[0]
	Dps_max = Dps[-1]
	n = PNSD2n(Dps, PNSD, 'e')
	m_BC = nBC + kBC * 1j
	m_shell = nShell
	n_BC = BCPNSD
	
	PNSD_fit = np.zeros(len(DBCps))
	n_fit = np.zeros(len(DBCps))
	n_noBC = np.zeros(len(DBCps))
	kext = 0
	
	for i in range(len(Dps)):
		if Dps[i]<DBCps_min or Dps[i]>DBCps_max:
			D_i, m_shell_RH_i = kappaKohler.RH2D(Dps[i], 0, m_BC, m_shell, kappa, RH, wl)
			Qext_noBCi = ps.MieQ(m_shell_RH_i, wl, D_i)[0]
			kext += Qext_noBCi * 1/4 * np.pi * (D_i*1e-9)**2 * n[i]*1e6
	
	for i in range(len(DBCps)):
		j = 0
		while j<len(Dps)-1:
			if Dps[j]<=DBCps[i] and Dps[j+1]>DBCps[i]:
				break
			j += 1
		if j<len(Dps)-1:
			PNSD_fit[i] = PNSD[j] + (PNSD[j+1]-PNSD[j]) * (np.log(DBCps[i])-np.log(Dps[j])) / (np.log(Dps[j+1])-np.log(Dps[j]))
			n_fit[i] = PNSD_fit[i] * cal_dlogDp(DBCps)
			for j in range(len(DBC)):
				if DBC[j]>DBCps[i]:
					n_BC[i,j] = 0
			if sum(n_BC[i])>n_fit[i]:
				n_BC[i] = n_BC[i] * n_fit[i] / sum(n_BC[i])
				# n_noBC[i] = 0
			else:
				n_noBC[i] = n_fit[i] - sum(n_BC[i])
			D_i, m_shell_RH_i = kappaKohler.RH2D(DBCps[i], 0, m_BC, m_shell, kappa, RH, wl)
			Qext_noBCi = ps.MieQ(m_shell_RH_i, wl, D_i)[0]
			kext += Qext_noBCi * 1/4 * np.pi * (D_i*1e-9)**2 * n_noBC[i]*1e6
		
		for j in range(len(DBC)):
			if DBC[j]<=DBCps[i]:
				D_ij, m_shell_RH_ij = kappaKohler.RH2D(DBCps[i], DBC[j], m_BC, m_shell, kappa, RH, wl)
				Qext_BCij = ps.MieQCoreShell(m_BC, m_shell_RH_ij, wl, DBC[j], D_ij*CT)[1] + ps.MieQCoreShell(m_BC, m_shell_RH_ij, wl, DBC[j], D_ij*CT)[2] * BCAE
				kext += Qext_BCij * 1/4 * np.pi * (D_ij*1e-9)**2 * n_BC[i,j]*1e6
	
	kext = kext * 1e6 # Mm^-1
	return kext

def cal_g_old(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, kappa, RH, wl, CT):
	'''
	This function is to use Dps, PNSD, RH and CT to calculate asymmetry factor g
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
	output:
		g       : asymmetry factor
	'''
	DBCps_min = DBCps[0]
	DBCps_max = DBCps[-1]
	Dps_min = Dps[0]
	Dps_max = Dps[-1]
	n = PNSD2n(Dps, PNSD, 'e')
	m_BC = nBC + kBC * 1j
	m_shell = nShell
	n_BC = BCPNSD
	
	PNSD_fit = np.zeros(len(DBCps))
	n_fit = np.zeros(len(DBCps))
	n_noBC = np.zeros(len(DBCps))
	g = 0
	
	for i in range(len(Dps)):
		if Dps[i]<DBCps_min or Dps[i]>DBCps_max:
			D_i, m_shell_RH_i = kappaKohler.RH2D(Dps[i], 0, m_BC, m_shell, kappa, RH, wl)
			g_noBCi = ps.MieQ(m_shell_RH_i, wl, D_i)[3]
			g += g_noBCi * n[i] / sum(n)
	
	for i in range(len(DBCps)):
		j = 0
		while j<len(Dps)-1:
			if Dps[j]<=DBCps[i] and Dps[j+1]>DBCps[i]:
				break
			j += 1
		if j<len(Dps)-1:
			PNSD_fit[i] = PNSD[j] + (PNSD[j+1]-PNSD[j]) * (np.log(DBCps[i])-np.log(Dps[j])) / (np.log(Dps[j+1])-np.log(Dps[j]))
			n_fit[i] = PNSD_fit[i] * cal_dlogDp(DBCps)
			for j in range(len(DBC)):
				if DBC[j]>DBCps[i]:
					n_BC[i,j] = 0
			if sum(n_BC[i])>n_fit[i]:
				n_BC[i] = n_BC[i] * n_fit[i] / sum(n_BC[i])
				# n_noBC[i] = 0
			else:
				n_noBC[i] = n_fit[i] - sum(n_BC[i])
			D_i, m_shell_RH_i = kappaKohler.RH2D(DBCps[i], 0, m_BC, m_shell, kappa, RH, wl)
			g_noBCi = ps.MieQ(m_shell_RH_i, wl, D_i)[3]
			g += g_noBCi * n_noBC[i] / np.nansum(n)
		
		for j in range(len(DBC)):
			if DBC[j]<=DBCps[i]:
				D_ij, m_shell_RH_ij = kappaKohler.RH2D(DBCps[i], DBC[j], m_BC, m_shell, kappa, RH, wl)
				g_BCij = ps.MieQCoreShell(m_BC, m_shell_RH_ij, wl, DBC[j], D_ij*CT)[1] + ps.MieQCoreShell(m_BC, m_shell_RH_ij, wl, DBC[j], D_ij*CT)[3]
				g += g_BCij * n_BC[i,j] / np.nansum(n)
	
	return g

def cal_ksca(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, kappa, RH, wl, CT):
	'''
	This function is to use Dps, PNSD, RH, CT and BCAE to calculate scattering coefficient
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
	output:
		ksca    : scattering coefficient, Mm-1
	'''
	m_BC = nBC + kBC * 1j
	m_shell = nShell
	n = PNSD2n(Dps, PNSD, 'd')
	BCPNSD = BCPNSD * cal_dlogDp(DBC)
	ksca = 0
	
	for i in range(len(Dps)):
		BCPNSD_i = DMASP2PNSD_Dp(DBCps, BCPNSD, Dps[i]) # array in shape of DBC
		n_noBC_i = n[i] - sum(BCPNSD_i)
		
		if n_noBC_i > 0: # have no BC particles
			D_noBC_i, m_shell_noBC_i = kappaKohler.RH2D(Dps[i], 0, m_BC, m_shell, kappa, RH, wl)
			Qsca_noBC_i = ps.MieQ(m_shell_noBC_i, wl, D_noBC_i)[1]
			ksca += Qsca_noBC_i * 1/4 * np.pi * (D_noBC_i*1e-9)**2 * n_noBC_i*1e6 # in m-1
		
		if sum(BCPNSD_i): # have BC particles
			for j in range(len(DBC)):
				D_BC_ij, m_shell_BC_ij = kappaKohler.RH2D(Dps[i], DBC[j], m_BC, m_shell, kappa, RH, wl)
				Qsca_BC_ij = ps.MieQCoreShell(m_BC, m_shell_BC_ij, wl, DBC[j], D_BC_ij*CT)[1]
				ksca += Qsca_BC_ij * 1/4 * np.pi * (D_BC_ij*1e-9)**2 * BCPNSD_i[j]*1e6
	
	ksca = ksca * 1e6 # in Mm-1
	return ksca

def cal_kext(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, kappa, RH, wl, CT, BCAE):
	'''
	This function is to use Dps, PNSD, RH, CT and BCAE to calculate extinct coefficient
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
		kext    : extinct coefficient, Mm-1
	'''
	m_BC = nBC + kBC * 1j
	m_shell = nShell
	n = PNSD2n(Dps, PNSD, 'd')
	BCPNSD = BCPNSD * cal_dlogDp(DBC)
	kext = 0
	
	for i in range(len(Dps)):
		BCPNSD_i = DMASP2PNSD_Dp(DBCps, BCPNSD, Dps[i]) # array in shape of DBC
		n_noBC_i = n[i] - sum(BCPNSD_i)
		
		if n_noBC_i > 0:
			D_noBC_i, m_shell_noBC_i = kappaKohler.RH2D(Dps[i], 0, m_BC, m_shell, kappa, RH, wl)
			Qext_noBC_i = ps.MieQ(m_shell_noBC_i, wl, D_noBC_i)[0]
			kext += Qext_noBC_i * 1/4 * np.pi * (D_noBC_i*1e-9)**2 * n_noBC_i*1e6 # in m-1
		
		if sum(BCPNSD_i):
			for j in range(len(DBC)):
				D_BC_ij, m_shell_BC_ij = kappaKohler.RH2D(Dps[i], DBC[j], m_BC, m_shell, kappa, RH, wl)
				Qext_BC_ij = ps.MieQCoreShell(m_BC, m_shell_BC_ij, wl, DBC[j], D_BC_ij*CT)[1] + ps.MieQCoreShell(m_BC, m_shell_BC_ij, wl, DBC[j], D_BC_ij*CT)[2] * BCAE
				kext += Qext_BC_ij * 1/4 * np.pi * (D_BC_ij*1e-9)**2 * BCPNSD_i[j]*1e6
	
	kext = kext * 1e6
	return kext

def cal_Mie(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, nMS, kappa, kappaMS, RH, wl, CT, BCAE):
	'''
	This function is to use Dps, PNSD, RH, CT and BCAE to calculate extinct coefficient, scattering coefficient and asymmetry factor g
	input:
		Dps     : diameter distribution, array, nm
		PNSD    : number distribution, array, dn/dlogDp
		DBCps   : BC particle diameter size, array, nm
		DBC     : BC core diameter size, array, nm
		BCPNSD  : BC particle number concentration, array, cm^-3
		kBC     : BC core imagine part of complex refractive index
		nBC     : BC core real part of complex refractive index
		nShell  : shell complex refractive index
		nMS     : n mixing state change rate by diameter
		kappa   : hygroscopicity parameter, float
		kappaMS : hygroscopicity parameter mixing state change rate by diameter
		RH      : relative humidity, float, percent
		wl      : wave length, nm
		CT      : BC core coating thickness adjust parameter
		BCAE    : BC core absorbing enhancement adjust parameter
	output:
		kext    : extinct coefficient, Mm-1
		ksca    : scattering coefficient, Mm-1
		g       : asymmetry factor
	'''
	m_BC = nBC + kBC * 1j
	#m_shell = nShell
	m_shell = np.zeros(len(Dps))
	kappa_shell = np.zeros(len(Dps))
	for i in range(len(Dps)):
		m_shell[i] = nShell + nMS * (i-round(len(Dps)/2)) / round(len(Dps)/2)
		kappa_shell[i] = kappa + kappaMS * (i-round(len(Dps)/2)) / round(len(Dps)/2)
	n = PNSD2n(Dps, PNSD, 'd')
	BCPNSD = BCPNSD * cal_dlogDp(DBC)
	kext = 0
	ksca = 0
	g = 0
	
	for i in range(len(Dps)):
		BCPNSD_i = DMASP2PNSD_Dp(DBCps, BCPNSD, Dps[i]) # array in shape of DBC
		n_noBC_i = n[i] - sum(BCPNSD_i)
		
		if n_noBC_i > 0:
			D_noBC_i, m_shell_noBC_i = kappaKohler.RH2D(Dps[i], 0, m_BC, m_shell[i], kappa_shell[i], RH, wl)
			Qsca_noBC_i = ps.MieQ(m_shell_noBC_i, wl, D_noBC_i)[1]
			ksca += Qsca_noBC_i * 1/4 * np.pi * (D_noBC_i*1e-9)**2 * n_noBC_i*1e6
			Qext_noBC_i = ps.MieQ(m_shell_noBC_i, wl, D_noBC_i)[0]
			kext += Qext_noBC_i * 1/4 * np.pi * (D_noBC_i*1e-9)**2 * n_noBC_i*1e6
			g_noBC_i = ps.MieQ(m_shell_noBC_i, wl, D_noBC_i)[3]
			g += g_noBC_i * n_noBC_i / sum(n)
		
		if sum(BCPNSD_i):
			for j in range(len(DBC)):
				D_BC_ij, m_shell_BC_ij = kappaKohler.RH2D(Dps[i], DBC[j], m_BC, m_shell[i], kappa_shell[i], RH, wl)
				Qsca_BC_ij = ps.MieQCoreShell(m_BC, m_shell_BC_ij, wl, DBC[j], D_BC_ij*CT)[1]
				ksca += Qsca_BC_ij * 1/4 * np.pi * (D_BC_ij*1e-9)**2 * BCPNSD_i[j]*1e6
				Qext_BC_ij = Qsca_BC_ij + ps.MieQCoreShell(m_BC, m_shell_BC_ij, wl, DBC[j], D_BC_ij*CT)[2] * BCAE
				kext += Qext_BC_ij * 1/4 * np.pi * (D_BC_ij*1e-9)**2 * BCPNSD_i[j]*1e6
				g_BC_ij = ps.MieQCoreShell(m_BC, m_shell_BC_ij, wl, DBC[j], D_BC_ij*CT)[3]
				g += g_BC_ij * BCPNSD_i[j] / sum(n)
	
	kext = kext * 1e6
	ksca = ksca * 1e6
	return kext, ksca, g

def cal_g(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, kappa, RH, wl, CT):
	'''
	This function is to use Dps, PNSD, RH, CT and BCAE to calculate asymmetry factor g
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
	output:
		kext    : extinct coefficient, Mm-1
	'''
	m_BC = nBC + kBC * 1j
	m_shell = nShell
	n = PNSD2n(Dps, PNSD, 'd')
	BCPNSD = BCPNSD * cal_dlogDp(DBC)
	g = 0
	
	for i in range(len(Dps)):
		BCPNSD_i = DMASP2PNSD_Dp(DBCps, BCPNSD, Dps[i]) # array in shape of DBC
		n_noBC_i = n[i] - sum(BCPNSD_i)
		
		if n_noBC_i > 0:
			D_noBC_i, m_shell_noBC_i = kappaKohler.RH2D(Dps[i], 0, m_BC, m_shell, kappa, RH, wl)
			g_noBC_i = ps.MieQ(m_shell_noBC_i, wl, D_noBC_i)[3]
			g += g_noBC_i * n_noBC_i / sum(n)
		
		if sum(BCPNSD_i):
			for j in range(len(DBC)):
				D_BC_ij, m_shell_BC_ij = kappaKohler.RH2D(Dps[i], DBC[j], m_BC, m_shell, kappa, RH, wl)
				g_BC_ij = ps.MieQCoreShell(m_BC, m_shell_BC_ij, wl, DBC[j], D_BC_ij*CT)[3]
				g += g_BC_ij * BCPNSD_i[j] / sum(n)
	
	return g

def cal_dtau(Dps, PNSD, wl, H, dH, m, **args):
	'''
	This function is to calculate delta aerosol optical depth in certern height and delta height
	input:
		Dps    : diameter distribution, array, nm
		PNSD   : number distribution, array, dn/dlogDp
		wl     : wave length, nm
		H      : height, m
		dH     : delta height, m
		m      : complex rafractive index
		**func : vertical distribution function type, A or B, default B
		**ddH  : step, default dH/100
	output:
		dtau   : delta aerosol optical depth
	'''
	if 'func' in args:
		func = args['func']
	else:
		func = 'B'
	if 'ddH' in args:
		ddH = args['ddH']
	else:
		ddH = dH / 100 # m
	
	n = round(dH/ddH)
	dtau = 0
	for i in range(n):
		PNSD_H = cal_vertical_distribution(Dps, PNSD, H+n*ddH, func=func)
		k_H = Dps2k(Dps, PNSD_H, wl, m) # Mm^-1
		dtau += k_H * ddH * 1e-6
	return dtau

def cal_VD(Dps, PNSD, z, VD):
	'''
	This function is to calculate vertical distribution
	type A: Na = N0, if H<H_PBL
		Na = N0 - k*(H-H_PBL), if H_PBL <= H < H_PBL+H_LD
		Na = (N0-k*H_LD)*exp(-(H-H_LD-H_PBL)/H_p), if H >= H_PBL + H_LD
		H_p=1800m, H_PBL=1000m, H_LD=200m, k=28cm^-3m^-1
	type B: Na = N0exp(-H/H_p)
		H_p=758m
	input:
		Dps    : diameter of particals, array, nm
		PNSD   : number distribution, array, cm^-3
		z      : height, m
		VD     : vertical distribution mixing, 0 for type B and 1 for type A
	output:
		ratio  : ratio of number distribution at z to the distribution on the ground
	'''
	n = PNSD2n(Dps, PNSD, 'd')
	H_p = 758 # m
	ratio_B = np.exp(-z/H_p)
	
	H_p = 1800 # m
	H_PBL = 1000 # m
	H_LD = 200 # m
	k = 28 # cm^-3m^-1, for N0=8000
	n_total = np.sum(n)
	k = k * n_total / 8000 # adapt coefficient k for n_total
	
	if z<H_PBL:
		n_total_H = n_total
	elif z<(H_PBL+H_LD):
		n_total_H = n_total - k * (z-H_PBL)
	else:
		n_total_H = (n_total-k*H_LD) * np.exp(-(z-H_LD-H_PBL)/H_p)
	ratio_A = n_total_H / n_total
	
	ratio = ratio_A * VD + ratio_B * (1-VD)
	
	return ratio

if __name__ == '__main__':
	import readTaizhou
	sp2 = readTaizhou.read_Taizhou('data/sp2/Taizhou.npy')
	Dps = sp2['Dps']
	PNSD = sp2['PNSD']
	PNSD = np.nanmean(PNSD, axis=0)
	PNSD_new, Dps_new = PNSD_Dps_adjust(PNSD, Dps, 0.5)
	import matplotlib.pyplot as plt
	plt.scatter(np.log(Dps), PNSD)
	plt.scatter(np.log(Dps_new), PNSD_new)
	plt.show()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
