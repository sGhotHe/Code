####################################################################################
# INTRODUCTION:
# This code is to do some BC related calculation
# Created by Hebs at 21/11/3/17:26
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import PyMieScatt as ps
import calPNSD
import readTaizhou
import phaseFunc
import kappaKohler

'''
some constant and calculation:
m_EC = 2.26+1.26j
m_air = 1
rho_EC = 1.8 # g/cm^3
rho_air = 1.293e-3 # g/cm^3
m_BC = m_EC * (1-V_air) + m_air * V_air
rho_BC = rho_BC * (1-V_air) + rho_air * V_air
'''
def rhoBC2mBC(rhoBC, **args):
	'''
	This function is to change rhoBC to mBC
	input:
		rhoBC         : BC density
		args:
			rhoEC : EC density, default 1.8 g/cm3
			nEC   : EC real part of complex refractive index, default 2.26
			mEC   : EC imagine part of complex refractive index, default 1.26
			nAir  : air real part of complex refractive index, default 1
	output:
		nBC           : BC real part of complex refractive index
		mBC           : BC imagine part of complex refractive index
	'''
	if 'rhoEC' in args:
		rhoEC = args['rhoEC']
	else:
		rhoEC = 1.8
	if 'nEC' in args:
		nEC = args['nEC']
	else:
		nEC = 2.26
	if 'kEC' in args:
		kEC = args['kEC']
	else:
		kEC = 1.26
	if 'nAir' in args:
		nAir = args['nAir']
	else:
		nAir = 1
	x = rhoBC / rhoEC
	kBC = kEC * x
	nBC = (nEC-nAir) * x + nAir
	return nBC, kBC

def DMASP22vBC2(Dps, DBCps, DBC, DMASP2):
	v_BC = np.zeros(len(DBC))
	DMASP2PNSD = calPNSD.DMASP22DMASP2PNSD(DBCps, DMASP2) # devide
	
	for i in range(len(Dps)):
		DMASP2PNSD_i = calPNSD.DMASP2PNSD_Dp(DBCps, DMASP2PNSD, Dps[i])
		DMASP2_i = calPNSD.PNSD2n(Dps, DMASP2PNSD_i, 'e') # multiply
		n_BC_i = calPNSD.PNSD2n(DBC, DMASP2_i, 'e')
		for j in range(len(DBC)):
			v_BC[j] += np.pi/6 * DBC[j]**3 * n_BC_i[j] # in nm^3/cm^3
	
	return v_BC*1e-9 # in um^3/cm^3

def DMASP22vBC(DBCps, DBC, DMASP2):
	'''
	This function is to use size-resolve DMA-sp2 data to calculate BCPNSD
	input:
		DBCps   : BC particle diameter size, array, nm
		DBC     : BC core diameter size, array, nm
		DMASP2  : BC size-resolved particle number size distribution, array, dn/dlogDBC
	output:
		v_BC    : BC particle volume concentration, array, um^3/cm^3
	'''
	v_BC = np.zeros(len(DBC))
	
	for i in range(len(DBCps)):
		n_BC_i = calPNSD.PNSD2n(DBC, DMASP2[i], 'e')
		for j in range(len(DBC)):
			v_BC[j] += np.pi/6 * DBC[j]**3 * n_BC_i[j] # in nm^3/cm^3
	
	return v_BC*1e-9 # in um^3/cm^3

def BCPNSD2vBC(DBC, BCPNSD):
	'''
	This function is to use BCPNSD data to calculate BCPVSD
	input:
		DBC     : BC core diameter size, array, nm
		BCPNSD  : BC particle number size distribution, array, dn/dlogDBC
	output:
		v_BC    : BC particle volume concentration, array, um^3/cm^3
	'''
	v_BC = np.zeros(len(DBC))
	n_BC = calPNSD.PNSD2n(DBC, BCPNSD, 'e')
	
	for i in range(len(DBC)):
		v_BC[i] = np.pi/6 * DBC[i]**3 * n_BC[i] # in nm^3/cm^3
	
	return v_BC*1e-9 # in um^3/cm^3

def BCPNSD2vBC2(Dps, DBC, BCPNSD):
	v_BC = 0
	
	for i in range(len(Dps)):
		BCPNSD_i = calPNSD.PNSD_Dp(DBC, BCPNSD, Dps[i])
		n_BC_i = calPNSD.PNSD2n(Dps, BCPNSD_i, 'e')
		v_BC += np.pi/6 * Dps[i]**3 * n_BC_i
	
	return v_BC*1e-9 # in um^3/cm^3

def cal_vertical_distribution_ratio(n_BC, H, **args):
	'''
	This function is to calculate aerosol vertical distribution
	input:
		n_BC   : BC core number distribution, array, cm^-3
		H      : height, m
		**func : vertical distribution function type, A or B, default B
			 type A: Na = N0, if H<H_PBL
			         Na = N0 - k*(H-H_PBL), if H_PBL <= H < H_PBL+H_LD
			         Na = (N0-k*H_LD)*exp(-(H-H_LD-H_PBL)/H_p), if H >= H_PBL + H_LD
			         H_p=1800m, H_PBL=1000m, H_LD=200m, k=28cm^-3m^-1
			 type B: Na = N0exp(-H/H_p)
			         H_p=758m
	output:
		ratio  : BC core number distribution in z to the distribution on the ground
	'''
	if 'func' in args:
		func = args['func']
	else:
		func = 'B'
	
	if func=='B':
		H_p = 758 # m
		ratio = np.exp(-H/H_p)
	else:
		H_p = 1800 # m
		H_PBL = 1000 # m
		H_LD = 200 # m
		k = 28 # cm^-3m^-1, for N0=8000
		n_total = np.sum(n_BC)
		k = k * n_total / 8000 # adapt coefficient k for n_total
		
		if H<H_PBL:
			n_total_H = n_total
		elif H<(H_PBL+H_LD):
			n_total_H = n_total - k * (H-H_PBL)
		else:
			n_total_H = (n_total-k*H_LD) * np.exp(-(H-H_LD-H_PBL)/H_p)
		ratio = n_total_H / n_total
	
	return ratio

def cal_dtau(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, wl, H, dH, **args):
	'''
	This function is to calculate delta aerosol optical depth in certern height and delta height
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
		dH      : delta height, m
		**func  : vertical distribution function type, A or B, default B
		**ddH   : step, default dH/100
	output:
		dtau    : delta aerosol optical depth
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
	kext = cal_kext(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, wl)
	H_integral = 0
	
	for i in range(n):
		H_integral += cal_vertical_distribution_ratio(n_BC, H+n*ddH, func=func) * ddH
	
	dtau = kext * 1e-6 * H_integral
	
	'''
	datu = 0
	for i in range(n):
		kext = cal_kext(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, wl)
		kext_H = cal_vertical_distribution_ratio(n_BC, H+n*ddH, func=func) * kext
		dtau += kext_H * ddH * 1e-6
	'''
	
	return dtau

def cal_SSA(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, wl, H, **args):
	'''
	This function is to calculate single scattering albedo in certern height
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
		**func  : vertical distribution function type, A or B, default B
	output:
		SSA     : single scattering albedo
	'''
	if 'func' in args:
		func = args['func']
	else:
		func = 'B'
	
	PNSD_H = PNSD * cal_vertical_distribution_ratio(n_BC, H, func='A')
	n_BC_H = n_BC * cal_vertical_distribution_ratio(n_BC, H, func='A')
	kext = cal_kext(Dps, PNSD_H, DBCps, DBC, n_BC_H, m_BC, m_shell, wl)
	ksca = cal_ksca(Dps, PNSD_H, DBCps, DBC, n_BC_H, m_BC, m_shell, wl)
	SSA = ksca / kext
	
	return SSA

def cal_ksca(Dps, PNSD, DBCps, DBC, DMASP2, m_BC, m_shell, wl):
	'''
	This function is to use DMA-SP2 and SMPS data to calculate scattering coefficient
	input:
		Dps     : particle diameter size, array, nm
		PNSD    : particle number size distribution, array, dn/dlogDp
		DBCps   : BC particle diameter size, array, nm
		DBC     : BC core diameter size, array, nm
		DMASP2  : BC size-resolved particle number size distribution, array, dn/dlogDBC
		m_BC    : BC core complex refractive index
		m_shell : shell complex refractive index
		wl      : wave length, float, nm
	output:
		k       : scattering coefficient, float, Mm^-1
	'''
	n = calPNSD.PNSD2n(Dps, PNSD, 'd')
	DMASP2PNSD = calPNSD.DMASP22DMASP2PNSD(DBCps, DMASP2)
	k = 0
	
	for i in range(len(Dps)):
		DMASP2PNSD_i = calPNSD.DMASP2PNSD_Dp(DBCps, DMASP2PNSD, Dps[i])
		DMASP2_i = calPNSD.PNSD2n(Dps, DMASP2PNSD_i, 'e')
		n_BC_i = calPNSD.PNSD2n(DBC, DMASP2_i, 'e')
		for j in range(len(DBC)):
			if (DBC[j]>=Dps[i]):
				k += ps.MieQ(m_BC, wl, DBC[j], asCrossSection=True)[1]*1e-18 * n_BC_i[j]*1e6 # in /m
			else:
				k += ps.MieQCoreShell(m_BC, m_shell, wl, DBC[j], Dps[i], asCrossSection=True)[1]*1e-18 * n_BC_i[j]*1e6
		if n[i]:
			n_noBC_i = n[i] - sum(n_BC_i)
			k += ps.MieQ(m_shell, wl, Dps[i], asCrossSection=True)[1]*1e-18 * n_noBC_i*1e6
	
	k *= 1e6 # in /Mm
	return k

def cal_kabs(Dps, PNSD, DBCps, DBC, DMASP2, m_BC, m_shell, wl):
	'''
	This function is to use DMA-SP2 and SMPS data to calculate absorbing coefficient
	input:
		Dps     : particle diameter size, array, nm
		PNSD    : particle number size distribution, array, dn/dlogDp
		DBCps   : BC particle diameter size, array, nm
		DBC     : BC core diameter size, array, nm
		DMASP2  : BC size-resolved particle number size distribution, array, dn/dlogDBC
		m_BC    : BC core complex refractive index
		m_shell : shell complex refractive index
		wl      : wave length, float, nm
	output:
		k       : absorbing coefficient, float, Mm^-1
	'''
	n = calPNSD.PNSD2n(Dps, PNSD, 'd')
	DMASP2PNSD = calPNSD.DMASP22DMASP2PNSD(DBCps, DMASP2)
	k = 0
	
	for i in range(len(Dps)):
		DMASP2PNSD_i = calPNSD.DMASP2PNSD_Dp(DBCps, DMASP2PNSD, Dps[i])
		DMASP2_i = calPNSD.PNSD2n(Dps, DMASP2PNSD_i, 'e')
		n_BC_i = calPNSD.PNSD2n(DBC, DMASP2_i, 'e') # in /cm^3
		for j in range(len(DBC)):
			if (DBC[j]>=Dps[i]):
				k += ps.MieQ(m_BC, wl, DBC[j], asCrossSection=True)[2]*1e-18 * n_BC_i[j]*1e6 # in /m
			else:
				k += ps.MieQCoreShell(m_BC, m_shell, wl, DBC[j], Dps[i], asCrossSection=True)[2]*1e-18 * n_BC_i[j]*1e6
		if n[i]:
			n_noBC_i = n[i] - sum(n_BC_i)
			k += ps.MieQ(m_shell, wl, Dps[i], asCrossSection=True)[2]*1e-18 * n_noBC_i*1e6
	
	k *= 1e6 # in /Mm
	return k

def cal_kext(Dps, PNSD, DBCps, DBC, DMASP2, m_BC, m_shell, wl):
	'''
	This function is to use DMA-SP2 and SMPS data to calculate extinction coefficient
	input:
		Dps     : particle diameter size, array, nm
		PNSD    : particle number size distribution, array, dn/dlogDp
		DBCps   : BC particle diameter size, array, nm
		DBC     : BC core diameter size, array, nm
		DMASP2  : BC size-resolved particle number size distribution, array, dn/dlogDBC
		m_BC    : BC core complex refractive index
		m_shell : shell complex refractive index
		wl      : wave length, nm
	output:
		k       : extinction coefficient, Mm^-1
	'''
	return cal_ksca(Dps, PNSD, DBCps, DBC, DMASP2, m_BC, m_shell, wl) + cal_kabs(Dps, PNSD, DBCps, DBC, DMASP2, m_BC, m_shell, wl)

def retrieve_from_k(Dps, PNSD, DBCps, DBC, DMASP2, ksca, wl_sca, kabs, wl_abs, **args):
	'''
	This function is to use PNSD, BCPNSD data, kabs and ksca data to retrieve m_BC, m_shell and rext
	input:
		Dps         : particle diameter size, array, nm
		PNSD        : number distribution, array, dn/dlogDp
		DBCps       : BC particle diameter size, array, nm
		DBC         : BC core diameter size, array, nm
		DMASP2      : BC size-resolved particle number size distribution, array, dn/dlogDBC
		ksca        : scattering coefficient, array, /Mm
		wl_sca      : scattering coefficient wavelength, array, nm
		kabs        : absorbing coefficient, array, /Mm
		wl_abs      : absorbing coefficient wavelength, array, nm
		**output    : whether to output calculating information, boolean, default True
		**rough_bin : rough calculating bin, float, default 0.1
		**fine_bin  : fine calculating bin, float, default 0.04
	output:
		m_BC        : BC core complex refractive index
		m_shell     : shell complex refractive index
		res         : residual of difference between calculation and observation, float
	'''
	if 'output' in args:
		output = args['output']
	else:
		output = True
	if 'rough_bin' in args:
		rough_bin = args['rough_bin']
	else:
		rough_bin = 0.1
	if 'fine_bin' in args:
		fine_bin = args['fine_bin']
	else:
		fine_bin = 0.02
	
	V_airs = np.arange(0, 1+1e-9, rough_bin)
	m_shells = np.arange(1.34, 1.46+1e-9, fine_bin)
	
	m_EC = 2.26+1.26j
	m_air = 1
	ksca_cals = np.zeros(len(wl_sca))
	kabs_cals = np.zeros(len(wl_abs))
	
	if output:
		print('ksca =', ksca, ', kabs =', kabs)
		print('abs rough calculating...')
	# rough calculate
	res_abs_min = 999999
	for i in range(len(V_airs)):
		m_BC_i = m_EC * (1-V_airs[i]) + m_air * V_airs[i]
		res_abs = 0
		for j in range(len(wl_abs)):
			kabs_cals[j] = cal_kabs(Dps, PNSD, DBCps, DBC, DMASP2, m_BC_i, 1.43+1e-7j, wl_abs[j])
			res_abs += (kabs_cals[j]-kabs[j])**2 / kabs[j]**2
		if np.isnan(res_abs):
			if output:
				print('invalid')
				print('ksca =', ksca, ', kabs =', kabs)
			return np.nan, np.nan, np.nan, np.nan, np.nan
		print(kabs_cals, kabs, res_abs)
		if res_abs<res_abs_min:
			res_abs_min = res_abs
			ii = i
			m_BC = m_BC_i
			kabs_cal = kabs_cals
		else:
			break
	
	if ii==0:
		V_airs = np.arange(V_airs[ii], V_airs[ii+1], fine_bin)
	elif ii==len(V_airs)-1:
		V_airs = np.arange(V_airs[ii-1], V_airs[ii], fine_bin)
	else:
		V_airs = np.arange(V_airs[ii-1], V_airs[ii+1], fine_bin)
	# fine calculate
	res_abs_min = 999999
	if output:
		print('done. abs fine calculating...')
	for i in range(len(V_airs)):
		m_BC_i = m_EC * (1-V_airs[i]) + m_air * V_airs[i]
		res_abs = 0
		for j in range(len(wl_abs)):
			kabs_cals[j] = cal_kabs(Dps, PNSD, DBCps, DBC, DMASP2, m_BC_i, 1.43+1e-7j, wl_abs[j])
			res_abs += (kabs_cals[j]-kabs[j])**2 / kabs[j]**2
		print(kabs_cals, kabs, res_abs)
		if res_abs<res_abs_min:
			res_abs_min = res_abs
			m_BC = m_BC_i
			kabs_cal = kabs_cals
		else:
			break
	
	if output:
		print('done, res_abs =', res_abs_min, 'm_BC =', m_BC)
		print('sca calculating...')	
	res_sca_min = 999999
	for i in range(len(m_shells)):
		res_sca = 0
		for j in range(len(wl_sca)):
			ksca_cals[j] = cal_ksca(Dps, PNSD, DBCps, DBC, DMASP2, m_BC, m_shells[i], wl_sca[j])
			res_sca += (ksca_cals[j]-ksca[j])**2 / ksca[j]**2
		if np.isnan(res_sca):
			if output:
				print('invalid')
				print('ksca =', ksca, ', kabs =', kabs)
			return np.nan, np.nan, np.nan, np.nan, np.nan
		print(ksca_cals, ksca, res_sca)
		if res_sca<res_sca_min:
			res_sca_min = res_sca
			m_shell = m_shells[i]
			ksca_cal = ksca_cals
		else:
			break
	
	res = res_sca_min + res_abs_min
	if output:
		print('done, res_sca =', res_sca_min, 'm_shell =', m_shell, ', res =', res)
	return m_BC, m_shell, ksca_cal, kabs_cal, res

def retrieves(Dps, PNSD, DBCps, DBC, DMASP2, ksca, wl_sca, kabs, wl_abs, **args):
	'''
	This function is to use PNSD, BCPNSD data, kabs and ksca data to retrieve m_BC, m_shell and rext, and save data into .npy
	input:
		Dps         : particle diameter size, array, nm
		PNSD        : number distribution, array, dn/dlogDp
		DBCps       : BC particle diameter size, array, nm
		DBC         : BC core diameter size, array, nm
		DMASP2      : BC size-resolved particle number size distribution, array, dn/dlogDBC
		ksca        : scattering coefficient, array, /Mm
		wl_sca      : scattering coefficient wavelength, array, nm
		kabs        : absorbing coefficient, array, /Mm
		wl_abs      : absorbing coefficient wavelength, array, nm
		**output    : whether to output calculating information, boolean, default True
		**rough_bin : rough calculating bin, float, default 0.1
		**fine_bin  : fine calculating bin, float, default 0.04
		**path      : data save path, string, default 'data/retrieve/retrieveFromK.npy'
	output:
		retrieveFromK.npy
	'''
	if 'output' in args:
		output = args['output']
	else:
		output = True
	if 'rough_bin' in args:
		rough_bin = args['rough_bin']
	else:
		rough_bin = 0.1
	if 'fine_bin' in args:
		fine_bin = args['fine_bin']
	else:
		fine_bin = 0.02
	if 'path' in args:
		path = args['path']
	else:
		path = 'data/retrieve/retrieveFromK.npy'
	
	m_BC = []
	m_shell = np.zeros(len(PNSD))
	ksca_cal = np.zeros((len(PNSD), len(wl_sca)))
	kabs_cal = np.zeros((len(PNSD), len(wl_abs)))
	res = np.zeros(len(PNSD))
	
	print('calculating...')
	for i in range(len(PNSD)):
		m_BC_i, m_shell[i], ksca_cal[i], kabs_cal[i], res[i] = retrieve_from_k(Dps, PNSD[i], DBCps, DBC, DMASP2[i], ksca[i], wl_sca, kabs[i], wl_abs, output=output, rough_bin=rough_bin, fine_bin=fine_bin, path=path)
		m_BC.append(m_BC_i)
		try:
			print(m_BC[i], round(m_shell[i]*1000)/1000, round(res[i]*100)/100)
		except:
			print('nan nan nan nan nan')
		print(round((i+1)/len(PNSD)*1000)/10, '% done...')
	m_BC = np.array(m_BC)
	print('done')
	
	data = dict(m_BC=m_BC, m_shell=m_shell, ksca_cal=ksca_cal, kabs_cal=kabs_cal, res=res)
	np.save(path, data)

def cal_ksca_bulk(Dps, PNSD, DBC, BCPNSD, m_BC, m_shell, wl, rext):
	'''
	This function is to use BC PNSD data to calculate scattering coefficient kabs
	input:
		Dps     : particle diameter size, array, nm
		PNSD    : number distribution, array, dn/dlogDp
		DBC     : BC core diameter size, array, nm
		BCPNSD  : BC particle number size distribution, array, dn/dlogDp
		m_BC    : BC core complex refractive index
		m_shell : shell complex refractive index
		wl      : wave length, nm
		rext    : external mixture rate, float, range from 0 to 1
	output:
		k       : scattering coefficient, Mm^-1
	'''
	k = 0
	n = calPNSD.PNSD2n(Dps, PNSD, 'd')
	
	for i in range(len(Dps)):
		BCPNSD_i = calPNSD.PNSD_Dp(DBC, BCPNSD, Dps[i])
		n_BC_i = calPNSD.PNSD2n(Dps, BCPNSD_i, 'e')
		n_ext_i = n_BC_i * rext # in /cm^3
		k += ps.MieQ(m_BC, wl, Dps[i], asCrossSection=True)[1]*1e-18 * n_ext_i*1e6
		if n[i]:
			n_int_i = n[i] - n_ext_i
			v_BC_int = n_BC_i * np.pi/6 * Dps[i]**3 * (1-rext) # in nm^3/cm^3
			D_BC_int = (6*v_BC_int/np.pi/n_int_i)**(1/3) # in nm
			if D_BC_int and not np.isnan(D_BC_int):
				k += ps.MieQCoreShell(m_BC, m_shell, wl, D_BC_int, Dps[i], asCrossSection=True)[1]*1e-18 * n_int_i*1e6
			else: # no BC core
				k += ps.MieQ(m_shell, wl, Dps[i], asCrossSection=True)[1]*1e-18 * n_int_i*1e6
	
	k *= 1e6 # in /Mm
	return k

def cal_kabs_bulk(Dps, PNSD, DBC, BCPNSD, m_BC, m_shell, wl, rext):
	'''
	This function is to use BC PVSD data to calculate absorbing coefficient kabs
	input:
		Dps     : particle diameter size, array, nm
		PNSD    : number distribution, array, dn/dlogDp
		DBC     : BC core diameter size, array, nm
		BCPNSD  : BC particle number size distribution, array, dn/dlogDp
		m_BC    : BC core complex refractive index
		m_shell : shell complex refractive index
		wl      : wave length, nm
		rext    : external mixture rate, float, range from 0 to 1
	output:
		k       : absorbing coefficient, Mm^-1
	'''
	k = 0
	n = calPNSD.PNSD2n(Dps, PNSD, 'd')
	
	for i in range(len(Dps)):
		BCPNSD_i = calPNSD.PNSD_Dp(DBC, BCPNSD, Dps[i])
		n_BC_i = calPNSD.PNSD2n(Dps, BCPNSD_i, 'e')
		n_ext_i = n_BC_i * rext # in /cm^3
		k += ps.MieQ(m_BC, wl, Dps[i], asCrossSection=True)[2]*1e-12 * n_ext_i
		if n[i]:
			n_int_i = n[i] - n_ext_i
			v_BC_int = n_BC_i * (1-rext) * np.pi/6 * Dps[i]**3 # in nm^3/cm^3
			D_BC_int = (6*v_BC_int/np.pi/n_int_i)**(1/3) # in nm
			if D_BC_int and not np.isnan(D_BC_int):
				k += ps.MieQCoreShell(m_BC, m_shell, wl, D_BC_int, Dps[i], asCrossSection=True)[2]*1e-12 * n_int_i
			else:
				k += ps.MieQ(m_shell, wl, Dps[i], asCrossSection=True)[2]*1e-18 * n_int_i*1e6
	
	k *= 1e6 # in /Mm
	return k

def retrieve_from_k_bulk(Dps, PNSD, DBC, BCPNSD, ksca, wl_sca, kabs, wl_abs, **args):
	'''
	This function is to use PNSD, BCPNSD data, kabs and ksca data to retrieve m_BC, m_shell and rext
	input:
		Dps         : particle diameter size, array, nm
		PNSD        : number distribution, array, dn/dlogDp
		DBC         : BC core diameter size, array, nm
		BCPNSD      : BC particle number size distribution, array, dn/dlogDp
		ksca        : scattering coefficient, array, /Mm
		wl_sca      : scattering coefficient wavelength, array, nm
		kabs        : absorbing coefficient, array, /Mm
		wl_abs      : absorbing coefficient wavelength, array, nm
		**output    : whether to output calculating information, boolean, default True
		**rough_bin : rough calculating bin, float, default 0.1
		**fine_bin  : fine calculating bin, float, default 0.04
	output:
		m_BC        : BC core complex refractive index
		m_shell     : shell complex refractive index
		rext        : external mixture rate, float, from 0 to 1
		res         : residual of difference between calculation and observation, float
	'''
	if 'output' in args:
		output = args['output']
	else:
		output = True
	if 'rough_bin' in args:
		rough_bin = args['rough_bin']
	else:
		rough_bin = 0.1
	if 'fine_bin' in args:
		fine_bin = args['fine_bin']
	else:
		fine_bin = 0.02
	
	rexts = np.arange(0, 1+1e-9, rough_bin)
	V_airs = np.arange(0, 1+1e-9, rough_bin)
	m_shells = np.arange(1.3, 1.5+1e-9, fine_bin)
	
	m_EC = 2.26+1.26j
	m_air = 1
	ksca_cals = np.zeros(len(wl_sca))
	kabs_cals = np.zeros(len(wl_abs))
	
	if output:
		print(ksca, kabs)
		print('rough calculating...')
	res_abs_min = 999999
	res_min = 999999
	for i in range(len(rexts)):
		for j in range(len(V_airs)):
			m_BC_j = m_EC * (1-V_airs[j]) + m_air * V_airs[j]
			res_abs = 0
			for k in range(len(wl_abs)):
				kabs_cals[k] = cal_kabs_bulk(Dps, PNSD, DBC, BCPNSD, m_BC_j, 1.43+1e-7j, wl_abs[k], rexts[i])
				res_abs += (kabs_cals[k]-kabs[k])**2 / kabs[k]**2
			if np.isnan(res_abs):
				if output:
					print('invalid')
					print('ksca =', ksca, ', kabs =', kabs)
				return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
			if res_abs<res_abs_min:
				print(kabs_cals, kabs, res_abs)
				res_abs_min = res_abs
				m_BC = m_BC_j
		
		for k in range(len(m_shells)):
			res_sca = 0
			for l in range(len(wl_sca)):
				ksca_cals[l] = cal_ksca_bulk(Dps, PNSD, DBC, BCPNSD, m_BC, m_shells[k], wl_sca[l], rexts[i])
				res_sca += (ksca_cals[l]-ksca[l])**2 / ksca[l]**2
			if np.isnan(res_sca):
				if output:
					print('invalid')
					print('ksca =', ksca, ', kabs =', kabs)
				return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
		
		res = res_abs + res_sca
		if res<res_min:
			print(ksca_cals, ksca, res)
			res_min = res
			ii = i
			jj = j
			rext = rexts[i]
			m_BC = m_BC_j
			m_shell = m_shells[k]
			ksca_cal = ksca_cals
			kabs_cal = kabs_cals
	
	if output:
		print('done. res =', res, ', fine calculating...')
	if ii==0:
		rexts = np.arange(rexts[ii], rexts[ii+1], fine_bin)
	elif ii==len(rexts)-1:
		rexts = np.arange(rexts[ii-1], rexts[ii], fine_bin)
	else:
		rexts = np.arange(rexts[ii-1], rexts[ii+1], fine_bin)
	if jj==0:
		V_airs = np.arange(V_airs[jj], V_airs[jj+1], fine_bin)
	elif jj==len(V_airs)-1:
		V_airs = np.arange(V_airs[jj-1], V_airs[jj], fine_bin)
	else:
		V_airs = np.arange(V_airs[jj-1], V_airs[jj+1], fine_bin)
	
	res_abs_min = 999999
	for i in range(len(rexts)):
		for j in range(len(V_airs)):
			m_BC_j = m_EC * (1-V_airs[j]) + m_air * V_airs[j]
			res_abs = 0
			for k in range(len(wl_abs)):
				kabs_cals[k] = cal_kabs_bulk(Dps, PNSD, DBC, BCPNSD, m_BC_j, 1.43+1e-7j, wl_abs[k], rexts[i])
				res_abs += (kabs_cals[k]-kabs[k])**2 / kabs[k]**2
			if res_abs<res_abs_min:
				res_abs_min = res_abs
				rext = rexts[i]
				m_BC = m_BC_j
				kabs_cal = kabs_cals
	
	if output:
		print('done, res_abs =', res_abs_min)
	return m_BC, m_shell, rext, ksca_cal, kabs_cal, res

def retrieves_bulk(Dps, PNSD, DBC, BCPNSD, ksca, wl_sca, kabs, wl_abs, **args):
	'''
	This function is to use PNSD, BCPNSD data, kabs and ksca data to retrieve m_BC, m_shell and rext, and save data into .npy
	input:
		Dps         : particle diameter size, array, nm
		PNSD        : number distribution, array, dn/dlogDp
		DBC         : BC core diameter size, array, nm
		BCPNSD      : BC particle number size distribution, array, dn/dlogDp
		ksca        : scattering coefficient, array, /Mm
		wl_sca      : scattering coefficient wavelength, array, nm
		kabs        : absorbing coefficient, array, /Mm
		wl_abs      : absorbing coefficient wavelength, array, nm
		**output    : whether to output calculating information, boolean, default True
		**rough_bin : rough calculating bin, float, default 0.1
		**fine_bin  : fine calculating bin, float, default 0.04
		**path      : data save path, string, default 'data/retrieve/retrieveFromKBulk.npy'
	output:
		retrieveFromK.npy
	'''
	if 'output' in args:
		output = args['output']
	else:
		output = True
	if 'rough_bin' in args:
		rough_bin = args['rough_bin']
	else:
		rough_bin = 0.1
	if 'fine_bin' in args:
		fine_bin = args['fine_bin']
	else:
		fine_bin = 0.02
	if 'path' in args:
		path = args['path']
	else:
		path = 'data/retrieve/retrieveFromKBulk.npy'
	
	m_BC = []
	m_shell = np.zeros(len(PNSD))
	rext = np.zeros(len(PNSD))
	ksca_cal = np.zeros((len(PNSD), len(wl_sca)))
	kabs_cal = np.zeros((len(PNSD), len(wl_abs)))
	res = np.zeros(len(PNSD))
	
	print('calculating...')
	for i in range(len(PNSD)):
		m_BC_i, m_shell[i], rext[i], ksca_cal[i], kabs_cal[i], res[i] = retrieve_from_k_bulk(Dps, PNSD[i], DBC, BCPNSD[i], ksca[i], wl_sca, kabs[i], wl_abs, output=output, rough_bin=rough_bin, fine_bin=fine_bin, path=path)
		m_BC.append(m_BC_i)
		try:
			print(m_BC[i], round(m_shell[i]*100)/100, round(rext[i]*100)/100, round(res[i]*100)/100)
		except:
			print('nan nan nan nan nan nan')
		print(round((i+1)/len(PNSD)*1000)/10, '% done...')
	m_BC = np.array(m_BC)
	print('done')
	
	data = dict(m_BC=m_BC, m_shell=m_shell, ksca_cal=ksca_cal, kabs_cal=kabs_cal, rext=rext, res=res)
	np.save(path, data)

def cal_BCI(DBC, DBCps, BCPNSD, **args):
	'''
	This function is to calculate BC mixing states entropy
	input:
		DBCps       : BC particle diameter size, array, nm
		DBC         : BC core diameter size, array, nm
		BCPNSD      : BC size-resolved particle number size distribution, array, dn/dlogDBC
		args:
			rho   : totaol particals average density, float, g/cm3, default 1
			rhoBC : BC core average density, float, g/cm3, default 0.95
	output:
		BCI         : BC mixing inhomogeneity
	'''
	if 'rho' in args:
		rho = args['rho']
	else:
		rho = 1
	if 'rhoBC' in args:
		rhoBC = args['rhoBC']
	else:
		rhoBC = 0.95
	
	BCPNSD = BCPNSD * calPNSD.cal_dlogDp(DBC) # in cm-3
	mBC = 0
	mTot = 0
	HAlpha = 0
	miBC = np.zeros(BCPNSD.shape)
	miTot = np.zeros(BCPNSD.shape)
	Hi = np.zeros(BCPNSD.shape)
	
	for i in range(len(DBCps)):
		for j in range(len(DBC)):
			miBC[i,j] = np.pi / 6 * rhoBC * (DBC[j]*1e-7)**3 * 1e12 * BCPNSD[i,j]
			if DBCps[i]>=DBC[j]:
				miTot[i,j] = miBC[i,j] + np.pi / 6 * rho * ((DBCps[i]*1e-7)**3-(DBC[j]*1e-7)**3) * 1e12 * BCPNSD[i,j]
			else:
				miTot[i,j] = miBC[i,j]
			if miBC[i,j]:
				piBC = miBC[i,j] / miTot[i,j]
			else:
				piBC = 0
			if piBC<1 and piBC>0:
				Hi[i,j] = -piBC * np.log(piBC) + (1-piBC) * np.log(1-piBC)
	
	for i in range(len(DBCps)):
		for j in range(len(DBC)):
			pi = miTot[i,j] / np.nansum(miTot)
			HAlpha += pi * Hi[i,j]
	
	DAlpha = np.exp(HAlpha)
	pBC = np.nansum(miBC) / np.nansum(miTot)
	HGamma = -pBC * np.log(pBC) + (1-pBC) * np.log(1-pBC)
	DGamma = np.exp(HGamma)
	BCI = (DAlpha-1) / (DGamma-1)
	
	return BCI

def cal_BCI_change_af(DBC, DBCps, BCPNSD, rate, **args):
	'''
	This function is to calculation when BCI change certain rate, how BC core diameter react. To satisfy other parameter stable, have to do this in conservation of total BC mass.
	input:
		DBCps        : BC particle diameter size, array, nm
		DBC          : BC core diameter size, array, nm
		BCPNSD       : BC size-resolved particle number size distribution, array, dn/dlogDBC
		rate         : BCI change rate
		args:
			rough bin : rough calculating bin, default 0.1
			fin bin   : fine calculating bin, default 0.01
	output:
		BCI_final    : BCI after adjust
		BCPNSD_final : BCPNSD after adjust
		af           : adjust factor for BCPNSD mixing parameter, new BCPNSD will be: 
				BCPNSD_origin * af + BCPNSD_maxetp * (1-af) or
				BCPNSD_origin * af + BCPNSD_minetp * (1-af)
				depends on rate > 1 or rate < 1
	'''
	if 'rough_bin' in args:
		rough_bin = args['rough_bin']
	else:
		rough_bin = 0.1
	if 'fine_bin' in args:
		fine_bin = args['fine_bin']
	else:
		fine_bin = 0.01
	
	# BCPNSD_maxetp: all BC mass distributed in every particals evenly
	# BCPNSD_minetp: all BC mass have same diameter as particals, thus DBC[i,j] = DBCps[i,-1]
	
	BCPNSD = BCPNSD * calPNSD.cal_dlogDp(DBC)
	BCI = cal_BCI(DBC, DBCps, BCPNSD)
	BCI = BCI * rate
	j_maxetp = round(len(DBC)/2)
	j_minetp = -1
	BCPNSD_maxetp = np.zeros(BCPNSD.shape)
	BCPNSD_minetp = np.zeros(BCPNSD.shape)
	
	for i in range(len(DBCps)):
		mBCi = 0
		for j in range(len(DBC)):
			mBCi += DBC[j]**3 * BCPNSD[i,j] # deleted unnecessary const
		BCPNSDij_maxetp = mBCi / DBC[j_maxetp]**3
		BCPNSD_maxetp[i,j_maxetp] = BCPNSDij_maxetp / calPNSD.cal_dlogDp(DBC)
		BCPNSDij_minetp = mBCi / DBC[j_minetp]**3
		BCPNSD_minetp[i,j_minetp] = BCPNSDij_minetp / calPNSD.cal_dlogDp(DBC)
		# except max entropy and min entropy DBC diameter, all is zero
	
	# calculate mixing parameter af
	# rough bin: 0.1, fine bin: 0.01
	afs = np.arange(0, 1+1e-9, rough_bin)
	res_min = 99999999
	ii = 0 # mark the minimum af index
	
	for i in range(len(afs)):
		if rate < 1:
			BCPNSD_new = BCPNSD * afs[i] + BCPNSD_maxetp * (1-afs[i])
		else:
			BCPNSD_new = BCPNSD * afs[i] + BCPNSD_minetp * (1-afs[i])
		BCI_new = cal_BCI(DBC, DBCps, BCPNSD_new)
		res = abs(BCI_new-BCI) / BCI
		if res < res_min:
			res_min = res
			ii = i
	
	if ii==0:
		afs = np.arange(afs[ii], afs[ii+1], fine_bin)
	elif ii==len(afs)-1:
		afs = np.arange(afs[ii-1], afs[ii], fine_bin)
	else:
		afs = np.arange(afs[ii-1], afs[ii+1], fine_bin)
	
	res_min = 9999999
	ii = 0
	
	for i in range(len(afs)):
		if rate < 1:
			BCPNSD_new = BCPNSD * afs[i] + BCPNSD_maxetp * (1-afs[i])
		else:
			BCPNSD_new = BCPNSD * afs[i] + BCPNSD_minetp * (1-afs[i])
		BCI_new = cal_BCI(DBC, DBCps, BCPNSD_new)
		res = abs(BCI_new-BCI) / BCI
		if res < res_min:
			res_min = res
			ii = i
			BCI_final = BCI_new
			BCPNSD_final = BCPNSD_new
	
	return BCI_final, BCPNSD_final, afs[ii]

if __name__ == '__main__':
	print('loading...')
	data = readTaizhou.read_Taizhou('data/sp2/Taizhou.npy')
	'''
	sp2 data:
	Dps     : partical diameter distribution, array, nm
	PNSD    : partical number concentration size distribution, array, cm^-3
	DBC     : BC core diameter distribution, array, nm
	DBCps   : BC core have partical diameter distribution, array, nm
	DMASP2  : BC size-resolved particle number size distribution, array, dn/dlogDBC
	BCPNSD  : BC partical number size distribution, array, dn/dlogDBC
	ksca    : scattering coefficient, Mm^-1
	wl_sca  : ksca wave length, nm
	kabs    : absorbing coefficient, Mm^-1
	wl_abs  : kabs wave length, nm
	'''
	
	Dps = data['Dps']
	PNSD = data['PNSD'][100:500]
	DBCps = data['DBCps']
	DBC = data['DBC']
	DMASP2 = data['DMASP2'][100:500]
	BCPNSD = data['BCPNSD'][100:500]
	BCPVSD = data['BCPVSD'][100:500]
	kabs = data['kabs'][100:500]
	ksca = data['ksca'][100:500]
	wl_abs = data['wl_abs']
	wl_sca = data['wl_sca']
	print('done')
	print('calculating...')
	'''
	# BC volume closure
	for i in range(len(PNSD)):
		v_BC1 = BCPNSD2vBC(DBC, BCPNSD[i])
		v_BC2 = calPNSD.PNSD2n(DBC, BCPVSD[i], 'e')
		v_BC3 = DMASP22vBC(DBCps, DBC, DMASP2[i])
		v_BC4 = DMASP22vBC2(Dps, DBCps, DBC, DMASP2[i])
		v_BC5 = BCPNSD2vBC2(Dps, DBC, BCPNSD[i])
		print(sum(v_BC1), sum(v_BC2), sum(v_BC3), sum(v_BC4), v_BC5)
	'''
	
	retrieves(Dps, PNSD, DBCps, DBC, DMASP2, ksca, wl_sca, kabs, wl_abs, path='data/retrieve/retrieveFromK_3wls.npy')
	'''
	retrieves(Dps, PNSD, DBCps, DBC, DMASP2, ksca, wl_sca, kabs[:,1:2], wl_abs[1:2], path='data/retrieve/retrieveFromK_532.npy')
	retrieves(Dps, PNSD, DBCps, DBC, DMASP2, ksca, wl_sca, kabs[:,:2], wl_abs[:2], path='data/retrieve/retrieveFromK_781&532.npy')
	retrieves(Dps, PNSD, DBCps, DBC, DMASP2, ksca, wl_sca, kabs[:,1:], wl_abs[1:], path='data/retrieve/retrieveFromK_532&405.npy')
	
	retrieves_bulk(Dps, PNSD, DBC, BCPNSD, ksca, wl_sca, kabs, wl_abs, path='data/retrieve/retrieveFromKBulk_3wls.npy')
	retrieves_bulk(Dps, PNSD, DBC, BCPNSD, ksca, wl_sca, kabs[:,1:2], wl_abs[1:2], path='data/retrieve/retrieveFromKBulk_532.npy')
	retrieves_bulk(Dps, PNSD, DBC, BCPNSD, ksca, wl_sca, kabs[:,:2], wl_abs[:2], path='data/retrieve/retrieveFromKBulk_781&532.npy')
	retrieves_bulk(Dps, PNSD, DBC, BCPNSD, ksca, wl_sca, kabs[:,1:], wl_abs[1:], path='data/retrieve/retrieveFromKBulk_532&405.npy')
	'''
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
