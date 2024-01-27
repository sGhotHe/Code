####################################################################################
# INTRODUCTION:
# This code is kappa-Kohler theory and some related calculation
# Created by Hebs at 21/11/9/14:46
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import PyMieScatt as ps
from scipy.optimize import fsolve
import calBC
import calPNSD
import readTaizhou
import sys

'''
How to retrieve kappa from f(RH)?
Here is kappa-Kohler theory:
S = a_w * exp(4*sigma_sa*M_w/R/T/rho_w/D)
1 / a_w = 1 + kappa * V_s / V_w
For core-shell structure, asume only one homogeneous composition coated, so there is:
V_w = a_w / (1-a_w) * kappa * V_s
There is
V_w = 1/6 * pi * (D**3-Dd**3)
V_s = 1/6 * pi * (Dd**3-DBC**3)
So here we have
S = 1 / (1+kappa*(Dd**3-DBC**3)/(D**3-Dd**3)) * exp(4*sigma_sa*M_w/R/T/rho_w/D)
Use RH and kappa data, we can calculate after hygroscopic growth particle diameter D, so as V_w
Use linear hypothesis, we can calculate complex refractive index for shell, by:
m_shell = m_shell_origin * V_s / (V_s+V_w) + 1 * V_w / (V_s+V_w)
Then we can use Mie theory to calculate particle scattering and absorbing cross section
Use PNSD and BCPNSD data, bulk scattring and absorbing coefficient can be calculated
Compare with f(RH) data, kappa can be asured
'''

def kappa_Kohler(D, Dd, DBC, kappa, **args):
	'''
	This function is kappa-Kohler theory
	input:
		D          : diameter, float, nm
		Dd         : dry diameter, float, nm
		DBC        : BC core diameter, float, nm
		kappa      : hygroscopicity parameter
		**sigma_sa : surface tension of the solution/air interface, float, J/m^2, default 0.072
		**M_w      : molecular weight of water, float, kg/mol, default 18/1000
		**R        : universal gas constant, float, default 8.31
		**T        : temperature, float, Kelvin, default 298.15
		**rho_w    : density of water, float, kg/m^3, default 1000
	output:
		S          : saturation ratio
	'''
	if 'sigma_sa' in args:
		sigma_sa = args['sigma_sa']
	else:
		sigma_sa = 0.072
	if 'M_w' in args:
		M_w = args['M_w']
	else:
		M_w = 18 / 1000
	if 'R' in args:
		R = args['R']
	else:
		R = 8.31
	if 'T' in args:
		T = args['T']
	else:
		T = 298.15
	if 'rho_w' in args:
		rho_w = args['rho_w']
	else:
		rho_w = 1000
	
	a = 1 / (1+kappa*(Dd**3-DBC**3)/(D**3-Dd**3))
	b = np.exp(4*sigma_sa*M_w/R/T/rho_w/(D*1e-9))
	S = a * b
	
	return S

def RH2D(Dd, DBC, m_BC, m_shell, kappa, RH, wl):
	'''
	This function is to calculate particle diameter and m_shell in certain RH
	input:
		Dd         : dry particle diameter, float, nm
		DBC        : BC core diameter, float, nm
		m_BC       : BC core complex refractive index, float
		m_shell    : shell complex refractive index, float
		kappa      : shell hygroscocipicity parameter, float
		RH         : relative humidity, float, percent
		wl         : wave length, float, nm
	output:
		D          : particle diameter after hygroscopic, float, nm
		m_shell_RH : m_shell after hygroscopic
	'''
	if DBC>Dd:
		return DBC, m_shell
	func_D = lambda D : kappa_Kohler(D, Dd, DBC, kappa) - RH / 100
	D = fsolve(func_D, Dd+1e-7)[0]
	V_s = 1/6 * np.pi * ((Dd*1e-9)**3-(DBC*1e-9)**3)
	V_w = 1/6 * np.pi * ((D*1e-9)**3-(Dd*1e-9)**3)
	m_shell_RH = m_shell * V_s / (V_s+V_w) + 1.3325 * V_w / (V_s+V_w)
	
	return D, m_shell_RH

def RH2ksca(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, kappa, RH, wl):
	'''
	This function is to use PNSD, BCPNSD, kappa and RH data to calculate bulk scattring coefficient ksca
	input:
		Dps     : particle diameter size, array, nm
		PNSD    : particle number size distribution, array, dn/dlogDp
		DBCps   : BC particle diameter size, array, nm
		DBC     : BC core diameter size, array, nm
		n_BC    : BC particle number concentration, array, cm^-3
		m_BC    : BC core complex refractive index
		m_shell : shell complex refractive index
		kappa   : hygroscopicity parameter, float
		RH      : relative humidity, float, percent
		wl      : wave length, float, nm
	output:
		ksca_RH : bulk scattering coefficient
	'''
	DBCps_min = DBCps[0]
	DBCps_max = DBCps[-1]
	Dps_min = Dps[0]
	Dps_max = Dps[-1]
	n = calPNSD.PNSD2n(Dps, PNSD, index='e')
	
	PNSD_fit = np.zeros(len(DBCps))
	n_fit = np.zeros(len(DBCps))
	n_noBC = np.zeros(len(DBCps))
	ksca_RH = 0
	
	for i in range(len(Dps)):
		if Dps[i]<DBCps_min or Dps[i]>DBCps_max:
			# use PNSD to calculate ksca
			D_i, m_shell_RH_i = RH2D(Dps[i], 0, m_BC, m_shell, kappa, RH, wl)
			Qsca_noBCi = ps.MieQ(m_shell_RH_i, wl, D_i)[1]
			ksca_RH += Qsca_noBCi * 1/4 * np.pi * (D_i*1e-9)**2 * n[i]*1e6
	
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
			n_fit[i] = PNSD_fit[i] * calPNSD.cal_dlogDp(DBCps)
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
			D_i, m_shell_RH_i = RH2D(DBCps[i], 0, m_BC, m_shell, kappa, RH, wl)
			Qsca_noBCi = ps.MieQ(m_shell_RH_i, wl, D_i)[1]
			ksca_RH += Qsca_noBCi * 1/4 * np.pi * (D_i*1e-9)**2 * n_noBC[i]*1e6
		# else:
		# DBCps[i] is either bigger than Dps_max or smaller than Dps_min
		# in this situation, n_fit[i] and n_noBC[i] is zero
		# just jump over
		
		for j in range(len(DBC)):
			if DBC[j]<=DBCps[i]:
				D_ij, m_shell_RH_ij = RH2D(DBCps[i], DBC[j], m_BC, m_shell, kappa, RH, wl)
				Qsca_BCij = ps.MieQCoreShell(m_BC, m_shell_RH_ij, wl, DBC[j], D_ij)[1]
				ksca_RH += Qsca_BCij * 1/4 * np.pi * (D_ij*1e-9)**2 * n_BC[i,j]*1e6
	
	ksca_RH = ksca_RH * 1e6 # Mm^-1
	return ksca_RH

def RH2kext(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, kappa, RH, wl):
	'''
	This function is to use PNSD, BCPNSD, kappa and RH data to calculate bulk extinction coefficient kext
	input:
		Dps     : particle diameter size, array, nm
		PNSD    : particle number size distribution, array, dn/dlogDp
		DBCps   : BC particle diameter size, array, nm
		DBC     : BC core diameter size, array, nm
		n_BC    : BC particle number concentration, array, cm^-3
		m_BC    : BC core complex refractive index
		m_shell : shell complex refractive index
		kappa   : hygroscopicity parameter, float
		RH      : relative humidity, float, percent
		wl      : wave length, float, nm
	output:
		kext_RH : bulk extinction coefficient
	'''
	DBCps_min = DBCps[0]
	DBCps_max = DBCps[-1]
	Dps_min = Dps[0]
	Dps_max = Dps[-1]
	n = calPNSD.PNSD2n(Dps, PNSD, index='e')
	
	PNSD_fit = np.zeros(len(DBCps))
	n_fit = np.zeros(len(DBCps))
	n_noBC = np.zeros(len(DBCps))
	kext_RH = 0
	
	for i in range(len(Dps)):
		if Dps[i]<DBCps_min or Dps[i]>DBCps_max:
			D_i, m_shell_RH_i = RH2D(Dps[i], 0, m_BC, m_shell, kappa, RH, wl)
			Qext_noBCi = ps.MieQ(m_shell_RH_i, wl, D_i)[0]
			kext_RH += Qext_noBCi * 1/4 * np.pi * (D_i*1e-9)**2 * n[i]*1e6
	
	for i in range(len(DBCps)):
		j = 0
		while j<len(Dps)-1:
			if Dps[j]<=DBCps[i] and Dps[j+1]>DBCps[i]:
				break
			j += 1
		if j<len(Dps)-1:
			PNSD_fit[i] = PNSD[j] + (PNSD[j+1]-PNSD[j]) * (np.log(DBCps[i])-np.log(Dps[j])) / (np.log(Dps[j+1])-np.log(Dps[j]))
			n_fit[i] = PNSD_fit[i] * calPNSD.cal_dlogDp(DBCps)
			for j in range(len(DBC)):
				if DBC[j]>DBCps[i]:
					n_BC[i,j] = 0
			if sum(n_BC[i])>n_fit[i]:
				n_BC[i] = n_BC[i] * n_fit[i] / sum(n_BC[i])
				# n_noBC[i] = 0
			else:
				n_noBC[i] = n_fit[i] - sum(n_BC[i])
			D_i, m_shell_RH_i = RH2D(DBCps[i], 0, m_BC, m_shell, kappa, RH, wl)
			Qext_noBCi = ps.MieQ(m_shell_RH_i, wl, D_i)[0]
			kext_RH += Qext_noBCi * 1/4 * np.pi * (D_i*1e-9)**2 * n_noBC[i]*1e6
		
		for j in range(len(DBC)):
			if DBC[j]<=DBCps[i]:
				D_ij, m_shell_RH_ij = RH2D(DBCps[i], DBC[j], m_BC, m_shell, kappa, RH, wl)
				Qext_BCij = ps.MieQCoreShell(m_BC, m_shell_RH_ij, wl, DBC[j], D_ij)[0]
				kext_RH += Qext_BCij * 1/4 * np.pi * (D_ij*1e-9)**2 * n_BC[i,j]*1e6
	
	kext_RH = kext_RH * 1e6 # Mm^-1
	return kext_RH

def RH2fRH(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, kappa, RH, wl):
	'''
	This function is to use PNSD, BCPNSD, kappa and RH data to calculate f(RH)
	input:
		Dps     : particle diameter size, array, nm
		PNSD    : particle number size distribution, array, dn/dlogDp
		DBCps   : BC particle diameter size, array, nm
		DBC     : BC core diameter size, array, nm
		n_BC    : BC particle number concentration, array, cm^-3
		m_BC    : BC core complex refractive index
		m_shell : shell complex refractive index
		kappa   : hygroscopicity parameter, float
		RH      : relative humidity, float, percent
		wl      : wave length, float, nm
	output:
		fRH    : f(RH)
	'''
	# calculate ksca_RH
	ksca_RH = RH2ksca(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, kappa, RH, wl)
	ksca_dry = calBC.cal_ksca(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, wl)
	fRH = ksca_RH / ksca_dry
	
	return fRH

def fRH2kappa(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, wl, RH, fRH):
	'''
	This function is to use f(RH) data to retrieve kappa
	input:
		Dps     : particle diameter size, array, nm
		PNSD    : particle number size distribution, array, dn/dlogDp
		DBCps   : BC particle diameter size, array, nm
		DBC     : BC core diameter size, array, nm
		n_BC    : BC particle number concentration, array, cm^-3
		m_BC    : BC core complex refractive index
		m_shell : shell complex refractive index
		wl      : wave length, float, nm
		RH      : relative humidity, float, percent
		fRH     : f(RH), float
	output:
		kappa      : hygroscopicity parameter, float
	'''
	func_kappa = lambda kappa : RH2fRH(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, kappa, RH, wl) - fRH
	kappa = fsolve(func_kappa, 1e-1)[0]
	
	return kappa

if __name__ == '__main__':
	data = readTaizhou.read_Taizhou('data/sp2/Taizhou.npy')
	Dps = data['Dps']
	PNSD = np.nanmean(data['PNSD'],axis=0)
	DBCps = data['DBCps']
	DBC = data['DBC']
	n_BC = np.nanmean(data['n_BC'],axis=0)
	m_BC = 1.67+0.67j
	m_shell = 1.5
	
	Dd = 100
	D_BC = 50
	m_BC = 1.67+0.67j
	m_shell = 1.5
	RH = 80
	kappa = 0.5
	wl = 525
	'''
	S = kappa_Kohler(Dd*1.1, Dd, D_BC, kappa)
	print(S)
	D, m_shell_RH = RH2D(Dd, D_BC, m_BC, m_shell, kappa, RH, wl)
	print(D, m_shell_RH)
	
	fRH = RH2fRH(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, kappa, RH, wl)
	print(fRH)
	kappa = fRH2kappa(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, wl, RH, fRH)
	print(kappa)
	'''
	ksca_RH = RH2ksca(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, kappa, RH, wl)
	print(ksca_RH)
	kext_RH = RH2kext(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, kappa, RH, wl)
	print(kext_RH)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
