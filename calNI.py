####################################################################################
# INTRODUCTION:
# This code is to run complex refractive index inhomogeneity's influence on Mie function
# Created by Hebs at 23/5/15/13:53
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import PyMieScatt as ps

import calPNSD
import kappaKohler
import readAtms

def cal_I2(par, binNum, I1, I2x, I2y):
	'''
	This function is to calculate complex refractive index's 2 level inhomogeneity
	input:
		par		: parameter, float
		binNum		: Dps bin number, int
		I1		: size resolved n change slope, float
		I2x		: peaks number, int
		I2y		: peaks separate, float
	output:
		par_new	: new parameter, float, array in shape [binNum,I2x]
	'''
	# define nI as a two dimensional parameter, 1st dimension: size resolve;
	# 2nd dimension: inhomogeneity in size bin.
	# the 1-d inhomogeneity use linear change to represent;
	# the 2-d inhomogeneity use several peaks to represent.
	# to upgrade calPNSD.cal_Mie() function.
	# use nI1 and nI2, which nI1 is size resolved n change slope,
	# nI2 is two peak's separate rate.
	# for x peak and y bin nI2, there is:
	# n[i] = n + n * y * (i-(x-1)/2)
	
	par_new = np.zeros((binNum,I2x))
	for i in range(binNum):
		par_new_i = par + I1 * (i-binNum/2) / (binNum/2)
		for j in range(I2x):
			par_new[i,j] = par_new_i + par_new_i * I2y * (j-(I2x-1)/2)
	return par_new

def cal_Mie(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, nI1, nI2x, nI2y, kappa, kappaI1, kappaI2x, kappaI2y, RH, wl, CT, BCAE):
	'''
	This function is to use Dps, PNSD, RH, CT and BCAE to calculate extinct coefficient, scattering coefficient and asymmetry factor g
	input:
		Dps      : diameter distribution, array, nm
		PNSD     : number distribution, array, dn/dlogDp
		DBCps    : BC particle diameter size, array, nm
		DBC      : BC core diameter size, array, nm
		BCPNSD   : BC particle number concentration, array, cm^-3
		kBC      : BC core imagine part of complex refractive index
		nBC      : BC core real part of complex refractive index
		nShell   : shell complex refractive index
		nI1      : n mixing state change rate by diameter
		nI2x     : n mixing state size bin peak number
		nI2y     : n mixing state size bin peak separate rate
		kappa    : hygroscopicity parameter, float
		kappaI   : hygroscopicity parameter mixing state change rate by diameter
		kappaI2x : hygroscopicity parameter mixing state size bin peak number
		kappaI2y : hygroscopicity parameter mixing state size bin peak separate rate
		RH       : relative humidity, float, percent
		wl       : wave length, nm
		CT       : BC core coating thickness adjust parameter
		BCAE     : BC core absorbing enhancement adjust parameter
	output:
		kext     : extinct coefficient, Mm-1
		ksca     : scattering coefficient, Mm-1
		g        : asymmetry factor
	'''
	m_BC = nBC + kBC * 1j
	#m_shell = nShell
	m_shell = cal_I2(nShell, len(Dps), nI1, nI2x, nI2y)
	kappa_shell = cal_I2(kappa, len(Dps), kappaI1, kappaI2x, kappaI2y)
	'''
	kappa_shell = np.zeros(len(Dps))
	for i in range(len(Dps)):
		kappa_shell[i] = kappa + kappaI * (i-len(Dps)/2) / (len(Dps)/2)
	'''
	n = calPNSD.PNSD2n(Dps, PNSD, 'd')
	BCPNSD = BCPNSD * calPNSD.cal_dlogDp(DBC)
	kext = 0
	ksca = 0
	g = 0
	
	for i in range(len(Dps)):
		BCPNSD_i = calPNSD.DMASP2PNSD_Dp(DBCps, BCPNSD, Dps[i]) # array in shape of DBC
		n_noBC_i = n[i] - sum(BCPNSD_i)
		
		if n_noBC_i > 0:
			for j in range(nI2x):
				for k in range(kappaI2x):
					D_noBC_i, m_shell_noBC_i = kappaKohler.RH2D(Dps[i], 0, m_BC, m_shell[i,j], kappa_shell[i,k], RH, wl)
					MieQ = ps.MieQ(m_shell_noBC_i, wl, D_noBC_i)
					Qsca_noBC_i = MieQ[1]
					ksca += Qsca_noBC_i * 1/4 * np.pi * (D_noBC_i*1e-9)**2 * n_noBC_i*1e6 / nI2x / kappaI2x # separate peak making n separate too
					Qext_noBC_i = MieQ[0]
					kext += Qext_noBC_i * 1/4 * np.pi * (D_noBC_i*1e-9)**2 * n_noBC_i*1e6 / nI2x / kappaI2x
					g_noBC_i = MieQ[3]
					g += g_noBC_i * n_noBC_i / sum(n) / nI2x
		
		if sum(BCPNSD_i):
			for j in range(len(DBC)):
				for k in range(nI2x):
					for l in range(kappaI2x):
						D_BC_ij, m_shell_BC_ij = kappaKohler.RH2D(Dps[i], DBC[j], m_BC, m_shell[i,k], kappa_shell[i,l], RH, wl)
						MieQCoreShell = ps.MieQCoreShell(m_BC, m_shell_BC_ij, wl, DBC[j], D_BC_ij*CT)
						Qsca_BC_ij = MieQCoreShell[1]
						ksca += Qsca_BC_ij * 1/4 * np.pi * (D_BC_ij*1e-9)**2 * BCPNSD_i[j]*1e6 / nI2x / kappaI2x
						Qext_BC_ij = Qsca_BC_ij + MieQCoreShell[2] * BCAE
						kext += Qext_BC_ij * 1/4 * np.pi * (D_BC_ij*1e-9)**2 * BCPNSD_i[j]*1e6 / nI2x / kappaI2x
						g_BC_ij = MieQCoreShell[3]
						g += g_BC_ij * BCPNSD_i[j] / sum(n) / nI2x
	
	kext = kext * 1e6
	ksca = ksca * 1e6
	return kext, ksca, g

def cal_total_AOD(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, nI1, nI2x, nI2y, kappa, kappaI1, kappaI2x, kappaI2y, CT, BCAE, VD, **args):
	'''
	This function is to calculate total atmospere air AOD
	input:
	output:
	'''
	if 'dz' in args:
		dz = args['dz']
	else:
		dz = z / 100
	
	# read atmosphere data from atms.dat
	atms = readAtms.read()
	z = atms['z'] * 1e3 # in m, from top to bottom
	dz = z[:-1] - z[1:] # dz in each two layers, in m
	wh = atms['wh']
	t = atms['t']
	RH = np.zeros(len(wh))
	for i in range(len(RH)):
		RH[i] = calRH.wh2RH(t[i], wh[i])
		if RH[i]<1:
			RH[i] = 1 # minimum RH set to 1
	
	AOD = 0
	for i in range(dz): # every level of atmosphere
		dh = dz * calPNSD.cal_VD(Dps, PNSD, h, VD)
		AOD += cal_Mie(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, nI1, nI2x, nI2y, kappa, kappaI1, kappaI2x, kappaI2y, RH, wl, CT, BCAE)[0] * 1e-6 * dh
	
	return AOD

if __name__ == '__main__':
	print('Hello!')
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
