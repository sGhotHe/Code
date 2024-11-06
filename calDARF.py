####################################################################################
# INTRODUCTION:
# This code is to do some DARF calculating
# Created by Hebs at 24/2/25/15:29
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import sys

import readAtms
import calNI
import calRH
import calPNSD

def DARF(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, nH1, kappa, kappaH1, MS, MSH1, MSH2x, MSH2y, VD, CT, CTH1, BCAE, wl, nn, moma, **args):
	'''
	This function is to use Mie modal to calculate aerosol parameters for SBDART
	input:
		Dps      : diameter distribution, array, nm
		PNSD     : number distribution, array, dn/dlogDp
		DBCps    : BC particle diameter size, array, nm
		DBC      : BC core diameter size, array, nm
		BCPNSD   : BC particle number concentration, array, cm^-3
		kBC      : BC core imagine part of complex refractive index, float
		nBC      : BC core real part of complex refractive index, float
		nShell   : shell complex refractive index, float
		nI       : n mixing state change rate by diameter, float
		nI2x     : n mixing state size bin peak number, int
		nI2y     : n mixing state size bin peak separate rate, flaot
		kappa    : hygroscopicity parameter, float
		kappaI   : hygroscopicity parameter mixing state change rate by diameter, float
		kappaI2x : hygroscopicity parameter mixing state size bin peak number, int
		kappaI2y : hygroscopicity parameter mixing state size bin peak separate rate, float
		MS       : mixing state, float
		VD       : vertical distribution parameter, float
		CT       : BC core coating thickness adjust parameter, float
		BCAE     : BC core absorbing enhancement adjust parameter, float
		wl       : wave length, nm, np.array, float
		nn       : aerosol layer level number, int
		moma     : moment of Legendre moments, int
		**args:
			debug				: debug flag, bool, default False
			angularResolution	: angular resolution of PyMieScatt phase function, degree, defualt 0.5
	output:
		aerosol optical parameters, from top to down
		dtau	: aerosol optical depth, np.array, in shape (len(wl),nn)
		waer	: single scattering albedo, np.array, in shape (len(wl),nn)
		pmom	: legendre moments of asymmetry factor, np.array, in shape (len(wl),nn,moma)
	'''
	if 'debug' in args:
		debug = args['debug']
	else:
		debug = False
	if 'angularResolution' in args:
		angularResolution = args['angularResolution']
	else:
		angularResolution = 0.5
	if 'atms_fn' in args:
		atms_fn = args['atms_fn']
	else:
		atms_fn = 'atms.dat'
	
	# read atms.dat data
	atms = readAtms.read(fn=atms_fn) # from EAR5 reanalysis data
	z = atms['z'] * 1e3 # in m, from top to bottom
	dz = np.zeros(len(z))
	dz[1:] = z[:-1] - z[1:] # dz in each two layers, in m
	dz[0] = dz[1]
	wh = atms['wh']
	t = atms['t']
	if debug:
		print('Max height :', z[0], 'm')
	# turn water vapor density in g/m3 to RH, in %
	RH = np.zeros(len(wh))
	for i in range(len(RH)):
		RH[i] = calRH.wh2RH(t[i], wh[i])
		if RH[i]<1:
			RH[i] = 1 # minimum RH set to 1
		elif RH[i]>=100:
			RH[i] = 99.9 # maximum RH set to 99.9
	if debug:
		print('RH : ', RH)
	
	# check the input data validity
	# if wl(k) < wl(k+1) doesn\'t meet, exit
	for i in range(len(wl)-1):
		if (wl[i+1]-wl[i])<0:
			print('wl(k) < wl(k+1) doesn\'t meet. Please check.')
			sys.exit()
	
	# calculate every level aerosol optical parameter
	
	if debug:
		print('calculating...')
	
	dtau = np.zeros((len(wl),nn))
	waer = np.zeros((len(wl),nn))
	g = np.zeros((len(wl),nn))
	pmom = np.zeros((len(wl),nn,moma))
	
	for i in range(len(wl)):
		for j in range(nn):
			VD_ratio = calPNSD.cal_VD(Dps, PNSD, z[atms['nn']-nn+j], VD)
			kext_ij, waer_ij, pmom_ij, g_ij = calNI.cal_Mie4(Dps, PNSD*VD_ratio, DBCps, DBC, BCPNSD*VD_ratio, kBC, nBC, nShell, nH1, kappa, kappaH1, MS, MSH1, MSH2x, MSH2y, RH[atms['nn']-nn+j], wl[i], moma, CT, CTH1, BCAE, angularResolution=angularResolution)
			#print('RH: ', RH[atms['nn']-nn+j])
			#print('height: ', z[atms['nn']-nn+j])
			#print('thickness: ', dz[atms['nn']-nn+j])
			#print('VD_ratio: ', VD_ratio)
			dtau[i,j] = kext_ij * 1e-6 * dz[atms['nn']-nn+j] # Mm^-1 to m^-1
			waer[i,j] = waer_ij
			g[i,j] = g_ij
			pmom[i,j] = pmom_ij
			if debug:
				print(round((j+i*nn+1)/(len(wl)*nn)*1000)/10, '% done...')
	
	if debug:
		print('done')
	
	return dtau, waer, pmom, g

if __name__ == '__main__':
	print('Hello!')
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
