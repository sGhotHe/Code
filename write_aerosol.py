####################################################################################
# INTRODUCTION:
# This code is to make aerosol.dat, including AOD, SSA and phase function, using aerosol size distribution file and parameterized vertical distribution
# Created by Hebs at 21/10/21/10:56
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import readSMPS
import readTaizhou
import phaseFunc
import readAtms
import calPNSD
import calBC
import calRH
import calMie
import sys
import re
import os
import datetime

def write_old(fn, nn, moma, wl, m, **args):
	'''
	This function is to use PNSD data to write aerosol.dat for SBDART
	input:
		fn          : PNSD data file name, smps.txt
		nn          : number of atmospheric levels for which aerosol information is specified
		momoa       : number of phase function moments
		wl          : the wavelength [ wl(k) < wl(k+1) ], array, nm
		m           : complex refractive index
		**func      : vertical distribution function type, A or B, default B
		**output    : output path, default ./aerosol.dat
	output:
		aerosol.dat : aerosol information for SBDART
		where 
                    nn          is the number of atmospheric levels for
                                which aerosol information is specified.

                    moma        number of phase function moments

                    wl(k)       is the wavelength [ wl(k) < wl(k+1) ]

                    dtau(i,k)   is the optical depth increment within
                                level i at wavelength k, information is
                                specified in top-down order. 

                    waer(i,k)   is the single scattering albedo

                    pmom(m,i,k) are legendre moments of the phase function.
                                Note that zeroeth moment is not read, it
                                is assumed to be 1.
	'''
	if 'func' in args:
		func = args['func']
	else:
		func = 'B'
	if 'output' in args:
		output = args['output']
	else:
		output = 'aerosol.dat'
	
	for i in range(len(wl)-1):
		if (wl[i+1]-wl[i])<0:
			print('wl(k) < wl(k+1) doesn\'t meet. Please check.')
			sys.exit()
	
	smps = readSMPS.readSMPS(fn)
	Dps = smps['Dps'][0]
	PNSD = np.mean(smps['PNSD'], 0) # average PNSD data
	
	atms = readAtms.read()
	z = atms['z'][::-1] * 1e3 # reverse, down to top, in m
	if nn>atms['nn']: # aerosol layers bigger than atmospheric layers
		print('Too many aerosol layers. Please check.')
		sys.exit()
	
	dz = z[1:] - z[:-1] # dz in each two layers, in m
	
	i = 0 # cut down from 20nm
	while Dps[i]<20:
		PNSD[i] = 0
		i += 1
	
	dtau = np.zeros((len(wl), nn))
	waer = np.zeros((len(wl), nn))
	pmom = np.zeros((len(wl), nn, moma))
	
	print('calculating...')
	
	for i in range(len(wl)):
		for j in range(nn):
			dtau[i,j] = calPNSD.cal_dtau(Dps, PNSD, wl[i], z[j], dz[j], m, func=func)
			waer[i,j] = calPNSD.cal_SSA(Dps, PNSD, wl[i], z[j], m, func=func)
			for k in range(6):
				pmom[i,j,k] = phaseFunc.cal_beta(Dps, PNSD, wl[i], z[j], m, k+1, func='A')
			print(round((j+i*nn+1)/(len(wl)*nn)*1000)/10, '% done...')
	
	print('writing...')
	
	with open(output, 'w') as f:
		f.write(str(nn)+'\t'+str(moma)+'\n')
		for i in range(len(wl)):
			f.write(str(wl[i])+'\n')
			for j in range(nn):
				f.write(str(dtau[i,j])+'\t'+str(waer[i,j])+'\t')
				for k in range(6):
					f.write(str(pmom[i,j,k])+'\t')
				f.write('\n')
	
	print('done')

def write_old(fn, nn, moma, wl, m_BC, **args):
	'''
	This function is to use DMA-SP2, SMPS and nephelometer data to write aerosol.dat for SBDART
	input:
		fn               : SP2, SMPS and nephelometer data file name
			for Taizhou.npy:
			Dps      : partical diameter distribution, array, nm
			PNSD     : partical number concentration distribution, array, dn/dlogDp
			DBCps    : partical diameter distribution, array, nm
			DBC      : BC core diameter distribution, array, nm
			n_BC     : BC core number concentration, array, cm^-3
			wl_sca   : nephelometer wave length, nm
			ksca     : nephelometer scattering coefficient, Mm^-1
		nn               : number of atmospheric levels for which aerosol information is specified
		momoa            : number of phase function moments
		wl               : the wavelength [ wl(k) < wl(k+1) ], array, nm
		m_BC             : BC core complex refractive index
		**func           : vertical distribution function type, A or B, default A
		**output         : output path, default ./aerosol.dat
		**m_BC_origin    : BC core origin complex refractive index, to judge if m_BC are changed in rate, default m_BC
		**m_shell_origin : shell complex refractive index, default read from wl_sca and ksca
		**write_g        : whether to write g information to g.dat, boolean, default False
		**write_mShell   : whether to write m_shell information to mShell.dat, boolean, default False
	output:
		aerosol.dat      : aerosol information for SBDART
		where 
                    nn          is the number of atmospheric levels for
                                which aerosol information is specified.

                    moma        number of phase function moments

                    wl(k)       is the wavelength [ wl(k) < wl(k+1) ]

                    dtau(i,k)   is the optical depth increment within
                                level i at wavelength k, information is
                                specified in top-down order. 

                    waer(i,k)   is the single scattering albedo

                    pmom(m,i,k) are legendre moments of the phase function.
                                Note that zeroeth moment is not read, it
                                is assumed to be 1.
	'''
	if 'func' in args:
		func = args['func']
	else:
		func = 'A'
	if 'output' in args:
		output = args['output']
	else:
		output = 'aerosol.dat'
	if 'm_BC_origin' in args:
		m_BC_origin = args['m_BC_origin']
	else:
		m_BC_origin = m_BC
	if 'm_shell_origin' in args:
		m_shell = args['m_shell_origin']
		cal_mshell = False
	else:
		cal_mshell = True
	if 'write_g' in args:
		write_g = args['write_g']
	else:
		write_g = False
	if 'write_mShell' in args:
		write_mShell = args['write_mShell']
	else:
		write_mShell = False
	
	# read data from SP2, SMPS and Nephelometer data
	data = readTaizhou.read_Taizhou(fn)
	Dps = data['Dps']
	PNSD = np.nanmean(data['PNSD'],axis=0)
	DBCps = data['DBCps']
	DBC = data['DBC']
	n_BC = np.nanmean(data['n_BC'],axis=0)
	
	for i in range(len(wl)-1):
		if (wl[i+1]-wl[i])<0:
			print('wl(k) < wl(k+1) doesn\'t meet. Please check.')
			sys.exit()
	
	atms = readAtms.read()
	z = atms['z'][::-1] * 1e3 # reverse, down to top, in m
	if nn>atms['nn']: # aerosol layers bigger than atmospheric layers
		print('Too many aerosol layers. Please check.')
		sys.exit()
	
	# calculate m_shell from ksca
	if cal_mshell:
		wl_sca = data['wl_sca']
		ksca = np.nanmean(data['ksca'],axis=0)
		m_shell = np.zeros(len(wl_sca))
		for i in range(len(wl_sca)):
			m_shell[i] = calBC.ksca2mshell(Dps, PNSD, DBCps, DBC, n_BC, m_BC_origin, ksca[i], wl_sca[i])
		m_shell = np.mean(m_shell)
	
	# write m_shell information to mShell.dat
	if write_mShell:
		with open('mShell.dat', 'w') as f:
			f.write(str(m_shell)+'\n')
	
	dz = z[1:] - z[:-1] # dz in each two layers, in m
	
	i = 0 # cut down from 20nm
	while Dps[i]<20:
		PNSD[i] = 0
		i += 1
	
	dtau = np.zeros((len(wl), nn))
	waer = np.zeros((len(wl), nn))
	pmom = np.zeros((len(wl), nn, moma))
	if write_g:
		g = np.zeros(len(wl))
	
	print('calculating...')
	
	for i in range(len(wl)):
		for j in range(nn):
			dtau[i,j] = calBC.cal_dtau(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, wl[i], z[j], dz[j], func='A')
			waer[i,j] = calBC.cal_SSA(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, wl[i], z[j], func='A')
			for k in range(6):
				pmom[i,j,k] = phaseFunc.cal_beta(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, wl[i], z[j], k+1, func='A')
			print(round((j+i*nn+1)/(len(wl)*nn)*1000)/10, '% done...')
		
		# calculate g data
		if write_g:
			P_bulk, theta = phaseFunc.cal_bulk_phase_function(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, wl[i], angularResolution=10)
			g[i] = phaseFunc.cal_g(P_bulk, theta)
	
	print('writing...')
	
	with open(output, 'w') as f:
		f.write(str(nn)+'\t'+str(moma)+'\n')
		for i in range(len(wl)):
			f.write(str(wl[i])+'\n')
			for j in range(nn):
				f.write(str(dtau[i,j])+'\t'+str(waer[i,j])+'\t')
				for k in range(6):
					f.write(str(pmom[i,j,k])+'\t')
				f.write('\n')
	if write_g:
		with open('g.dat', 'w') as f:
			for i in range(len(wl)):
				f.write(str(wl[i])+'\t'+str(g[i])+'\n')
	
	print('done')

def write(nn, moma, Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, kappa, CT, BCAE, wl, VD, **args):
	'''
	This function is to use PNSD data to write aerosol.dat for SBDART
	input:
		nn           : number of atmospheric levels for which aerosol information is specified
		momoa        : number of phase function moments
		Dps          : diameter distribution, array, nm
		PNSD         : number distribution, array, dn/dlogDp
		DBCps        : BC particle diameter size, array, nm
		DBC          : BC core diameter size, array, nm
		BCPNSD       : BC particle number size distribution, array, dn/dlogDBC
		nBC          : BC core real part of complex refractive index
		kBC          : BC core imagine part of complex refractive index
		nShell       : shell complex refractive index
		kappa        : hygroscopicity parameter, float
		CT           : BC coating thickness adjustment, float
		BCAE         : BC absorbing enhancement adjustment, float
		wl           : the wavelength [ wl(k) < wl(k+1) ], array, nm
		VD           : vertical distribution mixing, 0 for type A and 1 for type B
		**output     : aerosol.dat storage path, defualt 'aerosol.dat', string
		**debug      : debug flag, default False, bool
		**ddz        : AOD calculation step, default dz/10, float
		
	output:
		aerosol.dat      : aerosol information for SBDART
		where 
                    nn          is the number of atmospheric levels for
                                which aerosol information is specified.

                    moma        number of phase function moments

                    wl(k)       is the wavelength [ wl(k) < wl(k+1) ]

                    dtau(i,k)   is the optical depth increment within
                                level i at wavelength k, information is
                                specified in top-down order. 

                    waer(i,k)   is the single scattering albedo

                    pmom(m,i,k) are legendre moments of the phase function.
                                Note that zeroeth moment is not read, it
                                is assumed to be 1.
	'''
	if 'output' in args:
		output = args['output']
	else:
		output = 'aerosol.dat'
	if 'debug' in args:
		debug = args['debug']
	else:
		debug = False
	if 'ddz' in args:
		ddz = args['ddz']
	else:
		ddz = 10
	
	print('reading...')
	
	# read atmospere data from atms.dat
	atms = readAtms.read() # from EAR5 reanalysis data
	z = atms['z'] * 1e3 # in m, from top to bottom
	dz = z[:-1] - z[1:] # dz in each two layers, in m
	dz = np.insert(dz, 0, dz[0])
	wh = atms['wh']
	t = atms['t']
	
	if nn>atms['nn']: # aerosol layers bigger than atmospheric layers
		print('Too many aerosol layers. Please check.')
		sys.exit()
	
	# turn water vapor density in g/m3 to RH, in %
	RH = np.zeros(len(wh))
	for i in range(len(RH)):
		RH[i] = calRH.wh2RH(t[i], wh[i])
		if RH[i]<1:
			RH[i] = 1 # minimum RH set to 1
		if RH[i]>98:
			RH[i] = 98
	
	print('done\ncalculating...')
	
	dtau = np.zeros((len(wl), nn))
	waer = np.zeros((len(wl), nn))
	pmom = np.zeros((len(wl), nn, moma))
	
	for i in range(len(wl)):
		for j in range(nn):
			dtau[i,j] = calMie.cal_AOD(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, kappa, RH[atms['nn']-nn+j], wl[i], z[atms['nn']-nn+j], dz[atms['nn']-nn+j], CT, BCAE, VD, ddz=ddz)
			waer[i,j] = calMie.cal_SSA(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, kappa, RH[atms['nn']-nn+j], wl[i], CT, BCAE)
			for k in range(moma):
				pmom[i,j,k] = phaseFunc.cal_beta(Dps, PNSD, DBCps, DBC, BCPNSD, nBC+kBC*1j, nShell, kappa, RH[atms['nn']-nn+j], wl[i], z[atms['nn']-nn+j], VD, k+1, angularResolution=10)
			print(round((j+i*nn+1)/(len(wl)*nn)*1000)/10, '% done...')
	
	print('done\nwriting...')
	
	with open(output, 'w') as f:
		f.write(str(nn)+'\t'+str(moma)+'\n')
		for i in range(len(wl)):
			f.write(str(wl[i])+'\n')
			for j in range(nn):
				f.write(str(dtau[i,j])+'\t'+str(waer[i,j])+'\t')
				for k in range(moma):
					f.write(str(pmom[i,j,k])+'\t')
				f.write('\n')

def change_AOD(rate, **args):
	'''
	This function is to change aerosol.dat AOD in certain rate
	input:
		rate  : change rate, float
		**fn  : file path for albedo.dat, default ./aerosol.dat
	output:
		aerosol.dat for SBDART
	'''
	if 'fn' in args:
		fn = args['fn']
	else:
		fn = 'aerosol.dat'
	
	wl = []
	dtau = []
	waer = []
	pmom = []
	
	with open(fn, 'r') as f:
		line = f.readline()
		res = re.split('\t', line[:-1])
		nn = int(res[0])
		moma = int(res[1])
		data = f.readlines()
		wl_num = round(len(data)/(nn+1))
		data = np.array(data)
		data = data.reshape(wl_num, -1)
		for i in range(wl_num):
			wl.append(data[i,0][:-1])
			for j in range(nn):
				res = re.split('\t', data[i,j+1][:-2]) # delete \t and \n
				dtau.append(res[0])
				waer.append(res[1])
				pmom.append(res[2:])
	
	dtau = np.array(dtau, dtype=float)
	dtau = dtau.reshape(wl_num, -1)
	dtau = dtau * rate
	waer = np.array(waer, dtype=float)
	waer = waer.reshape(wl_num, -1)
	pmom = np.array(pmom, dtype=float)
	pmom = pmom.reshape(wl_num, nn, 6)
	
	if not os.path.exists(fn+'_origin'):
		os.system('mv '+fn+' '+fn+'_origin') # backup the origin file
	
	with open(fn, 'w') as f:
		f.write(str(nn)+'\t'+str(moma)+'\n')
		for i in range(len(wl)):
			f.write(str(wl[i])+'\n')
			for j in range(nn):
				f.write(str(dtau[i,j])+'\t'+str(waer[i,j])+'\t')
				for k in range(6):
					f.write(str(pmom[i,j,k])+'\t')
				f.write('\n')

def change_SSA(rate, **args):
	'''
	This function is to change aerosol.dat SSA in certain rate
	input:
		rate  : change rate, float
		**fn  : file path for albedo.dat, default ./aerosol.dat
	output:
		aerosol.dat for SBDART
	'''
	if 'fn' in args:
		fn = args['fn']
	else:
		fn = 'aerosol.dat'
	
	wl = []
	dtau = []
	waer = []
	pmom = []
	
	with open(fn, 'r') as f:
		line = f.readline()
		res = re.split('\t', line[:-1])
		nn = int(res[0])
		moma = int(res[1])
		data = f.readlines()
		wl_num = round(len(data)/(nn+1))
		data = np.array(data)
		data = data.reshape(wl_num, -1)
		for i in range(wl_num):
			wl.append(data[i,0][:-1])
			for j in range(nn):
				res = re.split('\t', data[i,j+1][:-2]) # delete \t and \n
				dtau.append(res[0])
				waer.append(res[1])
				pmom.append(res[2:])
	
	dtau = np.array(dtau, dtype=float)
	dtau = dtau.reshape(wl_num, -1)
	waer = np.array(waer, dtype=float)
	waer = waer.reshape(wl_num, -1)
	waer = waer * rate
	pmom = np.array(pmom, dtype=float)
	pmom = pmom.reshape(wl_num, nn, 6)
	
	if not os.path.exists(fn+'_origin'):
		os.system('mv '+fn+' '+fn+'_origin')
	
	with open(fn, 'w') as f:
		f.write(str(nn)+'\t'+str(moma)+'\n')
		for i in range(len(wl)):
			f.write(str(wl[i])+'\n')
			for j in range(nn):
				f.write(str(dtau[i,j])+'\t'+str(waer[i,j])+'\t')
				for k in range(6):
					f.write(str(pmom[i,j,k])+'\t')
				f.write('\n')

def change_back(**args):
	'''
	This function is to change albedo.dat back
	input:
		**fn  : file path for aerosol.dat, default ./aerosol.dat
	output:
		aerosol.dat for SBDART
	'''
	if 'fn' in args:
		fn = args['fn']
	else:
		fn = 'aerosol.dat'
	
	if os.path.exists(fn+'_origin'):
		if os.path.exists(fn):
			os.system('rm '+fn)
		os.system('mv '+fn+'_origin '+fn)
	else:
		print('No origin file aerosol.dat_origin. Please check.')

if __name__ == '__main__':
	sp2 = readTaizhou.read_Taizhou('data/sp2/Taizhou.npy')
	DBC = sp2['DBC']
	DBCps = sp2['DBCps']
	BCPNSD = np.zeros(sp2['DMASP2'][0].shape) # no BC core
	
	Dps = []
	PNSD = []
	
	with open('data/AOT/PNSD.txt', 'r') as f:
		infos = f.readlines()
		for line in infos:
			res = re.split('\t', line[:-1])
			Dps.append(float(res[0]))
			PNSD.append(float(res[1]))
	
	Dps = np.array(Dps)
	PNSD = np.array(PNSD)
	
	write(50, 6, Dps, PNSD, DBCps, DBC, BCPNSD, 1, 1, 1.4884+1e-3j, 1.24, 1, 1, [440,500,870,1640], 28/(28+39))
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
