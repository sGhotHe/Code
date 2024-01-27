####################################################################################
# INTRODUCTION:
# This code is to test f(RH) sensitivity of DARF in SBDART
# Created by Hebs at 21/11/9/14:46
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import kappaKohler
import readTaizhou
import calBC
import calPNSD
import phaseFunc
import readAtms
import calRH
import readFRH
import os
import write_aerosol
import read11

def cal_dtau_RH(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, kappa, RH, wl, H, dH, **args):
	'''
	This function is to calculate aerosol optical depth AOD in certain height, delta height and raletive humidity
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
	
	kext_RH = kappaKohler.RH2kext(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, kappa, RH, wl)
	H_integral = 0
	for i in range(n):
		H_integral += calBC.cal_vertical_distribution_ratio(n_BC, H+n*ddH, func=func) * ddH
	
	dtau = kext_RH * 1e-6 * H_integral
	
	return dtau

def cal_SSA_RH(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, kappa, RH, wl, H, **args):
	'''
	This function is to calculate single scattering albedo SSA in certain height and raletive humidity
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
		**func  : vertical distribution function type, A or B, default B
	output:
		SSA     : single scattering albedo, float
	'''
	if 'func' in args:
		func = args['func']
	else:
		func = 'B'
	
	PNSD_H = PNSD * calBC.cal_vertical_distribution_ratio(n_BC, H, func='A')
	n_BC_H = n_BC * calBC.cal_vertical_distribution_ratio(n_BC, H, func='A')
	kext_RH = kappaKohler.RH2kext(Dps, PNSD_H, DBCps, DBC, n_BC_H, m_BC, m_shell, kappa, RH, wl)
	ksca_RH = kappaKohler.RH2ksca(Dps, PNSD_H, DBCps, DBC, n_BC_H, m_BC, m_shell, kappa, RH, wl)
	SSA = ksca_RH / kext_RH
	
	return SSA

def cal_beta_RH(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, kappa, RH, wl, H, n, **args):
	'''
	This function is to calculate Legendre moments in certain height and raletive humidity
	input:
		Dps     : diameter distribution, array, nm
		PNSD    : number distribution, array, dn/dlogDp
		DBCps   : BC particle diameter size, array, nm
		DBC     : BC core diameter size, array, nm
		n_BC    : BC particle number concentration, array, cm^-3
		m_BC    : BC core complex refractive index
		m_shell : shell complex refractive index
		wl      : wave length, nm
		kappa   : hygroscopicity parameter, float
		RH      : relative humidity, float, percent
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
	P_bulk, theta = phaseFunc.cal_bulk_phase_function(Dps, PNSD*ratio, DBCps, DBC, n_BC*ratio, m_BC, m_shell, kappa, RH, wl, angularResolution=10)
	b = phaseFunc.beta(P_bulk, theta, n)
	
	return b

def write(fn, nn, moma, wl, m_BC, **args):
	'''
	This function is to write aerosol.dat in certain kappa
	input:
		fn               : SP2, SMPS and nephelometer data file name
			for Taizhou.npy:
			Dps    : partical diameter distribution, array, nm
			PNSD   : partical number concentration distribution, array, dn/dlogDp
			DBCps  : partical diameter distribution, array, nm
			DBC    : BC core diameter distribution, array, nm
			n_BC   : BC core number concentration, array, cm^-3
			wl_sca : nephelometer wave length, nm
			ksca   : nephelometer scattering coefficient, Mm^-1
		nn               : number of atmospheric levels for which aerosol information is specified
		momoa            : number of phase function moments
		wl               : the wavelength [ wl(k) < wl(k+1) ], array, nm
		m_BC             : BC core complex refractive index
		kappa            : hygroscopicity parameter, float
		**func           : vertical distribution function type, A or B, default A
		**output         : output path, default ./aerosol.dat
		**m_BC_origin    : BC core origin complex refractive index, to judge if m_BC are changed in rate, default m_BC
		**m_shell_origin : shell complex refractive index, default read from wl_sca and ksca
		**kappa_origin   : hygroscopicity parameter, float, default read from f(RH) data
		**fRH_fn         : f(RH) data file path, string, default data/fRH
		**write_g        : whether to write g information to g.dat, boolean, default False
		**write_mShell   : whether to write m_shell information to mShell.dat, boolean, default False
		**write_kappa    : whether to write kappa information to kappa.dat, boolean, default False
	output:
		aerosol.dat for SBDART
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
	if 'kappa_origin' in args:
		kappa = args['kappa_origin']
		cal_kappa = False
	else:
		cal_kappa = True
	if 'fRH_fn' in args:
		fRH_fn = args['fRH_fn']
	else:
		fRH_fn = 'data/fRH'
	if 'write_g' in args:
		write_g = args['write_g']
	else:
		write_g = False
	if 'write_mShell' in args:
		write_mShell = args['write_mShell']
	else:
		write_mShell = False
	if 'write_kappa' in args:
		write_kappa = args['write_kappa']
	else:
		write_kappa = False
	
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
	dz = z[1:] - z[:-1] # dz in each two layers, in m
	# need to read RH profile
	t = atms['t'][::-1]
	wh = atms['wh'][::-1]
	RH = np.zeros(len(z))
	for i in range(len(RH)):
		RH[i] = calRH.wh2RH(t[i], wh[i])
	
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
	
	# calculate kappa from f(RH)
	if cal_kappa:
		# read data from f(RH) data file path
		data = readFRH.readfRHs(fRH_fn)
		f80_635s = data['f80_635s']
		f80_525s = data['f80_525s']
		f80_450s = data['f80_450s']
		
		kappa = np.zeros((3,len(f80_635s)))
		
		print('kappa calculating...')
		for i in range(len(f80_635s)):
			kappa[0,i] = kappaKohler.fRH2kappa(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, 635, 80, f80_635s[i])
			kappa[1,i] = kappaKohler.fRH2kappa(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, 525, 80, f80_525s[i])
			kappa[2,i] = kappaKohler.fRH2kappa(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, 450, 80, f80_450s[i])
			print(round((i+1)/len(f80_635s)*1000)/10, '% done...')
		
		kappa = kappa.mean()
	
	if write_kappa:
		with open('kappa.dat', 'r') as f:
			f.write(str(kappa)+'\n')
	
	i = 0 # cut down from 20nm
	while Dps[i]<20:
		PNSD[i] = 0
		i += 1
	
	dtau = np.zeros((len(wl), nn))
	waer = np.zeros((len(wl), nn))
	pmom = np.zeros((len(wl), nn, moma))
	if write_g:
		g = np.zeros((len(wl), nn))
	
	print('calculating...')
	
	for i in range(len(wl)):
		for j in range(nn):
			dtau[i,j] = cal_dtau_RH(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, kappa, RH[j], wl[i], z[j], dz[j], func=func)
			waer[i,j] = cal_SSA_RH(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, kappa, RH[j], wl[i], z[j], func=func)
			for k in range(6):
				pmom[i,j,k] = cal_beta_RH(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, kappa, RH[j], wl[i], z[j], k+1, func=func)
			
			# calculate g data
			if write_g:
				P_bulk, theta = phaseFunc.cal_bulk_phase_function_RH(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, kappa, RH[j], wl[i], angularResolution=10)
				g[i,j] = phaseFunc.cal_g(P_bulk, theta)
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
	if write_g:
		with open('g.dat', 'w') as f:
			for i in range(len(wl)):
				f.write(str(wl[i])+'\n')
				for j in range(nn):
					f.write(str(g[i,j])+'\t')
				f.write('\n')
	
	print('done')

def run(fn, nn, moma, wl, m_BC, rate, **args):
	'''
	This function is to change kappa in certain rate
	input:
		fn              : SP2, SMPS and nephelometer data file name
			for Taizhou.npy:
			Dps    : partical diameter distribution, array, nm
			PNSD   : partical number concentration distribution, array, dn/dlogDp
			DBCps  : partical diameter distribution, array, nm
			DBC    : BC core diameter distribution, array, nm
			n_BC   : BC core number concentration, array, cm^-3
			wl_sca : nephelometer wave length, nm
			ksca   : nephelometer scattering coefficient, Mm^-1
		nn              : number of atmospheric levels for which aerosol information is specified
		momoa           : number of phase function moments
		wl              : the wavelength [ wl(k) < wl(k+1) ], array, nm
		m_BC            : BC core complex refractive index
		rate            : change rate, array
		**func          : vertical distribution function type, A or B, default A
		**output        : output path, default input/kappa/aerosol.dat_[rate]
		**m_shell_origin: shell complex refractive index, default read from wl_sca and ksca
		**kappa_origin  : hygroscopicity parameter, float, default read from f(RH) data
		**fRH_fn        : f(RH) data file path, string, default data/fRH
		**write_kappa   : whether to write kappa information to kappa.dat, boolean, default False
	output:
		output/kappa/iout11_[rate].txt
	'''
	if 'func' in args:
		func = args['func']
	else:
		func = 'A'
	if 'output' in args:
		output = args['output']
	else:
		output = 'input/kappa/aerosol.dat'
	if 'm_shell_origin' in args:
		m_shell = args['m_shell_origin']
		cal_mshell = False
	else:
		cal_mshell = True
	if 'kappa_origin' in args:
		kappa = args['kappa_origin']
		cal_kappa = False
	else:
		cal_kappa = True
	if 'fRH_fn' in args:
		fRH_fn = args['fRH_fn']
	else:
		fRH_fn = 'data/fRH'
	if 'write_kappa' in args:
		write_kappa = args['write_kappa']
	else:
		write_kappa = False
	
	if cal_mshell:
		wl_sca = data['wl_sca']
		ksca = np.nanmean(data['ksca'],axis=0)
		m_shell = np.zeros(len(wl_sca))
		for i in range(len(wl_sca)):
			m_shell[i] = calBC.ksca2mshell(Dps, PNSD, DBCps, DBC, n_BC, m_BC_origin, ksca[i], wl_sca[i])
		m_shell = np.mean(m_shell)
	
	if cal_kappa:
		data = readFRH.readfRHs(fRH_fn)
		f80_635s = data['f80_635s']
		f80_525s = data['f80_525s']
		f80_450s = data['f80_450s']
	
		kappa = np.zeros((3,len(f80_635s)))
		
		print('kappa calculating...')
		for i in range(len(f80_635s)):
			kappa[0,i] = kappaKohler.fRH2kappa(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, 635, 80, f80_635s[i])
			kappa[1,i] = kappaKohler.fRH2kappa(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, 525, 80, f80_525s[i])
			kappa[2,i] = kappaKohler.fRH2kappa(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, 450, 80, f80_450s[i])
			print(round((i+1)/len(f80_635s)*1000)/10, '% done...')
		
		kappa = kappa.mean()
	
	if write_kappa:
		with open('kappa.dat', 'r') as f:
			f.write(str(kappa)+'\n')
	
	for i in range(len(rate)):
		kappa_new = kappa * rate[i]
		write(fn, nn, moma, wl, m_BC, func=func, output=output, m_shell_origin=m_shell, kappa_origin=kappa_new)
		os.system('mv '+output+' '+output+'_'+str(round(rate[i]*100)/100))
	
	os.system('mv aerosol.dat aerosol.dat_origin') # backup origin file
	
	for i in range(len(rate)):
		os.system('cp '+output+'_'+str(round(rate[i]*100)/100)+' aerosol.dat')
		os.system('sbdart >output/kappa/iout11_'+str(round(rate[i]*100)/100)+'.txt')
	
	write_aerosol.change_back()

def read(rate):
	'''
	This function is to read sequence kappa change brought net radiative flux change at TOA
	input:
		rate    : kappa change rate, array, float
	output:
		nrf     : net radiative flux at TOA for each m_shell change rate, array, float
	'''
	nrf = np.zeros(len(rate)) # net radiative flux at TOA
	for i in range(len(rate)):
		filename = 'output/kappa/iout11_'+str(round(rate[i]*100)/100)+'.txt'
		iout11 = read11.read11(filename=filename)
		nrf[i] = iout11['fxdn'][0]-iout11['fxup'][0]
	return nrf

def read_parameter(rate):
	'''
	This function is to read sequence changed kappa
	input:
		rate    : kappa change rate, array, float
	output:
		kappa   : changed hygroscopicity parameter kappa, array in shape(len(rate))
	'''
	with open('kappa.dat', 'r') as f:
		kappa_origin = float(f.readline()[:-1])
	
	kappa = np.zeros(len(rate))
	for i in range(len(rate)):
		kappa[i] = kappa_origin * rate[i]
	
	return kappa

if __name__ == '__main__':
	# n change from -10% to 10%, bin in 1%
	rate = np.arange(-10,11) / 100 + 1
	'''
	write('data/sp2/Taizhou.npy', 50, 6, [440,500,870,1640], 1.8+0.54j, func='A', write_g=True, write_mShell=True, m_shell_origin=1.5833333333333333, kappa_origin=0.21537370768611558)
	run('data/sp2/Taizhou.npy', 50, 6, [440,500,870,1640], 1.8+0.54j, rate, func='A', m_shell_origin=1.5833333333333333, kappa_origin=0.21537370768611558)
	nrf = read(rate)
	print(nrf)
	write('data/sp2/Taizhou.npy', 50, 6, [440,500,870,1640], 1.8+0.54j, func='A', output='data/test/aerosol.dat_1e-9', m_shell_origin=1.5833333333333333, kappa_origin=1e-9)
	write_aerosol.write('data/sp2/Taizhou.npy', 50, 6, [440,500,870,1640], 1.8+0.54j, func='A', output='data/test/aerosol.dat_0', m_shell_origin=1.5833333333333333)
	'''
	write('data/sp2/Taizhou.npy', 50, 6, [440,500,870,1640], 1.8+0.54j, func='A', m_shell_origin=1.5833333333333333, kappa_origin=0.21537370768611558, output='output/location/zero/aerosol.dat')
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
