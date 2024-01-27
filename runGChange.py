####################################################################################
# INTRODUCTION:
# This code is to test asymmetry factor g sensitivity of DARF in SBDART
# Created by Hebs at 21/11/5/16:47
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import write_aerosol
import phaseFunc
import read11
import calBC
import readTaizhou
import readAtms
import os
import re
import runFRHChange
import calRH

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
	
	g_old = phaseFunc.cal_g(P, theta) * rate # target
	
	af = np.arange(af_min, af_max+af_bin_rough, af_bin_rough)
	ii = 0 # to mark the closest value subscript
	dif = 999999
	
	for i in range(len(af)):
		P_new = g_change(P, theta, rate, af=af[i])
		g_new = phaseFunc.cal_g(P_new, theta)
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
		g_new = phaseFunc.cal_g(P_new, theta)
		dg = abs(g_new-g_old)
		if dg<dif:
			dif = dg
			ii = i
	
	return af[ii]

def cal_beta_change(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, wl, H, n, rate, **args):
	'''
	This function is to use size distribution data and change rate to calculate Legendre moments change
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
	P_bulk, theta = phaseFunc.cal_bulk_phase_function(Dps, PNSD*ratio, DBCps, DBC, n_BC*ratio, m_BC, m_shell, wl, angularResolution=10)
	g = phaseFunc.cal_g(P_bulk,theta)
	af = cal_g_change_af(P_bulk, theta, rate)
	P_bulk_new = g_change(P_bulk, theta, rate, af=af)
	b = phaseFunc.beta(P_bulk_new, theta, n)
	
	return b

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
	P_bulk, theta = phaseFunc.cal_bulk_phase_function_RH(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, kappa, RH, wl, angularResolution=10)
	g = phaseFunc.cal_g(P_bulk,theta)
	af = cal_g_change_af(P_bulk, theta, rate)
	P_bulk_new = g_change(P_bulk, theta, rate, af=af)
	b = phaseFunc.beta(P_bulk_new, theta, n)
	
	return b

def write_old(fn, nn, moma, wl, m_BC, rate, **args):
	'''
	This function is to write aerosol.dat asymmetry factor g in certain rate
	input:
		fn          : SP2, SMPS and nephelometer data file name
			for Taizhou.npy:
			Dps    : partical diameter distribution, array, nm
			PNSD   : partical number concentration distribution, array, dn/dlogDp
			DBCps  : partical diameter distribution, array, nm
			DBC    : BC core diameter distribution, array, nm
			n_BC   : BC core number concentration, array, cm^-3
			wl_sca : nephelometer wave length, nm
			ksca   : nephelometer scattering coefficient, Mm^-1
		nn          : number of atmospheric levels for which aerosol information is specified
		momoa       : number of phase function moments
		wl          : the wavelength [ wl(k) < wl(k+1) ], array, nm
		m_BC        : BC core complex refractive index
		rate        : asymmetry factor change rate, float
		**func      : vertical distribution function type, A or B, default A
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
		func = 'A'
	if 'output' in args:
		output = args['output']
	else:
		output = 'aerosol.dat'
	
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
	
	wl_sca = data['wl_sca']
	ksca = np.nanmean(data['ksca'],axis=0)
	m_shell = np.zeros(len(wl_sca))
	for i in range(len(wl_sca)):
		m_shell[i] = calBC.ksca2mshell(Dps, PNSD, DBCps, DBC, n_BC, m_BC, ksca[i], wl_sca[i])
	m_shell = np.mean(m_shell)
	
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
			dtau[i,j] = calBC.cal_dtau(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, wl[i], z[j], dz[j], func='A')
			waer[i,j] = calBC.cal_SSA(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, wl[i], z[j], func='A')
			for k in range(6):
				pmom[i,j,k] = cal_beta_change(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, wl[i], z[j], k+1, rate, func='A')
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

def write(fn, nn, moma, wl, m_BC, rate, **args):
	'''
	This function is to write aerosol.dat in certain g change rate
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
		rate             : asymmetry factor change rate, float
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
		g = np.zeros((len(wl), len(RH)))
	
	print('calculating...')
	
	for i in range(len(wl)):
		for j in range(nn):
			dtau[i,j] = runFRHChange.cal_dtau_RH(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, kappa, RH[j], wl[i], z[j], dz[j], func='A')
			waer[i,j] = runFRHChange.cal_SSA_RH(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, kappa, RH[j], wl[i], z[j], func='A')
			for k in range(6):
				pmom[i,j,k] = cal_beta_change_RH(Dps, PNSD, DBCps, DBC, n_BC, m_BC, m_shell, kappa, RH[j], wl[i], z[j], k+1, rate, func='A')
			
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
				for j in range(len(RH)):
					f.write(str(g[i,j])+'\t')
				f.write('\n')
	
	print('done')

def run(fn, nn, moma, wl, m_BC, rate, **args):
	'''
	This function is to change aerosol.dat asymmetry factor g in certain rate
	input:
		fn              : data file name, string
		nn              : number of atmospheric levels for which aerosol information is specified
		momoa           : number of phase function moments
		wl              : the wavelength [ wl(k) < wl(k+1) ], array, nm
		m_BC            : BC complex refractive index
		rate            : change rate, array
		**func          : vertical distribution function type, A or B, default B
		**output        : output file path for albedo.dat, default input/g/aerosol.dat
		**m_shell_origin: shell complex refractive index, default read from wl_sca and ksca
		**kappa_origin  : hygroscopicity parameter, float, default read from f(RH) data
		**fRH_fn        : f(RH) data file path, string, default data/fRH
		**write_g       : whether to write g information to g.dat, boolean, default False
	output:
		output/g/iout11_[rate].txt
	'''
	if 'func' in args:
		func = args['func']
	else:
		func = 'B'
	if 'output' in args:
		output = args['output']
	else:
		output = 'input/g/aerosol.dat'
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
	
	for i in range(len(rate)):
		write(fn, nn, moma, wl, m_BC, rate=rate[i], func=func, output=output, m_BC_origin=m_BC_origin, m_shell_origin=m_shell, kappa_origin=kappa, write_g=True)
		os.system('mv '+output+' '+output+'_'+str(round(rate[i]*100)/100))
	
	os.system('mv aerosol.dat aerosol.dat_origin') # backup origin file
	
	for i in range(len(rate)):
		os.system('cp '+output+'_'+str(round(rate[i]*100)/100)+' aerosol.dat')
		os.system('sbdart >output/g/iout11_'+str(round(rate[i]*100)/100)+'.txt')
	
	write_aerosol.change_back()

def read(rate):
	'''
	This function is to read sequence g change brought net radiative flux change at TOA
	input:
		rate    : g change rate, array, float
	output:
		nrf     : net radiative flux at TOA for each n change rate, array, float
	'''
	nrf = np.zeros(len(rate)) # net radiative flux at TOA
	for i in range(len(rate)):
		filename = 'output/g/iout11_'+str(round(rate[i]*100)/100)+'.txt'
		iout11 = read11.read11(filename=filename)
		nrf[i] = iout11['fxdn'][0]-iout11['fxup'][0]
	return nrf

def read_parameter(rate):
	'''
	This function is to read sequence changed asymmetry factor g
	input:
		rate    : g change rate, array, float
	output:
		wl      : wave length, array, nm
		g       : changed g, array in shape(len(wl), len(rate), len(nn))
	'''
	with open('g.dat', 'r') as f:
		data = f.readlines()
	
	wl = np.zeros(round(len(data)/2))
	g = []
	
	for i in range(len(wl)):
		wl[i] = float(data[2*i][:-1])
		g_i_origin = np.array(re.split('\t', data[2*i+1][:-2]), dtype=float)
		g_i = []
		for j in range(len(rate)):
			g_ij = g_i_origin * rate[j]
			g_i.append(g_ij)
		g.append(g_i)
	
	g = np.array(g, dtype=float)
	
	return wl, g

if __name__ == '__main__':
	# n change from -10% to 10%, bin in 1%
	rate = np.arange(-10,11) / 100 + 1
	'''
	run('data/sp2/Taizhou.npy', 50, 6, [440,500,870,1640], 1.8+0.54j, rate, func='A', m_shell_origin=1.5833333333333333, kappa_origin=0.21537370768611558)
	nrf = read(rate)
	print(nrf)
	'''
	wl, g = read_parameter(rate)
	print(wl)
	print(g[0,0])
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
