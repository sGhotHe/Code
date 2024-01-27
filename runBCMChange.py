####################################################################################
# INTRODUCTION:
# This code is to test BC complex rafractive index sensitivity of DARF in SBDART
# Created by Hebs at 21/11/4/20:39
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import write_aerosol
import read11
import readTaizhou
import calBC
import os
import runFRHChange

def run_n(fn, nn, moma, wl, m_BC, rate, **args):
	'''
	This function is to change aerosol.dat BC complex refractive index real part n in certain rate
	input:
		fn              : data file name, string
		nn              : number of atmospheric levels for which aerosol information is specified
		momoa           : number of phase function moments
		wl              : the wavelength [ wl(k) < wl(k+1) ], array, nm
		m_BC            : BC complex refractive index
		rate            : change rate, array, float
		**func          : vertical distribution function type, A or B, default B
		**output        : output file path for albedo.dat, default input/n/aerosol.dat
		**m_BC_origin   : BC core origin complex refractive index, to judge if m_BC are changed in rate, default m_BC
		**m_shell_origin: shell complex refractive index, default read from wl_sca and ksca
		**kappa_origin  : hygroscopicity parameter, float, default read from f(RH) data
		**fRH_fn        : f(RH) data file path, string, default data/fRH
	output:
		output/n/iout11_[rate].txt
	'''
	if 'func' in args:
		func = args['func']
	else:
		func = 'B'
	if 'output' in args:
		output = args['output']
	else:
		output = 'input/n/aerosol.dat'
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
	
	for i in range(len(rate)):
		m_BC_new = m_BC.real * rate[i] + m_BC.imag * 1j
		runFRHChange.write(fn, nn, moma, wl, m_BC_new, func='A', output=output, m_BC_origin=m_BC_origin, m_shell_origin=m_shell, kappa_origin=kappa)
		os.system('mv '+output+' '+output+'_'+str(round(rate[i]*100)/100))
	
	os.system('mv aerosol.dat aerosol.dat_origin') # backup origin file
	
	for i in range(len(rate)):
		os.system('cp '+output+'_'+str(round(rate[i]*100)/100)+' aerosol.dat')
		os.system('sbdart >output/n/iout11_'+str(round(rate[i]*100)/100)+'.txt')
	
	write_aerosol.change_back()

def read_n(rate):
	'''
	This function is to read sequence n change brought net radiative flux change at TOA
	input:
		rate    : n change rate, array, float
	output:
		nrf     : net radiative flux at TOA for each n change rate, array, float
	'''
	nrf = np.zeros(len(rate)) # net radiative flux at TOA
	for i in range(len(rate)):
		filename = 'output/n/iout11_'+str(round(rate[i]*100)/100)+'.txt'
		iout11 = read11.read11(filename=filename)
		nrf[i] = iout11['fxdn'][0]-iout11['fxup'][0]
	return nrf

def read_n_parameter(rate):
	'''
	This function is to read sequence changed BC core complex refractive index real part n
	input:
		rate    : m_BC n change rate, array, float
	output:
		n       : changed m_BC n, array in shape(len(rate))
	'''
	m_BC = 1.67+0.67j
	n = np.zeros(len(rate))
	for i in range(len(rate)):
		n[i] = m_BC.real * rate[i]
	
	return n

def run_k(fn, nn, moma, wl, m_BC, rate, **args):
	'''
	This function is to change aerosol.dat BC complex refractive index imagine part k in certain rate
	input:
		fn              : data file name, string
		nn              : number of atmospheric levels for which aerosol information is specified
		momoa           : number of phase function moments
		wl              : the wavelength [ wl(k) < wl(k+1) ], array, nm
		m_BC            : BC complex refractive index
		rate            : change rate, array, float
		**func          : vertical distribution function type, A or B, default B
		**output        : output file path for albedo.dat, default input/k/aerosol.dat
		**m_BC_origin   : BC core origin complex refractive index, to judge if m_BC are changed in rate, default m_BC
		**m_shell_origin: shell complex refractive index, default read from wl_sca and ksca
		**kappa_origin  : hygroscopicity parameter, float, default read from f(RH) data
		**fRH_fn        : f(RH) data file path, string, default data/fRH
	output:
		output/k/iout11_[rate].txt
	'''
	if 'func' in args:
		func = args['func']
	else:
		func = 'B'
	if 'output' in args:
		output = args['output']
	else:
		output = 'input/k/aerosol.dat'
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
	
	for i in range(len(rate)):
		m_BC_new = m_BC.real + m_BC.imag * rate[i] * 1j
		runFRHChange.write(fn, nn, moma, wl, m_BC_new, func='A', output=output, m_BC_origin=m_BC_origin, m_shell_origin=m_shell, kappa_origin=kappa)
		os.system('mv '+output+' '+output+'_'+str(round(rate[i]*100)/100))
	
	os.system('mv aerosol.dat aerosol.dat_origin') # backup origin file
	
	for i in range(len(rate)):
		os.system('cp '+output+'_'+str(round(rate[i]*100)/100)+' aerosol.dat')
		os.system('sbdart >output/k/iout11_'+str(round(rate[i]*100)/100)+'.txt')
	
	write_aerosol.change_back()

def read_k(rate):
	'''
	This function is to read sequence k change brought net radiative flux change at TOA
	input:
		rate    : k change rate, array, float
	output:
		nrf     : net radiative flux at TOA for each k change rate, array, float
	'''
	nrf = np.zeros(len(rate)) # net radiative flux at TOA
	for i in range(len(rate)):
		filename = 'output/k/iout11_'+str(round(rate[i]*100)/100)+'.txt'
		iout11 = read11.read11(filename=filename)
		nrf[i] = iout11['fxdn'][0]-iout11['fxup'][0]
	return nrf

def read_k_parameter(rate):
	'''
	This function is to read sequence changed BC core complex refractive index imagine part k
	input:
		rate    : m_BC k change rate, array, float
	output:
		k       : changed m_BC k, array in shape(len(rate))
	'''
	m_BC = 1.67+0.67j
	k = np.zeros(len(rate))
	for i in range(len(rate)):
		k[i] = m_BC.imag * rate[i]
	
	return k

def run_mshell(fn, nn, moma, wl, m_BC, rate, **args):
	'''
	This function is to change aerosol.dat shell complex refractive index m_shell in certain rate
	input:
		fn              : data file name, string
		nn              : number of atmospheric levels for which aerosol information is specified
		momoa           : number of phase function moments
		wl              : the wavelength [ wl(k) < wl(k+1) ], array, nm
		m_BC            : BC complex refractive index
		rate            : change rate, float
		**func          : vertical distribution function type, A or B, default B
		**output        : output file path for albedo.dat, default input/mshell/aerosol.dat
		**m_BC_origin   : BC core origin complex refractive index, to judge if m_BC are changed in rate, default m_BC
		**m_shell_origin: shell complex refractive index, default read from wl_sca and ksca
		**kappa_origin  : hygroscopicity parameter, float, default read from f(RH) data
		**fRH_fn        : f(RH) data file path, string, default data/fRH
	output:
		output/mshell/iout11_[rate].txt
	'''
	if 'func' in args:
		func = args['func']
	else:
		func = 'B'
	if 'output' in args:
		output = args['output']
	else:
		output = 'input/mshell/aerosol.dat'
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
	
	for i in range(len(rate)):
		m_shell_new = m_shell * rate[i]
		runFRHChange.write(fn, nn, moma, wl, m_BC, func='A', output=output, m_BC_origin=m_BC_origin, m_shell_origin=m_shell_new, kappa_origin=kappa)
		os.system('mv '+output+' '+output+'_'+str(round(rate[i]*100)/100))
	
	os.system('mv aerosol.dat aerosol.dat_origin') # backup origin file
	
	for i in range(len(rate)):
		os.system('cp '+output+'_'+str(round(rate[i]*100)/100)+' aerosol.dat')
		os.system('sbdart >output/mshell/iout11_'+str(round(rate[i]*100)/100)+'.txt')
	
	write_aerosol.change_back()

def read_mshell(rate):
	'''
	This function is to read sequence m_shell change brought net radiative flux change at TOA
	input:
		rate    : m_shell change rate, array, float
	output:
		nrf     : net radiative flux at TOA for each m_shell change rate, array, float
	'''
	nrf = np.zeros(len(rate)) # net radiative flux at TOA
	for i in range(len(rate)):
		filename = 'output/mshell/iout11_'+str(round(rate[i]*100)/100)+'.txt'
		iout11 = read11.read11(filename=filename)
		nrf[i] = iout11['fxdn'][0]-iout11['fxup'][0]
	return nrf

def read_mshell_parameter(rate):
	'''
	This function is to read sequence changed shell complex refractive index real part n
	input:
		rate    : m_shell n change rate, array, float
	output:
		n       : changed m_shell n, array in shape(len(rate))
	'''
	with open('mShell.dat', 'r') as f:
		n_origin = float(f.readline()[:-1])
	
	n = np.zeros(len(rate))
	for i in range(len(rate)):
		n[i] = n_origin * rate[i]
	
	return n

if __name__ == '__main__':
	# n change from -10% to 10%, bin in 1%
	rate = np.arange(-10,11) / 100 + 1
	
	run_n('data/sp2/Taizhou.npy', 50, 6, [440,500,870,1640], 1.8+0.54j, rate, func='A', m_BC_origin=1.8+0.54j, m_shell_origin=1.5833333333333333, kappa_origin=0.21537370768611558)
	run_k('data/sp2/Taizhou.npy', 50, 6, [440,500,870,1640], 1.8+0.54j, rate, func='A', m_BC_origin=1.8+0.54j, m_shell_origin=1.5833333333333333, kappa_origin=0.21537370768611558)
	run_mshell('data/sp2/Taizhou.npy', 50, 6, [440,500,870,1640], 1.8+0.54j, rate, func='A', m_BC_origin=1.8+0.54j, m_shell_origin=1.5833333333333333, kappa_origin=0.21537370768611558)
	nrf = read_n(rate)
	print(nrf)
	nrf = read_k(rate)
	print(nrf)
	nrf = read_mshell(rate)
	print(nrf)
	'''
	n = read_k_parameter(rate)
	print(n)
	
	n = read_mshell_parameter(rate)
	print(n)
	'''
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
