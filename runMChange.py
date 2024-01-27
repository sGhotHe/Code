####################################################################################
# INTRODUCTION:
# This code is to test m sensitivity of DARF in SBDART
# Created by Hebs at 21/10/28/15:32
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import write_aerosol
import read11
import os

def run_n_inhomo():
	'''
	This function is to test how much n inhomogeneity can affect DARF sensitivity
	input:
	output:
	'''
	# how to express n inhomogeneity?
	# refer to black carbon, use two dimentional variable to represented,
	# black carbon use whole and single particle component ratio,
	# similary, n can use whole and single complex refractive index,
	# but firstly, use same way as black carbon, study how Mie parameters change
	# when black carbon inhomogeneity change and keep n constant,
	# according to Zhao 2020, n may have linear change by partical diameter,
	# study how Mie parameters change by slope rate change.

def run_n(fn, nn, moma, wl, m, rate, **args):
	'''
	This function is to change aerosol.dat complex refractive index real part n in certain rate
	input:
		fn          : PNSD data file name, smps.txt
		nn          : number of atmospheric levels for which aerosol information is specified
		momoa       : number of phase function moments
		wl          : the wavelength [ wl(k) < wl(k+1) ], array, nm
		m           : complex refractive index
		rate        : change rate, float
		**func      : vertical distribution function type, A or B, default B
		**output    : output file path for albedo.dat, default input/n/aerosol.dat
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
	
	for i in range(len(rate)):
		m_new = m.real * rate[i] + m.imag*1j
		write_aerosol.write(fn, nn, moma, wl, m_new, func=func, output=output)
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

def run_k(fn, nn, moma, wl, m, rate, **args):
	'''
	This function is to change aerosol.dat complex refractive index imagine part k in certain rate
	input:
		fn          : PNSD data file name, smps.txt
		nn          : number of atmospheric levels for which aerosol information is specified
		momoa       : number of phase function moments
		wl          : the wavelength [ wl(k) < wl(k+1) ], array, nm
		m           : complex refractive index
		rate        : change rate, float
		**func      : vertical distribution function type, A or B, default B
		**output    : output file path for albedo.dat, default input/n/aerosol.dat
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
	
	for i in range(len(rate)):
		m_new = m.real + m.imag * rate[i] * 1j
		write_aerosol.write(fn, nn, moma, wl, m_new, func=func, output=output)
		os.system('mv '+output+' '+output+'_'+str(round(rate[i]*100)/100))
	
	os.system('mv aerosol.dat aerosol.dat_origin') # backup origin file
	
	for i in range(len(rate)):
		os.system('cp '+output+'_'+str(round(rate[i]*100)/100)+' aerosol.dat')
		os.system('sbdart >output/k/iout11_'+str(round(rate[i]*100)/100)+'.txt')
	
	write_aerosol.change_back()

def read_k(rate):
	'''
	This function is to read sequence n change brought net radiative flux change at TOA
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

if __name__ == '__main__':
	#write_aerosol.write('data/smps/smps.txt', 50, 6, [440,500,870,1640], 1.53+0.1j, func='A')
	#n change from -10% to 10%, bin in 1%
	rate = np.arange(-10,11) / 100 + 1
	#run_n('data/smps/smps.txt', 50, 6, [440,500,870,1640], 1.53+0.1j, rate, func='A')
	print(read_n(rate))
	#run_k('data/smps/smps.txt', 50, 6, [440,500,870,1640], 1.53+0.1j, rate, func='A')
	print(read_k(rate))
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
