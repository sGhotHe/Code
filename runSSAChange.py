####################################################################################
# INTRODUCTION:
# This code is to test SSA sensitivity of DARF in SBDART
# Created by Hebs at 21/10/28/11:13
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import write_aerosol
import read11
import os
import re

def run(rate):
	'''
	This function is to run sequence aerosol.dat for SBDART
	input:
		rate    : SSA change rate, array, float
	output:
		output/SSA/iout11_[rate].txt
	'''
	for i in range(len(rate)):
		write_aerosol.change_SSA(rate[i])
		os.system('sbdart >output/SSA/iout11_'+str(round(rate[i]*100)/100)+'.txt')
		write_aerosol.change_back()
	
def read(rate):
	'''
	This function is to read sequence SSA change brought net radiative flux change at TOA
	input:
		rate    : SSA change rate, array, float
	output:
		nrf     : net radiative flux at TOA for each SSA change rate, array, float
	'''
	nrf = np.zeros(len(rate)) # net radiative flux at TOA
	for i in range(len(rate)):
		filename = 'output/SSA/iout11_'+str(round(rate[i]*100)/100)+'.txt'
		iout11 = read11.read11(filename=filename)
		nrf[i] = iout11['fxdn'][0]-iout11['fxup'][0]
	return nrf

def read_parameter(rate):
	'''
	This function is to read sequence changed single scattering albedo SSA
	input:
		rate    : ssa change rate, array
	output:
		wl      : wave length, array, nm
		ssa     : changed ssa, array in shape(len(wl), len(rate))
	'''
	with open('aerosol.dat', 'r') as f:
		res = re.split('\t', f.readline()[:-1])
		nn = int(res[0])
		data = np.array(f.readlines())
		data = data.reshape(round(len(data)/(nn+1)),-1)
		wl = np.zeros(len(data))
		ssa_origin = np.zeros(len(wl))
		for i in range(len(wl)):
			wl[i] = float(data[i,0][:-1])
			res = re.split('\t', data[i,1][:-1])
			ssa_origin[i] = float(res[1])
	
	ssa = np.zeros((len(wl), len(rate)))
	for i in range(len(wl)):
		for j in range(len(rate)):
			ssa[i,j] = ssa_origin[i] * rate[j]
	
	return wl, ssa

if __name__ == '__main__':
	#albedo change from -10% to 10%, bin in 1%
	rate = np.arange(-10,11) / 100 + 1
	
	run(rate)
	nrf = read(rate)
	
	print(nrf)
	'''
	wl, ssa = read_parameter(rate)
	print(wl)
	print(ssa[0])
	'''
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
