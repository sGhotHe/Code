####################################################################################
# INTRODUCTION:
# This code is to test AOD sensitivity of DARF in SBDART
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
		output/AOD/iout11_[rate].txt
	'''
	for i in range(len(rate)):
		write_aerosol.change_AOD(rate[i])
		os.system('sbdart >output/AOD/iout11_'+str(round(rate[i]*100)/100)+'.txt')
		write_aerosol.change_back()
	
def read(rate):
	'''
	This function is to read sequence AOD change brought net radiative flux change at TOA
	input:
		rate    : AOD change rate, array, float
	output:
		nrf     : net radiative flux at TOA for each AOD change rate, array, float
	'''
	nrf = np.zeros(len(rate)) # net radiative flux at TOA
	for i in range(len(rate)):
		filename = 'output/AOD/iout11_'+str(round(rate[i]*100)/100)+'.txt'
		iout11 = read11.read11(filename=filename)
		nrf[i] = iout11['fxdn'][0]-iout11['fxup'][0]
	return nrf

def read_parameter(rate):
	'''
	This function is to read sequence changed AOD
	input:
		rate    : albedo change rate, array, float
	output:
		wl      : wave length, array, nm
		aod     : changed total aerosol aod, array in shape(len(wl), len(rate))
	'''
	
	with open('aerosol.dat', 'r') as f:
		res = re.split('\t', f.readline()[:-1])
		nn = int(res[0])
		data = np.array(f.readlines())
		data = data.reshape(round(len(data)/(nn+1)),-1)
		wl = np.zeros(len(data))
		aod_origin = np.zeros(len(wl))
		aod = np.zeros((len(wl), len(rate)))
		for i in range(len(wl)):
			wl[i] = float(data[i,0][:-1])
			for j in range(nn):
				res = re.split('\t', data[i,j+1][:-1])
				aod_origin[i] += float(res[0])
	
	for i in range(len(wl)):
		for j in range(len(rate)):
			aod[i,j] = aod_origin[i] * rate[j]
	
	return wl, aod

if __name__ == '__main__':
	#AOD change from -10% to 10%, bin in 1%
	rate = np.arange(-10,11) / 100 + 1
	'''
	run(rate)
	nrf = read(rate)
	
	print(nrf)
	
	wl, aod = read_parameter(rate)
	print(wl)
	print(aod[0])
	'''
	name = ['Beijing', 'Hainan', 'Heilongjiang']
	
	for i in range(len(name)):
		path = 'input/albedo/' + name[i]
		fn = path + '/albedo.dat'
		for j in range(rate):
			write_aerosol.change_AOD(rate, fn=fn)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
