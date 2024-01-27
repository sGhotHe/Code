####################################################################################
# INTRODUCTION:
# This code is to test albedo sensitivity of DARF in SBDART
# Created by Hebs at 21/10/28/9:29
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import os
import write_albedo
import read11
import re

def run(rate, **args):
	'''
	This function is to run sequence albedo.dat SBDART
	input:
		rate     : albedo change rate, array, float
		**path   : file path, string, default ./albedo.dat
		**output : output path, string, default output/albedo
	output:
		output/albedo/iout11_[rate].txt
	'''
	if 'path' in args:
		path = args['path']
	else:
		path = 'albedo.dat'
	if 'output' in args:
		output = args['output']
	else:
		output = 'output/albedo'
	
	for i in range(len(rate)):
		write_albedo.change(rate[i], fn=path)
		os.system('sbdart >'+output+'/iout11_'+str(round(rate[i]*100)/100)+'.txt')
		write_albedo.change_back(fn=path)

def read(rate, **args):
	'''
	This function is to read sequence albedo change brought net radiative flux change at TOA
	input:
		rate    : albedo change rate, array, float output/albedo
		**path  : file path, string, default 
	output:
		nrf     : net radiative flux at TOA for each albedo change rate, array, float
	'''
	if 'path' in args:
		path = args['path']
	else:
		path = 'output/albedo'
	
	nrf = np.zeros(len(rate)) # net radiative flux at TOA
	for i in range(len(rate)):
		filename = path + '/iout11_'+str(round(rate[i]*100)/100)+'.txt'
		iout11 = read11.read11(filename=filename)
		nrf[i] = iout11['fxdn'][0]-iout11['fxup'][0]
	return nrf

def read_parameter(rate, **args):
	'''
	This function is to read sequence changed albedo
	input:
		rate    : albedo change rate, array, float
		**path  : albedo.dat file path, string, default ./albedo.dat
	output:
		wl      : wave length, array, nm
		albedo  : changed albedo, array in shape (len(wl), len(rate))
	'''
	if 'path' in args:
		path = args['path']
	else:
		path = 'albedo.dat'
	
	with open(path, 'r') as f:
		data = f.readlines()
		wl = np.zeros(len(data))
		albedo = np.zeros((len(data), len(rate)))
		albedo_origin = np.zeros(len(data))
		for i in range(len(data)):
			res = re.split('\t', data[i][:-1])
			wl[i] = res[0]
			albedo_origin[i] = res[1]
	
	for i in range(len(wl)):
		for j in range(len(rate)):
			albedo[i,j] = albedo_origin[i] * rate[j]
	
	return wl, albedo

def write_location():
	# albedo change from -10% to 10%, bin in 1%
	rate = np.arange(-10,11) / 100 + 1
	name = ['Beijing', 'Hainan', 'Heilongjiang']
	
	for i in range(len(name)):
		path = 'input/albedo/' + name[i]
		fn = path + '/albedo.dat'
		for j in range(12): # 12 months
			os.system('cp '+fn+'_'+str(j+1)+' '+fn)
			path_j = path + '/month_' + str(j+1)
			if not os.path.exists(path_j):
				os.system('mkdir '+path_j)
			for k in range(len(rate)):
				write_albedo.change(rate[k], fn=fn)
				os.system('mv '+fn+' '+path_j+'/albedo.dat_'+str(round(rate[k]*100)/100))
				write_albedo.change_back(fn=fn)
			os.system('rm '+fn)

if __name__ == '__main__':
	# albedo change from -10% to 10%, bin in 1%
	rate = np.arange(-10,11) / 100 + 1
	'''
	run(rate)
	nrf = read(rate)
	
	print(nrf)
	wl, albedo = read_parameter(rate)
	print(wl)
	print(albedo[0])
	'''
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
