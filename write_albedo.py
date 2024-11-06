####################################################################################
# INTRODUCTION:
# This code is to use MODIS data to make albedo.dat
# Created by Hebs at 21/5/24/18:06
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import readMODIS
import re
import os
import getfilesname
import sys

def write(lon, lat, year, **args):
	if 'month' in args:
		month = args['month']
	if 'day' in args:
		day = args['day']
	if 'jul' in args:
		jul = args['jul']
		times = jul2date(year,jul)
		month = times[1]
		day = times[2]
	if ('jul' not in args) and (('month' not in args) or ('day' not in args)) :
		print('The time variable is not appropriate, please check')
		sys.exit()
	if 'path' in args:
		path=args['path']
	else:
		path='/home/hebs/Code/SBDART/data/modis/albedo/'
	if 'satellite' in args:
		satellite=args['satellite']
	else:
		satellite=0
	if 'output' in args:
		output = args['output']
	else:
		output = 'albedo.dat'
	
	albedos = np.zeros(7) # 7: band num
	wls = np.zeros(7)
	for i in range(7):
		albedos[i], wls[i] = readMODIS.read_modis_albedo(lon, lat, year, month=month, day=day, band=i+1, path=path, satellite=satellite)
	
	albedos = albedos / 100 # in unit, not percent
	print('albedos:', albedos)
	print('wavelengths:', wls)
	with open(output, 'w') as f:
		for i in [2, 3, 0, 1, 4, 5, 6]: # from small to big
			f.write(str(wls[i]))
			f.write('\t')
			f.write(str(albedos[i]))
			f.write('\n')

def change(rate, **args):
	'''
	This function is to change albedo.dat
	input:
		rate  : change rate, float
		**fn  : file path for albedo.dat, default ./albedo.dat
	output:
		albedo.dat for SBDART
	'''
	if 'fn' in args:
		fn = args['fn']
	else:
		fn = 'albedo.dat'
	
	wl = []
	albedo = []
	
	with open(fn, 'r') as f:
		for line in f.readlines():
			res = re.split('\t', line[:-1])
			wl.append(res[0])
			albedo.append(res[1])
	
	if not os.path.exists(fn+'_origin'):
		os.system('mv '+fn+' '+fn+'_origin')
	
	albedo = np.array(albedo, dtype=float)
	albedo = albedo * rate
	
	with open(fn, 'w') as f:
		for i in range(len(wl)):
			f.write(wl[i]+'\t'+str(albedo[i])+'\n')

def change_back(**args):
	'''
	This function is to change albedo.dat back
	input:
		**fn  : file path for albedo.dat, default ./albedo.dat
	output:
		albedo.dat for SBDART
	'''
	if 'fn' in args:
		fn = args['fn']
	else:
		fn = 'albedo.dat'
	
	if os.path.exists(fn+'_origin'):
		if os.path.exists(fn):
			os.system('rm '+fn)
		os.system('mv '+fn+'_origin '+fn)
	else:
		print('No origin file albedo.dat_origin. Please check.')

def write_location():
	# Beijing: 116, 40
	# Hainan: 109, 18
	# Heilongjiang: 135, 53
	month = np.arange(1,13)
	for i in range(len(month)):
		write(116, 40, 2019, month=month[i], day=1, output='input/albedo/Beijing/albedo.dat')
		os.system('mv input/albedo/Beijing/albedo.dat input/albedo/Beijing/albedo.dat_'+str(month[i]))
			
		write(109, 18, 2019, month=month[i], day=1, output='input/albedo/Hainan/albedo.dat')
		os.system('mv input/albedo/Hainan/albedo.dat input/albedo/Hainan/albedo.dat_'+str(month[i]))
			
		write(135, 53, 2019, month=month[i], day=1, output='input/albedo/Heilongjiang/albedo.dat')
		os.system('mv input/albedo/Heilongjiang/albedo.dat input/albedo/Heilongjiang/albedo.dat_'+str(month[i]))

if __name__ == '__main__':
	write(116, 40, 2019, month=6, day=1)
	#write(126, 26, 2019, month=6, day=1)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
