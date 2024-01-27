####################################################################################
# INTRODUCTION:
# This code is to use EAR5 data to make atms.dat file for SBDART
# Created by Hebs at 21/10/26/16:33
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import readERA5
import re
import os

def write(fn, lat, lon, mon, **args):
	'''
	This function is to use EAR5 data to make atms.dat file for SBDART
	input:
		fn       : EAR5 data file name, string
		lat      : latitude, int
		lon      : longitude, int
		mon      : month, int
		**output : output path, string, default ./atms.dat
	output:
		atms.dat for SBDART
                           z  is the layer altitude in km 
                              (z must be monotonically decreasing)
                           p  is the pressure in millibars
                           t  is the temperature is Kelvin
                           wh is water vapor density g/m3
                           wo is ozone density g/m3
	'''
	if 'output' in args:
		output = args['output']
	else:
		output = 'atms.dat'
	
	era5 = readERA5.readERA5(fn, lat, lon)
	z = era5['z'][mon-1]
	p = era5['p']
	t = era5['t'][mon-1]
	wh = era5['wh'][mon-1]
	wo = era5['wo'][mon-1]
	
	with open(output, 'w') as f:
		f.write(str(len(z))+'\n')
		for i in range(len(z)):
			f.write(str(z[i])+'\t'+str(p[i])+'\t'+str(t[i])+'\t'+str(wh[i])+'\t'+str(wo[i])+'\n')

def double_layer(**args):
	'''
	This function is to use interpolate to double atms.dat layer, p/t/wh/wo for linear, z for log
	input:
		**path  : atms.dat path, default ./atms.dat
	output:
		layer number doubled atms.dat for SBDART
	'''
	if 'path' in args:
		path = args['path']
	else:
		path = 'atms.dat'
	
	z = []
	p = []
	t = []
	wh = []
	wo = []
	
	with open(path, 'r') as f:
		nn = int(f.readline()[:-1])
		for line in f.readlines():
			res = re.split('\t', line[:-1])
			z.append(res[0])
			p.append(res[1])
			t.append(res[2])
			wh.append(res[3])
			wo.append(res[4])
	
	z = np.array(z, dtype=float)
	p = np.array(p, dtype=float)
	t = np.array(t, dtype=float)
	wh = np.array(wh, dtype=float)
	wo = np.array(wo, dtype=float)
	
	z_new = []
	p_new = []
	t_new = []
	wh_new = []
	wo_new = []
	
	#attention! maxium atmospheric layer number is 65
	#if double layer number exceed 65, delete from top to down
	#the maxium new layer number is nn-(2*nn-65)=65-nn
	#if nn<65-nn, i.e.nn<32, layer number can be doubled
	
	if len(z)<32:
		for i in range(len(z)-1):
			z_new.append(z[i])
			z_new.append(np.exp(np.log(z[i+1]*z[i])/2))
			p_new.append(p[i])
			p_new.append((p[i]+p[i+1])/2)
			t_new.append(t[i])
			t_new.append((t[i]+t[i+1])/2)
			wh_new.append(wh[i])
			wh_new.append((wh[i]+wh[i+1])/2)
			wo_new.append(wo[i])
			wo_new.append((wo[i]+wo[i+1])/2)
		z_new.append(z[-1])
		p_new.append(p[-1])
		t_new.append(t[-1])
		wh_new.append(wh[-1])
		wo_new.append(wo[-1])
	else:
		for i in range(65-len(z)):
			z_new.append(z[i])
			z_new.append(np.exp(np.log(z[i+1]*z[i])/2))
			p_new.append(p[i])
			p_new.append((p[i]+p[i+1])/2)
			t_new.append(t[i])
			t_new.append((t[i]+t[i+1])/2)
			wh_new.append(wh[i])
			wh_new.append((wh[i]+wh[i+1])/2)
			wo_new.append(wo[i])
			wo_new.append((wo[i]+wo[i+1])/2)
		for i in range(65-len(z),len(z)):
			z_new.append(z[i])
			p_new.append(p[i])
			t_new.append(t[i])
			wh_new.append(wh[i])
			wo_new.append(wo[i])
	
	print(len(z_new))
	os.system('mv '+path+' '+path+'_origin') # backup the origin file
	
	with open(path, 'w') as f:
		f.write(str(len(z_new))+'\n')
		for i in range(len(z_new)):
			f.write(str(z_new[i])+'\t'+str(p_new[i])+'\t'+str(t_new[i])+'\t'+str(wh_new[i])+'\t'+str(wo_new[i])+'\n')

def change_back(**args):
	'''
	This function is to change atms.dat back
	input:
		**fn  : file path for atms.dat, default ./atms.dat
	output:
		atms.dat for SBDART
	'''
	if 'fn' in args:
		fn = args['fn']
	else:
		fn = 'atms.dat'
	
	if os.path.exists(fn+'_origin'):
		if os.path.exists(fn):
			os.system('rm '+fn)
		os.system('mv '+fn+'_origin '+fn)
	else:
		print('No origin file aerosol.dat_origin. Please check.')
		sys.exit()

def write_location():
	# Beijing: 116, 40
	# Hainan: 109, 18
	# Heilongjiang: 135, 53
	month = np.arange(1,13)
	for i in range(len(month)):
		write('data/era5/data.nc', 39, 116, month[i], output='input/atms/Beijing/atms.dat')
		double_layer(path='input/atms/Beijing/atms.dat')
		os.system('rm input/atms/Beijing/atms.dat_origin')
		os.system('mv input/atms/Beijing/atms.dat input/atms/Beijing/atms.dat_'+str(month[i]))
		write('data/era5/data.nc', 18, 109, month[i], output='input/atms/Hainan/atms.dat')
		double_layer(path='input/atms/Hainan/atms.dat')
		os.system('rm input/atms/Hainan/atms.dat_origin')
		os.system('mv input/atms/Hainan/atms.dat input/atms/Hainan/atms.dat_'+str(month[i]))
		write('data/era5/data.nc', 53, 135, month[i], output='input/atms/Heilongjiang/atms.dat')
		double_layer(path='input/atms/Heilongjiang/atms.dat')
		os.system('rm input/atms/Heilongjiang/atms.dat_origin')
		os.system('mv input/atms/Heilongjiang/atms.dat input/atms/Heilongjiang/atms.dat_'+str(month[i]))

if __name__ == '__main__':
	#write('data/era5/data.nc', 26, 126, 6)
	double_layer()
	#change_back()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
