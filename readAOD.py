####################################################################################
# INTRODUCTION:
# This code is to read AOD from aerosol.dat
# Created by Hebs at 23/6/5/14:59
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np

import re

def read(**args):
	'''
	This function is to read AOD from aerosol.dat
	input:
		**fn	: aerosol imformation file name, string, default 'aerosol.dat'
	output:
		wl	: wave length, array, int
		AOD	: aerosol optical depth, array, float
	'''
	if 'fn' in args:
		fn = args['fn']
	else:
		fn = 'aerosol.dat'
	
	with open(fn, 'r') as f:
		line = re.split('\t', f.readline()[:-1])
		nn = int(line[0])
		moma = int(line[1])
		infos = f.readlines()
		wl = []
		AOD = []
		i = 0
		while i < len(infos):
			wl.append(int(infos[i]))
			i += 1
			AOD_j = 0
			for j in range(nn):
				AOD_j += float(re.split('\t', infos[i][:-2])[0])
				i += 1
			AOD.append(AOD_j)
	
	return wl, AOD

if __name__ == '__main__':
	wl_367, AOD_367 = read(fn='1p367.dat')
	wl_49, AOD_49 = read(fn='1p49.dat')
	wl_51, AOD_51 = read(fn='1p51.dat')
	with open('AOD_infos.txt', 'w') as f:
		f.write('1.367:\n')
		f.write('wl: '+str(wl_367)+' , AOD: '+str(AOD_367)+'\n')
		f.write('1.49:\n')
		f.write('wl: '+str(wl_49)+' , AOD: '+str(AOD_49)+'\n')
		f.write('1.51:\n')
		f.write('wl: '+str(wl_51)+' , AOD: '+str(AOD_51)+'\n')
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
