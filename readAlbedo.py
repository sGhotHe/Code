####################################################################################
# INTRODUCTION:
# This code is to read albedo from albedo.dat
# Created by Hebs at 24/2/27/10:28
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import re

def read(**args):
	'''
	This function is to read albedo from albedo.dat
	input:
		**fn	: aerosol imformation file name, string, default 'aerosol.dat'
	output:
		wl		: wave length, array, nm, int
		albedo	: albedo, array, float
	'''
	if 'fn' in args:
		fn = args['fn']
	else:
		fn = 'albedo.dat'
	
	with open(fn, 'r') as f:
		lines = f.readlines()
		wl = []
		albedo = []
		i = 0
		while i < len(lines):
			infos = re.split('\t',lines[i])
			wl.append(round(float(infos[0])*1e3))
			albedo.append(float(infos[0]))
			i += 1
	
	return wl, albedo

if __name__ == '__main__':
	wl, albedo = read()
	print(wl)
	print(albedo)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
