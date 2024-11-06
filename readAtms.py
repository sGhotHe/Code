####################################################################################
# INTRODUCTION:
# This code is to read atms.dat data
# Created by Hebs at 21/10/25/10:20
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import re

def read(**args):
	'''
	This function is to read data from atms.dat
	input:
		**fn   : file name, string, default 'atms.dat'
	output:
		atms : atmosphere data, dictionary
			nn : number of atmosperic layers
			z  : layer altitude in km; p: pressure in millibars
			t  : temperature in Kelvin
			wh : water vapor density in g/m^3
			wo : ozone density in g/m^3
	'''
	if 'fn' in args:
		fn = args['fn']
	else:
		fn = 'atms.dat'
	
	with open(fn, 'r') as f:
		nn = int(f.readline())
		z = np.zeros(nn)
		p = np.zeros(nn)
		t = np.zeros(nn)
		wh = np.zeros(nn)
		wo = np.zeros(nn)
		
		i = 0
		for line in f.readlines():
			res = re.split('\t', line[:-1])
			z[i] = float(res[0])
			p[i] = float(res[1])
			t[i] = float(res[2])
			wh[i] = float(res[3])
			wo[i] = float(res[4])
			i += 1
	
	atms = dict(nn=nn, z=z, p=p, t=t, wh=wh, wo=wo)
	return atms

if __name__ == '__main__':
	atms = read()
	print(atms['z'])
	print(atms['z'][::-1])
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
