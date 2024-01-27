####################################################################################
# INTRODUCTION:
# This code is to read SBDART iout=11 output data
# Created by Hebs at 21/5/17/17:00
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import re

def clean(li):
	i = 0
	while i<len(li):
		if li[i]=='':
			li.pop(i)
		else:
			i = i + 1
	return li

def read11(**args):
	if 'filename' in args:
		filename = args['filename']
	else:
		filename = 'output/iout11.txt'
	
	lines = []
	with open(filename, 'r') as f:
		for res in f.readlines():
			res = res.rstrip('\n')
			res = res.lstrip(' ')
			res = re.split(' ', res)
			res = clean(res)
			lines.append(res)
	
	nz = lines[0][0] # number of atmospheric layers
	phidw = lines[0][1] # filter function equivalent width, um
	zz = [] # level altitudes, km
	pp = [] # level pressure, mb
	fxdn = [] # downward flux (direct+diffuse), W/m2
	fxup = [] # upward flux, W/m2
	fxdir = [] # downward flux, direct beam only, W/m2
	dfdz = [] # radiant energy flux divergence, mW/m3
	heat = [] # heating rate, K/day
	
	for i in range(1,len(lines)):
		line = lines[i]
		zz.append(line[0])
		pp.append(line[1])
		fxdn.append(line[2])
		fxup.append(line[3])
		fxdir.append(line[4])
		dfdz.append(line[5])
		heat.append(line[6])
	
	zz = np.array(zz, dtype=float)
	pp = np.array(pp, dtype=float)
	fxdn = np.array(fxdn, dtype=float)
	fxup = np.array(fxup, dtype=float)
	fxdir = np.array(fxdir, dtype=float)
	dfdz = np.array(dfdz, dtype=float)
	heat = np.array(heat, dtype=float)
	
	iout11 = dict(zz=zz, pp=pp, fxdn=fxdn, fxup=fxup, fxdir=fxdir, dfdz=dfdz, heat=heat)
	return iout11

if __name__ == '__main__':
	iout11 = read11(filename='output/albedo/iout11_0.9.txt')
	print(iout11['fxdn']-iout11['fxup'])
	print(iout11['zz'])
	iout11 = read11(filename='output/albedo/iout11_1.1.txt')
	print(iout11['fxdn']-iout11['fxup'])
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
