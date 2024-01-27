####################################################################################
# INTRODUCTION:
# This code is to read SBDART iout=7 output data
# Created by Hebs at 21/5/17/16:40
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

def read7(**args):
	if 'filename' in args:
		filename = args['filename']
	else:
		filename = 'iout7.txt'
	
	lines = []
	with open(filename, 'r') as f:
		for res in f.readlines():
			res = res.rstrip('\n')
			res = res.lstrip(' ')
			res = re.split(' ', res)
			res = clean(res)
			lines.append(res)
	
	nz = int(lines[2][0])
	line_num = len(lines)
	each_num = nz // 10 + 1
	each_num2 = 2 + 1 + (1+each_num)*5
	wl = []
	Z = [] # altitude, km
	fdird = [] # downward direct flux, w/m2/um
	fdifd = [] # downward diffuse flux, w/m2/um
	flxdn = [] # total downward flux, w/m2/um
	flxup = [] # total upward flux, w/m2/um
	
	i = 2
	while i<line_num-1:
		i += 3
		wl.append(float(lines[i][0]))
		
		i += 1
		zi = []
		for j in range(each_num):
			i += 1
			zi = zi + lines[i]
		Z.append(zi)
		
		i += 1
		fdirdi = []
		for j in range(each_num):
			i += 1
			fdirdi = fdirdi + lines[i]
		fdird.append(fdirdi)
		
		i += 1
		fdifdi = []
		for j in range(each_num):
			i += 1
			fdifdi = fdifdi + lines[i]
		fdifd.append(fdifdi)
		
		i += 1
		flxdni = []
		for j in range(each_num):
			i += 1
			flxdni = flxdni + lines[i]
		flxdn.append(flxdni)
		
		i += 1
		flxupi = []
		for j in range(each_num):
			i += 1
			flxupi = flxupi + lines[i]
		flxup.append(flxupi)
	
	Z = np.array(Z, dtype=float)
	wl = np.array(wl, dtype=float)
	fdird = np.array(fdird, dtype=float)
	fdifd = np.array(fdifd, dtype=float)
	flxdn = np.array(flxdn, dtype=float)
	flxup = np.array(flxup, dtype=float)
	
	iout7 = dict(Z=Z, wl=wl, fdird=fdird, fdifd=fdifd, flxdn=flxdn, flxup=flxup)
	return iout7

if __name__ == '__main__':
	iout7 = read7()
	print(iout7['Z'])
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
