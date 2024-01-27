####################################################################################
# INTRODUCTION:
# This code is to read SBDART iout=22 output data
# Created by Hebs at 21/5/18/13:08
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

def read22(**args):
	if 'filename' in args:
		filename = args['filename']
	else:
		filename = 'iout22.txt'
	
	lines = []
	with open(filename, 'r') as f:
		for res in f.readlines():
			res = res.rstrip('\n')
			res = res.lstrip(' ')
			res = re.split(' ', res)
			res = clean(res)
			lines.append(res)
	
	line = lines[0]
	nphi = int(line[0]) # number of user specified azimuth angles
	nzen = int(line[1]) # number of user specified zenith angles
	nz = int(line[2]) # number of atmospheric levels
	ffew = float(line[3]) # filter function equivalent width, um
	
	i = 0
	phi_line_num = nphi // 10 + 1
	phi = [] # user specified anizmuth angles, degrees
	for j in range(phi_line_num):
		i += 1
		phi += lines[i]
	
	uzen_line_num = nzen // 10 + 1
	uzen = [] # user specified zenith angles, degrees
	for j in range(uzen_line_num):
		i += 1
		uzen += lines[i]
	
	z_line_num = nz // 10 + 1
	z = [] # altitudes of atmospheric layers, km
	for j in range(z_line_num):
		i += 1
		z += lines[i]
	
	fxdn = [] # downward flux (direct+diffuse), W/m2
	for j in range(z_line_num):
		i += 1
		fxdn += lines[i]
	
	fxup = [] # upward flux, W/m2
	for j in range(z_line_num):
		i += 1
		fxup += lines[i]
	
	fxdir = [] # downward flux, direct beam only, W/m2
	for j in range(z_line_num):
		i += 1
		fxdir += lines[i]
	
	phi = np.array(phi, dtype=float)
	uzen = np.array(uzen, dtype=float)
	z = np.array(z, dtype=float)
	fxdn = np.array(fxdn, dtype=float)
	fxup = np.array(fxup, dtype=float)
	fxdir = np.array(fxdir, dtype=float)
	
	i += 1
	uurl_infos = lines[i:]
	uurl = np.zeros((nz, nzen, nphi)) # radiance at each layer, W/m2/str
	for i in range(nz):
		for j in range(nzen):
			for k in range(nphi):
				ii = (i+j+k) // 10
				jj = (i+j+k) % 10
				uurl[i,j,k] = uurl_infos[ii][jj]
	uurl = np.array(uurl, dtype=float)
	
	iout22 = dict(nphi=nphi, nzen=nzen, nz=nz, ffew=ffew, phi=phi, uzen=uzen, z=z, fxdn=fxdn, fxup=fxup, fxdir=fxdir, uurl=uurl)
	return iout22

if __name__ == '__main__':
	iout22 = read22()
	print(iout22['uurl'][0,0])
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
