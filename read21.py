####################################################################################
# INTRODUCTION:
# This code is to read SBDART iout=21 output data
# Created by Hebs at 21/5/18/9:54
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

def read21(**args):
	if 'filename' in args:
		filename = args['filename']
	else:
		filename = 'iout21.txt'
	
	lines = []
	with open(filename, 'r') as f:
		for res in f.readlines():
			res = res.rstrip('\n')
			res = res.lstrip(' ')
			res = re.split(' ', res)
			res = clean(res)
			lines.append(res) 
	
	line = lines[0]
	wlinf = float(line[0])
	wlsup = float(line[1])
	ffew = float(line[2]) # filter function equivalent width, um
	topdn = float(line[3]) # total downward flux at ZOUT(2) km, w/m2
	topup = float(line[4]) # total upward flux at ZOUT(2) km, w/m2
	topdir = float(line[5]) # direct downward flux at ZOUT(2) km, w/m2
	botdn = float(line[6]) # total downward flux at ZOUT(1) km, w/m2
	botup = float(line[7]) # total upward flux at ZOUT(1) km, w/m2
	botdir = float(line[8]) # direct downward flux at ZOUT(1) km, w/m2
	line = lines[1]
	nphi = int(line[0]) # number of user azimuth angles
	nzen = int(line[1]) # number of user zenith angles
	
	phi_line_num = nphi // 10 + 1
	uzen_line_num = nzen // 10 + 1
	phi = [] # user relative azimuth angles (nphi values)
	uzen = [] # user zenith angles (nzen values)
	i = 1
	for j in range(phi_line_num):
		i += 1
		phi = phi + lines[i]
	for j in range(uzen_line_num):
		i += 1
		uzen = uzen + lines[i]
	phi = np.array(phi, dtype=float)
	uzen = np.array(uzen, dtype=float)
	
	r = np.zeros((nzen, nphi))
	for j in range(nzen):
		i += 1
		r[j] = lines[i]
	
	iout20 = dict(wlinf=wlinf, wlsup=wlsup, ffew=ffew, topdn=topdn, topup=topup, topdir=topdir, botdn=botdn, botdir=botdir, nphi=nphi, nzen=nzen, phi=phi, uzen=uzen, r=r)
	return iout20

if __name__ == '__main__':
	iout21 = read21()
	print(iout21['r'][10])
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
