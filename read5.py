####################################################################################
# INTRODUCTION:
# This code is to read SBDART iout=5 output data
# Created by Hebs at 21/5/18/14:44
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

def read5(**args):
	if 'filename' in args:
		filename = args['filename']
	else:
		filename = 'iout5.txt'
	
	lines = []
	with open(filename, 'r') as f:
		for res in f.readlines():
			res = res.rstrip('\n')
			res = res.lstrip(' ')
			res = re.split(' ', res)
			res = clean(res)
			lines.append(res)
	
	nw = lines[2][0] # number of wavelength
	line_num = len(lines)
	wl = [] # wavelength, um
	ffv = [] # filter function value
	topdn = [] # total downward flux at ZOUT(2) km, w/m2
	topup = [] # total upward flux at ZOUT(2) km, w/m2
	topdir = [] # direct downward flux at ZOUT(2) km, w/m2
	botdn = [] # total downward flux at ZOUT(1) km, w/m2
	botup = [] # total upward flux at ZOUT(1) km, w/m2
	botdir = [] # direct downward flux at ZOUT(1) km, w/m2
	nphi = [] # number of user azimuth angles
	nzen = [] # number of user zenith angles
	phi = [] # user relative azimuth angles (nphi values)
	uzen = [] # user zenith angles (nzen values)
	uurs = [] # radiance at user angles at altitude ZOUT(2), w/m2/um/str
	
	i = 2
	while i<line_num-1:
		i += 1
		line = lines[i]
		wl.append(line[0])
		ffv.append(line[1])
		topdn.append(line[2])
		topup.append(line[3])
		topdir.append(line[4])
		botdn.append(line[5])
		botup.append(line[6])
		botdir.append(line[7])
		
		i += 1
		nphii = int(lines[i][0])
		nzeni = int(lines[i][1])
		nphi.append(nphii)
		nzen.append(nzeni)
		
		phii_line_num = nphii // 10 + 1
		phii = []
		for j in range(phii_line_num):
			i += 1
			phii += lines[i]
		phi.append(phii)
		
		uzeni_line_num = nzeni // 10 + 1
		uzeni = []
		for j in range(uzeni_line_num):
			i += 1
			uzeni += lines[i]
		uzen.append(uzeni)
		
		uursi = np.zeros((nzeni, nphii))
		for j in range(nzeni):
			uursii = []
			for k in range(phii_line_num):
				i += 1
				uursii += lines[i]
			uursi[j] = uursii
		uurs.append(uursi)
	
	wl = np.array(wl, dtype=float)
	ffv = np.array(ffv, dtype=float)
	topdn = np.array(topdn, dtype=float)
	topup = np.array(topup, dtype=float)
	topdir = np.array(topdir, dtype=float)
	botdn = np.array(botdn, dtype=float)
	botup = np.array(botup, dtype=float)
	botdir = np.array(botdir, dtype=float)
	nphi = np.array(nphi, dtype=int)
	nzen = np.array(nzen, dtype=int)
	phi = np.array(phi, dtype=float)
	uzen = np.array(uzen, dtype=float)
	uurs = np.array(uurs, dtype=float)
	
	iout5 = dict(wl=wl, ffv=ffv, topdn=topdn, topup=topup, topdir=topdir, botdn=botdn, botup=botup, botdir=botdir, nphi=nphi, nzen=nzen, phi=phi, uzen=uzen, uurs=uurs)
	return iout5

if __name__ == '__main__':
	iout5 = read5()
	print(iout5['wl'])
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
