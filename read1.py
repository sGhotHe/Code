####################################################################################
# INTRODUCTION:
# This code is to read SBDART iout=1 output data
# Created by Hebs at 21/5/18/15:50
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

def read1(**args):
	if 'filename' in args:
		filename = args['filename']
	else:
		filename = 'output/iout1.txt'
	
	lines = []
	with open(filename, 'r') as f:
		for res in f.readlines():
			res = res.rstrip('\n')
			res = res.lstrip(' ')
			res = re.split(' ', res)
			res = clean(res)
			lines.append(res)
	
	nw = int(lines[2][0])
	lines = lines[3:]
	wl = np.zeros(nw)
	ffv = np.zeros(nw)
	topdn = np.zeros(nw)
	topup = np.zeros(nw)
	topdir = np.zeros(nw)
	botdn = np.zeros(nw)
	botup = np.zeros(nw)
	botdir = np.zeros(nw)
	
	for i in range(nw):
		line = lines[i]
		wl[i] = float(line[0])
		ffv[i] = float(line[1])
		topdn[i] = float(line[2])
		topup[i] = float(line[3])
		topdir[i] = float(line[4])
		botdn[i] = float(line[5])
		botup[i] = float(line[6])
		botdir[i] = float(line[7])
	
	iout1 = dict(wl=wl, ffv=ffv, topdn=topdn, topup=topup, topdir=topdir, botdn=botdn, botup=botup, botdir=botdir)
	return iout1

if __name__ == '__main__':
	iout1 = read1()
	print(iout1['wl'])
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
