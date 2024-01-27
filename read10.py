####################################################################################
# INTRODUCTION:
# This code is to read SBDART iout=10 output data
# Created by Hebs at 21/5/18/14:27
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

def read10(**args):
	if 'filename' in args:
		filename = args['filename']
	else:
		filename = 'iout10.txt'
	
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
	ffew = float(line[2])
	topdn = float(line[3])
	topup = float(line[4])
	topdir = float(line[5])
	botdn = float(line[6])
	botup = float(line[7])
	botdir = float(line[8])
	
	iout10 = dict(wlinf=wlinf, wlsup=wlsup, ffew=ffew, topdn=topdn, topup=topup, topdir=topdir, botdn=botdn, botup=botup, botdir=botdir)
	return iout10

if __name__ == '__main__':
	iout10 = read10()
	print(iout10['topdn'])
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
