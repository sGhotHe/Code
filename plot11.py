####################################################################################
# INTRODUCTION:
# This code is to read SBDART iout=11 output data
# Created by Hebs at 21/5/17/17:11
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import matplotlib.pyplot as plt
import read11

if __name__ == '__main__':
	iout11 = read11.read11()
	zz = iout11['zz']
	pp = iout11['pp']
	fxdn = iout11['fxdn']
	fxup = iout11['fxup']
	fxdir = iout11['fxdir']
	dfdz = iout11['dfdz']
	heat = iout11['heat']
	
	fig = plt.figure()
	plt.plot(heat, zz)
	plt.show()
	plt.close()
