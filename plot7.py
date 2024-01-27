####################################################################################
# INTRODUCTION:
# This code is to plot SBDART iout=7 output data
# Created by Hebs at 21/5/17/16:45
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import matplotlib.pyplot as plt
import read7

if __name__ == '__main__':
	iout7 = read7.read7()
	Z = iout7['Z']
	wl = iout7['wl']
	fdird = iout7['fdird']
	fdifd = iout7['fdifd']
	flxdn = iout7['flxdn']
	flxup = iout7['flxup']
	
	fig = plt.figure()
	plt.plot(flxdn[29]-flxup[29], Z[29], label=str(wl[29])+' um')
	plt.title('Net downward flux')
	plt.ylim(0)
	plt.legend()
	plt.grid(which='major', lw=1.5, c='k', linestyle=':', alpha=0.8)
	plt.show()
	plt.close()
