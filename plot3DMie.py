####################################################################################
# INTRODUCTION:
# This code is to use pymiescatt to plot 3D Mie function results,
# beta (k*a), n and Qext for x, y, z separately
# Created by Hebs at 24/6/20/9:08
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import matplotlib.pyplot as plt
import PyMieScatt as ps

def plot(wl):
	'''
	This function is main plotting function
	input:
		wl	: wavelength, nm, float
	output:
		3D Mie function results figure
	'''
	beta_log = np.arange(-1,1,0.01)
	beta = 10**beta_log
	n = np.arange(1.3,1.5,0.01)
	
	mBC = 1.67+0.67j
	
	Qext = np.zeros((len(beta),len(n)))
	Qsca = np.zeros((len(beta),len(n)))
	g = np.zeros((len(beta),len(n)))
	
	for i in range(len(beta)):
		for j in range(len(n)):
			MieQCoreShell = ps.MieQCoreShell(mBC, n[j]+1e-3j, wl, beta[i]*wl*0.5, beta[i]*wl)
			Qext[i,j] = MieQCoreShell[0]
			Qsca[i,j] = MieQCoreShell[1]
			g[i,j] = MieQCoreShell[3]
	
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	X, Y = np.meshgrid(n, beta_log)
	ax.plot_surface(X, Y, Qext)
	ax.set_xlabel('n')
	ax.set_ylabel('beta')
	ax.set_zlabel('Qext')
	plt.show()

if __name__ == '__main__':
	plot(550)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
