####################################################################################
# INTRODUCTION:
# This code is Planck Function
# Created by Hebs at 21/5/18/17:18
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np

def B(lamda, T):
	h = 6.6262e-34 # J s
	c = 2.99793e8 # m/s
	k = 1.3806e-23 # J/K
	lamda = lamda * 1e-6 #um
	B = 2*h*c**2/lamda**5/(np.exp(h*c/lamda/k/T)-1)
	return B*1e-6 # W/m2/um

if __name__ == '__main__':
	print(B(10, 280))
