####################################################################################
# INTRODUCTION:
# This code is to do FAST
# FAST is Fourier Amplitude Sensitivity Test, a way to do sensitivity analyse
# Created by Hebs at 23/5/8/13:00
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import random
from scipy.stats import norm
import matplotlib.pyplot as plt

def EFAST_norm(num, loc, scale, freq):
	'''
	This function is extended Fourier Amplitude Sensitivity Test method code for normal distribution
	input:
		num     : sample number, int
		loc     : normal distribution location parameter, float
		scale   : normal distribution scale parameter, float
		freq    : sin signal frequency, float
	output:
		x       : parameter list, array
	'''
	x = np.arange(num*1.0)
	
	phi = random.uniform(-np.pi, np.pi)
	for i in range(len(x)):
		omega = 2 * np.pi * freq / num
		x[i] = 0.5 + 1 / np.pi * np.arcsin(np.sin(omega*x[i]+phi))
		x[i] = norm(loc=loc, scale=scale).ppf(x[i])
	return x

def cal_FAST():
	'''
	This function is to calculate FAST
	input:
	output:
	'''
	# to do FAST, need to do: 1st, use search curve to get parameters; 2nd, run function;
	# 3rd, use output and parameters' frequence to calculate FAST
	x = EFAST_norm(100, 1, 1, 30)
	print(x)
	plt.plot(np.arange(100),x)
	plt.show()

if __name__ == '__main__':
	cal_FAST()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
