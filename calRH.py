####################################################################################
# INTRODUCTION:
# This code is to do some calculation of water vapor in air
# Created by Hebs at 21/10/26/16:08
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np

def es(T):
	'''
	C-C function
	input:
		T   : temperature, float, Kelvin
	output:
		es  : saturated vapor pressure in T, float, hPa
	'''
	Lv = 40.8e3 # J/mol
	R = 8.31
	T0 = 273.15
	es0 = 6.082e2 # Pa
	
	es = es0 * np.exp(Lv/R*(1/T0-1/T)) / 100
	return es

def dpt2wh(t, dpt):
	'''
	This function is to use temperature and dew point temperature to calculate water mass concentration in air
	input:
		t      : temperature, K
		dpt    : dew point temperature, K
	output:
		wh     : water mass concentration, g/m^3
	'''
	M = 18 # g/mol
	R = 8.31
	e = es(dpt) # water vapor pressure, hPa
	wh = e*1e2 * M / R / t # g/m^3
	return wh

def wh2RH(t, wh):
	'''
	This function is to use temperature and water mass concentration to calculate raletive humidity RH
	input:
		t    : temperature, float, Kelvin
		wh   : water mass concentration, float, g/m^3
	output:
		RH   : raletive humidity, float, percent
	'''
	M = 18 # g/mol
	R = 8.31
	e = wh * R * t / M / 100 # hPa
	RH = e / es(t) * 100
	return RH

def q2rho(t, q, p):
	'''
	This function is to use water mixing ratio q to calculate density of water rho
	input:
		t    : temperature, float, Kelvin
		q    : water mixing ratio, float, kg/kg
		p    : pressure, float, hPa
	output:
		rho  : density of water, g/m^3
	'''
	M = 28.959 # g/mol, for air
	R = 8.31
	rho = p*1e2 * M * q / R / t
	return rho

if __name__ == '__main__':
	wh = dpt2wh(3.9+273.15, -7.3+273.15)
	print(wh)
	RH = wh2RH(3.9+273.15, wh)
	print(RH)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
