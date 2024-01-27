####################################################################################
# INTRODUCTION:
# This code is to read AERONET monthly average AOD data (.lev20)
# Created by Hebs at 21/5/24/13:34
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import re
import phaseFunc

def clean(li):
	i = 0
	while i<len(li):
		if li[i]=='':
			li.pop(i)
		else:
			i = i + 1
	return li

def readfile(fn):
	lines = []
	with open(fn, 'r') as f:
		for res in f.readlines():
			res = res.rstrip('\n')
			res = re.split(',', res)
			res = clean(res)
			lines.append(res)
	return lines

def readAERONET(path):
	fn_aod = path + '20190101_20191231_Beijing.aod'
	fn_ssa = path + '20190101_20191231_Beijing.ssa'
	fn_asy = path + '20190101_20191231_Beijing.asy'
	fn_pfn = path + '20190101_20191231_Beijing.pfn'
	aod_infos = readfile(fn_aod)[7:19] # total 12 month data
	ssa_infos = readfile(fn_ssa)[7:19]
	asy_infos = readfile(fn_asy)[7:19]
	pfn_infos = readfile(fn_pfn)[7:19]
	
	month = np.arange(1,13)
	wl = np.array([440, 675, 870, 1020])
	angles = np.array([180.000000, 178.290000, 176.070000, 173.840000, 171.610000, 169.370000, 167.140000, 164.900000, 162.670000, 160.430000, 158.200000, 155.960000, 153.720000, 151.490000, 149.250000, 147.020000, 144.780000, 142.550000, 140.310000, 138.070000, 135.840000, 133.600000, 131.370000, 129.130000, 126.890000, 124.660000, 122.420000, 120.190000, 117.950000, 115.710000, 113.480000, 111.240000, 109.010000, 106.770000, 104.530000, 102.300000, 100.060000, 97.830000, 95.590000, 93.350000, 91.120000, 90.000000, 88.880000, 86.650000, 84.410000, 82.170000, 79.940000, 77.700000, 75.470000, 73.230000, 70.990000, 68.760000, 66.520000, 64.290000, 62.050000, 59.810000, 57.580000, 55.340000, 53.110000, 50.870000, 48.630000, 46.400000, 44.160000, 41.930000, 39.690000, 37.450000, 35.220000, 32.980000, 30.750000, 28.510000, 26.280000, 24.040000, 21.800000, 19.570000, 17.330000, 15.100000, 12.860000, 10.630000, 8.390000, 6.160000, 3.930000, 1.710000, 0.000000])
	
	aod = np.zeros((len(month),len(wl)))
	ssa = np.zeros((len(month),len(wl)))
	asy = np.zeros((len(month),len(wl)))
	for i in range(len(month)):
		for j in range(len(wl)):
			aod[i,j] = float(aod_infos[i][j+1])
			ssa[i,j] = float(ssa_infos[i][j+1])
			asy[i,j] = float(asy_infos[i][j+1])
	pfn = np.zeros((len(month),len(wl),len(angles)))
	for i in range(len(month)):
		for j in range(len(wl)):
			for k in range(len(angles)):
				pfn[i,j,k] = float(pfn_infos[i][j*len(wl)+k+1])
	
	aeronet = dict(month=month, wl=wl, aod=aod, ssa=ssa, asy=asy, angles=angles, pfn=pfn)
	return aeronet

if __name__ == '__main__':
	aeronet = readAERONET('data/aeronet/')
	aod = aeronet['aod']
	print(aod)
	angles = aeronet['angles']
	f = aeronet['pfn']
	thetas = phaseFunc.angle2theta(angles)
	nstr = 6
	beta = np.zeros((3,4,nstr))
	for i in range(3):
		for j in range(4):
			for k in range(nstr):
				beta[i,j,k] = phaseFunc.beta(f[i,j], thetas, k+1)
	print(beta)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
