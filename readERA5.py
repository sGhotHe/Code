####################################################################################
# INTRODUCTION:
# This code is to read ERA5 data
# Created by Hebs at 21/5/24/15:45
# Contact: hebishuo@pku.edu.cn
####################################################################################

import netCDF4 as nc
import numpy as np
import calRH

def readERA5(fn, lat, lon):
	'''
	This function is to read data from ERA5
	input:
		fn    : file path, string
		lat   : latitude, int
		lon   : longitude, int
	output:
		era5  : ERA5 data dictionary
			z  is the layer altitude in km (z must be monotonically decreasing)
			p  is the pressure in millibars
			t  is the temperature is Kelvin
			wh is water vapor density g/m3
			wo is ozone density g/m3
	'''
	data = nc.Dataset(fn)
	p = data.variables['level'][:] # hPa
	z = data.variables['z'] # (month, level, latitude, longitude), m^2/s^2, g=9.80665m/s^2
	o3 = data.variables['o3']
	q = data.variables['q']
	t = data.variables['t']
	lats = data.variables['latitude'] # 721
	lons = data.variables['longitude'] # 1140
	
	# lat_BJ = 41
	# lat_BJ_sub = 196
	# lon_BJ = 116
	# lon_BJ_sub = 464
	lat_sub = round((len(lats)-1) / 2 - lat / 90 * (len(lats)-1)/2)
	lon_sub = round(lon / 360 * len(lons))
	
	month = np.arange(1,13)
	z = z[:,:,lat_sub,lon_sub] / 9.80665 / 1000 # in km
	o3 = o3[:,:,lat_sub,lon_sub]
	q = q[:,:,lat_sub,lon_sub]
	t = t[:,:,lat_sub,lon_sub]
	
	wh = np.zeros((12,p.shape[0]))
	wo = np.zeros((12,p.shape[0]))
	for i in range(12):
		for j in range(p.shape[0]):
			wh[i,j] = calRH.q2rho(t[i,j], q[i,j], p[j])
			wo[i,j] = calRH.q2rho(t[i,j], o3[i,j], p[j])
	
	era5 = dict(month=month, z=z, p=p, t=t, wh=wh, wo=wo)
	return era5

if __name__ == '__main__':
	era5 = readERA5('data/era5/data.nc', 26, 126)
	print(era5['z'][0])
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
