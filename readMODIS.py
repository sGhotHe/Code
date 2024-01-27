####################################################################################
# INTRODUCTION:
# This code is to read MODIS data, including MCD43C3 for surface albedo and MOD04_L2 (Terra) or MYD04_L2 (Aqua) for AOD
# Created by Hebs at 21/5/24/18:06
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
from pyhdf.SD import SD
from date2jul import date2jul,jul2date
import sys, os

def read_modis_albedo(lon, lat, year, **args):
	'''
	This function is to read MODIS data from Terra/Aqua for MCD43C3
	input:
		lon    : longitude
		lat    : latitude
		year   : year
		**args :
			month    : month
			day      : day of month
			jul      : julian day
			path     : file path
			band     : specific spectral band, default 4 for 0.55um
			strcase  : to determine use white sky albedo (WSA) or black sky albedo (BSA), default is BSA
	output:
		albedo : surface albedo, %
		wv     : wavelength, um
	'''
	if 'month' in args:
		month = args['month']
	if 'day' in args:
		day = args['day']
	if 'jul' in args:
		jul = args['jul']
		times = jul2date(year,jul)
		month = times[1]
		day = times[2]
	if ('jul' not in args) and (('month' not in args) or ('day' not in args)) :
		print('The time variable is not appropriate, please check')
		sys.exit()
	if 'path' in args:
		path=args['path']
	else:
		path='/home/hebs/Code/SBDART/data/modis/albedo/'
	if 'band' in args:
		band = args['band']
	else:
		band = 4
	if 'white' in args:
		strcase = 'WSA'
	else:
		strcase = 'BSA'
	
	strbands = [
		"Band1",
		"Band2",
		"Band3",
		"Band4",
		"Band5",
		"Band6",
		"Band7",
		"vis",
		"nir",
		"shortwave"
	]
	bandwvs = [
		"0.67",
		"0.86",
		"0.47",
		"0.55",
		"1.24",
		"1.64",
		"2.1",
		"0.3-0.7",
		"0.7-5",
		"0.3-5",
	]
	jul = round(date2jul(year, month, day))
	strjul = str(jul).rjust(3, '0') # in form of 001, 002, ..., 999
	strband = strbands[band-1]
	varname = 'Albedo_' + strcase + "_" + strband
	flist=os.listdir(path)
	
	def fileselect(ls, start):
		ls_new = []
		for i in range(len(ls)):
			if ls[i].startswith(start):
				ls_new.append(ls[i])
		return ls_new
	
	flist = fileselect(flist, 'MCD43C3.A'+str(year)+strjul)
	albedo = 32767
	for fn in flist:
		try:
			hdf = SD(path+fn)
		except:
			print(('Corresponding file does not Exist'))
			return np.nan, np.nan
		albedos = hdf.select(varname).get()
		lat_index = round(abs(lat - 90.0) * 20.0)
		lon_index = round((lon + 180.0) * 20.0)
		albedo = albedos[lat_index][lon_index]
		if albedo<32767:
			break
	
	if albedo==32767: # Fill value
		return np.nan, bandwvs[band-1]
	return albedo*0.1, bandwvs[band-1]

def read_modis_aod(lon, lat, year, **args):
	'''
	This function is to read MODIS data from Terra for MOD04_L2 or Aqua for MCD04_L2
	input:
		lon    : longitude
		lat    : latitude
		year   : year
		**args :
			month     : month
			day       : day of month
			jul       : julian day
			path      : file path
			satellite : 0 for Terra, 1 for Aqua, default 0
			form      : 0 for no distance output, 1 for both aod and distance output, default 0
	output:
		aod      : aerosol optical depth
		distance : distance from data point to given location
	'''
	if 'month' in args:
		month = args['month']
	if 'day' in args:
		day = args['day']
	if 'jul' in args:
		jul = args['jul']
		times = jul2date(year,jul)
		month = times[1]
		day = times[2]
	if ('jul' not in args) and (('month' not in args) or ('day' not in args)) :
		print('The time variable is not appropriate, please check')
		sys.exit()
	if 'path' in args:
		path=args['path']
	else:
		path='/home/hebs/Desktop/SBDART/data/modis/aod/'
	if 'satellite' in args:
		satellite=args['satellite']
	else:
		satellite=0
	if 'form' in args:
		form = args['form']
	else:
		form = 0
	
	products = ['MOD04', 'MYD04']
	satellites = ['Terra', 'Aqua']
	strmodis = products[satellite]
	strsatellite = satellites[satellite]
	jul = round(date2jul(year, month, day))
	strjul = str(jul).rjust(3, '0') # in form of 001, 002, ..., 999
	varname = 'Deep_Blue_Aerosol_Optical_Depth_550_Land_Best_Estimate'
	flist=os.listdir(path)
	
	def fileselect(ls, start):
		ls_new = []
		for i in range(len(ls)):
			if ls[i].startswith(start):
				ls_new.append(ls[i])
		return ls_new
	
	flist = fileselect(flist, strmodis+'_L2.A'+str(year)+strjul)
	
	delta_lon = 0.5 # 需要的经纬度和数据经纬度的最小差
	delta_lat = 0.5
	lond = 6371 * np.cos(np.pi/180*lat) * (np.pi/180) # longitude distance per degree
	latd = 6371 * np.pi / 180 # latitude distance per degree
	
	exist = False
	for fn in flist:
		try:
			hdf=SD(path+fn)
		except:
			print(('Corresponding '+strsatellite+' File does not Exist'))
			if form:
				return np.nan, np.nan
			return np.nan
		
		lons = hdf.select('Longitude').get()
		lats = hdf.select('Latitude').get()
		lonr = [lons.min(),lons.max()] # data longitude range
		latr = [lats.min(),lats.max()] # data latitude range
		
		if lonr[0]<=lon and lonr[1]>=lon and latr[0]<=lat and latr[1]>=lat:
			exist = True
			aods = hdf.select(varname).get()
			aod = aods[(abs(lons-lon)<=delta_lon)&(abs(lats-lat)<=delta_lat)] * 0.001
			if len(aod)==0: # no data matched
				aod = np.nan
				distance = np.nan
			else:
				lons1 = lons[(abs(lons-lon)<=delta_lon)&(abs(lats-lat)<=delta_lat)] - lon
				lats1 = lats[(abs(lons-lon)<=delta_lon)&(abs(lats-lat)<=delta_lat)] - lat
				distances = np.sqrt((lons1*lond)**2+(lats1*latd)**2)
				aod_sorted = aod[distances.argsort()]
				distances=distances[distances.argsort()]
				if len(distances)==1:
					distance = distances[0]
					aod = aod_sorted[0]
					if aod<0:
						aod = -9999.
				else:
					for distance,aod in zip(distances,aod_sorted):
						if aod>0:
							break
					if aod<0:
						aod = -9999.
						distance = distances[0]
	
	if exist==False or aod<0:
		aod = np.nan
		distance = np.nan
	
	if form:
		return aod, distance
	return aod

if __name__ == '__main__':
	#aod = read_modis_aod(116, 39, 2019, month=1, day=2, satellite=0, form=0)
	#print(aod)
	for i in range(1,8):
		albedo = read_modis_albedo(116, 39, 2019, month=1, day=1, band=i)
		print(albedo)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
