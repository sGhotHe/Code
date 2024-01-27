####################################################################################
# INTRODUCTION:
# This code is to test how different location change DARF in SBDART
# Created by Hebs at 21/11/8/15:28
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import write_INPUT
import write_aerosol
import write_albedo
import write_atms
import read11
import os
import date2jul

def run(location, iday, **args):
	'''
	This function is to run different location annual average DARF at 12:00 for SBDART
	input:
		location    : location for SBDART, dictionary
			name : location name, string
			alat : latitude of point on earth's surface, int
			alon : east longitude of point on earth's surface, int
		iday        : day in standard year, int, 1 to 365
		**fn        : INPUT path, string, default ./INPUT
		**fn_albedo : albedo.dat path, string, default ./albedo.dat
		**fn_atms   : atms.dat path, string, default ./atms.dat
		**output    : output path, string, default output/location/[name]/iout11_[iday].txt
	output:
		iout11.txt for each iday
	'''
	if 'fn' in args:
		fn = args['fn']
	else:
		fn = 'INPUT'
	if 'fn_albedo' in args:
		fn_albedo = args['fn_albedo']
		change_albedo_path = True
	else:
		change_albedo_path = False
	if 'fn_atms' in args:
		fn_atms = args['fn_atms']
		change_atms_path = True
	else:
		change_atms_path = False
	if 'output' in args:
		output = args['output']
	else:
		output = 'output/location'
	
	name = location['name']
	alat = location['alat']
	alon = location['alon']
	
	time = 12 - 8 # Beijing time to UTC time
	
	if os.path.exists(fn):
		os.system('mv '+fn+' '+fn+'_origin') # back up
	
	if change_albedo_path:
		os.system('mv albedo.dat albedo.dat_origin') # back up
		os.system('cp '+fn_albedo+' albedo.dat')
	
	if change_atms_path:
		os.system('mv atms.dat atms.dat_origin') # back up
		os.system('cp '+fn_atms+' atms.dat')
	
	path = output + '/' + name
	if not os.path.exists(path):
		os.system('mkdir '+path)
	write_INPUT.write(alat=alat, alon=alon, iday=iday, time=time)
	os.system('sbdart >'+path+'/iout11_'+str(iday)+'.txt')
	
	write_INPUT.change_back()
	if change_albedo_path:
		write_albedo.change_back()
	if change_atms_path:
		write_atms.change_back()

def read(location, **args):
	'''
	This function is to read different location annual average DARF at 12:00
	input:
		location : location for SBDART, dictionary
			name : location name, string
			alat : latitude of point on earth's surface, int
			alon : east longitude of point on earth's surface, int
		**path   : file path, string, default output/location
	output:
		nrf      : annual average net radiative flux, array in shape (len(name)), W/m^2
	'''
	if 'path' in args:
		path = args['path']
	else:
		path = 'output/location'
	
	name = location['name']
	
	iday = np.arange(1,366)
	nrf = np.zeros((len(name),len(iday)))
	
	for i in range(len(name)):
		for j in range(len(iday)):
			fn = output+'/'+name[i]+'/iout11_'+str(iday[j])+'.txt'
			iout11 = read11.read11(filename=fn)
			nrf[i,j] = iout11['fxdn'][0]-iout11['fxup'][0]
	
	nrf_avg = np.nanmean(nrf,axis=1)
	
	return nrf_avg

def run_albedo(rate, location):
	'''
	This function is to run albedo change in different location and annual average at 12:00
	input:
		rate     : albedo change rate, array
		location : location for SBDART, dictionary
			name : location name, string
			alat : latitude of point on earth's surface, int
			alon : east longitude of point on earth's surface, int
	output:
		output/location/albedo/[rate]/[name]/iout11_[iday].txt
	'''
	name = location['name']
	alat = location['alat']
	alon = location['alon']
	
	loc = []
	for i in range(len(name)):
		loc.append(dict(name=name[i],alat=alat[i],alon=alon[i]))
	iday = np.arange(1,366)
	
	print('running...')
	
	for i in range(len(loc)):
		name = loc[i]['name']
		path_albedo = 'input/albedo/' + name
		path_atms = 'input/atms/' + name
		for j in range(len(rate)):
			output = 'output/location/albedo/' + str(round(rate[j]*100)/100)
			if not os.path.exists(output):
				os.system('mkdir '+output)
			for k in range(len(iday)):
				fn_albedo = path_albedo + '/month_' + str(date2jul.iday2month(iday[k])) + '/albedo.dat'
				os.system('cp '+fn_albedo+'_'+str(round(rate[j]*100)/100)+' '+fn_albedo)
				fn_atms = path_atms + '/atms.dat'
				os.system('cp '+fn_atms+'_'+str(date2jul.iday2month(iday[k]))+' '+fn_atms)
				run(loc[i], iday[k], fn_albedo=fn_albedo, fn_atms=fn_atms, output=output)
				os.system('rm '+fn_albedo)
				os.system('rm '+fn_atms)
			print(round((j+1+i*len(rate))/len(rate)/len(loc)*1000)/10, '% done...')
	
	print('done')

def read_albedo(rate, location):
	'''
	This function is to read albedo change leaded DARF change in different location and annual average at 12:00
	input:
		rate     : albedo change rate, array
		location : location for SBDART, dictionary
			name : location name, string
			alat : latitude of point on earth's surface, int
			alon : east longitude of point on earth's surface, int
	output:
		nrf      : annual average net radiative flux, array in shape (len(name), len(rate)), W/m^2
	'''
	name = location['name']
	iday = np.arange(1,366)
	
	nrf = np.zeros((len(name),len(rate),len(iday)))
	
	for i in range(len(name)):
		for j in range(len(rate)):
			path = 'output/location/albedo/' + str(round(rate[j]*100)/100) + '/' + name[i]
			for k in range(len(iday)):
				fn = path + '/iout11_'+str(iday[k])+'.txt'
				iout11 = read11.read11(filename=fn)
				nrf[i,j,k] = iout11['fxdn'][0]-iout11['fxup'][0]
	
	nrf_avg = np.nanmean(nrf,axis=2)
	
	return nrf_avg

def read_albedo_parameter(rate, location, **args):
	'''
	This function is to read sequence changed albedo
	input:
		rate     : albedo change rate, array
		location : location for SBDART, dictionary
			name : location name, string
			alat : latitude of point on earth's surface, int
			alon : east longitude of point on earth's surface, int
		**path   : file path, string, default input/albedo
	output:
		wls      : wave length, array in shape(len(location), len(wl)), nm
		albedos  : changed albedo, array in shape (len(location), len(wl), len(rate))
	'''
	if 'path' in args:
		path = args['path']
	else:
		path = 'input/albedo'
	
	import runAlbedoChange
	name = location['name']
	wls = []
	albedos = []
	
	for i in range(len(name)):
		path_i = path + '/' + name[i]
		albedos_i = []
		for j in range(12): # 12 months
			fn = path_i + '/albedo.dat_' + str(j+1)
			wl, albedo = runAlbedoChange.read_parameter(rate, path=fn)
			albedos_i.append(albedo)
		wls.append(wl)
		albedos.append(np.nanmean(albedos_i,0))
	
	wls = np.array(wls)
	albedos = np.array(albedos)
	
	return wls, albedos

def run_AOD(rate, location):
	'''
	This function is to run AOD change in different location and annual average at 12:00
	input:
		rate     : AOD change rate, array
		location : location for SBDART, dictionary
			name : location name, string
			alat : latitude of point on earth's surface, int
			alon : east longitude of point on earth's surface, int
	output:
		output/location/AOD/[rate]/[name]/iout11_[iday].txt
	'''
	name = location['name']
	alat = location['alat']
	alon = location['alon']
	
	loc = []
	for i in range(len(name)):
		loc.append(dict(name=name[i],alat=alat[i],alon=alon[i]))
	
	iday = np.arange(1,366)
	time = 12 - 8 # Beijing time to UTC time
	
	print('calculating...')
	
	for i in range(len(loc)):
		name = loc[i]['name']
		path_albedo = 'input/albedo/' + name
		path_atms = 'input/atms/' + name
		for j in range(len(rate)):
			output = 'output/location/AOD/' + str(round(rate[j]*100)/100)
			if not os.path.exists(output):
				os.system('mkdir '+output)
			#change AOD
			write_aerosol.change_AOD(rate[j])
			for k in range(len(iday)):
				fn_albedo = path_albedo + '/month_' + str(date2jul.iday2month(iday[k])) + '/albedo.dat'
				os.system('cp '+fn_albedo+'_'+str(round(rate[j]*100)/100)+' '+fn_albedo)
				fn_atms = path_atms + '/atms.dat'
				os.system('cp '+fn_atms+'_'+str(date2jul.iday2month(iday[k]))+' '+fn_atms)
				run(loc[i], iday[k], fn_albedo=fn_albedo, fn_atms=fn_atms, output=output)
				os.system('rm '+fn_albedo)
				os.system('rm '+fn_atms)
			#change AOD back
			write_aerosol.change_back()
			print(round((j+1+i*len(rate))/len(rate)/len(loc)*1000)/10, '% done...')
	
	print('done')

def read_AOD(rate, location):
	'''
	This function is to read AOD change leaded DARF change in different location and annual average at 12:00
	input:
		rate     : AOD change rate, array
		location : location for SBDART, dictionary
			name : location name, string
			alat : latitude of point on earth's surface, int
			alon : east longitude of point on earth's surface, int
	output:
		nrf      : annual average net radiative flux, array in shape (len(name), len(rate)), W/m^2
	'''
	name = location['name']
	iday = np.arange(1,366)
	
	nrf = np.zeros((len(name),len(rate),len(iday)))
	
	for i in range(len(name)):
		for j in range(len(rate)):
			path = 'output/location/AOD/' + str(round(rate[j]*100)/100) + '/' + name[i]
			for k in range(len(iday)):
				fn = path + '/iout11_'+str(iday[k])+'.txt'
				iout11 = read11.read11(filename=fn)
				nrf[i,j,k] = iout11['fxdn'][0]-iout11['fxup'][0]
	
	nrf_avg = np.nanmean(nrf,axis=2)
	
	return nrf_avg

def run_g(rate, location):
	'''
	This function is to run g change in different location and annual average at 12:00
	input:
		rate     : g change rate, array
		location : location for SBDART, dictionary
			name : location name, string
			alat : latitude of point on earth's surface, int
			alon : east longitude of point on earth's surface, int
	output:
		output/location/g/[rate]/[name]/iout11_[iday].txt
	'''
	name = location['name']
	alat = location['alat']
	alon = location['alon']
	
	loc = []
	for i in range(len(name)):
		loc.append(dict(name=name[i],alat=alat[i],alon=alon[i]))
	
	iday = np.arange(1,366)
	time = 12 - 8 # Beijing time to UTC time
	
	print('calculating...')
	
	for i in range(len(loc)):
		name = loc[i]['name']
		path_albedo = 'input/albedo/' + name
		path_atms = 'input/atms/' + name
		for j in range(len(rate)):
			output = 'output/location/g/' + str(round(rate[j]*100)/100)
			if not os.path.exists(output):
				os.system('mkdir '+output)
			#change g
			os.system('mv aerosol.dat aerosol.dat_origin')
			os.system('cp input/g/aerosol.dat_'+str(round(rate[j]*100)/100)+' aerosol.dat')
			for k in range(len(iday)):
				fn_albedo = path_albedo + '/month_' + str(date2jul.iday2month(iday[k])) + '/albedo.dat'
				os.system('cp '+fn_albedo+'_'+str(round(rate[j]*100)/100)+' '+fn_albedo)
				fn_atms = path_atms + '/atms.dat'
				os.system('cp '+fn_atms+'_'+str(date2jul.iday2month(iday[k]))+' '+fn_atms)
				run(loc[i], iday[k], fn_albedo=fn_albedo, fn_atms=fn_atms, output=output)
				os.system('rm '+fn_albedo)
				os.system('rm '+fn_atms)
			#change g back
			write_aerosol.change_back()
			print(round((j+1+i*len(rate))/len(rate)/len(loc)*1000)/10, '% done...')
	
	print('done')

def read_g(rate, location):
	'''
	This function is to read g change leaded DARF change in different location and annual average at 12:00
	input:
		rate     : g change rate, array
		location : location for SBDART, dictionary
			name : location name, string
			alat : latitude of point on earth's surface, int
			alon : east longitude of point on earth's surface, int
	output:
		nrf      : annual average net radiative flux, array in shape (len(name), len(rate)), W/m^2
	'''
	name = location['name']
	iday = np.arange(1,366)
	
	nrf = np.zeros((len(name),len(rate),len(iday)))
	
	for i in range(len(name)):
		for j in range(len(rate)):
			path = 'output/location/g/' + str(round(rate[j]*100)/100) + '/' + name[i]
			for k in range(len(iday)):
				fn = path + '/iout11_'+str(iday[k])+'.txt'
				iout11 = read11.read11(filename=fn)
				nrf[i,j,k] = iout11['fxdn'][0]-iout11['fxup'][0]
	
	nrf_avg = np.nanmean(nrf,axis=2)
	
	return nrf_avg

def run_k(rate, location):
	'''
	This function is to run k change in different location and annual average at 12:00
	input:
		rate     : k change rate, array
		location : location for SBDART, dictionary
			name : location name, string
			alat : latitude of point on earth's surface, int
			alon : east longitude of point on earth's surface, int
	output:
		output/location/k/[rate]/[name]/iout11_[iday].txt
	'''
	name = location['name']
	alat = location['alat']
	alon = location['alon']
	
	loc = []
	for i in range(len(name)):
		loc.append(dict(name=name[i],alat=alat[i],alon=alon[i]))
	
	iday = np.arange(1,366)
	time = 12 - 8 # Beijing time to UTC time
	
	print('calculating...')
	
	for i in range(len(loc)):
		name = loc[i]['name']
		path_albedo = 'input/albedo/' + name
		path_atms = 'input/atms/' + name
		for j in range(len(rate)):
			output = 'output/location/k/' + str(round(rate[j]*100)/100)
			if not os.path.exists(output):
				os.system('mkdir '+output)
			#change k
			os.system('mv aerosol.dat aerosol.dat_origin')
			os.system('cp input/k/aerosol.dat_'+str(round(rate[j]*100)/100)+' aerosol.dat')
			for k in range(len(iday)):
				fn_albedo = path_albedo + '/month_' + str(date2jul.iday2month(iday[k])) + '/albedo.dat'
				os.system('cp '+fn_albedo+'_'+str(round(rate[j]*100)/100)+' '+fn_albedo)
				fn_atms = path_atms + '/atms.dat'
				os.system('cp '+fn_atms+'_'+str(date2jul.iday2month(iday[k]))+' '+fn_atms)
				run(loc[i], iday[k], fn_albedo=fn_albedo, fn_atms=fn_atms, output=output)
				os.system('rm '+fn_albedo)
				os.system('rm '+fn_atms)
			#change k back
			write_aerosol.change_back()
			print(round((j+1+i*len(rate))/len(rate)/len(loc)*1000)/10, '% done...')
	
	print('done')

def read_k(rate, location):
	'''
	This function is to read g change leaded DARF change in different location and annual average at 12:00
	input:
		rate     : g change rate, array
		location : location for SBDART, dictionary
			name : location name, string
			alat : latitude of point on earth's surface, int
			alon : east longitude of point on earth's surface, int
	output:
		nrf      : annual average net radiative flux, array in shape (len(name), len(rate)), W/m^2
	'''
	name = location['name']
	iday = np.arange(1,366)
	
	nrf = np.zeros((len(name),len(rate),len(iday)))
	
	for i in range(len(name)):
		for j in range(len(rate)):
			path = 'output/location/k/' + str(round(rate[j]*100)/100) + '/' + name[i]
			for k in range(len(iday)):
				fn = path + '/iout11_'+str(iday[k])+'.txt'
				iout11 = read11.read11(filename=fn)
				nrf[i,j,k] = iout11['fxdn'][0]-iout11['fxup'][0]
	
	nrf_avg = np.nanmean(nrf,axis=2)
	
	return nrf_avg

def run_kappa(rate, location):
	'''
	This function is to run kappa change in different location and annual average at 12:00
	input:
		rate     : kappa change rate, array
		location : location for SBDART, dictionary
			name : location name, string
			alat : latitude of point on earth's surface, int
			alon : east longitude of point on earth's surface, int
	output:
		output/location/kappa/[rate]/[name]/iout11_[iday].txt
	'''
	name = location['name']
	alat = location['alat']
	alon = location['alon']
	
	loc = []
	for i in range(len(name)):
		loc.append(dict(name=name[i],alat=alat[i],alon=alon[i]))
	
	iday = np.arange(1,366)
	time = 12 - 8 # Beijing time to UTC time
	
	print('calculating...')
	
	for i in range(len(loc)):
		name = loc[i]['name']
		path_albedo = 'input/albedo/' + name
		path_atms = 'input/atms/' + name
		for j in range(len(rate)):
			output = 'output/location/kappa/' + str(round(rate[j]*100)/100)
			if not os.path.exists(output):
				os.system('mkdir '+output)
			#change kappa
			os.system('mv aerosol.dat aerosol.dat_origin')
			os.system('cp input/kappa/aerosol.dat_'+str(round(rate[j]*100)/100)+' aerosol.dat')
			for k in range(len(iday)):
				fn_albedo = path_albedo + '/month_' + str(date2jul.iday2month(iday[k])) + '/albedo.dat'
				os.system('cp '+fn_albedo+'_'+str(round(rate[j]*100)/100)+' '+fn_albedo)
				fn_atms = path_atms + '/atms.dat'
				os.system('cp '+fn_atms+'_'+str(date2jul.iday2month(iday[k]))+' '+fn_atms)
				run(loc[i], iday[k], fn_albedo=fn_albedo, fn_atms=fn_atms, output=output)
				os.system('rm '+fn_albedo)
				os.system('rm '+fn_atms)
			#change kappa back
			write_aerosol.change_back()
			print(round((j+1+i*len(rate))/len(rate)/len(loc)*1000)/10, '% done...')
	
	print('done')

def read_kappa(rate, location):
	'''
	This function is to read kappa change leaded DARF change in different location and annual average at 12:00
	input:
		rate     : kappa change rate, array
		location : location for SBDART, dictionary
			name : location name, string
			alat : latitude of point on earth's surface, int
			alon : east longitude of point on earth's surface, int
	output:
		nrf      : annual average net radiative flux, array in shape (len(name), len(rate)), W/m^2
	'''
	name = location['name']
	iday = np.arange(1,366)
	
	nrf = np.zeros((len(name),len(rate),len(iday)))
	
	for i in range(len(name)):
		for j in range(len(rate)):
			path = 'output/location/kappa/' + str(round(rate[j]*100)/100) + '/' + name[i]
			for k in range(len(iday)):
				fn = path + '/iout11_'+str(iday[k])+'.txt'
				iout11 = read11.read11(filename=fn)
				nrf[i,j,k] = iout11['fxdn'][0]-iout11['fxup'][0]
	
	nrf_avg = np.nanmean(nrf,axis=2)
	
	return nrf_avg

def run_mshell(rate, location):
	'''
	This function is to run mshell change in different location and annual average at 12:00
	input:
		rate     : mshell change rate, array
		location : location for SBDART, dictionary
			name : location name, string
			alat : latitude of point on earth's surface, int
			alon : east longitude of point on earth's surface, int
	output:
		output/location/mshell/[rate]/[name]/iout11_[iday].txt
	'''
	name = location['name']
	alat = location['alat']
	alon = location['alon']
	
	loc = []
	for i in range(len(name)):
		loc.append(dict(name=name[i],alat=alat[i],alon=alon[i]))
	
	iday = np.arange(1,366)
	time = 12 - 8 # Beijing time to UTC time
	
	print('calculating...')
	
	for i in range(len(loc)):
		name = loc[i]['name']
		path_albedo = 'input/albedo/' + name
		path_atms = 'input/atms/' + name
		for j in range(len(rate)):
			output = 'output/location/mshell/' + str(round(rate[j]*100)/100)
			if not os.path.exists(output):
				os.system('mkdir '+output)
			#change mshell
			os.system('mv aerosol.dat aerosol.dat_origin')
			os.system('cp input/mshell/aerosol.dat_'+str(round(rate[j]*100)/100)+' aerosol.dat')
			for k in range(len(iday)):
				fn_albedo = path_albedo + '/month_' + str(date2jul.iday2month(iday[k])) + '/albedo.dat'
				os.system('cp '+fn_albedo+'_'+str(round(rate[j]*100)/100)+' '+fn_albedo)
				fn_atms = path_atms + '/atms.dat'
				os.system('cp '+fn_atms+'_'+str(date2jul.iday2month(iday[k]))+' '+fn_atms)
				run(loc[i], iday[k], fn_albedo=fn_albedo, fn_atms=fn_atms, output=output)
				os.system('rm '+fn_albedo)
				os.system('rm '+fn_atms)
			#change mshell back
			write_aerosol.change_back()
			print(round((j+1+i*len(rate))/len(rate)/len(loc)*1000)/10, '% done...')
	
	print('done')

def read_mshell(rate, location):
	'''
	This function is to read mshell change leaded DARF change in different location and annual average at 12:00
	input:
		rate     : mshell change rate, array
		location : location for SBDART, dictionary
			name : location name, string
			alat : latitude of point on earth's surface, int
			alon : east longitude of point on earth's surface, int
	output:
		nrf      : annual average net radiative flux, array in shape (len(name), len(rate)), W/m^2
	'''
	name = location['name']
	iday = np.arange(1,366)
	
	nrf = np.zeros((len(name),len(rate),len(iday)))
	
	for i in range(len(name)):
		for j in range(len(rate)):
			path = 'output/location/mshell/' + str(round(rate[j]*100)/100) + '/' + name[i]
			for k in range(len(iday)):
				fn = path + '/iout11_'+str(iday[k])+'.txt'
				iout11 = read11.read11(filename=fn)
				nrf[i,j,k] = iout11['fxdn'][0]-iout11['fxup'][0]
	
	nrf_avg = np.nanmean(nrf,axis=2)
	
	return nrf_avg

def run_n(rate, location):
	'''
	This function is to run n change in different location and annual average at 12:00
	input:
		rate     : n change rate, array
		location : location for SBDART, dictionary
			name : location name, string
			alat : latitude of point on earth's surface, int
			alon : east longitude of point on earth's surface, int
	output:
		output/location/n/[rate]/[name]/iout11_[iday].txt
	'''
	name = location['name']
	alat = location['alat']
	alon = location['alon']
	
	loc = []
	for i in range(len(name)):
		loc.append(dict(name=name[i],alat=alat[i],alon=alon[i]))
	
	iday = np.arange(1,366)
	time = 12 - 8 # Beijing time to UTC time
	
	print('calculating...')
	
	for i in range(len(loc)):
		name = loc[i]['name']
		path_albedo = 'input/albedo/' + name
		path_atms = 'input/atms/' + name
		for j in range(len(rate)):
			output = 'output/location/n/' + str(round(rate[j]*100)/100)
			if not os.path.exists(output):
				os.system('mkdir '+output)
			#change n
			os.system('mv aerosol.dat aerosol.dat_origin')
			os.system('cp input/n/aerosol.dat_'+str(round(rate[j]*100)/100)+' aerosol.dat')
			for k in range(len(iday)):
				fn_albedo = path_albedo + '/month_' + str(date2jul.iday2month(iday[k])) + '/albedo.dat'
				os.system('cp '+fn_albedo+'_'+str(round(rate[j]*100)/100)+' '+fn_albedo)
				fn_atms = path_atms + '/atms.dat'
				os.system('cp '+fn_atms+'_'+str(date2jul.iday2month(iday[k]))+' '+fn_atms)
				run(loc[i], iday[k], fn_albedo=fn_albedo, fn_atms=fn_atms, output=output)
				os.system('rm '+fn_albedo)
				os.system('rm '+fn_atms)
			#change n back
			write_aerosol.change_back()
			print(round((j+1+i*len(rate))/len(rate)/len(loc)*1000)/10, '% done...')
	
	print('done')

def read_n(rate, location):
	'''
	This function is to read n change leaded DARF change in different location and annual average at 12:00
	input:
		rate     : n change rate, array
		location : location for SBDART, dictionary
			name : location name, string
			alat : latitude of point on earth's surface, int
			alon : east longitude of point on earth's surface, int
	output:
		nrf      : annual average net radiative flux, array in shape (len(name), len(rate)), W/m^2
	'''
	name = location['name']
	iday = np.arange(1,366)
	
	nrf = np.zeros((len(name),len(rate),len(iday)))
	
	for i in range(len(name)):
		for j in range(len(rate)):
			path = 'output/location/n/' + str(round(rate[j]*100)/100) + '/' + name[i]
			for k in range(len(iday)):
				fn = path + '/iout11_'+str(iday[k])+'.txt'
				iout11 = read11.read11(filename=fn)
				nrf[i,j,k] = iout11['fxdn'][0]-iout11['fxup'][0]
	
	nrf_avg = np.nanmean(nrf,axis=2)
	
	return nrf_avg

def run_SSA(rate, location):
	'''
	This function is to run SSA change in different location and annual average at 12:00
	input:
		rate     : SSA change rate, array
		location : location for SBDART, dictionary
			name : location name, string
			alat : latitude of point on earth's surface, int
			alon : east longitude of point on earth's surface, int
	output:
		output/location/SSA/[rate]/[name]/iout11_[iday].txt
	'''
	name = location['name']
	alat = location['alat']
	alon = location['alon']
	
	loc = []
	for i in range(len(name)):
		loc.append(dict(name=name[i],alat=alat[i],alon=alon[i]))
	
	iday = np.arange(1,366)
	time = 12 - 8 # Beijing time to UTC time
	
	print('calculating...')
	
	for i in range(len(loc)):
		name = loc[i]['name']
		path_albedo = 'input/albedo/' + name
		path_atms = 'input/atms/' + name
		for j in range(len(rate)):
			output = 'output/location/SSA/' + str(round(rate[j]*100)/100)
			if not os.path.exists(output):
				os.system('mkdir '+output)
			#change SSA
			write_aerosol.change_SSA(rate[j])
			for k in range(len(iday)):
				fn_albedo = path_albedo + '/month_' + str(date2jul.iday2month(iday[k])) + '/albedo.dat'
				os.system('cp '+fn_albedo+'_'+str(round(rate[j]*100)/100)+' '+fn_albedo)
				fn_atms = path_atms + '/atms.dat'
				os.system('cp '+fn_atms+'_'+str(date2jul.iday2month(iday[k]))+' '+fn_atms)
				run(loc[i], iday[k], fn_albedo=fn_albedo, fn_atms=fn_atms, output=output)
				os.system('rm '+fn_albedo)
				os.system('rm '+fn_atms)
			#change SSA back
			write_aerosol.change_back()
			print(round((j+1+i*len(rate))/len(rate)/len(loc)*1000)/10, '% done...')
	
	print('done')

def read_SSA(rate, location):
	'''
	This function is to read SSA change leaded DARF change in different location and annual average at 12:00
	input:
		rate     : SSA change rate, array
		location : location for SBDART, dictionary
			name : location name, string
			alat : latitude of point on earth's surface, int
			alon : east longitude of point on earth's surface, int
	output:
		nrf      : annual average net radiative flux, array in shape (len(name), len(rate)), W/m^2
	'''
	name = location['name']
	iday = np.arange(1,366)
	
	nrf = np.zeros((len(name),len(rate),len(iday)))
	
	for i in range(len(name)):
		for j in range(len(rate)):
			path = 'output/location/SSA/' + str(round(rate[j]*100)/100) + '/' + name[i]
			for k in range(len(iday)):
				fn = path + '/iout11_'+str(iday[k])+'.txt'
				iout11 = read11.read11(filename=fn)
				nrf[i,j,k] = iout11['fxdn'][0]-iout11['fxup'][0]
	
	nrf_avg = np.nanmean(nrf,axis=2)
	
	return nrf_avg

def run_zero(location):
	'''
	This function is to run zero in different location and annual average at 12:00
	input:
		location : location for SBDART, dictionary
			name : location name, string
			alat : latitude of point on earth's surface, int
			alon : east longitude of point on earth's surface, int
	output:
		output/location/zero/[name]/iout11_[month].txt
	'''
	name = location['name']
	alat = location['alat']
	alon = location['alon']
	
	loc = []
	for i in range(len(name)):
		loc.append(dict(name=name[i],alat=alat[i],alon=alon[i]))
	iday = np.arange(1,366)
	
	os.system('mv aerosol.dat aerosol.dat_origin')
	os.system('cp input/zero/aerosol.dat aerosol.dat')
	
	for i in range(len(loc)):
		name = loc[i]['name']
		path_albedo = 'input/albedo/' + name
		path_atms = 'input/atms/' + name
		output = 'output/location/zero'
		if not os.path.exists(output):
			os.system('mkdir '+output)
		for j in range(len(iday)):
			fn_albedo = path_albedo + '/albedo.dat'
			os.system('cp '+fn_albedo+'_'+str(date2jul.iday2month(iday[j]))+' '+fn_albedo)
			fn_atms = path_atms + '/atms.dat'
			os.system('cp '+fn_atms+'_'+str(date2jul.iday2month(iday[j]))+' '+fn_atms)
			run(loc[i], iday[j], fn_albedo=fn_albedo, fn_atms=fn_atms, output=output)
			os.system('rm '+fn_albedo)
			os.system('rm '+fn_atms)
	
	write_aerosol.change_back()

def read_zero(location):
	'''
	This function is to read zero DARF in different location and annual average at 12:00
	input:
		location : location for SBDART, dictionary
			name : location name, string
			alat : latitude of point on earth's surface, int
			alon : east longitude of point on earth's surface, int
	output:
		nrf      : annual average net radiative flux, array in shape (len(name)), W/m^2
	'''
	name = location['name']
	iday = np.arange(1,366)
	
	nrf = np.zeros((len(name), len(iday)))
	
	for i in range(len(name)):
		path = 'output/location/zero/' + name[i]
		for j in range(len(iday)):
			fn = path + '/iout11_'+str(iday[j])+'.txt'
			iout11 = read11.read11(filename=fn)
			nrf[i,j] = iout11['fxdn'][0]-iout11['fxup'][0]
	
	nrf_avg = np.nanmean(nrf,axis=1)
	
	return nrf_avg

if __name__ == '__main__':
	name = ['Beijing']
	alat = [40]
	alon = [116]
	location = dict(name=name, alat=alat, alon=alon)
	
	# albedo change from -10% to 10%, bin in 1%
	rate = np.arange(-10,11) / 100 + 1
	
	run_albedo(rate, location)
	run_AOD(rate, location)
	run_g(rate, location)
	run_k(rate, location)
	run_kappa(rate, location)
	run_mshell(rate, location)
	run_n(rate, location)
	run_SSA(rate, location)
	'''
	nrf = read_albedo(rate, location)
	print(nrf[0]) # for Beijing
	albedo_wls, albedos = read_albedo_parameter(rate, location)
	print(albedo_wls[0])
	print(albedos[0,0])
	nrf = read_AOD(rate, location)
	print(nrf[0])
	nrf = read_g(rate, location)
	print(nrf[0])
	nrf = read_k(rate, location)
	print(nrf[0])
	nrf = read_kappa(rate, location)
	print(nrf[0])
	nrf = read_mshell(rate, location)
	print(nrf[0])
	nrf = read_n(rate, location)
	print(nrf[0])
	nrf = read_SSA(rate, location)
	print(nrf[0])
	
	run_zero(location)
	zero_nrf = read_zero(location)
	print(zero_nrf)
	'''
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
