####################################################################################
# INTRODUCTION:
# This code is to read fRH data
# Created by Hebs at 21/11/9/10:58
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import re
import date2jul
import calRH
import getfilesname
from scipy.optimize import curve_fit
import sys

def readfRH(fn):
	'''
	This function is to read one cycle of data points of humidified nephelometer system
	input:
		fn       :   the filename
	output:
		fRH      :   a dictionary which stores data points of this cycle,
			fRH = {year, jul, date, time, wv, SCA, BSCA, state, ST, Tset,
		                    Ttrue, T1, RH1, T2, RH2, RH, angstrom}
	'''
	f = open(fn)
	lines = f.readlines()
	NP = len(lines)
	juls = np.zeros(NP)
	dates = []
	times = []
	SCAs = np.zeros((NP, 3, 2))
	BSCAs = np.zeros((NP, 3, 2))
	wls = np.array((635, 525, 450))
	STs = np.zeros(NP)
	states = []
	Tsets = np.zeros(NP)
	Ttrues = np.zeros(NP)
	T1s = np.zeros(NP)
	RH1s = np.zeros(NP)
	T2s = np.zeros(NP)
	RH2s = np.zeros(NP)
	RHs = np.zeros(NP) # calculated from ST, T1, RH1, T2, RH2
    
	for line, i in zip(lines, list(range(NP))):
		infos = re.split(",|;", line.strip())
		
		time_infos = re.split(":|_|\.", infos[0])
		year = int(time_infos[0])
		month = int(time_infos[1])
		day = int(time_infos[2])
		hour = int(time_infos[3])
		minute = int(time_infos[4])
		second = int(time_infos[5])
		jul = date2jul.date2jul(year, month, day, hour=hour, minute=minute, second=second)
		juls[i] = jul
		
		time_infos = re.split("_", infos[0])
		dates.append(time_infos[0])
		times.append(time_infos[1])
		
		SCAs_dry = np.array([float(x) for x in infos[1:4]])
		BSCAs_dry = np.array([float(x) for x in infos[4:7]])
		SCAs_wet = np.array([float(x) for x in infos[10:13]])
		BSCAs_wet = np.array([float(x) for x in infos[13:16]])
		SCAs[i, :, 0] = SCAs_dry
		SCAs[i, :, 1] = SCAs_wet
		BSCAs[i, :, 0] = BSCAs_dry
		BSCAs[i, :, 1] = BSCAs_wet
		
		STs[i] = infos[16]
		states.append(infos[20])
		
		Tset_infos = re.split(" ", infos[21].strip())
		Tsets[i] = Tset_infos[2]
		Ttrue_infos = re.split(" ", infos[22].strip())
		try:
			Ttrues[i] = Ttrue_infos[2]
		except:
			Ttrues[i] = -1
		
		T1_infos = re.split(" ", infos[23].strip())
		T1s[i] = T1_infos[2]
		RH1_infos = re.split(" ", infos[24].strip())
		RH1s[i] = RH1_infos[2]
		T2_infos = re.split(" ", infos[25].strip())
		T2s[i] = T2_infos[2]
		RH2_infos = re.split(" ", infos[26].strip())
		RH2s[i] = RH2_infos[2]
		
		H_avgs = 1/2 * calRH.es(T1s+273.15) * RH1s/100 + 1/2 * calRH.es(T2s+273.15) * RH2s/100 # calculate RH
		RHs = H_avgs / calRH.es(STs+273.15) * 100
		
	if NP == 0:
		fRH = {}
	else:
		fRH = {
			"year"     : year, # int
			"jul"      : juls, # array float
			"date"     : dates, # array string
			"time"     : times, # array string
			"wl"       : wls, # [635, 525, 450]
			"SCA"      : SCAs, # array [[dry_sca_635, wet_sca_635], [...], [...]]
			"BSCA"     : BSCAs, # array
			"state"    : states, # array string, N presents normal
			"ST"       : STs, # sample T
			"Tset"     : Tsets,
			"Ttrue"    : Ttrues,
			"T1"       : T1s,
			"RH1"      : RH1s,
			"T2"       : T2s,
			"RH2"      : RH2s,
			"RH"       : RHs, # sample RH
		}
	return fRH

def Kotchenruther_fit(RH, a, b):
	"""
	This function is used to fit fRH measurements and proposed by Robert A. Kotchenruther in the paper published in 1999: "Humidification factors for atmospheric aerosols off the mid-Atlantic coast of the United States"
	---Keywords-----------------------------
		RH       :  Relative humidify (%)
		a,b      :  Two fitting parameters
	---Return-------------------------------
		fRH      :  aerosol light scattering enhancement factor at RH 
	"""
	fRH = 1 + a * (RH / 100.0) ** b
	return fRH

def cal_fRH_fit_parameter(fRHs, RHs, **args):
	'''
	This function is to fit fRH curve
	input:
		fRHs    : fRH data, array
		RHs     : RH data, array, percent
		**func  : fit function selector, int, default 1
	output:
		paras   : fit parameters
	'''
	if 'func' in args:
		func = args['func']
	else:
		func = 1
	
	def Kotchenruther_fit(RH, a, b):
		'''
		This function is used to fit fRH measurements and proposed by Robert A. Kotchenruther in the paper published in 1999: "Humidification factors for atmospheric aerosols off the mid-Atlantic coast of the United States"
		input:
			RH       :  Relative humidify (%)
			a,b      :  Two fitting parameters
		output:
			fRH      :  aerosol light scattering enhancement factor at RH 
		'''
		fRH = 1 + a * (RH/100)**b
		return fRH
	
	def Brock_fit(RH, ksca):
		'''
		This function is used to fit fRH measurements and porposed by Charles A.Brock in the paper published in 2016 : "Aerosol optical properties in the southeastern United States in summer: Part 1: Hygroscopic growth"
	    	input:
			RH      :  Relative humidity (%)
			ksca    :  The fitting parameter
		output:
			fRH     :  aerosol light scattering enhancement factor at RH
		'''
		fRH = 1 + ksca * RH / (100-RH)
		return fRH
	
	if func==1:
		paras, pcov = curve_fit(Kotchenruther_fit, RHs, fRHs)
	elif func==2:
		paras, pcov = curve_fit(Brock_fit, RHs, fRHs)
	else:
		print('Unknown fit function. Please check.')
		sys.exit()
	
	return paras

def fRH_fit(fRHs, RHs, **args):
	'''
	This function is to calculate fitted fRH curve function
	input:
		fRHs      : fRH data, array
		RHs       : RH data, array, percent
		**func    : fit function selector, int, default 1
	output:
		fRH_curve : fitted fRH curve, function
			input:
				RH    : relative humidity, float, percent
			output:
				fRH   : f(RH) data, float
	'''
	if 'func' in args:
		func = args['func']
	else:
		func = 1
	
	if func==1:
		a, b = cal_fRH_fit_parameter(fRHs, RHs)
		def fRH_curve(RH):
			fRH_curve = 1 + a * (RH/100)**b
			return fRH_curve
	elif func==2:
		ksca = cal_fRH_fit_parameter(fRHs, RHs, func=2)
		def fRH_curve(RH):
			fRH_curve = 1 + ksca * RH / (100-RH)
			return fRH_curve
	else:
		print('Unknown fit function. Please check.')
		sys.exit()
	
	return fRH_curve

def fRH_revise(fRHs, RHs, **args):
	'''
	This function is to revise fRH, delete abnormal data and downward data
	input:
		fRHs     : raw fRH data, array
		RHs      : raw RH data, array, percent
		**minNum : minimum data point number, int, default 0
	output:
		fRHs     : revised fRH data, array
		RHs      : revised RH data, array, percent
	'''
	if 'minNum' in args:
		minNum = args['minNum']
	else:
		minNum = 0
	
	if len(RHs)<=minNum: # delete smaller than minNum data cycle
		return [], []
	
	index = np.where(RHs==RHs.min())[0][0]
	RHs = RHs[index:]
	fRHs = fRHs[index:]
	# donwward RH data have been deleted
	
	index = np.where(fRHs==fRHs.min())[0][0]
	RHs = RHs[index:]
	fRHs = fRHs[index:]
	# donwward fRH data have been deleted
	
	if len(RHs)<=10: # delete small quantity data
		return [], []
	
	return fRHs, RHs

def cal_angstrom(wls, kscas):
	'''
	This function is to calculate angstrom index from fRH data
	input:
		wls      : wave length, array, nm
		kscas    : scattering coefficients, array, Mm^-1
	output:
		angstrom : angstrom index, if fit failed, return -99999
	'''
	def angstrom_fit(wl, A, C):
		logksca = -A*np.log(wl)+C
		return logksca
	
	try:
		paras, pcov = curve_fit(angstrom_fit, wls, np.log(kscas))
	except:
		return -99999
	
	return paras[0]

def readfRHs(fn, **args):
	'''
	This function is to read total fRH data files
	input:
		fn          : file name
		**minNum    : minimum data point number, int, default 0
		**minDyrSca : minimum dry scattering coefficient to delete, float, default 0 Mm-1
		**revise    : whether to revise f(RH) data, boolean, default True
	output:
		cycles      : fRH cycles
			cycles: years, dates, times, juls, wls, drySca_635s, drySca_525s, drySca_450s, fRH_635s, fRH_525s, fRH_450s, f80_635s, f80_525s, f80_450s, ksca_635s, ksca_525s, ksca_450s, RHs, angstroms
	'''
	if 'minNum' in args:
		minNum = args['minNum']
	else:
		minNum = 0
	if 'minDrySca' in args:
		minDrySca = args['minDrySca']
	else:
		minDrySca = 0
	if 'revise' in args:
		revise = args['revise']
	else:
		revise = True
	
	print("loading...")
	
	root, files = getfilesname.file_name(fn)
	files.sort() # sort file names in time sequence
	file_num = len(files)
	years = []
	dates = []
	times = []
	juls = []
	wls = []
	drySca_635s = []
	drySca_525s = []
	drySca_450s = []
	fRH_635s = []
	fRH_525s = []
	fRH_450s = []
	ksca_635s = []
	ksca_525s = []
	ksca_450s = []
	f80_635s = []
	f80_525s = []
	f80_450s = []
	RHs = []
	angstroms = []
	num = 0 # count valuable data number
	
	for i in range(file_num): # load every cycle data
		file_name = root+'/'+files[i]
		cycle = readfRH(file_name)
		'''
		"year" : year, # int
		"jul"  : juls, # array float
		"date" : dates, # array string
		"time" : times, # array string
		"wv"   : wvs, # [635, 525, 450]
		"SCA"  : SCAs, # array [[dry_sca_635, wet_sca_635], [...], [...]]
		"BSCA" : BSCAs, # array
		"state": states, # array string, N presents normal
		"ST"   : STs, # sample T
		"Tset" : Tsets,
		"Ttrue": Ttrues,
		"T1"   : T1s,
		"RH1"  : RH1s,
		"T2"   : T2s,
		"RH2"  : RH2s,
		"RH"   : RHs # sample RH
		'''
		if cycle=={}:
			continue # delete empty data
		
		state = cycle['state']
		if state[0]!='N': # delete parallel data
			continue
		
		SCA_dry = cycle['SCA'][:,:,0]
		if SCA_dry.mean()<minDrySca: # delete too clean data
			continue
		
		RH = cycle['RH']
		SCA = cycle['SCA']
		fRH_635 = SCA[:,0,1] / SCA[:,0,0]
		fRH_525 = SCA[:,1,1] / SCA[:,1,0]
		fRH_450 = SCA[:,2,1] / SCA[:,2,0]
		if revise:
			fRH_635, RH_635 = fRH_revise(fRH_635,RH,minNum=minNum)
			fRH_525, RH_525 = fRH_revise(fRH_525,RH,minNum=minNum)
			fRH_450, RH_450 = fRH_revise(fRH_450,RH,minNum=minNum)
		else:
			RH_635 = RH
			RH_525 = RH
			RH_450 = RH
		if len(fRH_635)<minNum or len(fRH_525)<minNum or len(fRH_450)<minNum: # delete too little data
			continue
		
		num = num +1
		
		ksca_635 = cal_fRH_fit_parameter(fRH_635, RH_635, func=2)
		ksca_525 = cal_fRH_fit_parameter(fRH_525, RH_525, func=2)
		ksca_450 = cal_fRH_fit_parameter(fRH_450, RH_450, func=2)
		
		fRH_635 = fRH_fit(fRH_635, RH_635)(RH)
		fRH_525 = fRH_fit(fRH_525, RH_525)(RH)
		fRH_450 = fRH_fit(fRH_450, RH_450)(RH)
		
		f80_635 = fRH_fit(fRH_635, RH)(80)
		f80_525 = fRH_fit(fRH_525, RH)(80)
		f80_450 = fRH_fit(fRH_450, RH)(80)
		
		drySca_635 = SCA[:,0,0]
		drySca_525 = SCA[:,1,0]
		drySca_450 = SCA[:,2,0]
		wl = cycle['wl']
		angstrom = np.zeros(len(drySca_635))
		for i in range(len(angstrom)):
			angstrom[i] = cal_angstrom(wl, [drySca_635[i], drySca_525[i], drySca_450[i]])
		
		years.append(cycle['year'])
		dates.append(cycle['date'])
		times.append(cycle['time'])
		juls.append(cycle['jul'])
		wls.append(wl)
		drySca_635s.append(drySca_635)
		drySca_525s.append(drySca_525)
		drySca_450s.append(drySca_450)
		fRH_635s.append(fRH_635)
		fRH_525s.append(fRH_525)
		fRH_450s.append(fRH_450)
		f80_635s.append(f80_635)
		f80_525s.append(f80_525)
		f80_450s.append(f80_450)
		ksca_635s.append(ksca_635)
		ksca_525s.append(ksca_525)
		ksca_450s.append(ksca_450)
		RHs.append(RH)
		angstroms.append(angstrom)
		
	print("done")
	print("total", num, "valid fRH data from ", dates[0][0], times[0][0], "to", dates[-1][-1], times[-1][-1])
	
	cycles = dict(years=years, dates=dates, times=times, juls=juls, wls=wls, drySca_635s=drySca_635s, drySca_525s=drySca_525s, drySca_450s=drySca_450s, fRH_635s=fRH_635s, fRH_525s=fRH_525s, fRH_450s=fRH_450s, f80_635s=f80_635s, f80_525s=f80_525s, f80_450s=f80_450s, ksca_635s=ksca_635s, ksca_525s=ksca_525s, ksca_450s=ksca_450s, RHs=RHs, angstroms=angstroms)
	
	return cycles

if __name__ == '__main__':
	fRH = readfRH('data/fRH/HNeph_2017-10-31_21-36-04_NCycle.txt')
	print(fRH['RH'])
	fRHs = readfRHs('data/test')
	print(fRHs['RHs'][0])
	print(fRHs['fRH_635s'][0])
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
