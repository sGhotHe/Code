####################################################################################
# INTRODUCTION:
# This function is to read SMPS data
# Created by Hebs at 20/12/1/17:37
####################################################################################

import numpy as np
import re
import date2jul
import getfilesname
import matplotlib.pyplot as plt

def readSMPS(fn):
	'''
	This function is to read SMPS data
	input:
		fn   : file name, string
	output:
		Dps  : diameter distribution, array, nm
		PNSD : number distribution, array, dn/dlogDp
		n    : total concentration, #/cm^3
		date : date
		time : time
		jul  : julian day
	'''
	begin_sub = 4 # diameter or concentration begin subscript
	end_sub = -26
	n_sub = -2
	date_sub = 1
	time_sub = 2
	num = -1 # store data number, minus introduction line
	bin = 0 # store data bin number
	
	with open(fn,"r") as f:
		for line in f.readlines():
			res = re.split('\t', line)
			if len(res)<10:
				continue
			num = num + 1
			bin = len(res)-begin_sub+end_sub
	
	Dps = np.zeros((num,bin))
	PNSD = np.zeros((num,bin))
	n = np.zeros(num)
	dates = []
	times = []
	juls = np.zeros(num)
	i = 0 # store line number
	
	with open(fn,"r") as f:
		for line in f.readlines():
			res = re.split('\t', line)
			if len(res)<10:
				continue
			if i==0:
				Dp = res[begin_sub:end_sub]
			else:
				PNSD[i-1] = res[begin_sub:end_sub]
				n[i-1] = res[n_sub]
				dates.append(res[date_sub])
				date_infos = re.split('/', res[date_sub])
				year = int(date_infos[2]) + 2000
				month = int(date_infos[0])
				day = int(date_infos[1])
				times.append(res[time_sub])
				time_infos = re.split(':', res[time_sub])
				hour = int(time_infos[0])
				minute = int(time_infos[1])
				second = int(time_infos[2])
				juls[i-1] = date2jul.date2jul(year, month, day, hour=hour, minute=minute, second=second)
			i = i + 1
	
	###########################################################################
	PNSD = PNSD * 2.5 # 0.4Lpm from DMA, 0.6Lpm from clean air
	n = n * 2.5
	###########################################################################
	
	for i in range(num):
		Dps[i] = Dp
	
	smps = {
		"Dps": Dps,
		"PNSD": PNSD,
		"n": n,
		"date": dates,
		"time": times,
		"jul": juls
	}
	return smps

def jul2year(jul):
	'''
	This function is to judge jul is 2020 or 2021
	input:
		jul  : julian day
	output:
		year : year
	'''
	if jul<300:
		return 2021
	return 2020

def readSMPSs(fn):
	'''
	This function is to read multi SMPS datas
	'''
	print("loading...")
	
	root, files = getfilesname.file_name(fn)
	files.sort() # sort file names in time sequence
	file_num = len(files)
	
	Dps = []
	PNSD = []
	n = []
	year = []
	date = []
	time = []
	jul = []
	jul_last = 0 # to record last jul
	
	for i in range(file_num): # load every cycle data
		file_name = root+'/'+files[i]
		smps = readSMPS(file_name)
		tmp = smps['jul']
		if i==0:
			Dps = smps['Dps'][0] # all data have same Dps bin
		# need to do some splicing work here with two files
		# 1. find last jul and next jul
		# 2. build missing jul every 5 min in this gap
		# 3. append to date, time and year, set PNSD to -1
		# 4. attention: should careful when calculate total volumn
		else:
			jul_next = tmp[0] # to record next jul
			jul_interval = jul_next - jul_last
			nd = int(jul_interval)
			hour = int((jul_interval - nd) * 24.0)
			minute = int((jul_interval - nd - hour / 24.0) * 24.0 * 60)
			# attention: minute could be X0, X5 and X4, X9
			jul_num = hour * 12 + round(minute/5.) - 1
			for j in range(jul_num): # append nan data
				jul_append = jul_last + 5 * (j+1) / 24. / 60.
				year_append = jul2year(jul_append)
				date_append = date2jul.jul2strDate(year_append, jul_append, form='M/D/Y')
				time_append = date2jul.jul2strDate(year_append, jul_append, form='h:m:s')
				PNSD_append = np.zeros(109)
				n_append = 0
				
				jul.append(jul_append)
				date.append(date_append)
				time.append(time_append)
				year.append(year_append)
				PNSD.append(PNSD_append)
				n.append(n_append)
		
		for j in range(len(smps['Dps'])):
			PNSD.append(smps['PNSD'][j])
			n.append(smps['n'][j])
			date.append(smps['date'][j])
			time.append(smps['time'][j])
			year.append(jul2year(tmp[j]))
			jul.append(tmp[j])
		jul_last = jul[-1]
	
	print("done")
	print("SMPS data from "+date[0]+" "+time[0]+" to "+date[-1]+" "+time[-1])
	smps = dict(Dps=Dps, PNSD=PNSD, n=n, year=year, date=date, time=time, jul=jul)
	return smps

if __name__ == '__main__':
	smps = readSMPSs('data/smps')
	print(smps['n'])
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
