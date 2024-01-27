####################################################################################
# INTRODUCTION:
# This code is to read vertical distribution data from data.cma.cn
# Created by Hebs at 21/10/26/12:43
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import re

def read(fn):
	'''
	This function is to read vertical distribution data from data.cma.cn
	input:
		fn    : data file name, string
	output:
		data  : virtical distribution data, dictionary, year: year; mon: month; day: day; hour: hour; nn: layer number; z: altitude, m; p: pressure, hPa; t: temperature, Celsius; dpt: dew point temperature, Celsius
		if t or dpt equals to 999999, data are not detected
	'''
	year = []
	mon = []
	day = []
	hour = []
	z = []
	p = []
	t = []
	dpt = []
	
	with open(fn, 'r') as f:
		f.readline()
		for line in f.readlines():
			res = re.split(' ', line[:-2]) # delete space and \n
			year.append(int(res[1]))
			mon.append(int(res[2]))
			day.append(int(res[3]))
			hour.append(int(res[4]))
			z.append(float(res[7]))
			p.append(float(res[6]))
			if float(res[8])>99999:
				t.append(np.nan)
			else:
				t.append(float(res[8]))
			if float(res[9])>99999:
				dpt.append(np.nan)
			else:
				dpt.append(float(res[9]))
	
	year = year[0]
	mon = mon[0]
	day = day[0]
	hour = np.array(hour)
	z = np.array(z)
	p = np.array(p)
	t = np.array(t)
	dpt = np.array(dpt)
	
	data = dict(nn=len(z), year=year, mon=mon, day=day, hour=hour, z=z, p=p, t=t, dpt=dpt)
	return data

if __name__ == '__main__':
	data = read('data/vertical/data.txt')
	print(data['dpt'])
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
