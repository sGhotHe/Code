####################################################################################
# INTRODUCTION:
# This code is to plot SBDART iout=7 output data
# Created by Hebs at 21/5/17/16:45
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import matplotlib.pyplot as plt
import read1
import planck

if __name__ == '__main__':
	iout1 = read1.read1()
	botdn = iout1['botdn']
	topup = iout1['topup']
	wl = iout1['wl']
	'''
	iout1_2 = read1.read1(filename='output/iout1_2.txt')
	botdn_2 = iout1_2['botdn']
	topup_2 = iout1_2['topup']
	iout1_3 = read1.read1(filename='output/iout1_3.txt')
	botdn_3 = iout1_3['botdn']
	topup_3 = iout1_3['topup']
	'''
	
	B_280 = np.zeros(len(wl))
	B_260 = np.zeros(len(wl))
	B_240 = np.zeros(len(wl))
	B_220 = np.zeros(len(wl))
	for i in range(len(wl)):
		B_280[i] = planck.B(wl[i], 280)
		B_260[i] = planck.B(wl[i], 260)
		B_240[i] = planck.B(wl[i], 240)
		B_220[i] = planck.B(wl[i], 220)
	
	plt.rcParams['ytick.direction'] = 'in'
	plt.rcParams['xtick.direction'] = 'in'
	
	fig = plt.figure(figsize=(8,4))
	plt.subplots_adjust(bottom=0.15)
	ax = plt.subplot(111)
	ax.tick_params(top='on', right='on', which='both')
	ax.set_title('Downward Surface Irradiance')
	plt.plot(wl, B_280*np.pi, c='k', linestyle='dashdot', alpha=0.7, lw=1)
	plt.text(20.4, 9.3, 'T = 280K', fontsize=12)
	plt.plot(wl, B_260*np.pi, c='k', linestyle='dashdot', alpha=0.7, lw=1)
	plt.text(20.4, 7.4, 'T = 260K', fontsize=12)
	plt.plot(wl, B_240*np.pi, c='k', linestyle='dashdot', alpha=0.7, lw=1)
	plt.text(20.4, 5.7, 'T = 240K', fontsize=12)
	plt.plot(wl, B_220*np.pi, c='k', linestyle='dashdot', alpha=0.7, lw=1)
	plt.text(20.4, 4, 'T = 220K', fontsize=12)
	#plt.plot(wl, botdn_3, c='b', lw=2, alpha=0.7, label='0')
	#plt.plot(wl, botdn_2, c='g', lw=2, alpha=0.7, label='1')
	plt.plot(wl, botdn, c='r', lw=2, alpha=0.7, label='0')
	plt.text(0.9, 0.9, 'Tcloud', fontsize=12, transform=ax.transAxes)
	plt.ylim(0)
	plt.yticks(fontsize=12)
	plt.ylabel('$W/m^2/um$', fontsize=12)
	plt.xlim(4,22)
	plt.xticks([5,10,15,20], [5,10,15,20], fontsize=12)
	plt.xlabel('wavelength, um', fontsize=12)
	plt.minorticks_on()
	plt.legend(bbox_to_anchor=(1,0.9), frameon=False)
	#plt.savefig('figure/plot1.jpg')
	plt.show()
	plt.close()
	'''
	plt.rcParams['ytick.direction'] = 'in'
	plt.rcParams['xtick.direction'] = 'in'
	
	fig = plt.figure(figsize=(8,4))
	plt.subplots_adjust(bottom=0.15)
	ax = plt.subplot(111)
	ax.tick_params(top='on', right='on', which='both')
	ax.set_title('Upward TOA Irradiance')
	plt.plot(wl, B_280*np.pi, c='k', linestyle='dashdot', alpha=0.7, lw=1)
	plt.text(20.4, 9.3, 'T = 280K', fontsize=12)
	plt.plot(wl, B_260*np.pi, c='k', linestyle='dashdot', alpha=0.7, lw=1)
	plt.text(20.4, 7.4, 'T = 260K', fontsize=12)
	plt.plot(wl, B_240*np.pi, c='k', linestyle='dashdot', alpha=0.7, lw=1)
	plt.text(20.4, 5.7, 'T = 240K', fontsize=12)
	plt.plot(wl, B_220*np.pi, c='k', linestyle='dashdot', alpha=0.7, lw=1)
	plt.text(20.4, 4, 'T = 220K', fontsize=12)
	plt.plot(wl, topup_3, c='b', lw=2, alpha=0.7, label='0')
	plt.plot(wl, topup_2, c='g', lw=2, alpha=0.7, label='1')
	plt.plot(wl, topup, c='r', lw=2, alpha=0.7, label='5')
	plt.text(0.9, 0.9, 'Tcloud', fontsize=12, transform=ax.transAxes)
	plt.ylim(0)
	plt.yticks(fontsize=12)
	plt.ylabel('$W/m^2/um$', fontsize=12)
	plt.xlim(4,22)
	plt.xticks([5,10,15,20], [5,10,15,20], fontsize=12)
	plt.xlabel('wavelength, um', fontsize=12)
	plt.minorticks_on()
	plt.legend(bbox_to_anchor=(1,0.9), frameon=False)
	#plt.savefig('figure/plot1_2.jpg')
	plt.show()
	plt.close()
	'''
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
