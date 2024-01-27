####################################################################################
# INTRODUCTION:
# This code is to plot all output result to compare parameter's sensitivity in different location, including Beijing, Hainan and Heilongjiang
# Created by Hebs at 21/11/19/9:59
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import matplotlib.pyplot as plt
import runLocationChange
import runAODChange
import runGChange
import runBCMChange
import runFRHChange
import runSSAChange
import os

def cal_xticks(ls, label):
	ticks = np.zeros(5)
	ticks[0] = ls[0]
	ticks[1] = ls[5]
	ticks[2] = ls[10]
	ticks[3] = ls[15]
	ticks[4] = ls[20]
	
	if label:
		for i in range(len(ticks)):
			ticks[i] = round(ticks[i]*100)/100
	
	return ticks

def cal_yticks(ls, label):
	ticks = np.zeros(5)
	ticks[0] = ls[0]
	ticks[1] = 0.75 * ls[0] + 0.25 * ls[20]
	ticks[2] = 0.5 * ls[0] + 0.5 * ls[20]
	ticks[3] = 0.25 * ls[0] + 0.75 * ls[20]
	ticks[4] = ls[20]
	
	if label:
		for i in range(len(ticks)):
			ticks[i] = round(ticks[i])
	
	return ticks

def plot(rate, location):
	'''
	Plot main function
	input:
		rate     : SSA change rate, array
		location : location for SBDART, dictionary
			name : location name, string
			alat : latitude of point on earth's surface, int
			alon : east longitude of point on earth's surface, int
	'''
	print('load data...')
	
	if os.path.exists('output/location/data.npy'):
		data = np.load('output/location/data.npy', allow_pickle=True).item()
		wl_albedo = data['wl_albedo']
		albedo = data['albedo']
		albedo_nrf = data['albedo_nrf']
		wl_aod = data['wl_aod']
		aod = data['aod']
		aod_nrf = data['aod_nrf']
		wl_ssa = data['wl_ssa']
		ssa = data['ssa']
		ssa_nrf = data['ssa_nrf']
		wl_g = data['wl_g']
		g = data['g']
		g_nrf = data['g_nrf']
		m_BC_k = data['m_BC_k']
		m_BC_k_nrf = data['m_BC_k_nrf']
		m_BC_n = data['m_BC_n']
		m_BC_n_nrf = data['m_BC_n_nrf']
		m_shell = data['m_shell']
		m_shell_nrf = data['m_shell_nrf']
		kappa = data['kappa']
		kappa_nrf = data['kappa_nrf']
		zero_nrf = data['zero_nrf']
	else:
		wl_albedo, albedo = runLocationChange.read_albedo_parameter(rate, location)
		albedo_nrf = runLocationChange.read_albedo(rate, location)
		wl_aod, aod = runAODChange.read_parameter(rate)
		aod_nrf = runLocationChange.read_AOD(rate, location)
		wl_g, g = runGChange.read_parameter(rate)
		g_nrf = runLocationChange.read_g(rate, location)
		m_BC_k = runBCMChange.read_k_parameter(rate)
		m_BC_k_nrf = runLocationChange.read_k(rate, location)
		kappa = runFRHChange.read_parameter(rate)
		kappa_nrf = runLocationChange.read_kappa(rate, location)
		m_shell = runBCMChange.read_mshell_parameter(rate)
		m_shell_nrf = runLocationChange.read_mshell(rate, location)
		m_BC_n = runBCMChange.read_n_parameter(rate)
		m_BC_n_nrf = runLocationChange.read_n(rate, location)
		wl_ssa, ssa = runSSAChange.read_parameter(rate)
		ssa_nrf = runLocationChange.read_SSA(rate, location)
		zero_nrf = runLocationChange.read_zero(location)
		data = dict(wl_albedo=wl_albedo, albedo=albedo, albedo_nrf=albedo_nrf, wl_aod=wl_aod, aod=aod, aod_nrf=aod_nrf, wl_ssa=wl_ssa, ssa=ssa, ssa_nrf=ssa_nrf, wl_g=wl_g, g=g, g_nrf=g_nrf, m_BC_k=m_BC_k, m_BC_k_nrf=m_BC_k_nrf, m_BC_n=m_BC_n, m_BC_n_nrf=m_BC_n_nrf, m_shell=m_shell, m_shell_nrf=m_shell_nrf, kappa=kappa, kappa_nrf=kappa_nrf, zero_nrf=zero_nrf)
		np.save('output/location/data.npy', data)
	
	print('done')
	
	name = location['name']
	alat = location['alat']
	for i in range(len(name)):
		print(name[i], ': latitude:', alat[i], ', AOD average:', np.mean(aod[i]), ', albedo average:', np.mean(albedo[i]))
	
	print('plotting...')
	
	g = np.nanmean(g, axis=2)
	albedo_i = 3 # 0.86um
	aod_i = 2 # 870nm
	g_i = 2 # 870nm
	ssa_i = 2 # 870nm
	c = ['r', 'b', 'orange']
	yticks = [-80, -60, -40, -20, 0, 20]
	fontsize = 16
	lw = 5
	savefig = False
	
	print('plot DARF')
	
	plt.figure(figsize=(12,12))
	plt.subplots_adjust(left=0.08, right=0.98, top=0.98, bottom=0.06, wspace=0.18, hspace=0.33)
	
	plt.subplot(331)
	plt.plot(m_BC_n, m_BC_n_nrf[0]-zero_nrf[0], c='r', label='Beijing', lw=lw)
	plt.plot(m_BC_n, m_BC_n_nrf[1]-zero_nrf[1], c='b', label='Hainan', lw=lw)
	plt.plot(m_BC_n, m_BC_n_nrf[2]-zero_nrf[2], c='orange', label='Heilongjiang', lw=lw)
	ax = plt.gca()
	plt.xticks(cal_xticks(m_BC_n,False), cal_xticks(m_BC_n,True), fontsize=fontsize, fontweight='bold')
	plt.xlabel('$m_{BC}$ n', fontsize=fontsize, fontweight='bold')
	plt.yticks(yticks, yticks, fontsize=fontsize, fontweight='bold')
	plt.ylabel('DARF, $W/m^2$', fontsize=fontsize, fontweight='bold')
	plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
	plt.legend(prop={'size':'large', 'weight':'bold'})
	
	plt.subplot(332)
	for i in range(len(name)):
		plt.plot(m_BC_k, m_BC_k_nrf[i]-zero_nrf[i], c=c[i], lw=lw)
	ax = plt.gca()
	plt.xticks(cal_xticks(m_BC_k,False), cal_xticks(m_BC_k,True), fontsize=fontsize, fontweight='bold')
	plt.xlabel('$m_{BC}$ k', fontsize=fontsize, fontweight='bold')
	plt.yticks(yticks, yticks, fontsize=fontsize, fontweight='bold')
	plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
	
	plt.subplot(333)
	for i in range(len(name)):
		plt.plot(m_shell, m_shell_nrf[i]-zero_nrf[i], c=c[i], lw=lw)
	ax = plt.gca()
	plt.xticks(cal_xticks(m_shell,False), cal_xticks(m_shell,True), fontsize=fontsize, fontweight='bold')
	plt.xlabel('$m_{shell}$ n', fontsize=fontsize, fontweight='bold')
	plt.yticks(yticks, yticks, fontsize=fontsize, fontweight='bold')
	plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
	
	plt.subplot(334)
	for i in range(len(name)):
		plt.plot(aod[aod_i], aod_nrf[i]-zero_nrf[i], c=c[i], lw=lw)
	ax = plt.gca()
	plt.xticks(cal_xticks(aod[aod_i],False), cal_xticks(aod[aod_i],True), fontsize=fontsize, fontweight='bold')
	plt.xlabel('AOD', fontsize=fontsize, fontweight='bold')
	plt.yticks(yticks, yticks, fontsize=fontsize, fontweight='bold')
	plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
	
	plt.subplot(335)
	for i in range(len(name)):
		plt.plot(g[g_i], g_nrf[i]-zero_nrf[i], c=c[i], lw=lw)
	ax = plt.gca()
	plt.xticks(cal_xticks(g[g_i],False), cal_xticks(g[g_i],True), fontsize=fontsize, fontweight='bold')
	plt.xlabel('g', fontsize=fontsize, fontweight='bold')
	plt.yticks(yticks, yticks, fontsize=fontsize, fontweight='bold')
	plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
	
	plt.subplot(336)
	for i in range(len(name)):
		plt.plot(rate, albedo_nrf[i]-zero_nrf[i], c=c[i], lw=lw)
	ax = plt.gca()
	plt.xticks(cal_xticks(rate,False), cal_xticks(rate,True), fontsize=fontsize, fontweight='bold')
	plt.xlabel('albedo, change rate', fontsize=fontsize, fontweight='bold')
	plt.yticks(yticks, yticks, fontsize=fontsize, fontweight='bold')
	plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
	
	plt.subplot(325)
	for i in range(len(name)):
		plt.plot(ssa[ssa_i], ssa_nrf[i]-zero_nrf[i], c=c[i], lw=lw)
	ax = plt.gca()
	plt.xticks(cal_xticks(ssa[ssa_i],False), cal_xticks(ssa[ssa_i],True), fontsize=fontsize, fontweight='bold')
	plt.xlabel('SSA', fontsize=fontsize, fontweight='bold')
	plt.yticks(yticks, yticks, fontsize=fontsize, fontweight='bold')
	plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
	
	plt.subplot(326)
	for i in range(len(name)):
		plt.plot(kappa, kappa_nrf[i]-zero_nrf[i], c=c[i], lw=lw)
	ax = plt.gca()
	plt.xticks(cal_xticks(kappa,False), cal_xticks(kappa,True), fontsize=fontsize, fontweight='bold')
	plt.xlabel('kappa', fontsize=fontsize, fontweight='bold')
	plt.yticks(yticks, yticks, fontsize=fontsize, fontweight='bold')
	plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
	
	if savefig:
		plt.savefig("figure/plotLocation1.pdf")
	plt.show()
	plt.close()
	
	plt.figure(figsize=(12,12))
	plt.subplots_adjust(left=0.1, right=0.99, top=0.96, bottom=0.05, hspace=0.4)
	for i in range(len(name)):
		plt.subplot(311+i)
		plt.title(name[i], fontsize=20, fontweight='bold')
		label = ['AOD', 'SSA', 'g', 'Albedo', 'kappa', '$m_{shell}$ n', '$m_{BC}$ k', '$m_{BC}$ n']
		value = [aod_nrf[i,-1]-aod_nrf[i,10], ssa_nrf[i,-1]-ssa_nrf[i,10], g_nrf[i,-1]-g_nrf[i,10], albedo_nrf[i,-1]-albedo_nrf[i,10], kappa_nrf[i,-1]-kappa_nrf[i,10], m_shell_nrf[i,-1]-m_shell_nrf[i,10], m_BC_k_nrf[i,-1]-m_BC_k_nrf[i,10], m_BC_n_nrf[i,-1]-m_BC_n_nrf[i,10]]
		plt.bar(label, value, 0.5, color=c[i])
		plt.axhline(y=0, c='gray', linewidth=1.5)
		ax = plt.gca()
		if i==0:
			plt.text(0.04, 0.8, '$\Delta$=10%', fontsize=20, fontweight='bold', fontstyle='italic', transform=ax.transAxes)
			plt.ylabel('$\Delta$DARF, $W/m^2$', fontsize=20, fontweight='bold')
		plt.xticks(fontsize=20, fontweight='bold')
		plt.yticks(fontsize=20, fontweight='bold')
	
	if savefig:
		plt.savefig("figure/plotLocation2.pdf")
	plt.show()
	plt.close()
	
	plt.figure(figsize=(12,6))
	plt.subplots_adjust(left=0.13, right=0.99, top=0.93, bottom=0.07, hspace=0.8)
	for i in range(len(name)):
		plt.subplot(311+i)
		plt.title(name[i], fontsize=20, fontweight='bold')
		label = ['AOD', 'SSA', 'g', 'Albedo', 'kappa', '$m_{shell}$ n', '$m_{BC}$ k', '$m_{BC}$ n']
		value = [aod_nrf[i,11]-aod_nrf[i,10], ssa_nrf[i,11]-ssa_nrf[i,10], g_nrf[i,11]-g_nrf[i,10], albedo_nrf[i,11]-albedo_nrf[i,10], kappa_nrf[i,11]-kappa_nrf[i,10], m_shell_nrf[i,11]-m_shell_nrf[i,10], m_BC_k_nrf[i,11]-m_BC_k_nrf[i,10], m_BC_n_nrf[i,11]-m_BC_n_nrf[i,10]]
		p = [aod[aod_i,10], ssa[ssa_i,10], g[g_i,10], albedo[i,albedo_i,10], kappa[10], m_shell[10], m_BC_k[10], m_BC_n[10]]
		s = np.array(value, dtype=float) / (np.array(p, dtype=float)*0.01)
		plt.bar(label, s, 0.5, color=c[i])
		plt.axhline(y=0, c='gray', linewidth=1.5)
		ax = plt.gca()
		if i==1:
			plt.ylabel('$S_i$($F_S$), $W/m^2$', fontsize=20, fontweight='bold')
		plt.xticks(fontsize=20, fontweight='bold')
		plt.yticks(fontsize=20, fontweight='bold')
	
	if savefig:
		plt.savefig("figure/plotLocation3.pdf")
	plt.show()
	plt.close()
	
	plt.figure(figsize=(12,6))
	plt.subplots_adjust(left=0.13, right=0.99, top=0.92, bottom=0.07, hspace=0.8)
	for i in range(len(name)):
		plt.subplot(311+i)
		plt.title(name[i], fontsize=20, fontweight='bold')
		label = ['AOD', 'SSA', 'g', 'Albedo', 'kappa', '$m_{shell}$ n', '$m_{BC}$ k', '$m_{BC}$ n']
		value = [aod_nrf[i,11]-aod_nrf[i,10], ssa_nrf[i,11]-ssa_nrf[i,10], g_nrf[i,11]-g_nrf[i,10], albedo_nrf[i,11]-albedo_nrf[i,10], kappa_nrf[i,11]-kappa_nrf[i,10], m_shell_nrf[i,11]-m_shell_nrf[i,10], m_BC_k_nrf[i,11]-m_BC_k_nrf[i,10], m_BC_n_nrf[i,11]-m_BC_n_nrf[i,10]]
		p = [aod[aod_i,10], ssa[ssa_i,10], g[g_i,10], albedo[i,albedo_i,10], kappa[10], m_shell[10], m_BC_k[10], m_BC_n[10]]
		s = np.array(value, dtype=float) / (np.array(p, dtype=float)*0.01)
		delta_p = np.array([0.01, 0.03, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1])
		s = abs(s) * delta_p
		plt.bar(label, s, 0.5, color=c[i])
		plt.axhline(y=0, c='gray', linewidth=1.5)
		ax = plt.gca()
		if i==1:
			plt.ylabel('$\Delta$DARF, $W/m^2$', fontsize=20, fontweight='bold')
		plt.xticks(fontsize=20, fontweight='bold')
		plt.yticks(fontsize=20, fontweight='bold')
	
	if savefig:
		plt.savefig("figure/plotLocation4.pdf")
	plt.show()
	plt.close()

if __name__ == '__main__':
	# parameter change from -10% to 10%, bin in 1%
	rate = np.arange(-10,11) / 100 + 1
	
	name = ['Beijing', 'Hainan', 'Heilongjiang']
	alat = [40, 18, 53]
	alon = [116, 109, 135]
	location = dict(name=name, alat=alat, alon=alon)
	
	plot(rate, location)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
