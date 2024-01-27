####################################################################################
# INTRODUCTION:
# This code is to plot all output result to compare parameter's sensitivity
# Created by Hebs at 21/11/6/13:55
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import matplotlib.pyplot as plt
import runAlbedoChange
import runAODChange
import runBCMChange
import runGChange
import runSSAChange
import read11

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

def plot():
	'''
	Main plot function
	plot parameter:
		albedo
		AOD
		SSA
		g
		BC core complex refractive index m_BC, both real part and imagine part
		shell complex refractive index m_shell
	'''
	#albedo change from -10% to 10%, bin in 1%
	rate = np.arange(-10,11) / 100 + 1
	
	wl_albedo, albedo = runAlbedoChange.read_parameter(rate)
	albedo_nrf = runAlbedoChange.read(rate)
	wl_aod, aod = runAODChange.read_parameter(rate)
	aod_nrf = runAODChange.read(rate)
	m_BC_n = runBCMChange.read_n_parameter(rate)
	m_BC_n_nrf = runBCMChange.read_n(rate)
	m_BC_k = runBCMChange.read_k_parameter(rate)
	m_BC_k_nrf = runBCMChange.read_k(rate)
	m_shell = runBCMChange.read_mshell_parameter(rate)
	m_shell_nrf = runBCMChange.read_mshell(rate)
	wl_g, g = runGChange.read_parameter(rate)
	g_nrf = runGChange.read(rate)
	wl_ssa, ssa = runSSAChange.read_parameter(rate)
	ssa_nrf = runSSAChange.read(rate)
	
	albedo_i = 3 # 860nm
	aod_i = 2 # 870nm
	g_i = 2 # 870nm
	ssa_i = 2 # 870nm
	
	iout11 = read11.read11(fn='output/zero/iout11.txt')
	zero_nrf = iout11['fxdn'][0]-iout11['fxup'][0]
	
	plt.figure(figsize=(12,12))
	plt.subplots_adjust(left=0.08, right=0.99, top=0.99, bottom=0.07, wspace=0.25)
	
	plt.subplot(331)
	plt.plot(m_BC_n, m_BC_n_nrf-zero_nrf)
	ax = plt.gca()
	plt.text(0.075, 0.075, '(a) $m_{BC}$ n', fontsize=12, fontweight='bold', fontstyle='italic', transform=ax.transAxes)
	plt.xticks(cal_xticks(m_BC_n,False), cal_xticks(m_BC_n,True), fontsize=12, fontweight='bold')
	plt.yticks(cal_yticks(m_BC_n_nrf-zero_nrf,False), cal_yticks(m_BC_n_nrf-zero_nrf,True), fontsize=12, fontweight='bold')
	plt.ylabel('Net Radiative Flux, $W/m^2$', fontsize=12, fontweight='bold')
	plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
	
	plt.subplot(332)
	plt.plot(m_BC_k, m_BC_k_nrf-zero_nrf)
	ax = plt.gca()
	plt.text(0.075, 0.075, '(b) $m_{BC}$ k', fontsize=12, fontweight='bold', fontstyle='italic', transform=ax.transAxes)
	plt.xticks(cal_xticks(m_BC_k,False), cal_xticks(m_BC_k,True), fontsize=12, fontweight='bold')
	plt.yticks(cal_yticks(m_BC_k_nrf-zero_nrf,False), cal_xticks(m_BC_k_nrf-zero_nrf,True), fontsize=12, fontweight='bold')
	plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
	
	plt.subplot(333)
	plt.plot(m_shell, m_shell_nrf-zero_nrf)
	ax = plt.gca()
	plt.text(0.075, 0.075, '(c) $m_{shell}$ n', fontsize=12, fontweight='bold', fontstyle='italic', transform=ax.transAxes)
	plt.xticks(cal_xticks(m_shell,False), cal_xticks(m_shell,True), fontsize=12, fontweight='bold')
	plt.yticks(cal_yticks(m_shell_nrf-zero_nrf,False), cal_yticks(m_shell_nrf-zero_nrf,True), fontsize=12, fontweight='bold')
	plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
	
	plt.subplot(334)
	plt.plot(albedo[albedo_i], albedo_nrf-zero_nrf)
	ax = plt.gca()
	plt.text(0.075, 0.075, '(d) Albedo, '+str(wl_albedo[albedo_i]*1e3)+' nm', fontsize=12, fontweight='bold', fontstyle='italic', transform=ax.transAxes)
	plt.xticks(cal_xticks(albedo[albedo_i],False), cal_xticks(albedo[albedo_i],True), fontsize=12, fontweight='bold')
	plt.yticks(cal_yticks(albedo_nrf-zero_nrf,False), cal_yticks(albedo_nrf-zero_nrf,True), fontsize=12, fontweight='bold')
	plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
	
	plt.subplot(335)
	plt.plot(aod[aod_i], aod_nrf-zero_nrf)
	ax = plt.gca()
	plt.text(0.075, 0.075, '(e) AOD, '+str(wl_aod[aod_i])+' nm', fontsize=12, fontweight='bold', fontstyle='italic', transform=ax.transAxes)
	plt.xticks(cal_xticks(aod[aod_i],False), cal_xticks(aod[aod_i],True), fontsize=12, fontweight='bold')
	plt.yticks(cal_yticks(aod_nrf-zero_nrf,False), cal_yticks(aod_nrf-zero_nrf,True), fontsize=12, fontweight='bold')
	plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
	
	plt.subplot(337)
	plt.plot(g[g_i], g_nrf-zero_nrf)
	ax = plt.gca()
	plt.text(0.075, 0.075, '(f) g, '+str(wl_g[g_i])+' nm', fontsize=12, fontweight='bold', fontstyle='italic', transform=ax.transAxes)
	plt.xticks(cal_xticks(g[g_i],False), cal_xticks(g[g_i],True), fontsize=12, fontweight='bold')
	plt.yticks(cal_yticks(g_nrf-zero_nrf,False), cal_yticks(g_nrf-zero_nrf,True), fontsize=12, fontweight='bold')
	plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
	
	plt.subplot(338)
	plt.plot(ssa[ssa_i], ssa_nrf-zero_nrf)
	ax = plt.gca()
	plt.text(0.075, 0.075, '(g) SSA, '+str(wl_ssa[ssa_i])+' nm', fontsize=12, fontweight='bold', fontstyle='italic', transform=ax.transAxes)
	plt.xticks(cal_xticks(ssa[ssa_i],False), cal_xticks(ssa[ssa_i],True), fontsize=12, fontweight='bold')
	plt.yticks(cal_yticks(ssa_nrf-zero_nrf,False), cal_yticks(ssa_nrf-zero_nrf,True), fontsize=12, fontweight='bold')
	plt.grid(which='major', linewidth=1.5, alpha=0.8, color='gray', linestyle=':')
	
	grid = plt.GridSpec(3,3,wspace=0.5,hspace=0.5)
	plt.subplot(grid[1:,2])
	label = ['SSA', 'g', 'AOD', 'Albedo', '$m_{shell}$ n', '$m_{BC}$ k', '$m_{BC}$ n']
	value = [ssa_nrf[-1]-ssa_nrf[10], g_nrf[-1]-g_nrf[10], aod_nrf[-1]-aod_nrf[10], albedo_nrf[-1]-albedo_nrf[10], m_shell_nrf[-1]-m_shell_nrf[10], m_BC_k_nrf[-1]-m_BC_k_nrf[10], m_BC_n_nrf[-1]-m_BC_n_nrf[10]]
	plt.barh(label, value, 0.5)
	plt.axvline(x=0, c='k', linewidth=1.5)
	ax = plt.gca()
	plt.text(0.075, 0.93, '$\Delta$=10%', fontsize=12, fontweight='bold', fontstyle='italic', transform=ax.transAxes)
	plt.xticks(fontsize=12, fontweight='bold')
	plt.xlabel('$\Delta$DARF, $W/m^2$', fontsize=12, fontweight='bold')
	plt.yticks(fontsize=12, fontweight='bold')
	
	#plt.savefig("figure/plotAll.pdf")
	plt.show()
	plt.close()

if __name__ == '__main__':
	plot()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
