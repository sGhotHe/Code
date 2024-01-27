####################################################################################
# INTRODUCTION:
# This code is to plot optical closure, data from Taizhou.npy
# Created by Hebs at 21/12/22/10:46
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import scipy.stats as stats
import readTaizhou
import calBC
import calPNSD
import date2jul

def calculate():
	print('loading...')
	data = readTaizhou.read_Taizhou('data/sp2/Taizhou.npy')
	Dps = data['Dps']
	PNSD = data['PNSD'][1000:1200]
	DBCps = data['DBCps']
	DBC = data['DBC']
	DMASP2 = data['DMASP2'][1000:1200]
	BCPNSD = data['BCPNSD'][1000:1200]
	kabs = data['kabs'][1000:1200, 1:]
	ksca = data['ksca'][1000:1200]
	wl_abs = data['wl_abs'][1:]
	wl_sca = data['wl_sca']
	
	data = np.load('data/retrieve/retrieveFromKBulk_2wls.npy', allow_pickle=True).item()
	m_BC_bulk = data['m_BC']
	m_shell_bulk = data['m_shell']
	rext = data['rext']
	
	print('done')
	print('calculating...')
	'''
	ksca_cal_bulk = np.zeros((len(PNSD),len(wl_sca)))
	kabs_cal_bulk = np.zeros((len(PNSD),len(wl_abs)))
	
	for i in range(len(PNSD)):
		for j in range(len(wl_sca)):
			ksca_cal_bulk[i,j] = calBC.cal_ksca_bulk(Dps, PNSD[i], DBC, BCPNSD[i], m_BC_bulk[i], m_shell_bulk[i], wl_sca[j], rext[i])
		for j in range(len(wl_abs)):
			kabs_cal_bulk[i,j] = calBC.cal_kabs_bulk(Dps, PNSD[i], DBC, BCPNSD[i], m_BC_bulk[i], m_shell_bulk[i], wl_abs[j], rext[i])
		print(round((i+1)/len(PNSD)*10000)/100, '% done...')
	
	data = dict(ksca_cal_bulk=ksca_cal_bulk, kabs_cal_bulk=kabs_cal_bulk)
	np.save('data/k/kBulk_2wls.npy', data)
	'''
	data = np.load('data/retrieve/retrieveFromKBulk_2wls.npy', allow_pickle=True).item()
	m_BC = data['m_BC']
	m_shell = data['m_shell']
	rext = data['rext']
	res = data['res']
	data = np.load('data/k/kBulk_2wls.npy', allow_pickle=True).item()
	ksca_cal = data['ksca_cal_bulk']
	kabs_cal = data['kabs_cal_bulk']
	data = dict(m_BC=m_BC, m_shell=m_shell, rext=rext, ksca_cal=ksca_cal, kabs_cal=kabs_cal, res=res)
	np.save('data/retrieve/retrieveFromKBulk_2wls_2.npy', data)
	
	print('done')

def plot():
	print('loading...')
	data = readTaizhou.read_Taizhou('data/sp2/Taizhou.npy')
	Dps = data['Dps']
	PNSD = data['PNSD'][100:500]
	DBCps = data['DBCps']
	DBC = data['DBC']
	DMASP2 = data['DMASP2'][100:500]
	BCPNSD = data['BCPNSD'][100:500]
	kabs = data['kabs'][100:500]
	ksca = data['ksca'][100:500]
	wl_abs = data['wl_abs']
	wl_sca = data['wl_sca']
	jul = data['jul']
	
	data = np.load('data/retrieve/retrieveFromK_3wls.npy', allow_pickle=True).item()
	m_BC = data['m_BC']
	m_shell = data['m_shell']
	ksca_cal = data['ksca_cal']
	kabs_cal = data['kabs_cal']
	
	print('done\ncalculating...')
	BCPVSD = np.zeros((len(PNSD), len(DBC)))
	for i in range(len(PNSD)):
		v_BC = calBC.DMASP22vBC(DBCps, DBC, DMASP2[i])
		BCPVSD[i] = calPNSD.n2PNSD(DBC, v_BC, 'e')
	
	x = np.arange(len(PNSD))
	x_ticks = [0, 50, 100, 150, 200, 250, 300, 350, 399]
	x_juls = jul[x_ticks]
	x_dates = []
	for i in range(len(x_juls)):
		x_dates.append(date2jul.jul2strDate(2018, x_juls[i], form='h:m'))
	
	print('done\nplotting...')
	left = 0.064
	bottom = 0.05
	right = 0.979
	top = 0.96
	wspace = 0.145
	hspace = 0.11
	dw = 0.146
	dh = 0.03
	dy = 0.02
	lw = 3
	fs = 12
	fw = 'bold'
	font = dict(weight=fw, size=fs)
	
	fig = plt.figure(figsize=(16,12))
	plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
	
	# PNSD
	#plt.subplot(611)
	plt.axes([left, bottom+(top-bottom)/6*5+dy, right-left, (top-bottom)/6-dh])
	ax = plt.pcolormesh(x, Dps, np.transpose(PNSD)/1e4, cmap='jet', vmin=0, vmax=4)
	plt.yscale('log')
	plt.ylim(Dps[0], 1000)
	plt.yticks([20, 50, 100, 200, 500, 1000], [20, 50, 100, 200, 500, 1000], fontsize=fs, fontweight=fw)
	plt.ylabel('Dps, nm', fontsize=fs, fontweight=fw)
	plt.xticks([0, 50, 100, 150, 200, 250, 300, 350, 399], ['', '', '', '', '', '', '', '', ''], fontsize=fs, fontweight=fw)
	plt.grid(which='major', linewidth=1.5, color='black', linestyle=':')
	cb = plt.colorbar(ax, extend='both', pad=0.01)
	cb.set_label('dN/dlog(Dp), $\\times 10^4 /cm^3$', fontdict=font)
	cb.ax.tick_params(labelsize=fs)
	
	# ksca
	plt.axes([left, bottom+(top-bottom)/6*4+dy, right-left-dw, (top-bottom)/6-dh])
	plt.plot(x, ksca[:,1], lw=lw)
	plt.xlim(0, len(PNSD)-1)
	plt.xticks([0, 50, 100, 150, 200, 250, 300, 350, 399], ['', '', '', '', '', '', '', '', ''], fontsize=fs, fontweight=fw)
	plt.yticks(fontsize=fs, fontweight=fw)
	plt.ylabel('ksca 525nm', fontsize=fs, fontweight=fw)
	plt.grid(which='major', linewidth=1.5, color='gray', linestyle=':', alpha=0.5)
	
	# BCPVSD
	#plt.subplot(613)
	plt.axes([left, bottom+(top-bottom)/6*3+dy, right-left, (top-bottom)/6-dh])
	ax = plt.pcolormesh(x, DBC, np.transpose(BCPVSD), cmap='jet', vmin=0, vmax=4)
	plt.yscale('log')
	plt.ylim(DBC[5])
	plt.yticks([100, 200, 300, 400], [100, 200, 300, 400], fontsize=fs, fontweight=fw)
	plt.ylabel('DBC, nm', fontsize=fs, fontweight=fw)
	plt.xticks([0, 50, 100, 150, 200, 250, 300, 350, 399], ['', '', '', '', '', '', '', '', ''], fontsize=fs, fontweight=fw)
	plt.grid(which='major', linewidth=1.5, color='black', linestyle=':')
	cb = plt.colorbar(ax, extend='both', pad=0.01)
	cb.set_label('dV/dlog(DBC), $\\times 10^{-12}$', fontdict=font)
	cb.ax.tick_params(labelsize=fs)
	
	# kabs
	plt.axes([left, bottom+(top-bottom)/6*2+dy, right-left-dw, (top-bottom)/6-dh])
	plt.plot(x, kabs[:,1], lw=lw)
	plt.xlim(0, len(PNSD)-1)
	plt.xticks([0, 50, 100, 150, 200, 250, 300, 350, 399], ['', '', '', '', '', '', '', '', ''], fontsize=fs, fontweight=fw)
	plt.yticks(fontsize=fs, fontweight=fw)
	plt.ylabel('kabs 535nm', fontsize=fs, fontweight=fw)
	plt.grid(which='major', linewidth=1.5, color='gray', linestyle=':', alpha=0.5)
	
	# m_BC
	plt.axes([left, bottom+(top-bottom)/6+dy, right-left-dw, (top-bottom)/6-dh])
	plt.plot(x, m_BC.real, lw=lw)
	plt.xlim(0, len(PNSD)-1)
	plt.xticks([0, 50, 100, 150, 200, 250, 300, 350, 399], ['', '', '', '', '', '', '', '', ''], fontsize=fs, fontweight=fw)
	plt.ylabel('$m_{BC}$', fontsize=fs, fontweight=fw)
	plt.yticks(fontsize=fs, fontweight=fw)
	plt.grid(which='major', linewidth=1.5, color='gray', linestyle=':', alpha=0.5)
	
	# m_shell
	plt.axes([left, bottom+dy, right-left-dw, (top-bottom)/6-dh])
	plt.plot(x, m_shell, lw=lw)
	plt.xlim(0, len(PNSD)-1)
	plt.xlabel('Time', fontsize=fs, fontweight=fw)
	plt.xticks([0, 50, 100, 150, 200, 250, 300, 350, 399], x_dates, fontsize=fs, fontweight=fw)
	plt.ylim(1.3, 1.5)
	plt.ylabel('$m_{shell}$', fontsize=fs, fontweight=fw)
	plt.yticks(fontsize=fs, fontweight=fw)
	plt.grid(which='major', linewidth=1.5, color='gray', linestyle=':', alpha=0.5)
	
	plt.show()
	plt.close()
	
	# m_BC and m_shell distribution statistics
	m_shell = m_shell[~np.isnan(m_shell)]
	for i in range(len(m_BC)):
		if np.isnan(m_BC[i].real) or np.isnan(m_BC[i].imag):
			m_BC[i] = 1.5+0.5j
	m_BC_n = m_BC.real
	m_BC_k = m_BC.imag
	
	# three part: polluted, less polluted and clean
	data1 = [m_shell[330:], m_BC_n[330:], []]
	data2 = [[], [], m_BC_k[330:]]
	c = ['r', 'g', 'b']
	alpha = 0.7
	boxprops = dict(linewidth=2) # setting figure style
	medianprops = dict(linewidth=3)
	meanpointprops = dict(marker='+', markersize=15, markeredgecolor='k', markeredgewidth=2)
	capprops = dict(linewidth=2)
	whiskerprops = dict(linewidth=2)
	
	fig = plt.figure(figsize=(8,6))
	plt.subplots_adjust(left=0.12, bottom=0.074, right=0.874, top=0.94, wspace=0.2, hspace=0.22)
	
	ax1 = fig.add_subplot(111)
	bp1 = ax1.boxplot(data1, widths=0.2, showmeans=True, showfliers=False, medianprops=medianprops, meanprops=meanpointprops, boxprops=boxprops, capprops=capprops, whiskerprops=whiskerprops)
	for i in range(len(data1)):
		plt.setp(bp1['boxes'][i], c=c[i], alpha=alpha)
		plt.setp(bp1['whiskers'][2*i], c=c[i], alpha=alpha)
		plt.setp(bp1['whiskers'][2*i+1], c=c[i], alpha=alpha)
		plt.setp(bp1['medians'][i], c=c[i], alpha=alpha)
		plt.setp(bp1['caps'][2*i], c=c[i], alpha=alpha)
		plt.setp(bp1['caps'][2*i+1], c=c[i], alpha=alpha)
	ax1.set_ylabel('n', font=font)
	ax1.set_ylim(1, 1.7)
	plt.yticks(fontweight=fw, fontsize=fs)
	plt.xticks([1, 2, 3], ['$m_{shell}$', '$m_{BC}$ n', '$m_{BC}$ k'], fontweight=fw, fontsize=fs)
	
	ax2 = ax1.twinx()
	bp2 = ax2.boxplot(data2, widths=0.2, showmeans=True, showfliers=False, medianprops=medianprops, meanprops=meanpointprops, boxprops=boxprops, capprops=capprops, whiskerprops=whiskerprops)
	for i in range(len(data2)):
		plt.setp(bp2['boxes'][i], c=c[i], alpha=alpha)
		plt.setp(bp2['whiskers'][2*i], c=c[i], alpha=alpha)
		plt.setp(bp2['whiskers'][2*i+1], c=c[i], alpha=alpha)
		plt.setp(bp2['medians'][i], c=c[i], alpha=alpha)
		plt.setp(bp2['caps'][2*i], c=c[i], alpha=alpha)
		plt.setp(bp2['caps'][2*i+1], c=c[i], alpha=alpha)
	ax2.set_ylabel('k', font=font)
	ax2.set_ylim(0, 0.7)
	plt.yticks(fontweight=fw, fontsize=fs)
	plt.xticks([1, 2, 3], ['$m_{shell}$', '$m_{BC}$ n', '$m_{BC}$ k'], fontweight=fw, fontsize=fs)
	plt.grid(which='major', linewidth=1.5, alpha=0.8, linestyle=':')
	
	def statistics(m_shell, m_BC):
		m_shell_mean = np.nanmean(m_shell)
		m_shell_std = np.std(m_shell)
		m_BC_n_mean = np.nanmean(m_BC.real)
		m_BC_n_std = np.std(m_BC.real)
		m_BC_k_mean = np.nanmean(m_BC.imag)
		m_BC_k_std = np.std(m_BC.imag)
		str_m_shell = str(round(m_shell_mean*1000)/1000) + ' \xB1 ' + str(round(m_shell_std*1000)/1000)
		str_m_BC_n = str(round(m_BC_n_mean*1000)/1000) + ' \xB1 ' + str(round(m_BC_n_std*1000)/1000)
		str_m_BC_k = str(round(m_BC_k_mean*1000)/1000) + ' \xB1 ' + str(round(m_BC_k_std*1000)/1000)
		return str_m_shell, str_m_BC_n, str_m_BC_k
	str_m_shell, str_m_BC_n, str_m_BC_k = statistics(m_shell[330:], m_BC[330:])
	plt.text(0.02, 0.4, str_m_shell, fontsize=fs, fontweight=fw, fontstyle='italic', transform=ax1.transAxes)
	plt.text(0.34, 0.9, str_m_BC_n, fontsize=fs, fontweight=fw, fontstyle='italic', transform=ax1.transAxes)
	plt.text(0.68, 0.9, str_m_BC_k, fontsize=fs, fontweight=fw, fontstyle='italic', transform=ax1.transAxes)
	plt.text(0.1, 0.9, 'Clean', fontsize=21, fontweight=fw, fontstyle='italic', transform=ax1.transAxes)
	
	#plt.show()
	plt.close()
	
	# optical closure results
	def fit(x, y):
		for i in range(len(x)):
			if np.isnan(x[i]) or np.isnan(y[i]):
				x[i] = 0
				y[i] = 0
		A = optimize.curve_fit(lambda X,a:a*X, x, y)[0][0]
		R = stats.pearsonr(x, y)[0]
		return A, R
	
	plt.figure(figsize=(12,4))
	plt.subplots_adjust(left=0.08, bottom=0.15, right=0.979, top=0.94, wspace=0.2, hspace=0.22)
	for i in range(3):
		plt.subplot(131+i)
		plt.scatter(ksca[:,i], ksca_cal[:,i]) #####################
		A, R = fit(ksca[:,i], ksca_cal[:,i]) #####################
		plt.plot(np.arange(300+i*100), A*np.arange(300+i*100), ls='--', lw=2, label='1:1 line', c='orange')
		ax = plt.gca()
		text = 'y = ' + str(round(A*100)/100) + ' * x, $R^2$ = ' + str(round(R**2*1e3)/1e3)
		plt.text(0.2, 0.1, text, fontsize=fs, fontweight=fw, fontstyle='italic', transform=ax.transAxes)
		plt.text(0.1, 0.9, str(wl_sca[i])+' nm', fontsize=fs, fontweight=fw, fontstyle='italic', transform=ax.transAxes)
		plt.grid(which='major', lw=1.5, alpha=0.8, c='gray', ls=':')
		plt.xlim(0, 300+i*100)
		plt.xticks(fontsize=15, fontweight='bold')
		plt.ylim(0, 300+i*100)
		plt.yticks(fontsize=15, fontweight='bold')
		if i==0:
			plt.ylabel('ksca_cal, /Mm', fontsize=fs, fontweight=fw)
			plt.xlabel('Neph, /Mm', fontsize=fs, fontweight=fw)
	#plt.show()
	plt.close()
	
	plt.figure(figsize=(12,4))
	plt.subplots_adjust(left=0.08, bottom=0.15, right=0.979, top=0.94, wspace=0.2, hspace=0.22)
	for i in range(3):
		plt.subplot(131+i)
		plt.scatter(kabs[:,i], kabs_cal[:,i]) #####################
		A, R = fit(kabs[:,i], kabs_cal[:,i]) #####################
		plt.plot(np.arange(60+i*60), A*np.arange(60+i*60), ls='--', lw=2, label='1:1 line', c='orange')
		ax = plt.gca()
		text = 'y = ' + str(round(A*100)/100) + ' * x, $R^2$ = ' + str(round(R**2*1e3)/1e3)
		plt.text(0.25, 0.15, text, fontsize=fs, fontweight=fw, fontstyle='italic', transform=ax.transAxes)
		plt.text(0.1, 0.9, str(wl_abs[i])+' nm', fontsize=fs, fontweight=fw, fontstyle='italic', transform=ax.transAxes)
		plt.grid(which='major', lw=1.5, alpha=0.8, c='gray', ls=':')
		plt.xlim(0, 60+i*60)
		plt.xticks(fontsize=15, fontweight='bold')
		plt.ylim(0, 80+i*10)
		plt.yticks(fontsize=15, fontweight='bold')
		if i==0:
			plt.ylabel('kabs_cal, /Mm', fontsize=fs, fontweight=fw)
			plt.xlabel('PASS, /Mm', fontsize=fs, fontweight=fw)
	#plt.show()
	plt.close()
	print('done')

if __name__ == '__main__':
	#calculate()
	plot()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
