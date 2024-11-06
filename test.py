import numpy as np
import readTaizhou
import calPNSD
import PyMieScatt as ps
import re
import matplotlib.pyplot as plt
import time

if __name__ == '__main__':
	sp2 = readTaizhou.read_Taizhou('data/sp2/Taizhou.npy')
	Dps = sp2['Dps']
	DBC = sp2['DBC']
	DBCps = sp2['DBCps']
	PNSD = sp2['PNSD']
	DMASP2 = sp2['DMASP2'] # dn/dlogDBC
	
	# turn time sequent data to mean data
	PNSD = np.nanmean(PNSD, axis=0)
	DMASP2 = np.nanmean(DMASP2, axis=0)
	# turn DMASP2 to BCPNSD,
	# due to DBCps has different dlogDBCps by bin, have to do special treatment
	BCPNSD = np.zeros(DMASP2.shape) # turn dn/dlogDBC to dn/dlogDBCps/dlogDBC
	for i in range(len(DBCps)):
		dlogDBC = calPNSD.cal_dlnDp(DBC)
		if i<len(DBCps)-1:
			dlogDBCps_i = np.log10(DBCps[i+1]/DBCps[i])
		else:
			dlogDBCps_i = np.log10(DBCps[i]/DBCps[i-1])
		BCPNSD[i] = DMASP2[i] / dlogDBCps_i
	
	clean_PNSD = PNSD / 10 # to fit Beijing aerosol number distribution level
	dirty_PNSD = PNSD / 2
	clean_BCPNSD = BCPNSD / 10
	dirty_BCPNSD = BCPNSD / 2
	
	PNSD = 0.5 * clean_PNSD + 0.5 * dirty_PNSD
	BCPNSD = 0.5 * clean_BCPNSD + 0.5 * dirty_BCPNSD
	
	Dps_m = []
	PNSD_m = []
	
	with open('data/AOT/PNSD.txt', 'r') as f:
		infos = f.readlines()
		for line in infos:
			res = re.split('\t', line[:-1])
			Dps_m.append(float(res[0]))
			PNSD_m.append(float(res[1]))
	
	#Dps = np.array(Dps_m)
	#PNSD = np.array(PNSD_m)
	
	n = calPNSD.PNSD2n(Dps, PNSD, 'd')
	
	y1 = np.zeros(len(Dps))
	y2 = np.zeros(len(Dps))
	ksca1 = 0
	ksca2 = 0
	x_trans = 0
	
	for i in range(len(Dps)):
		Qsca1 = ps.MieQ(1.484+0.05j, 500, Dps[i])[1]
		Qsca2 = ps.MieQ(1.53+0.05j, 500, Dps[i])[1]
		ksca1 += Qsca1 * 1/4 * np.pi * (Dps[i]*1e-9)**2 * n[i] * 1e6
		ksca2 += Qsca2 * 1/4 * np.pi * (Dps[i]*1e-9)**2 * n[i] * 1e6
		y1[i] = Qsca2 - Qsca1
		if y1[i]<0 and x_trans==0:
			x_trans = (Dps[i]+Dps[i-1]) / 2
		y2[i] = y1[i] * 1/4 * np.pi * (Dps[i]*1e-9)**2 * n[i] * 1e6 * 1e6
	
	print(ksca1*1e6, ksca2*1e6)
	print(x_trans)
	
	fig = plt.figure()
	plt.subplots_adjust(left=0.15, right=0.85, top=0.9, bottom=0.15, wspace=0.05, hspace=0.1)
	ax = fig.add_subplot(111)
	ax.plot(Dps, PNSD, c='k')
	ax.set_xlabel('Dps, nm', fontweight='bold')
	ax.set_ylabel('PNSD, $cm^{-3}$', fontweight='bold')
	plt.xscale('log')
	plt.yscale('log')
	ax.set_ylim(1e-1,1e4)
	ax.set_xlim(1e1,5e3)
	ax2 = ax.twinx()
	ax2.axhline(c='gray')
	ax2.plot(Dps, y2, c='b')
	ax2.text(x_trans-5e2, -7e-2, 'Dps = '+str(round(x_trans*100)/100))
	ax2.set_ylabel('dksca', fontweight='bold')
	ax2.scatter(x_trans, 0, c='r')
	ax2.fill_between(Dps, 0, y2, y2>0, color='b', alpha=0.25)
	ax2.fill_between(Dps, y2, 0, y2<0, color='r', alpha=0.25)
	plt.text(0.05, 0.5, 'RI = 1.484, ksca = '+str(round(ksca1*1e6*100)/100), transform=ax.transAxes)
	plt.text(0.05, 0.45, 'RI = 1.53, ksca = '+str(round(ksca2*1e6*100)/100), transform=ax.transAxes)
	plt.text(0.05, 0.4, '\u0394ksca =  '+str(round((ksca2-ksca1)/ksca1*10000)/100)+'%', transform=ax.transAxes)
	save_path = 'figure/test/'
	#plt.savefig(save_path+time.strftime('%Y%m%d%H%M%S', time.localtime())+'_OATvsCF.pdf')
	plt.show()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
