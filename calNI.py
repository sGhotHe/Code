####################################################################################
# INTRODUCTION:
# This code is to run complex refractive index inhomogeneity's influence on Mie function
# Created by Hebs at 23/5/15/13:53
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import PyMieScatt as ps

import calPNSD
import kappaKohler
import readAtms
import phaseFunc

def nomal(x, miu, sigma):
	'''
	nomal distribution function
	f(x) = 1/sqrt(2*pi)/sigma*exp(-(x-miu)**2/2/sigma**2)
	input:
		x		: x
		miu		: miu
		sigma	: sigma
	output:
		y		: f(x)
	'''
	y = 1 / np.sqrt(2*np.pi) / sigma * np.exp(-(x-miu)**2/2/sigma**2)
	return y

def cal_H1(val, binNum, H1):
	'''
	This function is to calculate parameter's 1 level heterogeneity
	input:
		val			: parameter value, float
		binNum		: Dps bin number, int
		H1			: size resolved value change rate, float
	output:
		val_new		: new parameter, float, array in shape [binNum]
	'''
	val_new = np.empty(binNum,dtype=float)
	for i in range(binNum):
		val_new[i] = val + val * H1 * (i-binNum/2) / binNum
	
	return val_new

def cal_H2(val, binNum, H1, H2x, H2y):
	'''
	This function is to calculate parameter's 2 level heterogeneity
	input:
		val			: parameter value, float
		binNum		: Dps bin number, int
		H1			: size resolved value change rate, float
		H2x			: peaks number, int
		H2y			: peaks separate rate, float
	output:
		val_new		: new parameter, float, array in shape [binNum,H2x]
		val_p		: proportion of each peak, float, array in shape [H2x]
	'''
	# define nI as a two dimensional parameter, 1st dimension: size resolve;
	# 2nd dimension: inhomogeneity in size bin.
	# the 1-d inhomogeneity use linear change to represent;
	# the 2-d inhomogeneity use nomal change to represent. update in 24/5/29
	# to upgrade calPNSD.cal_Mie() function.
	# use H1 and H2, which H1 is size resolved n change slope,
	# H2y is two peak's separate rate, H2x is peak's number,
	# for nomal distribution function f(x,miu,sigma), 
	# the proportion is f(H2y[j],(H2x-1)/2,sigma), need to be nomalized.
	
	val_new = np.empty((binNum,H2x),dtype=float)
	val_p = np.empty(H2x,dtype=float)
	
	for j in range(H2x):
		val_p[j] = nomal(j,(H2x-1)/2,1)
	val_p_total = np.sum(val_p)
	val_p = val_p / val_p_total
	
	for i in range(binNum):
		val_new_i = val + val * H1 * (i-binNum/2) / binNum
		for j in range(H2x):
			val_new[i,j] = val_new_i + val * H2y * (j-(H2x-1)/2) / H2x
		
	return val_new, val_p

def cal_Mie5(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, nShellH1, kappa, kappaH1, MS, MSH1, MSH2x, MSH2y, RH, wl, moma, CT, CTH1, BCAE, **args):
	'''
	Same as cal_Mie4, but output particle size information
	input:
		Dps      : diameter distribution, array, nm
		PNSD     : number distribution, array, dn/dlogDps, cm^-3
		DBCps    : BC particle diameter size, array, nm
		DBC      : BC core diameter size, array, nm
		BCPNSD   : BC particle number concentration, array, dn/dlogDBCps/dlogDBC, cm^-3
		kBC      : BC core imagine part of complex refractive index, float
		nBC      : BC core real part of complex refractive index, float
		nShell   : shell complex refractive index, float
		nShellH1 : n mixing state change rate by diameter, float
		kappa    : hygroscopicity parameter, float
		kappaH1  : hygroscopicity parameter mixing state change rate by diameter, float
		MS       : mixing state, float
		MSH1     : mixing state change rate by diameter, float
		MSH2x    : mixing state bin peak number, int
		MSH2y    : mixing state bin peak saperate rate, float
		RH       : relative humidity, float, percent
		wl       : wave length, nm, float
		CT       : BC core coating thickness adjust parameter, float
		CTH1     : BC core coating thickness change rate by diameter, float
		BCAE     : BC core absorbing enhancement adjust parameter, float
		**args:
			angularResolution : angular resolution of PyMieScatt phase function, degree, defualt 0.5
	output:
		kext     : extinct coefficient, Mm-1, float
		waer     : single scattering coefficient, float
		pmom     : Legendre moments of asymmetry factor, np.array, in shape (moma)
		g        : asymmetry factor, float
	'''
	if 'angularResolution' in args:
		angularResolution = args['angularResolution']
	else:
		angularResolution = 0.5
	
	m_BC = nBC + kBC * 1j
	n = calPNSD.PNSD2n(Dps, PNSD, 'd') # in cm^-3
	dlogDps = calPNSD.cal_dlnDp(Dps)
	dlogDBC = calPNSD.cal_dlnDp(DBC)
	
	m_shell = cal_H1(nShell, len(Dps), nShellH1)
	#print('m : ', m_shell)
	kappa_shell = cal_H1(kappa, len(Dps), kappaH1)
	#print('kappa : ', kappa_shell)
	MS_shell, MS_p = cal_H2(MS, len(Dps), MSH1, MSH2x, MSH2y)
	#print('MS : ', MS_shell)
	#print('MS proportion : ', MS_p)
	CT_shell = cal_H1(CT, len(Dps), CTH1)
	#print('CT : ', CT_shell)
	
	kext = np.zeros(len(Dps))
	ksca = np.zeros(len(Dps))
	g = np.zeros(len(Dps))
	
	for i in range(len(Dps)):
		for k in range(MSH2x):
			DMASP2_i = calPNSD.DMASP2PNSD_Dp(DBCps, BCPNSD*dlogDBC, Dps[i])# * dlogDps # array in shape of DBC, in dn/dlogDBC
			nBC_i = DMASP2_i# * dlogDBC # in cm^-3
			# three parts: BC and coated material, only BC and only coated material, use MS to decide how much
			n_noBC_i = n[i] - sum(nBC_i) * MS_shell[i,k] # no BC part and external part
			n_onlyBC_i = nBC_i * (1-MS_shell[i,k]) # external part, in shape of (len(DBC))
			n_int_i = nBC_i * MS_shell[i,k] # internal part, in shape of (len(DBC))
			#print(round(Dps[i]), n_noBC_i, sum(n_onlyBC_i), sum(n_int_i))
			
			# calculate external part
			if n_noBC_i > 0:
				D_noBC_i, m_shell_noBC_i = kappaKohler.RH2D(Dps[i], 0, m_BC, m_shell[i], kappa_shell[i], RH, wl)
				#print(D_noBC_i, m_shell_noBC_i, Dps[i], m_shell[i])
				#print('D_noBC_i: ', D_noBC_i)
				MieQ = ps.MieQ(m_shell_noBC_i, wl, D_noBC_i)
				ScatteringFunction = ps.ScatteringFunction(m_shell_noBC_i, wl, D_noBC_i, angularResolution=angularResolution, normalization='t')
				
				Qsca = MieQ[1]
				ksca[i] += Qsca * 1/4 * np.pi * (D_noBC_i*1e-9)**2 * n_noBC_i*1e6 * MS_p[k]
				# separate peak making n separate too
				
				Qext = MieQ[0]
				kext[i] += Qext * 1/4 * np.pi * (D_noBC_i*1e-9)**2 * n_noBC_i*1e6 * MS_p[k]
				#print('noBC: ',Qsca,' ',Qext)
				#print('noBC: ',Qsca/Qext)
				
				g[i] += MieQ[3] * n_noBC_i / sum(n) * MS_p[k]
			
			if sum(n_onlyBC_i) > 0:
				for j in range(len(DBC)):
					# no hygrosocpic groth
					#print('DBC[j]: ', DBC[j])
					MieQ = ps.MieQ(m_BC, wl, DBC[j])
					ScatteringFunction = ps.ScatteringFunction(m_BC, wl, DBC[j], angularResolution=angularResolution, normalization='t')
					
					Qsca = MieQ[1]
					ksca[i] += Qsca * 1/4 * np.pi * (DBC[j]*1e-9)**2 * n_onlyBC_i[j]*1e6 * MS_p[k]
					
					Qext = MieQ[0]
					kext[i] += Qext * 1/4 * np.pi * (DBC[j]*1e-9)**2 * n_onlyBC_i[j]*1e6 * MS_p[k]
					#print('onlyBC: ',Qsca,' ',Qext)
					#print('onlyBC: ',Qsca/Qext)
					
					g[i] += MieQ[3] * n_onlyBC_i[j] / sum(n) * MS_p[k]
			
			# calculate internal mix part
			if sum(n_int_i) > 0:
				for j in range(len(DBC)):
					D_BC_ij, m_shell_BC_ij = kappaKohler.RH2D(Dps[i], DBC[j], m_BC, m_shell[i], kappa_shell[i], RH, wl)
					#print('D_BC_ij, Dps[i]: ', D_BC_ij, ' ', Dps[i])
					MieQCoreShell = ps.MieQCoreShell(m_BC, m_shell_BC_ij, wl, DBC[j], D_BC_ij*CT_shell[i])
					CoreShellScatteringFunction = ps.CoreShellScatteringFunction(m_BC, m_shell_BC_ij, wl, DBC[j], D_BC_ij*CT_shell[i], angularResolution=180/(180/angularResolution+1), normed=True)
					
					Qsca = MieQCoreShell[1]
					ksca[i] += Qsca * 1/4 * np.pi * (D_BC_ij*1e-9)**2 * n_int_i[j]*1e6 * MS_p[k]
					
					Qext = MieQCoreShell[1] + MieQCoreShell[2] * BCAE
					kext[i] += Qext * 1/4 * np.pi * (D_BC_ij*1e-9)**2 * n_int_i[j]*1e6 * MS_p[k]
					#print('internal: ',Qsca,' ',Qext)
					#print('internal: ',Qsca/Qext)
					
					g[i] += MieQCoreShell[3] * n_int_i[j] / sum(n) * MS_p[k]
	
	waer = ksca / kext
	kext = kext * 1e6 # m^-1 to Mm^-1
	#print(sum(PNSD))
	#print(ksca * 1e6)
	#print(kext, waer, g)
	
	return kext, waer, g

def cal_Mie4(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, nShellH1, kappa, kappaH1, MS, MSH1, MSH2x, MSH2y, RH, wl, moma, CT, CTH1, BCAE, **args):
	'''
	This function is to use Dps, PNSD, RH, CT and BCAE to calculate extinct coefficient, 
	scattering coefficient and Legendre moments of asymmetry factor.
	Add new heterogeneity parameter: MS, CT, nShell and kappa
	input:
		Dps      : diameter distribution, array, nm
		PNSD     : number distribution, array, dn/dlogDps, cm^-3
		DBCps    : BC particle diameter size, array, nm
		DBC      : BC core diameter size, array, nm
		BCPNSD   : BC particle number concentration, array, dn/dlogDBCps/dlogDBC, cm^-3
		kBC      : BC core imagine part of complex refractive index, float
		nBC      : BC core real part of complex refractive index, float
		nShell   : shell complex refractive index, float
		nShellH1 : n mixing state change rate by diameter, float
		kappa    : hygroscopicity parameter, float
		kappaH1  : hygroscopicity parameter mixing state change rate by diameter, float
		MS       : mixing state, float
		MSH1     : mixing state change rate by diameter, float
		MSH2x    : mixing state bin peak number, int
		MSH2y    : mixing state bin peak saperate rate, float
		RH       : relative humidity, float, percent
		wl       : wave length, nm, float
		CT       : BC core coating thickness adjust parameter, float
		CTH1     : BC core coating thickness change rate by diameter, float
		BCAE     : BC core absorbing enhancement adjust parameter, float
		**args:
			angularResolution : angular resolution of PyMieScatt phase function, degree, defualt 0.5
	output:
		kext     : extinct coefficient, Mm-1, float
		waer     : single scattering coefficient, float
		pmom     : Legendre moments of asymmetry factor, np.array, in shape (moma)
		g        : asymmetry factor, float
	'''
	if 'angularResolution' in args:
		angularResolution = args['angularResolution']
	else:
		angularResolution = 0.5
	
	m_BC = nBC + kBC * 1j
	n = calPNSD.PNSD2n(Dps, PNSD, 'd') # in cm^-3
	dlogDps = calPNSD.cal_dlnDp(Dps)
	dlogDBC = calPNSD.cal_dlnDp(DBC)
	
	m_shell = cal_H1(nShell, len(Dps), nShellH1)
	#print('m : ', m_shell)
	kappa_shell = cal_H1(kappa, len(Dps), kappaH1)
	#print('kappa : ', kappa_shell)
	MS_shell, MS_p = cal_H2(MS, len(Dps), MSH1, MSH2x, MSH2y)
	#print('MS : ', MS_shell)
	#print('MS proportion : ', MS_p)
	CT_shell = cal_H1(CT, len(Dps), CTH1)
	#print('CT : ', CT_shell)
	
	kext = 0
	ksca = 0
	g = 0
	theta = np.zeros(round(180/angularResolution+1))
	SU_bulk = np.zeros(len(theta))
	pmom = np.zeros(moma)
	
	for i in range(len(Dps)):
		for k in range(MSH2x):
			DMASP2_i = calPNSD.DMASP2PNSD_Dp(DBCps, BCPNSD*dlogDBC, Dps[i])# * dlogDps # array in shape of DBC, in dn/dlogDBC
			nBC_i = DMASP2_i# * dlogDBC # in cm^-3
			# three parts: BC and coated material, only BC and only coated material, use MS to decide how much
			n_noBC_i = n[i] - sum(nBC_i) * MS_shell[i,k] # no BC part and external part
			n_onlyBC_i = nBC_i * (1-MS_shell[i,k]) # external part, in shape of (len(DBC))
			n_int_i = nBC_i * MS_shell[i,k] # internal part, in shape of (len(DBC))
			#print(round(Dps[i]), n_noBC_i, sum(n_onlyBC_i), sum(n_int_i))
			
			# calculate external part
			if n_noBC_i > 0:
				D_noBC_i, m_shell_noBC_i = kappaKohler.RH2D(Dps[i], 0, m_BC, m_shell[i], kappa_shell[i], RH, wl)
				#print(D_noBC_i, m_shell_noBC_i, Dps[i], m_shell[i])
				#print('D_noBC_i: ', D_noBC_i)
				MieQ = ps.MieQ(m_shell_noBC_i, wl, D_noBC_i)
				ScatteringFunction = ps.ScatteringFunction(m_shell_noBC_i, wl, D_noBC_i, angularResolution=angularResolution, normalization='t')
				
				Qsca = MieQ[1]
				ksca += Qsca * 1/4 * np.pi * (D_noBC_i*1e-9)**2 * n_noBC_i*1e6 * MS_p[k]
				# separate peak making n separate too
				
				Qext = MieQ[0]
				kext += Qext * 1/4 * np.pi * (D_noBC_i*1e-9)**2 * n_noBC_i*1e6 * MS_p[k]
				#print('noBC: ',Qsca,' ',Qext)
				#print('noBC: ',Qsca/Qext)
				
				g += MieQ[3] * n_noBC_i / sum(n) * MS_p[k]
				
				theta = ScatteringFunction[0] # array of the angles used in calculations
				SU = ScatteringFunction[3] / sum(ScatteringFunction[3]) # array of scattered intensity of unpolarized light
				#print('pmom: ')
				#for i in range(moma):
				#	print(phaseFunc.beta(SU, theta, i+1))
				SU_bulk += SU * Qsca * 1/4 * np.pi * (D_noBC_i*1e-9)**2 * n_noBC_i*1e6  * MS_p[k]
				#SU += SU * n_noBC_i / n[i] / nI2x / kappaI2x
				# use scattering coefficient weighted to calculate bulk pmom
				# to normolize the result, need to devide total scattering coefficient at last
			
			if sum(n_onlyBC_i) > 0:
				for j in range(len(DBC)):
					# no hygrosocpic groth
					#print('DBC[j]: ', DBC[j])
					MieQ = ps.MieQ(m_BC, wl, DBC[j])
					ScatteringFunction = ps.ScatteringFunction(m_BC, wl, DBC[j], angularResolution=angularResolution, normalization='t')
					
					Qsca = MieQ[1]
					ksca += Qsca * 1/4 * np.pi * (DBC[j]*1e-9)**2 * n_onlyBC_i[j]*1e6 * MS_p[k]
					
					Qext = MieQ[0]
					kext += Qext * 1/4 * np.pi * (DBC[j]*1e-9)**2 * n_onlyBC_i[j]*1e6 * MS_p[k]
					#print('onlyBC: ',Qsca,' ',Qext)
					#print('onlyBC: ',Qsca/Qext)
					
					g += MieQ[3] * n_onlyBC_i[j] / sum(n) * MS_p[k]
					
					theta = ScatteringFunction[0]
					SU = ScatteringFunction[3] / sum(ScatteringFunction[3])
					#print('pmom: ')
					#for i in range(moma):
					#	print(phaseFunc.beta(SU, theta, i+1))
					SU_bulk += SU * Qsca * 1/4 * np.pi * (DBC[j]*1e-9)**2 * n_onlyBC_i[j]*1e6 * MS_p[k]
					#SU += SU * n_onlyBC_i[j] / n[i] / nI2x / kappaI2x
			
			# calculate internal mix part
			if sum(n_int_i) > 0:
				for j in range(len(DBC)):
					D_BC_ij, m_shell_BC_ij = kappaKohler.RH2D(Dps[i], DBC[j], m_BC, m_shell[i], kappa_shell[i], RH, wl)
					#print('D_BC_ij, Dps[i]: ', D_BC_ij, ' ', Dps[i])
					MieQCoreShell = ps.MieQCoreShell(m_BC, m_shell_BC_ij, wl, DBC[j], D_BC_ij*CT_shell[i])
					CoreShellScatteringFunction = ps.CoreShellScatteringFunction(m_BC, m_shell_BC_ij, wl, DBC[j], D_BC_ij*CT_shell[i], angularResolution=180/(180/angularResolution+1), normed=True)
					
					Qsca = MieQCoreShell[1]
					ksca += Qsca * 1/4 * np.pi * (D_BC_ij*1e-9)**2 * n_int_i[j]*1e6 * MS_p[k]
					
					Qext = MieQCoreShell[1] + MieQCoreShell[2] * BCAE
					kext += Qext * 1/4 * np.pi * (D_BC_ij*1e-9)**2 * n_int_i[j]*1e6 * MS_p[k]
					#print('internal: ',Qsca,' ',Qext)
					#print('internal: ',Qsca/Qext)
					
					g += MieQCoreShell[3] * n_int_i[j] / sum(n) * MS_p[k]
					
					theta = CoreShellScatteringFunction[0]
					SU = CoreShellScatteringFunction[3] / sum(CoreShellScatteringFunction[3])
					#print('pmom: ')
					#for i in range(moma):
					#	print(phaseFunc.beta(SU, theta, i+1))
					SU_bulk += SU * Qsca * 1/4 * np.pi * (D_BC_ij*1e-9)**2 * n_int_i[j]*1e6 * MS_p[k]
					#SU += SU * n_int_i[j] / n[i] / nI2x / kappaI2x
	
	SU = SU_bulk / ksca
	for i in range(moma):
		pmom[i] = phaseFunc.beta(SU, theta, i+1)
	waer = ksca / kext
	kext = kext * 1e6 # m^-1 to Mm^-1
	#print(sum(PNSD))
	#print(ksca * 1e6)
	#print(kext, waer, g)
	
	return kext, waer, pmom, g

def cal_Mie3(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, nShellH1, kappa, kappaH1, MS, MSH1, MSH2x, MSH2y, RH, wl, CT, CTH1, BCAE):
	'''
	This function is to use Dps, PNSD, RH, CT and BCAE to calculate extinct coefficient, 
	scattering coefficient and asymmetry factor.
	Add new heterogeneity parameter: MS, CT, nShell and kappa
	input:
		Dps      : diameter distribution, array, nm
		PNSD     : number distribution, array, dn/dlogDps, cm^-3
		DBCps    : BC particle diameter size, array, nm
		DBC      : BC core diameter size, array, nm
		BCPNSD   : BC particle number concentration, array, dn/dlogDBCps/dlogDBC, cm^-3
		kBC      : BC core imagine part of complex refractive index, float
		nBC      : BC core real part of complex refractive index, float
		nShell   : shell complex refractive index, float
		nShellH1 : n mixing state change rate by diameter, float
		kappa    : hygroscopicity parameter, float
		kappaH1  : hygroscopicity parameter mixing state change rate by diameter, float
		MS       : mixing state, float
		MSH1     : mixing state change rate by diameter, float
		MSH2x    : mixing state bin peak number, int
		MSH2y    : mixing state bin peak saperate rate, float
		RH       : relative humidity, float, percent
		wl       : wave length, nm, float
		CT       : BC core coating thickness adjust parameter, float
		CTH1     : BC core coating thickness change rate by diameter, float
		BCAE     : BC core absorbing enhancement adjust parameter, float
	output:
		kext     : extinct coefficient, Mm-1
		ksca     : scattering coefficient, Mm-1
		g        : asymmetry factor
	'''
	
	m_BC = nBC + kBC * 1j
	n = calPNSD.PNSD2n(Dps, PNSD, 'd') # in cm^-3
	dlogDps = calPNSD.cal_dlnDp(Dps)
	dlogDBC = calPNSD.cal_dlnDp(DBC)
	
	m_shell = cal_H1(nShell, len(Dps), nShellH1)
	#print('m : ', m_shell)
	kappa_shell = cal_H1(kappa, len(Dps), kappaH1)
	#print('kappa : ', kappa_shell)
	MS_shell, MS_p = cal_H2(MS, len(Dps), MSH1, MSH2x, MSH2y)
	#print('MS : ', MS_shell)
	CT_shell = cal_H1(CT, len(Dps), CTH1)
	#print('CT : ', CT_shell)
	
	kext = 0
	ksca = 0
	g = 0
	
	for i in range(len(Dps)):
		for k in range(MSH2x):
			DMASP2_i = calPNSD.DMASP2PNSD_Dp(DBCps, BCPNSD*dlogDBC, Dps[i])# * dlogDps # array in shape of DBC, in dn/dlogDBC
			nBC_i = DMASP2_i# * dlogDBC # in cm^-3
			# three parts: BC and coated material, only BC and only coated material, use MS to decide how much
			n_noBC_i = n[i] - sum(nBC_i) * MS_shell[i,k] # no BC part and external part
			n_onlyBC_i = nBC_i * (1-MS_shell[i,k]) # external part, in shape of (len(DBC))
			n_int_i = nBC_i * MS_shell[i,k] # internal part, in shape of (len(DBC))
			#print(round(Dps[i]), n_noBC_i, sum(n_onlyBC_i), sum(n_int_i))
			
			# calculate external part
			if n_noBC_i > 0:
				D_noBC_i, m_shell_noBC_i = kappaKohler.RH2D(Dps[i], 0, m_BC, m_shell[i], kappa_shell[i], RH, wl)
				#print(D_noBC_i, m_shell_noBC_i, Dps[i], m_shell[i])
				#print('D_noBC_i: ', D_noBC_i)
				MieQ = ps.MieQ(m_shell_noBC_i, wl, D_noBC_i)
				
				Qsca = MieQ[1]
				ksca += Qsca * 1/4 * np.pi * (D_noBC_i*1e-9)**2 * n_noBC_i*1e6 * MS_p[k]
				# separate peak making n separate too
				
				Qext = MieQ[0]
				kext += Qext * 1/4 * np.pi * (D_noBC_i*1e-9)**2 * n_noBC_i*1e6 * MS_p[k]
				#print('noBC: ',Qsca,' ',Qext)
				#print('noBC: ',Qsca/Qext)
				
				g += MieQ[3] * n_noBC_i / sum(n) * MS_p[k]
			
			if sum(n_onlyBC_i) > 0:
				for j in range(len(DBC)):
					# no hygrosocpic groth
					#print('DBC[j]: ', DBC[j])
					MieQ = ps.MieQ(m_BC, wl, DBC[j])
					
					Qsca = MieQ[1]
					ksca += Qsca * 1/4 * np.pi * (DBC[j]*1e-9)**2 * n_onlyBC_i[j]*1e6 * MS_p[k]
					
					Qext = MieQ[0]
					kext += Qext * 1/4 * np.pi * (DBC[j]*1e-9)**2 * n_onlyBC_i[j]*1e6 * MS_p[k]
					#print('onlyBC: ',Qsca,' ',Qext)
					#print('onlyBC: ',Qsca/Qext)
					
					g += MieQ[3] * n_onlyBC_i[j] / sum(n) * MS_p[k]
			
			# calculate internal mix part
			if sum(n_int_i) > 0:
				for j in range(len(DBC)):
					D_BC_ij, m_shell_BC_ij = kappaKohler.RH2D(Dps[i], DBC[j], m_BC, m_shell[i], kappa_shell[i], RH, wl)
					#print('D_BC_ij, Dps[i]: ', D_BC_ij, ' ', Dps[i])
					MieQCoreShell = ps.MieQCoreShell(m_BC, m_shell_BC_ij, wl, DBC[j], D_BC_ij*CT_shell[i])
					
					Qsca = MieQCoreShell[1]
					ksca += Qsca * 1/4 * np.pi * (D_BC_ij*1e-9)**2 * n_int_i[j]*1e6 * MS_p[k]
					
					Qext = MieQCoreShell[1] + MieQCoreShell[2] * BCAE
					kext += Qext * 1/4 * np.pi * (D_BC_ij*1e-9)**2 * n_int_i[j]*1e6 * MS_p[k]
					#print('internal: ',Qsca,' ',Qext)
					#print('internal: ',Qsca/Qext)
					
					g += MieQCoreShell[3] * n_int_i[j] / sum(n) * MS_p[k]
	
	kext = kext * 1e6 # m^-1 to Mm^-1
	ksca = ksca * 1e6
	print(kext, ksca, g)
	
	return kext, ksca, g

def cal_Mie2(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, nI1, nI2x, nI2y, kappa, kappaI1, kappaI2x, kappaI2y, MS, RH, wl, moma, CT, BCAE, **args):
	'''
	This function is to use Dps, PNSD, RH, CT and BCAE to calculate extinct coefficient, single scattering coefficient
	and Legendre moments of asymmetry factor
	input:
		Dps      : diameter distribution, array, nm
		PNSD     : number distribution, array, dn/dlogDps, cm^-3
		DBCps    : BC particle diameter size, array, nm
		DBC      : BC core diameter size, array, nm
		BCPNSD   : BC particle number concentration, array, dn/dlogDBCps/dlogDBC, cm^-3
		kBC      : BC core imagine part of complex refractive index, float
		nBC      : BC core real part of complex refractive index, float
		nShell   : shell complex refractive index, float
		nI1      : n mixing state change rate by diameter, float
		nI2x     : n mixing state size bin peak number, int
		nI2y     : n mixing state size bin peak separate rate, flaot
		kappa    : hygroscopicity parameter, float
		kappaI   : hygroscopicity parameter mixing state change rate by diameter, float
		kappaI2x : hygroscopicity parameter mixing state size bin peak number, int
		kappaI2y : hygroscopicity parameter mixing state size bin peak separate rate, float
		MS       : mixing state, float
		RH       : relative humidity, float, percent
		wl       : wave length, nm, float
		moma     : moment of Legendre moments, int
		CT       : BC core coating thickness adjust parameter, float
		BCAE     : BC core absorbing enhancement adjust parameter, float
		**args:
			angularResolution : angular resolution of PyMieScatt phase function, degree, defualt 0.5
	output:
		kext     : extinct coefficient, Mm-1, float
		waer     : single scattering coefficient, float
		pmom     : Legendre moments of asymmetry factor, np.array, in shape (moma)
		g        : asymmetry factor, float
	'''
	if 'angularResolution' in args:
		angularResolution = args['angularResolution']
	else:
		angularResolution = 0.5
	
	m_BC = nBC + kBC * 1j
	m_shell = cal_I2(nShell, len(Dps), nI1, nI2x, nI2y)
	kappa_shell = cal_I2(kappa, len(Dps), kappaI1, kappaI2x, kappaI2y)
	n = calPNSD.PNSD2n(Dps, PNSD, 'd') # in cm^-3
	dlogDps = calPNSD.cal_dlnDp(Dps)
	dlogDBC = calPNSD.cal_dlnDp(DBC)
	
	kext = 0
	ksca = 0
	g = 0
	theta = np.zeros(round(180/angularResolution+1))
	SU_bulk = np.zeros(len(theta))
	pmom = np.zeros(moma)
	
	for i in range(len(Dps)):
		DMASP2_i = calPNSD.DMASP2PNSD_Dp(DBCps, BCPNSD*dlogDBC, Dps[i])# * dlogDps # array in shape of DBC, in dn/dlogDBC
		nBC_i = DMASP2_i# * dlogDBC # in cm^-3
		# three parts: BC and coated material, only BC and only coated material, use MS to decide how much
		n_noBC_i = n[i] - sum(nBC_i) * MS # no BC part and external part
		n_onlyBC_i = nBC_i * (1-MS) # external part, in shape of (len(DBC))
		n_int_i = nBC_i * MS # internal part, in shape of (len(DBC))
		#print(n_noBC_i, sum(n_onlyBC_i), sum(n_int_i))
		
		# calculate external part
		if n_noBC_i > 0:
			for j in range(nI2x):
				for k in range(kappaI2x):
					D_noBC_i, m_shell_noBC_i = kappaKohler.RH2D(Dps[i], 0, m_BC, m_shell[i,j], kappa_shell[i,k], RH, wl)
					#print(D_noBC_i, m_shell_noBC_i, Dps[i], m_shell[i,j])
					#print('D_noBC_i: ', D_noBC_i)
					MieQ = ps.MieQ(m_shell_noBC_i, wl, D_noBC_i)
					ScatteringFunction = ps.ScatteringFunction(m_shell_noBC_i, wl, D_noBC_i, angularResolution=angularResolution, normalization='t')
					
					Qsca = MieQ[1]
					ksca += Qsca * 1/4 * np.pi * (D_noBC_i*1e-9)**2 * n_noBC_i*1e6 / nI2x / kappaI2x
					# separate peak making n separate too
					
					Qext = MieQ[0]
					kext += Qext * 1/4 * np.pi * (D_noBC_i*1e-9)**2 * n_noBC_i*1e6 / nI2x / kappaI2x
					#print('noBC: ',Qsca,' ',Qext)
					#print('noBC: ',Qsca/Qext)
					
					g += MieQ[3] * n_noBC_i / sum(n) / nI2x / kappaI2x
					
					theta = ScatteringFunction[0] # array of the angles used in calculations
					SU = ScatteringFunction[3] / sum(ScatteringFunction[3]) # array of scattered intensity of unpolarized light
					#print('pmom: ')
					#for i in range(moma):
					#	print(phaseFunc.beta(SU, theta, i+1))
					SU_bulk += SU * Qsca * 1/4 * np.pi * (D_noBC_i*1e-9)**2 * n_noBC_i*1e6 / nI2x / kappaI2x
					#SU += SU * n_noBC_i / n[i] / nI2x / kappaI2x
					# use scattering coefficient weighted to calculate bulk pmom
					# to normolize the result, need to devide total scattering coefficient at last
		
		if sum(n_onlyBC_i) > 0:
			for j in range(len(DBC)):
				for k in range(nI2x):
					for l in range(kappaI2x):
						# no hygrosocpic groth
						#print('DBC[j]: ', DBC[j])
						MieQ = ps.MieQ(m_BC, wl, DBC[j])
						ScatteringFunction = ps.ScatteringFunction(m_BC, wl, DBC[j], angularResolution=angularResolution, normalization='t')
						
						Qsca = MieQ[1]
						ksca += Qsca * 1/4 * np.pi * (DBC[j]*1e-9)**2 * n_onlyBC_i[j]*1e6 / nI2x / kappaI2x
						
						Qext = MieQ[0]
						kext += Qext * 1/4 * np.pi * (DBC[j]*1e-9)**2 * n_onlyBC_i[j]*1e6 / nI2x / kappaI2x
						#print('onlyBC: ',Qsca,' ',Qext)
						#print('onlyBC: ',Qsca/Qext)
						
						g += MieQ[3] * n_onlyBC_i[j] / sum(n) / nI2x / kappaI2x
						
						theta = ScatteringFunction[0]
						SU = ScatteringFunction[3] / sum(ScatteringFunction[3])
						#print('pmom: ')
						#for i in range(moma):
						#	print(phaseFunc.beta(SU, theta, i+1))
						SU_bulk += SU * Qsca * 1/4 * np.pi * (DBC[j]*1e-9)**2 * n_onlyBC_i[j]*1e6 / nI2x / kappaI2x
						#SU += SU * n_onlyBC_i[j] / n[i] / nI2x / kappaI2x
		
		# calculate internal mix part
		if sum(n_int_i) > 0:
			for j in range(len(DBC)):
				for k in range(nI2x):
					for l in range(kappaI2x):
						D_BC_ij, m_shell_BC_ij = kappaKohler.RH2D(Dps[i], DBC[j], m_BC, m_shell[i,k], kappa_shell[i,l], RH, wl)
						#print('D_BC_ij, DBC[j]: ', D_BC_ij, ' ', DBC[j])
						MieQCoreShell = ps.MieQCoreShell(m_BC, m_shell_BC_ij, wl, DBC[j], D_BC_ij*CT)
						CoreShellScatteringFunction = ps.CoreShellScatteringFunction(m_BC, m_shell_BC_ij, wl, DBC[j], D_BC_ij*CT, angularResolution=180/(180/angularResolution+1), normed=True)
						
						Qsca = MieQCoreShell[1]
						ksca += Qsca * 1/4 * np.pi * (D_BC_ij*1e-9)**2 * n_int_i[j]*1e6 / nI2x / kappaI2x
						
						Qext = MieQCoreShell[1] + MieQCoreShell[2] * BCAE
						kext += Qext * 1/4 * np.pi * (D_BC_ij*1e-9)**2 * n_int_i[j]*1e6 / nI2x / kappaI2x
						#print('internal: ',Qsca,' ',Qext)
						#print('internal: ',Qsca/Qext)
						
						g += MieQCoreShell[3] * n_int_i[j] / sum(n) / nI2x / kappaI2x
						
						theta = CoreShellScatteringFunction[0]
						SU = CoreShellScatteringFunction[3] / sum(CoreShellScatteringFunction[3])
						#print('pmom: ')
						#for i in range(moma):
						#	print(phaseFunc.beta(SU, theta, i+1))
						SU_bulk += SU * Qsca * 1/4 * np.pi * (D_BC_ij*1e-9)**2 * n_int_i[j]*1e6 / nI2x / kappaI2x
						#SU += SU * n_int_i[j] / n[i] / nI2x / kappaI2x
	
	SU = SU_bulk / ksca
	for i in range(moma):
		pmom[i] = phaseFunc.beta(SU, theta, i+1)
	waer = ksca / kext
	kext = kext * 1e6 # m^-1 to Mm^-1
	#print(sum(PNSD))
	#print(ksca * 1e6)
	
	return kext, waer, pmom, g

def cal_Mie(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, nI1, nI2x, nI2y, kappa, kappaI1, kappaI2x, kappaI2y, RH, wl, CT, BCAE):
	'''
	This function is to use Dps, PNSD, RH, CT and BCAE to calculate extinct coefficient, scattering coefficient and asymmetry factor g
	input:
		Dps      : diameter distribution, array, nm
		PNSD     : number distribution, array, dn/dlogDp
		DBCps    : BC particle diameter size, array, nm
		DBC      : BC core diameter size, array, nm
		BCPNSD   : BC particle number concentration, array, cm^-3
		kBC      : BC core imagine part of complex refractive index
		nBC      : BC core real part of complex refractive index
		nShell   : shell complex refractive index
		nI1      : n mixing state change rate by diameter
		nI2x     : n mixing state size bin peak number
		nI2y     : n mixing state size bin peak separate rate
		kappa    : hygroscopicity parameter, float
		kappaI   : hygroscopicity parameter mixing state change rate by diameter
		kappaI2x : hygroscopicity parameter mixing state size bin peak number
		kappaI2y : hygroscopicity parameter mixing state size bin peak separate rate
		RH       : relative humidity, float, percent
		wl       : wave length, nm
		CT       : BC core coating thickness adjust parameter
		BCAE     : BC core absorbing enhancement adjust parameter
	output:
		kext     : extinct coefficient, Mm-1
		ksca     : scattering coefficient, Mm-1
		g        : asymmetry factor
	'''
	m_BC = nBC + kBC * 1j
	#m_shell = nShell
	m_shell = cal_I2(nShell, len(Dps), nI1, nI2x, nI2y)
	kappa_shell = cal_I2(kappa, len(Dps), kappaI1, kappaI2x, kappaI2y)
	'''
	kappa_shell = np.zeros(len(Dps))
	for i in range(len(Dps)):
		kappa_shell[i] = kappa + kappaI * (i-len(Dps)/2) / (len(Dps)/2)
	'''
	n = calPNSD.PNSD2n(Dps, PNSD, 'd')
	BCPNSD = BCPNSD * calPNSD.cal_dlogDp(DBC)
	kext = 0
	ksca = 0
	g = 0
	
	for i in range(len(Dps)):
		BCPNSD_i = calPNSD.DMASP2PNSD_Dp(DBCps, BCPNSD, Dps[i]) # array in shape of DBC
		n_noBC_i = n[i] - sum(BCPNSD_i)
		
		if n_noBC_i > 0:
			for j in range(nI2x):
				for k in range(kappaI2x):
					D_noBC_i, m_shell_noBC_i = kappaKohler.RH2D(Dps[i], 0, m_BC, m_shell[i,j], kappa_shell[i,k], RH, wl)
					MieQ = ps.MieQ(m_shell_noBC_i, wl, D_noBC_i)
					Qsca_noBC_i = MieQ[1]
					ksca += Qsca_noBC_i * 1/4 * np.pi * (D_noBC_i*1e-9)**2 * n_noBC_i*1e6 / nI2x / kappaI2x # separate peak making n separate too
					Qext_noBC_i = MieQ[0]
					kext += Qext_noBC_i * 1/4 * np.pi * (D_noBC_i*1e-9)**2 * n_noBC_i*1e6 / nI2x / kappaI2x
					g_noBC_i = MieQ[3]
					g += g_noBC_i * n_noBC_i / sum(n) / nI2x / kappaI2x
		
		if sum(BCPNSD_i):
			for j in range(len(DBC)):
				for k in range(nI2x):
					for l in range(kappaI2x):
						D_BC_ij, m_shell_BC_ij = kappaKohler.RH2D(Dps[i], DBC[j], m_BC, m_shell[i,k], kappa_shell[i,l], RH, wl)
						MieQCoreShell = ps.MieQCoreShell(m_BC, m_shell_BC_ij, wl, DBC[j], D_BC_ij*CT)
						Qsca_BC_ij = MieQCoreShell[1]
						ksca += Qsca_BC_ij * 1/4 * np.pi * (D_BC_ij*1e-9)**2 * BCPNSD_i[j]*1e6 / nI2x / kappaI2x
						Qext_BC_ij = Qsca_BC_ij + MieQCoreShell[2] * BCAE
						kext += Qext_BC_ij * 1/4 * np.pi * (D_BC_ij*1e-9)**2 * BCPNSD_i[j]*1e6 / nI2x / kappaI2x
						g_BC_ij = MieQCoreShell[3]
						g += g_BC_ij * BCPNSD_i[j] / sum(n) / nI2x / kappaI2x
	
	kext = kext * 1e6
	ksca = ksca * 1e6
	return kext, ksca, g

def cal_total_AOD(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, nI1, nI2x, nI2y, kappa, kappaI1, kappaI2x, kappaI2y, CT, BCAE, VD, **args):
	'''
	This function is to calculate total atmospere air AOD
	input:
	output:
	'''
	if 'dz' in args:
		dz = args['dz']
	else:
		dz = z / 100
	
	# read atmosphere data from atms.dat
	atms = readAtms.read()
	z = atms['z'] * 1e3 # in m, from top to bottom
	dz = z[:-1] - z[1:] # dz in each two layers, in m
	wh = atms['wh']
	t = atms['t']
	RH = np.zeros(len(wh))
	for i in range(len(RH)):
		RH[i] = calRH.wh2RH(t[i], wh[i])
		if RH[i]<1:
			RH[i] = 1 # minimum RH set to 1
	
	AOD = 0
	for i in range(dz): # every level of atmosphere
		dh = dz * calPNSD.cal_VD(Dps, PNSD, h, VD)
		AOD += cal_Mie(Dps, PNSD, DBCps, DBC, BCPNSD, kBC, nBC, nShell, nI1, nI2x, nI2y, kappa, kappaI1, kappaI2x, kappaI2y, RH, wl, CT, BCAE)[0] * 1e-6 * dh
	
	return AOD

if __name__ == '__main__':
	'''
	# for Qiujie
	import readTaizhou
	import re
	sp2 = readTaizhou.read_Taizhou('data/sp2/Taizhou.npy')
	DBC = sp2['DBC']
	DBCps = sp2['DBCps']
	BCPNSD = np.zeros(sp2['DMASP2'][0].shape) # no BC core
	
	Dps = []
	PNSD = []
	
	with open('data/AOT/PNSD.txt', 'r') as f:
		infos = f.readlines()
		for line in infos:
			res = re.split('\t', line[:-1])
			Dps.append(float(res[0]))
			PNSD.append(float(res[1]))
	
	Dps = np.array(Dps)
	PNSD = np.array(PNSD)
	
	nShell = [1.484, 1.51, 1.54]
	kappa = 1.2
	RH = 40
	wl = 500
	
	ksca = np.zeros(len(nShell))
	g = np.zeros(len(nShell))
	
	for i in range(len(nShell)):
		kext, ksca[i], g[i] = cal_Mie(Dps, PNSD, DBCps, DBC, BCPNSD, 1, 1, nShell[i], 1, 1, 1, kappa, 1, 1, 1, RH, wl, 1, 1)
		
	with open('result.txt', 'w') as f:
		f.write('RI\tksca\tg\n')
		for i in range(len(nShell)):
			f.write(str(round(nShell[i]))+'\t'+str(ksca[i])+'\t'+str(g[i])+'\n')
	'''
	val_new, val_p = cal_H2(5,10,0.1,5,0.1)
	print(val_p)
	print(val_new)
	val_new = cal_H1(5,10,0.1)
	print(val_new)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
