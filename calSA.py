####################################################################################
# INTRODUCTION:
# This code is to calculate runDARF outputs' sensitivity
# Created by Hebs at 24/3/10/10:49
# Contact: hebishuo@pku.edu.cn
####################################################################################

default_range = 0.3

run_paras = [] # default parameters list
run_paras.append('n')
run_paras.append('nH1')
run_paras.append('nBC')
run_paras.append('kBC')
run_paras.append('PNSD')
run_paras.append('MS')
run_paras.append('MSH1')
#run_paras.append('MSH2')
run_paras.append('VD')
run_paras.append('CT')
run_paras.append('CTH1')
run_paras.append('kappa')
run_paras.append('kappaH1')
run_paras.append('rhoBC')
run_paras.append('BCPNSD')
run_paras.append('BCAE')
run_paras.append('amb')
run_paras.append('albedo')

default_paras = [] # default parameters list
default_paras.append('n')
default_paras.append('PNSD')
default_paras.append('kappa')

default_paras.append('nBC')
default_paras.append('kBC')
default_paras.append('rhoBC')
default_paras.append('BCPNSD')

default_paras.append('MS')
default_paras.append('CT')
default_paras.append('BCAE')

default_paras.append('VD')
default_paras.append('amb')
default_paras.append('albedo')

default_paras.append('nH1')
default_paras.append('kappaH1')
default_paras.append('CTH1')
default_paras.append('MSH1')
#default_paras.append('MSH2')

default_outputs = [] # default outputs list
default_outputs.append('RF_top')
default_outputs.append('RF_bot')
default_outputs.append('dtau_bot')
default_outputs.append('waer_bot')
default_outputs.append('g_bot')

import numpy as np
import os
import runDARF
import readTaizhou
import calPNSD

def run_OAT(dtime, par_range, **args):
	'''
	This function is to calculate OAT sensitivity of runDARF
	input:
		dtime		: runDARF.py output time, string, example: '230101'
		par_range	: parameter change rate range, float
		**path		: runDARF.py output save path, string, default 'output/DARF/'
		**debug		: debug flat, bool, default False
	output:
		runDARF.py output
	'''
	if 'path' in args:
		path = args['path']
	else:
		path = 'output/DARF/'
	if 'debug' in args:
		debug = args['debug']
	else:
		debug = False
	
	path = path + dtime + '/OAT/'
	
	# read data
	
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
	
	clean_PNSD = PNSD / 1000 # to fit Beijing aerosol number distribution level
	dirty_PNSD = PNSD / 20
	clean_BCPNSD = BCPNSD / 50
	dirty_BCPNSD = BCPNSD / 10
	
	# set parameters running list
	
	parameters = run_paras
	
	par_num = len(default_paras)
	paras_rate = np.zeros((par_num, par_num, 2))
	
	for i in range(par_num):
		for j in range(par_num):
			if i==j:
				paras_rate[i,j] = np.array([1-par_range, 1+par_range])
			else:
				paras_rate[i,j] = np.ones(2)
	
	# make directory for data saving
	
	os.system('mkdir '+path)
	os.system('mkdir '+path+'all/')
	for i in range(par_num):
		fn = path + parameters[i]
		os.system('mkdir '+fn)
	
	# running
	
	if debug:
		print('done\ncalculating...')
	
	for i in range(par_num):
		infos = []
		output_path = path + parameters[i] + '/'
		for j in range(2):
			output_names = ['0_'+str(j)+'.txt', '1_'+str(j)+'.txt', 'atms_'+str(j)+'.dat', 'aerosol_'+str(j)+'.dat', 'albedo_'+str(j)+'.dat']
			info, paras = runDARF.run([525], 15, 6,
			n_rate			= paras_rate[i,0,j], 
			nH1_rate		= paras_rate[i,1,j], 
			nBC_rate		= paras_rate[i,2,j], 
			kBC_rate		= paras_rate[i,3,j], 
			PNSD_rate		= paras_rate[i,4,j], 
			MS_rate			= paras_rate[i,5,j], 
			MSH1_rate		= paras_rate[i,6,j], 
			VD_rate			= paras_rate[i,7,j], 
			CT_rate			= paras_rate[i,8,j], 
			CTH1_rate		= paras_rate[i,9,j], 
			kappa_rate		= paras_rate[i,10,j], 
			kappaH1_rate	= paras_rate[i,11,j], 
			rhoBC_rate		= paras_rate[i,12,j], 
			BCPNSD_rate		= paras_rate[i,13,j], 
			BCAE_rate		= paras_rate[i,14,j], 
			amb_rate		= paras_rate[i,15,j], 
			albedo_rate		= paras_rate[i,16,j], 
			Dps=Dps, DBC=DBC, DBCps=DBCps, clean_PNSD=clean_PNSD, dirty_PNSD=dirty_PNSD, clean_BCPNSD=clean_BCPNSD, dirty_BCPNSD=dirty_BCPNSD, angularResolution=30, debug=debug, output_path=output_path, output_names=output_names)
			infos.append(info)
			np.save(path+parameters[i]+'/infos.npy', infos)
			if i==0 and j==0:
				np.save(path+'all/paras.npy', dict(paras=paras,paras_rate=paras_rate))
			if debug:
				print(round((i*2+j+1)/par_num/2*100), '% done...')
				
	if debug:
		print('done')

def read_OAT(dtime, **args):
	'''
	This function is to read OAT running results
	input:
		dtime		: runDARF.py output time, string, example: '230101'
		**path		: runDARF.py output save path, string, default 'output/DARF/'
		**paras		: parameters list, array, string, default default_paras
		**outputs	: outputs list, array, string, default default_outputs
	output:
		in dictionary
		paras	: parameters list, array, string
		outputs	: outputs list, array, string
		OAT		: OAT reasults, array, in shape(len(outputs),len(paras))
	'''
	if 'path' in args:
		path = args['path']
	else:
		path = 'output/DARF/'
	if 'paras' in args:
		paras = args['paras']
	else:
		paras = default_paras
	if 'outputs' in args:
		outputs = args['outputs']
	else:
		outputs = default_outputs
	
	path = path + dtime + '/OAT/'
	
	# read data
	
	paras_data = np.load(path+'all/paras.npy', allow_pickle=True).item()
	paras_rate = paras_data['paras_rate']
	par_num = len(default_paras)
	out_num = len(outputs)
	
	OAT = np.zeros((out_num,par_num))
	
	for i in range(par_num):
		if default_paras[i] in paras:
			data_path = path + default_paras[i] + '/'
			data = np.load(data_path+'infos.npy', allow_pickle=True)
			
			for j in range(out_num):
				OAT[j,i] = abs((data[1][outputs[j]]-data[0][outputs[j]])/(paras_rate[i,i,1]-paras_rate[i,i,0])) + 1e-9
	
	OAT = OAT[np.where(OAT!=0)].reshape(out_num,-1)
	infos = dict(paras=paras, outputs=outputs, OAT=OAT)
	return infos

def cal_CP3(dtime, par_range, **args):
	'''
	This function is to calculate Constraint Parameter sensitivity of runDARF,
	but for list and wavelength divided outputs, including dtau, waer, g
	input:
		dtime		: runDARF.py output time, string, example: '230101'
		par_range	: parameter change rate range, float
		**path		: runDARF.py output save path, string, default 'output/DARF/'
		**paras		: parameters list, array, string, default default_paras
		**outputs	: outputs list, array, string, default ['dtau','waer','g']
	output:
		in dictionary
		paras		: parameters list, array, string
		run_num		: running number, int
		CP			: Constraint Parameter Sensitivity Analysis result, array, float,
				  	  in shape (outputs_num, paras_num), from top to bottom
	'''
	if 'path' in args:
		path = args['path']
	else:
		path = 'output/DARF/'
	if 'paras' in args:
		paras = args['paras']
	else:
		paras = default_paras
	if 'outputs' in args:
		outputs = args['outputs']
	else:
		outputs = ['dtau', 'waer', 'g']
	
	path = path + dtime + '/'
	
	# read data
	data1_path = path + 'all/'
	data1 = np.load(data1_path+'infos.npy', allow_pickle=True)
	run_num = len(data1)
	par_num = len(default_paras)
	out_num = len(outputs)
	nn = len(data1[0][outputs[0]][0]) # the level number of SBDART running outputs
	
	out_datas1 = np.zeros((out_num,nn,run_num))
	
	for i in range(run_num):
		for j in range(nn):
			for k in range(out_num):
				out_datas1[k,j,i] = data1[i][outputs[k]][0][j]
	
	# read factor values
	x = np.load(data1_path+'paras.npy', allow_pickle=True).item()
	
	CP = np.zeros((out_num,nn,par_num))
	
	for i in range(par_num):
		if default_paras[i] in paras:
			data2_path = path + default_paras[i] + '/'
			data2 = np.load(data2_path+'infos.npy', allow_pickle=True)
			
			out_datas2 = np.zeros((out_num,nn,run_num))
			
			for j in range(run_num):
				for k in range(nn):
					for l in range(out_num):
						out_datas2[l,k,j] = data2[j][outputs[l]][0][k]
			
			for j in range(out_num):
				for k in range(nn):
					CP[j,k,i] = np.sqrt(abs((np.nanvar(out_datas1[j,k])-np.nanvar(out_datas2[j,k]))/0.33**2/(default_range**2-par_range**2)))
		else:
			CP[:,:,i] = np.ones((out_num,nn)) * 99999 # mask value
	
	CP = CP[np.where(CP!=99999)].reshape(out_num,nn,-1)
	infos = dict(run_num=run_num, paras=paras, outputs=outputs, CP=CP, nn=nn)
	return infos

def cal_CP2(dtime, par_range, **args):
	'''
	This function is to calculate Constraint Parameter sensitivity of runDARF,
	but for list outputs, including RF, heat
	input:
		dtime		: runDARF.py output time, string, example: '230101'
		par_range	: parameter change rate range, float
		**path		: runDARF.py output save path, string, default 'output/DARF/'
		**paras		: parameters list, array, string, default default_paras
		**outputs	: outputs list, array, string, default ['RF','heat1','heat0']
	output:
		in dictionary
		paras		: parameters list, array, string
		run_num		: running number, int
		CP			: Constraint Parameter Sensitivity Analysis result, array, float,
				  	  in shape (outputs_num, paras_num)
	'''
	if 'path' in args:
		path = args['path']
	else:
		path = 'output/DARF/'
	if 'paras' in args:
		paras = args['paras']
	else:
		paras = default_paras
	if 'outputs' in args:
		outputs = args['outputs']
	else:
		outputs = ['RF', 'heat1', 'heat0']
	
	path = path + dtime + '/'
	
	# read data
	data1_path = path + 'all/'
	data1 = np.load(data1_path+'infos.npy', allow_pickle=True)
	run_num = len(data1)
	par_num = len(default_paras)
	out_num = len(outputs)
	nn = len(data1[0][outputs[0]]) # the level number of SBDART running outputs
	
	out_datas1 = np.zeros((out_num,nn,run_num))
	
	for i in range(run_num):
		for j in range(nn):
			for k in range(out_num):
				out_datas1[k,j,i] = data1[i][outputs[k]][j]
	
	# read factor values
	x = np.load(data1_path+'paras.npy', allow_pickle=True).item()
	
	CP = np.zeros((out_num,nn,par_num))
	
	for i in range(par_num):
		if default_paras[i] in paras:
			data2_path = path + default_paras[i] + '/'
			data2 = np.load(data2_path+'infos.npy', allow_pickle=True)
			
			out_datas2 = np.zeros((out_num,nn,run_num))
			
			for j in range(run_num):
				for k in range(nn):
					for l in range(out_num):
						out_datas2[l,k,j] = data2[j][outputs[l]][k]
			
			for j in range(out_num):
				for k in range(nn):
					CP[j,k,i] = np.sqrt(abs((np.nanvar(out_datas1[j,k])-np.nanvar(out_datas2[j,k]))/0.33**2/(default_range**2-par_range**2)))
		else:
			CP[:,:,i] = np.ones((out_num,nn)) * 99999 # mask value
	
	CP = CP[np.where(CP!=99999)].reshape(out_num,nn,-1)
	infos = dict(run_num=run_num, paras=paras, outputs=outputs, CP=CP, nn=nn)
	return infos

def cal_CP(dtime, par_range, **args):
	'''
	This function is to calculate Constraint Parameter sensitivity of runDARF
	input:
		dtime		: runDARF.py output time, string, example: '230101'
		par_range	: parameter change rate range, float
		**path		: runDARF.py output save path, string, default 'output/DARF/'
		**paras		: parameters list, array, string, default default_paras
		**outputs	: outputs list, array, string, default default_outputs
	output:
		in dictionary
		paras		: parameters list, array, string
		run_num		: running number, int
		CP			: Constraint Parameter Sensitivity Analysis result, array, float,
				  	  in shape (outputs_num, paras_num)
	'''
	if 'path' in args:
		path = args['path']
	else:
		path = 'output/DARF/'
	if 'paras' in args:
		paras = args['paras']
	else:
		paras = default_paras
	if 'outputs' in args:
		outputs = args['outputs']
	else:
		outputs = default_outputs
	
	path = path + dtime + '/'
	
	# read data
	data1_path = path + 'all/'
	data1 = np.load(data1_path+'infos.npy', allow_pickle=True)
	run_num = len(data1)
	par_num = len(default_paras)
	out_num = len(outputs)
	
	out_datas1 = np.zeros((out_num,run_num))
	
	for i in range(run_num):
		for j in range(out_num):
			out_datas1[j,i] = data1[i][outputs[j]]
	
	# read factor values
	x = np.load(data1_path+'paras.npy', allow_pickle=True).item()
	
	CP = np.zeros((out_num,par_num))
	
	for i in range(par_num):
		if default_paras[i] in paras:
			data2_path = path + default_paras[i] + '/'
			data2 = np.load(data2_path+'infos.npy', allow_pickle=True)
			
			out_datas2 = np.zeros((out_num,run_num))
			
			for j in range(run_num):
				for k in range(out_num):
					out_datas2[k,j] = data2[j][outputs[k]]
			
			for j in range(out_num):
				CP[j,i] = np.sqrt(abs((np.nanvar(out_datas1[j])-np.nanvar(out_datas2[j]))/0.33**2/(default_range**2-par_range**2)))
		else:
			CP[:,i] = np.ones(out_num) * 99999 # mask value
	
	CP = CP[np.where(CP!=99999)].reshape(out_num,-1)
	infos = dict(run_num=run_num, paras=paras, outputs=outputs, CP=CP)
	return infos

if __name__ == '__main__':
	
	run_OAT('241010', 0.1, debug=True)
	infos = read_OAT('241010')
	print(infos['OAT'][0])
	'''
	infos = cal_CP('241012', 0.1)
	print(default_paras)
	print(infos['outputs'])
	print(infos['CP'][0])
	
	infos = cal_CP('241009', 0.05)
	print(infos['CP'][0])
	infos = read_OAT('241009')
	print(infos['OAT'][0])
	
	infos = cal_CP2('240310', 0.05)
	print(infos['outputs'])
	print(infos['CP'][:,0,:]) # 0 for top, -1 for bottom
	infos = cal_CP3('240310', 0.05)
	print(infos['outputs'])
	print(infos['CP'][:,0,:])
	print(default_paras)
	'''
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
