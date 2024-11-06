####################################################################################
# INTRODUCTION:
# This code is to make INPUT for SBDART
# Created by Hebs at 21/11/8/15:30
# Contact: hebishuo@pku.edu.cn
####################################################################################

import numpy as np
import os
import sys

def write(**args):
	'''
	This function is to write INPUT for SBDART
	input:
		**fn       : output file path, string, default ./INPUT
		**wlinf    : lower wavelength limit when ISAT=0, > 0.250 microns, float, default 0.47
		**wlsup    : upper wavelength limit when ISAT=0, < 100 microns, float, default 2.1
		**wlinc    : this parameter specifies the spectral resolution of the SBDART run, float, default 0.01
		**ngrid    : sets the number of grid points, int, default 50
		**zgrid1   : controls the resolution near the bottom of the grid, default 0.1
		**zgrid2   : sets the maximum permissible step size (at the top of the grid), default 1
		**nstr     : number of computational zenith angles used, int and must be diviseble by 2, default 4
		**iday     : the number of days into a standard year, int, default 0
		**time     : UTC time (Grenwich) in decimal hours, int, default 0
		**alat     : latitude of point on earth's surface, int, default 0
		**alon     : east longitude of point on earth's surface, int, default 0
			time, alat and alon are ignored if iday equals to 0
		**csza     : cosine of solar zenith angle, float, default -1
		**sca      : solar zenith angle, float, default 0
			sza is ignored if csza is non-negative or iday is non-zero
		**saza     : solar azimuth angle, float, default 0
			saza is ignored if iday is non-zero
		
		**idatm    : atmospheric profile, int, default 0
		                             default          ozone(atm-cm)
                ------------------      -------------------  -------------
                0 User Specified
                1 TROPICAL                     4.117         0.253   .0216
                2 MID-LATITUDE SUMMER          2.924         0.324   .0325
                3 MID-LATITUDE WINTER          0.854         0.403   .0336
                4 SUB-ARCTIC SUMMER            2.085         0.350   .0346
                5 SUB-ARCTIC WINTER            0.418         0.486   .0340
                6 US62                         1.418         0.349   .0252
                **iaer    : boundary layer aerosol type selector, int, default -1
                -1-read aerosol optical depth and scattering parameters
                  from aerosol.dat.
                0-no boundary layer aerosols (all BLA parameters ignored)
                1-rural       
                2-urban       
                3-oceanic     
                4-tropospheric
                5-user defined spectral dependence of BLA 
		
		**isalb   : surface albedo feature, int, default -1
                -1 -spectral surface albedo read from "albedo.dat"
                 0 -user specified, spectrally uniform albedo set with ALBCON
                 1 -snow              
                 2 -clear water       
                 3 -lake water         
                 4 -sea  water
                 5 -sand           (data range 0.4 - 2.3um)
                 6 -vegetation     (data range 0.4 - 2.6um)
                 7 -ocean water brdf, includes bio-pigments, foam, and sunglint
                    additional input parameters provided in SC.
                 8 -Hapke analytic brdf model. additional input parameters
                    provided in SC. 
                 9 -Ross-thick Li-sparse brdf model. additional input
                    parameters in SC.
		
		**iout     : standard output selector, int, default 11
		iout=0:
			no standard output is produced, DISORT subroutine is not called, but diagnostics selected by idb in gas absorption or aerosol subroutines are active.
		iout=1:
			one output record for each wavelength.
			wl: wavelength, um
			ffv: filter function value
			topdn: total downward flux at ZOUT(2) km, w/m2
			topup: total upward flux at ZOUT(2) km, w/m2
			topdir: direct downward flux at ZOUT(2) km, w/m2
			botdn: total downward flux at ZOUT(1) km, w/m2
			botup: total upward flux at ZOUT(1) km, w/m2
			botdir: direct downward flux at ZOUT(1) km, w/m2
		iout=5:
			nzen+3 records for each wavelength
			wl: wavelength, um
			ffv: filter function value
			topdn: total downward flux at ZOUT(2) km, w/m2
			topup: total upward flux at ZOUT(2) km, w/m2
			topdir: direct downward flux at ZOUT(2) km, w/m2
			botdn: total downward flux at ZOUT(1) km, w/m2
			botup: total upward flux at ZOUT(1) km, w/m2
			botdir: direct downward flux at ZOUT(1) km, w/m2
			nphi: number of user azimuth angles
			nzen: number of user zenith angles
			phi: user specified azimuth angles, degree
			uzen: user specified zenith angles, degree
			vzen: user specified nadir angles, degrees
			uurs: radiance at user angles at altitude ZOUT(2) , w/m2/um/str
		iout=6:
			same as iout=5 except radiance is for ZOUT(1) altitude
		iout=7:
			radiative flux at each layer for each wavelength.  This output option can produce a huge amount of output if many wavelength sample points are used.
			radiative flux at each layer for each wavelength
			fzw: block id
			nz: number of z levels
			nw: number of wavelengths
			z: altitude, km
			fdird: downward direct flux, w/m2/um
			fdifd: downward diffuse flux, w/m2/um
			flxdn: total downward flux, w/m2/um
			flxup: total upward flux, w/m2/um
		iout=10:
			out output record per run, integrated over wavelength.
			wlinf, wlsup, ffew, topdn, topup, topdir, botdn, botup, botdir
			ffew: filter function equivalent width, um
			topdn: total downward flux at ZOUT(2) km, w/m2
			topup: total upward flux at ZOUT(2) km, w/m2
			topdir: direct downward flux at ZOUT(2) km, w/m2
			botdn: total downward flux at ZOUT(1) km, w/m2
			botup: total upward flux at ZOUT(1) km, w/m2
			botdir: direct downward flux at ZOUT(1) km, w/m2
		iout=11:
			radiant fluxes at each atmospheric layer integrated over wavelength.
			nz: number of atmospheric layers
			phidw: filter function equivalent width, um
			zz: level altitudes, km
			pp: level pressure, mb
			fxdn: downward flux (direct+diffuse), W/m2
			fxup: upward flux, W/m2
			fxdir: downward flux, direct beam only, W/m2
			dfdz: radiant energy flux divergence, mW/m3
			heat: heating rate, K/day
		iout=20:
			radiance output at ZOUT(2) km.
			wlinf, wlsup, ffew, topdn, topup, topdir, botdn, botup, botdir
			ffew: filter function equivalent width, um
			topdn: total downward flux at ZOUT(2) km, w/m2
			topup: total upward flux at ZOUT(2) km, w/m2
			topdir: direct downward flux at ZOUT(2) km, w/m2
			botdn: total downward flux at ZOUT(1) km, w/m2
			botup: total upward flux at ZOUT(1) km, w/m2
			botdir: direct downward flux at ZOUT(1) km, w/m2
			nphi: number of user azimuth angles
			nzen: number of user zenith angles
			phi: user relative azimuth angles (nphi values)
			uzen: user zenith angles (nzen values)
			r: radiance array (nphi,nzen), W/m2/sr
		iout=21:
			same as iout=20 except radiance output at ZOUT(1) km
		iout=22:
			radiance and flux at each atmospheric layer integrated over wavelength.
			nphi: number of user specified azimuth angles
			nzen: number of user specified zenith angles
			nz: number of atmospheric levels
			ffew: filter function equivalent width, um
			phi: user specified anizmuth angles, degrees
			uzen: user specified zenith angles, degrees
			z: altitudes of atmospheric layers, km
			fxdn: downward flux (direct+diffuse), W/m2
			fxup: upward flux, W/m2
			fxdir: downward flux, direct beam only, W/m2
			uurl: radiance at each layer, W/m2/str
	output:
		INPUT for SBDART
	'''
	if 'fn' in args:
		fn = args['fn']
	else:
		fn = 'INPUT'
	if 'wlinf' in args:
		wlinf = args['wlinf']
	else:
		wlinf = 0.47
	if 'wlsup' in args:
		wlsup = args['wlsup']
	else:
		wlsup = 2.1
	if 'wlinc' in args:
		wlinc = args['wlinc']
	else:
		wlinc = 0.01
	if 'ngrid' in args:
		ngrid = args['ngrid']
	else:
		ngrid = 15
	if 'zgrid1' in args:
		zgrid1 = args['zgrid1']
	else:
		zgrid1 = 0.2
	if 'zgrid2' in args:
		zgrid2 = args['zgrid2']
	else:
		zgrid2 = 10
	if 'nstr' in args:
		nstr = args['nstr']
	else:
		nstr = 4
	if 'iday' in args:
		iday = args['iday']
	else:
		iday = 0
	if 'time' in args:
		time = args['time']
	else:
		time = 0
	if 'alat' in args:
		alat = args['alat']
	else:
		alat = 0
	if 'alon' in args:
		alon = args['alon']
	else:
		alon = 0
	if 'csza' in args:
		csza = args['csza']
	else:
		csza = -1
	if 'sza' in args:
		sza = args['sza']
	else:
		sza = 0
	if 'saza' in args:
		saza = args['saza']
	else:
		saza = 0
	if 'idatm' in args:
		idatm = args['idatm']
	else:
		idatm = 0
	if 'iaer' in args:
		iaer = args['iaer']
	else:
		iaer = -1
	if 'isalb' in args:
		isalb = args['isalb']
	else:
		isalb = -1
	if 'iout' in args:
		iout = args['iout']
	else:
		iout = 11
	
	with open(fn, 'w') as f:
		f.write('&INPUT\n\t')
		
		f.write('wlinf = '+str(wlinf)+'\n\t')
		f.write('wlsup = '+str(wlsup)+'\n\t')
		f.write('wlinc = '+str(wlinc)+'\n\t')
		
		f.write('ngrid = '+str(ngrid)+'\n\t')
		f.write('zgrid1 = '+str(zgrid1)+'\n\t')
		f.write('zgrid2 = '+str(zgrid2)+'\n\t')
		
		f.write('nstr = '+str(nstr)+'\n\t')
		
		f.write('iday = '+str(iday)+'\n\t')
		f.write('time = '+str(time)+'\n\t')
		f.write('alat = '+str(alat)+'\n\t')
		f.write('alon = '+str(alon)+'\n\t')
		
		f.write('csza = '+str(csza)+'\n\t')
		f.write('sza = '+str(sza)+'\n\t')
		f.write('saza = '+str(saza)+'\n\t')
		
		f.write('idatm = '+str(idatm)+'\n\t')
		f.write('iaer = '+str(iaer)+'\n\t')
		f.write('isalb = '+str(isalb)+'\n\t')
		f.write('iout = '+str(iout)+'\n')
		
		f.write('/\n')

def change_back(**args):
	'''
	This function is to change INPUT back
	input:
		**fn  : file path for INPUT, default ./INPUT
	output:
		albedo.dat for SBDART
	'''
	if 'fn' in args:
		fn = args['fn']
	else:
		fn = 'INPUT'
	
	if os.path.exists(fn+'_origin'):
		if os.path.exists(fn):
			os.system('rm '+fn)
		os.system('mv '+fn+'_origin '+fn)
	else:
		print('No origin file INPUT_origin. Please check.')
		sys.exit()

if __name__ == '__main__':
	write()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
