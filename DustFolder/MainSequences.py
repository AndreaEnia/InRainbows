#!/usr/bin/python
# -*- coding: latin-1 -*-

import numpy as np

def RenziniPeng(Mstar):
	'''
	Main Sequence from Kennicutt+98
	I/O: Mstar (in Msol), SFR (in Msol yr^-1)
	'''
	return 10**(-7.64 + 0.76*Mstar)

def Kennicutt98(SigmaStar):
	'''
	Main Sequence from Kennicutt+98
	I/O: SigmaStar (in Msol), SigmaSFR (in Msol yr^-1)
	'''
	return 10**(-10.41 + 0.99*SigmaStar)

def Kennicutt98_4x(SigmaStar):
	'''
	Main Sequence from Kennicutt+98, 4 times elevated
	I/O: SigmaStar (in Msol), SigmaSFR-4x (in Msol yr^-1)
	'''
	return 10**(-10.41 + 0.99*SigmaStar + 0.6021)

def Kennicutt98_Nx(SigmaStar, N):
	'''
	Main Sequence from Kennicutt+98, N times elevated
	I/O: SigmaStar (in Msol), SigmaSFR-Nx (in Msol yr^-1)
	'''
	if N > 0: return 10**(-10.41 + 0.99*SigmaMStar + np.log10(N))
	elif N == 0: return 10**(-10.41 + 0.99*SigmaMStar)
	elif N < 0: return 10**(-10.41 + 0.99*SigmaMStar - np.log10(abs(N)))

def CD16(SigmaStar):
	'''
	Main Sequence from Cano-Diaz+16
	I/O: SigmaStar (in Msol), SigmaSFR (in Msol yr^-1)
	'''
	return 10**(-7.95 + 0.72*SigmaStar)

def CD16_4x(SigmaStar):
	'''
	Main Sequence from Cano-Diaz+16, 4 times elevated
	I/O: SigmaStar (in Msol), SigmaSFR-4x (in Msol yr^-1)
	'''
	return 10**(-7.95 + 0.72*SigmaStar + 0.6021)

def REM(SigmaMStar, m, q):
	'''
	Main Sequence from Rodighiero-Enia-Morselli
	I/O: SigmaStar (in Msol), SigmaSFR (in Msol yr^-1)
	'''
	return 10**(q + m*SigmaMStar)

def REM_4x(SigmaMStar, m, q):
	'''
	Main Sequence from Rodighiero-Enia-Morselli, 4 times elevated
	I/O: SigmaStar (in Msol), SigmaSFR-4x (in Msol yr^-1)
	'''
	return 10**(q + m*SigmaMStar + 0.6021)

def REM_Nx(SigmaMStar, m, q, N):
	'''
	Main Sequence from Rodighiero-Enia-Morselli, N times elevated
	I/O: SigmaStar (in Msol), SigmaSFR-Nx (in Msol yr^-1)
	'''
	if N > 0: return 10**(q + m*SigmaMStar + np.log10(N))
	elif N == 0: return 10**(q + m*SigmaMStar)
	elif N < 0: return 10**(q + m*SigmaMStar - np.log10(abs(N)))

def Custom_Nx(SigmaMStar, m, q, N):
	return 10**(q + m*SigmaMStar + np.log10(N))