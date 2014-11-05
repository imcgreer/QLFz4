#!/usr/bin/env python

import numpy as np

def sim_DRW_lightcurve(t,SFinf,tau,mean_mag):
	'''Simulate a DRW lightcurve for a given time series t, with parameters
	   (SFinf, tau), and mean magnitude.
	   Uses equations A4 and A5 in KBS09 (see also MacLeod+10 sec 2.2).'''
	mu = mean_mag
	mag = np.zeros(len(t),dtype=np.float32)
	mag[0] = mean_mag
	dt = np.diff(t)
	for i in range(1,len(t)):
		loc = np.exp(-dt[i-1]/tau)*mag[i-1] + mu*(1-np.exp(-dt[i-1]/tau))
		var = 0.5 * SFinf**2 * (1-np.exp(-2*dt[i-1]/tau))
		mag[i] = np.random.normal(loc=loc,scale=np.sqrt(var))
	return mag

