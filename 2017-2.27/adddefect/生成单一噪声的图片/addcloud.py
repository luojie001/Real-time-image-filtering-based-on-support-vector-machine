#!/usr/bin/python
# -*- coding: utf-8 -*-
from astropy.io import fits
import numpy as np 
import glob,os,math
#####################
def randn(M,N):
	c_real=np.random.normal(0,1,(M,N))
	c_imag=np.random.normal(0,1,(M,N))
	cc=np.complex_(c_real)
	cc.imag=c_imag
	return cc
##################################################################################
def ft_phase_screen(r0,N,delta,L0,l0):
	del_f = 1.0/(N*delta)
	fx = np.array(range(-N/2,N/2)) * del_f
	fx,fy = np.meshgrid(fx,fx)
	th=np.arctan2(fy,fx)
	f = np.hypot(fx,fy)
	fm = 5.92/l0/(2*math.pi)
	f0 = 1.0/L0
	PSD_phi = 0.023*(r0**(-5.0/3))*(f**(-11.0/3))
	PSD_phi[N/2,N/2] = 0
	cn =math.sqrt(2)*randn(N,N)*np.sqrt(PSD_phi)*del_f
	phz = ift2(cn,1)
	return phz
##################################################################################
def ift2(G,delta_l):
	Gpd=np.fft.ifftshift(G)
	Goverpad=np.fft.ifft2(Gpd)
	g=np.fft.ifftshift(Goverpad)*(N*delta_l)**2
	return g
############################  adcloud ###################################
def adcloud(f):
	phz_hi = ft_phase_screen(r0, N, delta, L0, l0)
	phz=phz_hi
	f=gray*0
	for i in range(1,N):
		for j in range(1,N):
			if phz[i][j]>1:
				f[i,j]=20*phz[i,j]+300
			else:
				f[i,j]=gray[i,j]
	return f
##########################################################	
PN=1.0
D = 1.0
r0 = 0.1
N = 1024
L0 = 100
l0 = 0.01
delta = D/N
path='/home/awen/Awen/training_set/1000/*.fit*'
if __name__== '__main__':
    rec = glob.glob(path)
    rec.sort() #按读入顺序排列
for k in range(0,1):
	hh=fits.open(rec[k]) 
	gray=hh[0].data
	phz_hi = ft_phase_screen(r0, N, delta, L0, l0)
	phz=phz_hi
	f=gray*0
	for i in range(1,N):
		for j in range(1,N):
			if phz[i][j]>1:
				f[i,j]=20*phz[i,j]+300
			else:
				f[i,j]=gray[i,j]
	show = fits.PrimaryHDU(f)
	showlist = fits.HDUList([show])
	showlist.writeto('/home/awen/Awen/training_set/yun1.fits')   
