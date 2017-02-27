#!/usr/bin/python
# -*- coding: utf-8 -*-

import os,io,types,math,glob
from astropy.io import fits
import numpy as np
from snrmy import *
from AstroImgNoise import *

rec = []
path='/home/awen/Awen/training_set/1000/*.fit*'
if __name__== '__main__':
    rec = glob.glob(path)
    rec.sort() #按读入顺序排列
length=len(rec)
list = []
for k in range(0,1):
	hh=fits.open(rec[k]) 
	f=hh[0].data
	print rec[k]
	Source=[10,2.0e11,f,15,3]
	CCD=[1,100,0.99999,100,51,1,10000,0.99*2]
	Exposuretime=0.0001
	OutImg=AstroImgNoise1(Source,CCD,Exposuretime)
	snr=snr(OutImg)
	print "snr=",snr
	print "SNR=",SNR
	#a=str(k)
	#show = fits.PrimaryHDU(OutImg)
	#showlist = fits.HDUList([show])
	#showlist.writeto('/home/awen/Awen/training_set/badsnr'+a+'.fits')     	
	

