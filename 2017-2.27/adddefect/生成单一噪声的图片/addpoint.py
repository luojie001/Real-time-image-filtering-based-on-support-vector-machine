#!/usr/bin/python
# -*- coding: utf-8 -*-

import os,random,glob
from astropy.io import fits
import numpy as np

rec = []
path=r'/home/awen/Awen/1000/*'
if __name__== '__main__':
	rec = glob.glob(path)
	rec.sort() #按读入顺序排列
length=len(rec)
print length
list = []
for i in range(0,200):
	hh=fits.open(rec[i]) #read the data of fits
	f=hh[0].data 
	fmax=np.max(f)
	Imgw=len(f)
	Imgl=len(f[0])
	for n in range(5):
		dx=random.randint(0,1023)
		dy=random.randint(0,1023)
		f[dx,dy]=fmax
		print i,dx,dy
	a=str(i+1)
	show = fits.PrimaryHDU(f)
	showlist = fits.HDUList([show])
	showlist.writeto('/home/awen/Awen/training_set/cosmic_rays/point'+a+'.fits')     
