#!/usr/bin/python
# -*- coding: utf-8 -*-

import os,random,glob
from astropy.io import fits
import numpy as np

rec = []
path=r'/home/awen/Awen/training_set/100/*'
if __name__== '__main__':
	rec = glob.glob(path)
	rec.sort() #按读入顺序排列
length=len(rec)
print length
list = []
for i in range(0,200):
	hh=fits.open(rec[i]) #read the data of fits
	f=hh[0].data 
	n=np.random.normal(100,60,(1024,1024))
	f=f+n
	a=str(i+1)
	show = fits.PrimaryHDU(f)
	showlist = fits.HDUList([show])
	showlist.writeto('/home/awen/Awen/training_set/noice/n'+a+'.fits')     
