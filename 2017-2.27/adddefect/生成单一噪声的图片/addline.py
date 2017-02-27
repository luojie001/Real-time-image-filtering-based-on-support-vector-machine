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
for i in range(350,400):
	hh=fits.open(rec[i]) #read the data of fits
	f=hh[0].data 
	fmax=np.max(f)
	Imgw=len(f)
	Imgl=len(f[0])
	d=1
	for n in range(3):
		x=random.randint(0,Imgw)
		y=random.randint(0,Imgl)
		f[x:x+d,y:y+d]=fmax
		for j in range(50):
			dx=random.randint(0,1)
			dy=random.randint(-1,1)
			x=x+dx
			y=y+dy
			f[x:x+d,y:y+d]=fmax
	a=str(i+1)
	show = fits.PrimaryHDU(f)
	showlist = fits.HDUList([show])
	showlist.writeto('/home/awen/Awen/training_set/line/line'+a+'.fits')     
