#!/usr/bin/python
# -*- coding: utf-8 -*-

import os,io,types,math,glob  #调用函数库
from astropy.io import fits
import numpy as np

def touyin(f):
	Imgw=len(f)
	Imgl=len(f[0])
	fmax=np.max(f)
	Imgl=len(f[0])
	#投影到x,y轴
	xline=[]
	yline=[]
	for j in range(0,Imgl):
		xmax=np.max(f[:,j])
		ymax=np.max(f[j,:])
		if xmax>=0.9*fmax:
			xline.append(j)
		if ymax>=0.9*fmax:
			yline.append(j)	
	#旋转45°
	Imgw2=2*Imgw-1
	cc=Imgw-1
	g=np.zeros([Imgw2,Imgw2])
	for i in range(Imgw):
		for j in range(Imgl):
			g[i+j][i+j+cc]=f[i][j]
		cc=cc-2
	#投影到45°,135°
	x2line=[]
	y2line=[]
	for j in range(0,Imgw2):
		x2max=np.max(g[:,j])
		y2max=np.max(g[j,:])
		if x2max>=0.9*fmax:
			x2line.append(j)
		if y2max>=0.9*fmax:
			y2line.append(j)
	a=len(xline)
	b=len(yline)
	c=len(x2line)
	d=len(y2line)
	#print a,b,c,d	
	panduan=[a,b,c,d]
	wlk=0.0
	if np.min(panduan)==1:
		wlk=(np.max(panduan)-1)/(np.min(panduan)*1.0)
	else:
		wlk=np.max(panduan)/(np.min(panduan)*1.0)
	return float(wlk)
