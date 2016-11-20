#!/usr/bin/python 
# -*- coding: utf-8 -*-

import os,random,math
from astropy.io import fits
import numpy as np

def WL0(f):
	w=len(f)
	l=len(f[0])
	fmin=np.min(f)
	ones=f*0+65537-65536 #将数据类型由f的uint16转化为uint32
	g=f-ones*np.min(f)
	f1=g*(w-1)/np.max(g)
	f1=f1.astype(np.int16)
	#print f1,np.max(f1),np.min(f1)
	t=f*0
	#print f1,type(f1)
	for i in range(w):
		for j in range(l-1):
			t[f1[i,j+1],f1[i,j]] +=1
			t[f1[i,j],f1[i,j+1]] +=1
	#print t
	ASM=0.0
	CON=0.0
	IDF=0.0
	Hxy=0.0
	s=sum(sum(t))+0.0
	#print "S=",s
	for x in range(w):
		for y in range(l):
			ASM=ASM+pow(t[x,y],2)/pow(s,2)
			IDF=IDF+t[x,y]/s/(1+pow((x-y),2))
			if t[x,y]!=0:
				Hxy=Hxy-t[x,y]/s*(np.log(t[x,y]/s)/np.log(2))
			if 0<=y+x & y+x<=w-1 :
				CON=CON+pow(x,2)*t[y,y+x]
			if 0<=y-x & y-x<=w-1:
				CON=CON+pow(x,2)*t[y,y-x]
	ASM=ASM*10000	
	CON=np.log(CON)
	IDF=IDF*100
	#print "ASM=","%0.3f"%ASM,"CON=","%0.1f"%CON,"IDF=",IDF,"Hxy=","%0.1f"%Hxy
	return float("%0.1f"%ASM),float("%0.1f"%CON),float("%0.1f"%IDF),float("%0.1f"%Hxy)
