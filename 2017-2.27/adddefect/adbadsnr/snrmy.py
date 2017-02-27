#!/usr/bin/python
# -*- coding: utf-8 -*-
import math  #调用函数库
import numpy as np

def snr(f):
	stdmax=0            #stdmax 初始化
	stdmin=1000          #stdmin初始化
	Imgw=len(f)          #获得数组的w，l
	Imgl=len(f[0])
	modesizew=9         #给定模板的大小9
	modesizel=9
	for i in range(0,Imgw-modesizew):       #对每一副图算snr
		for j in range(0,Imgl-modesizel):
			a1=np.zeros((modesizew,modesizel),dtype=int)
			a1=f[i:i+modesizew,j:j+modesizel]
			std=np.std(a1)
			if std > stdmax:
				stdmax=std
			if std < stdmin:
				stdmin=std
	snrmy=10*math.log(stdmax/stdmin)
	return snrmy
