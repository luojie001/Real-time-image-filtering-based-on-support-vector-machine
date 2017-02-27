#!/usr/bin/python
# -*- coding: utf-8 -*-

import os,io,types,math
from astropy.io import fits
import numpy as np
import scipy.stats as stats

def AstroImgNoise1(Source,CCD,Exposuretime):
	#源图像参数传入
	Mag=Source[0]
	Zeropoint=Source[1]
	Img=Source[2]
	SkyB=Source[3]
	TeleA=Source[4]
	#CCD参数传入
	Type=CCD[0]
	RONs=CCD[1]
	CTE=CCD[2]
	DC=CCD[3]
	Gain=CCD[4]
	Bining=CCD[5]
	FWell=CCD[6]
	QE=CCD[7]
	Imgw=len(Img)
	Imgl=len(Img[0])
	if Type==1:
		Gain=Gain*np.ones([Imgw,Imgl])
	elif Type==2:
		RandAmp=0.01
		Gain =Gain*(np.ones([Imgw,Imgl])+RandAmp*np.random.random([Imgw,Imgl]))
	else:
		print "Wrong Type of CCD"
	#III 光度学计算 图像有效灰度值转化为ADU
	Nphoton=Zeropoint*10**(-0.4*Mag)*Exposuretime*np.pi*TeleA**2
	#背景光子数计算
	Bphoton=Zeropoint*10**(-0.4*SkyB)*Exposuretime*np.pi*TeleA**2
	#产生像素分布的概率密度函数pdf
	Imgpdf=np.round(Img*Nphoton/(np.sum(np.sum(Img))))
	#产生入射图像光电子矩阵
	Img=Imgpdf*QE+stats.poisson.rvs(mu=Bphoton/Imgw/Imgl,size=([Imgw,Imgl]))
	#IV 产生读出噪声矩阵RON (包含Bias和电路噪声,Possion分布)
	RONMatrix=stats.poisson.rvs(mu=RONs,size=([Imgw,Imgl]));
	#产生暗电流噪声矩阵 
	DCMatrix=date=stats.poisson.rvs(mu=DC*Exposuretime,size=([Imgw,Imgl]))
	#根据电荷转移率产生电荷损失矩阵
	CTI=(Imgw+Imgl)*(1-CTE)*(1+0.01*np.random.random([Imgw,Imgl]))*Img
	#根据增益和电荷矩阵产生读出电子图像矩阵
	Img=(Img-CTI)*Gain+RONMatrix+DCMatrix
	#根据CCD井深产生最终图像
	Img[Img>FWell]=FWell
	#IIV 根据是否bining确定输出矩阵大小
	Outimg=np.zeros([Imgw/Bining,Imgl/Bining])
	for i in range(0,Imgw/Bining):
		for j in range(0,Imgl/Bining):
			for k in range(0,Bining):
		    		Outimg[i,j]=Outimg[i,j]+Img[(i-1)*Bining+k,(j-1)*Bining+k]
	
	SignalNoiseRatio=20*np.log(Nphoton/(Imgw*Imgl*(RONs**2+DC)+Bphoton+Nphoton)**(0.5))
	#IIIV 最终输出图像
	return Outimg
