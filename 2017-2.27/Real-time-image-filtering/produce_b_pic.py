#!/usr/bin/python
# -*- coding: utf-8 -*-
#产生带5种带噪声的图片数据
import os,random,glob,math
from astropy.io import fits
import numpy as np
import scipy.stats as stats
from time import time

start = time()
##############cosmic_rays######################
def cosmic_rays(f,n):
	fmax=np.max(f)
	for n in range(n):
		dx=random.randint(0,len(f)-1)
		dy=random.randint(0,len(f[0])-1)
		f[dx,dy]=fmax
	return f
#################linear#######################
def linear(f,n):
	d=1
	fmax=np.max(f)
	for n in range(n):
		x=random.randint(0,len(f)-1)
		y=random.randint(0,len(f[0])-1)
		f[x:x+d,y:y+d]=fmax
		for j in range(50):
			dx=random.randint(0,1)
			dy=random.randint(-1,1)
			x=x+dx
			y=y+dy
			f[x:x+d,y:y+d]=fmax
	return f
########################  adcloud ###################
def randn(M,N):
	c_real=np.random.normal(0,1,(M,N))
	c_imag=np.random.normal(0,1,(M,N))
	cc=np.complex_(c_real)
	cc.imag=c_imag
	return cc
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
def ift2(G,delta_l):
	Gpd=np.fft.ifftshift(G)
	Goverpad=np.fft.ifft2(Gpd)
	g=np.fft.ifftshift(Goverpad)*(N*delta_l)**2
	return g
def adcloud(f,n):
	phz_hi = ft_phase_screen(r0, N, delta, L0, l0)
	phz=phz_hi.real
	g=f*0
	for i in range(1,N):
		for j in range(1,N):
			if phz[i][j]>-n:
				g[i,j]=20*phz[i,j]+300
			else:
				g[i,j]=f[i,j]
	return g

PN=1.0
D = 1.0
r0 = 0.1
N = 1024
L0 = 100
l0 = 0.01
delta = D/N
#################################AstroImgNoise##########################
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
#############################################
def adnoise(f):
	Source=[10,2.0e11,f,15,3]
	CCD=[1,100,0.99999,100,51,1,10000,0.99*2]
	Exposuretime=0.0001
	OutImg=AstroImgNoise1(Source,CCD,Exposuretime)
	return OutImg
rec = []
path=r'/home/awen/Awen/training_set/100/*.fit*'
if __name__== '__main__':
	rec = glob.glob(path)
	rec.sort() 
print "There are",len(rec),"pairs of images under this folder"
##########
dl=0 
ul=100
##########
d=str(dl)
u=str(ul)
for zz in range(5):
	if zz==0:
		for k in range(dl,ul):
			hh=fits.open(rec[k]) 
			f=hh[0].data 
			f=linear(f,3)
			f=adcloud(f,0)#[-10,10]
			f=adnoise(f)
			a=str(k)
			b=str(zz)
			show = fits.PrimaryHDU(f)
			showlist = fits.HDUList([show])
			showlist.writeto('/home/awen/Awen/training_set/100bad/'+b+'test'+a+'.fits') 
	elif zz==1:
		for k in range(dl,ul):
			hh=fits.open(rec[k]) #read the data of fits
			f=hh[0].data 
			f=cosmic_rays(f,5)
			f=adcloud(f,0)#[-10,10]
			f=adnoise(f)
			a=str(k)
			b=str(zz)
			show = fits.PrimaryHDU(f)
			showlist = fits.HDUList([show])
			showlist.writeto('/home/awen/Awen/training_set/100bad/'+b+'test'+a+'.fits') 
	elif zz==2:
		for k in range(dl,ul):
			hh=fits.open(rec[k]) #read the data of fits
			f=hh[0].data 
			f=cosmic_rays(f,5)
			f=linear(f,3)
			f=adnoise(f)
			a=str(k)
			b=str(zz)
			show = fits.PrimaryHDU(f)
			showlist = fits.HDUList([show])
			showlist.writeto('/home/awen/Awen/training_set/100bad/'+b+'test'+a+'.fits') 
	elif zz==3:
		for k in range(dl,ul):
			hh=fits.open(rec[k]) #read the data of fits
			f=hh[0].data 
			f=cosmic_rays(f,5)
			f=linear(f,3)
			f=adnoise(f)
			a=str(k)
			b=str(zz)
			show = fits.PrimaryHDU(f)
			showlist = fits.HDUList([show])
			showlist.writeto('/home/awen/Awen/training_set/100bad/'+b+'test'+a+'.fits') 
	else :    
		zz==4
		for k in range(dl,ul):
			hh=fits.open(rec[k]) #read the data of fits
			f=hh[0].data 
			f=cosmic_rays(f,5)
			f=linear(f,3)
			f=adcloud(f,0)#[-10,10]
			f=adnoise(f)
			a=str(k)
			b=str(zz)
			show = fits.PrimaryHDU(f)
			showlist = fits.HDUList([show])
			showlist.writeto('/home/awen/Awen/training_set/100bad/'+b+'test'+a+'.fits') 
stop = time()
print(str(stop-start) + "秒")
