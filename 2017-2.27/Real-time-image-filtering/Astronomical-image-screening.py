#!/usr/bin/python
# -*- coding: utf-8 -*-

import os,io,types,math,glob,shutil  #调用函数库
from astropy.io import fits
import numpy as np
import pycuda.autoinit
import pycuda.driver as drv
from pycuda.compiler import SourceModule
from svmutil import *
from svm import *
from wenli0 import *
from touyin import *
import time
 
#start = time()
#print("Start: " + str(start))

mod = SourceModule("""
__global__ void Fun_snr(float *S_td, float *f_picture,int *N,int *Modesize)
{
	const int Imgw = N[0,0];
	const int modesize = Modesize[0,0];
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	float s_td=0.0;
	float sum=0.0;   
	if(blockIdx.x <4060&& !((i+1)%Imgw/(Imgw-modesize+1)))
		{
		for(int p=0 ; p<modesize ; p++)
			for(int q=0 ; q<modesize ; q++)
			{
			sum  =sum + f_picture[ (blockIdx.x+Imgw/blockDim.x*p) * blockDim.x + threadIdx.x + q];	
			}
		sum = sum /pow(float(modesize),2);
		for(int p=0 ; p<modesize ; p++)
			for(int q=0 ; q<modesize  ; q++)
			{
			s_td =s_td + pow((f_picture[ (blockIdx.x+Imgw/blockDim.x*p)  * blockDim.x + threadIdx.x + q]-sum),2);
	 		}
		S_td[i]= sqrt(s_td / pow(float(modesize),2));
		}
	else
		{
			S_td[i]=0;
		}
}
__global__ void Fun_ROAD(float *ROAD, float *f_picture,int *N)
{	
	const int Imgw = N[0,0];
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if(blockIdx.x <4088&& !((i+1)%Imgw/(Imgw-2)))
		{
		float a[9];
		int a_i=1;
		for(int p=0 ; p<3 ; p++)
			for(int q=0 ; q<3; q++)
			{
			a[a_i]=f_picture[(blockIdx.x+Imgw/blockDim.x*p) * blockDim.x + threadIdx.x + q]-f_picture[(blockIdx.x+Imgw/blockDim.x*1) * blockDim.x + threadIdx.x + 1];
			a[a_i]=abs(a[a_i]);
			a_i+=1;
	 		}
	 		{
				int ii,j;
				float t;
				for(ii=0;ii<8;ii++)
					for(j=ii+1;j<9;j++)
					{
						if(a[ii]>a[j])
						{
						t=a[ii];
						a[ii]=a[j];
						a[j]=t;
						}
					}
			}
		float a_sum=0;
		for(int kk=0;kk<5;kk++)
			a_sum=a_sum+a[kk];
		ROAD[i]= a_sum;
		}
	else
		{
			ROAD[i]=0;
		}
		__syncthreads();
}
""")
Fun_snr= mod.get_function("Fun_snr")
Fun_ROAD= mod.get_function("Fun_ROAD")

#数据传入,参数设定
rec = []
#输入图像的路径
##########
dl=0 
ul=1000
##########
d=str(dl)
u=str(ul)
path = r'/home/awen/Awen/tdate/test'+d+'.'+u+'/*.fit*'
#产生的数据保存
#str1="/home/awen/Awen/printdate.txt"
#file1=open(str1,'w+')  

y, x = svm_read_problem('/home/awen/Test/svm实时图像处理/Real-time-image-filtering/test')
prob  = svm_problem(y, x)
#c与g的值有ezsy.py对train训练所得
param = svm_parameter('-t 0 -c 8 -g 0.125 -q') #线性核 cost参数(1)   是否估算正确概率(0)
model = svm_train(prob,param)   #训练好的SVM模型
ww=0
while ww<3:
#按读入顺序排列
	if __name__== '__main__':
	    rec = glob.glob(path)
	    rec.sort() 
	length=len(rec)
	if length>0:
		print "该",path,"下现在有",length,"副图片"
		for k in range(0,length):
			hh=fits.open(rec[k]) #read the data of fits
			print k,rec[k]
			f=hh[0].data
			f=f*256.0/np.max(f)
			Imgw=len(f) 
			Imgl=len(f[0])
			modesize=9  
			fw=f.astype(np.float)
			wenli0=WL0(fw)
			wlk=touyin(f)*10

			#把Imgw作为数组传入cuda
			N=np.zeros([1,1])
			N[0,0]=Imgw
			N= np.int32(N)  

			#把modesize传入cuda
			Modesize=np.zeros([1,1])
			Modesize[0,0]=modesize
			Modesize=np.int32(Modesize)  

			#将f一维化传入cuda
			f=f.astype(np.float32)
			f_picture=f.flatten()  #一维化
			S_td=np.zeros_like(f_picture).astype(np.float32)
			ROAD=np.zeros_like(f_picture).astype(np.float32)

			# GPU run
			nTheads = 256   
			nBlocks = Imgw*Imgl/nTheads
			Fun_snr(drv.Out(S_td),drv.In(f_picture),drv.In(N),drv.In(Modesize), block=( nTheads, 1, 1 ), grid=( nBlocks, 1 ,1) )
			Fun_ROAD(drv.Out(ROAD), drv.In(f_picture),drv.In(N), block=(nTheads, 1, 1 ), grid=( nBlocks, 1 ,1) )
			ROAD.sort()
			ROAD=ROAD*0.1
			S_td=S_td[S_td>0]
			snrmy=10*math.log(max(S_td)/min(S_td))
			#wenli='0 1:'+str(float("%0.1f"%snrmy))+' 2:'+str(float("%0.1f"%np.max(ROAD)))+' 3:'+str(float("%0.1f"%wlk))+' 4:'+str("%0.1f"wenli0[0]*0.1)+' 5:'+str(wenli0[1])+' 6:'+str(wenli0[2])+' 7:'+str(wenli0[3]*10)
			yt=[1]
			xt =[{1: float("%0.1f"%snrmy), 2: float("%0.1f"%np.max(ROAD)), 3: float("%0.1f"%wlk), 4: float("%0.1f"%(wenli0[0]*0.1)), 5: wenli0[1], 6: wenli0[2], 7: wenli0[3]}]
			print xt
			p_label, p_acc, p_val = svm_predict(yt,xt, model)
			
			if p_acc[0]<1:
				shutil.move('/home/awen/Awen/tdate/test'+d+'.'+u+'/'+os.path.basename(rec[k]),'/home/awen/Awen/tdate/test'+d+'.'+u+'bad/'+os.path.basename(rec[k]))
			else:
				shutil.move('/home/awen/Awen/tdate/test'+d+'.'+u+'/'+os.path.basename(rec[k]),'/home/awen/Awen/tdate/test'+d+'.'+u+'god/'+os.path.basename(rec[k]))
				
		#file1.write("%s\n"%wenli)
		ww=0
	else:
		time.sleep(10)
		ww+=1
#stop = time()
#print("Stop: " + str(stop))
#print(str(stop-start) + "秒")
