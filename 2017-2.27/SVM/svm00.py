#!/usr/bin/python
# -*- coding: utf-8 -*-
from svmutil import *
from svm import *

y, x = svm_read_problem('/home/awen/Awen/training_set/test.txt')
yt, xt = svm_read_problem('/home/awen/Awen/training_set/train.txt')
prob  = svm_problem(y, x)
xc=12
xg=1
param = svm_parameter('-t 0 -c '+str(xc*0.01)+' -g '+str(xg*0.01)+' -p 0.1 -b 0') #线性核 cost参数(1)   是否估算正确概率(0)
model = svm_train(prob,param)   #训练好的SVM模型
print('test:')
p_label, p_acc, p_val = svm_predict(yt, xt, model)
"""
Acc_max=p_acc[0]
XP=0
for xp in range(1,100,1):
	param = svm_parameter('-t 0 -c 0.12 -g 0.01 -p '+str(xp*0.01)+' -b 0') 
	model = svm_train(prob,param)   #训练好的SVM模型
	print('test:')
	p_label, p_acc, p_val = svm_predict(yt, xt, model)
	if p_acc[0]>Acc_max:
		Acc_max=p_acc[0]
		XP=xp
print XP,Acc_max
"""
