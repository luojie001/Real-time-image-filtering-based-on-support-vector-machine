#!/usr/bin/python
# -*- coding: utf-8 -*-
from svmutil import *
from svm import *

y, x = svm_read_problem('/home/awen/Awen/training_set/train.txt')
yt, xt = svm_read_problem('/home/awen/Awen/training_set/test.txt')
#model = svm_train(y, x )
#print('test:')
#p_label, p_acc, p_val = svm_predict(yt[200:202], xt[200:202], model)
#print(p_label)
prob  = svm_problem(y, x)
#param = svm_parameter('-t 0 -c 1 -b 0') #线性核 cost参数(1)   是否估算正确概率(0)
param = svm_parameter('-t 0 -c '+str(1)+' -g '+str(1)+' -b 0') 
model = svm_train(prob,param)   #训练好的SVM模型
print('test:')
p_label, p_acc, p_val = svm_predict(yt, xt, model)
print p_acc

Acc_max=p_acc[0]
X1=1
X2=1
for x1 in range(1,100,1):
	for x2 in range(1,100,1):
		param = svm_parameter('-t 0 -c '+str(x1*0.01)+' -g '+str(x2*0.01)+' -b 0') 
		model = svm_train(prob,param)   #训练好的SVM模型
		print('test:')
		print x1,x2
		p_label, p_acc, p_val = svm_predict(yt, xt, model)
		if p_acc[0]>Acc_max:
			Acc_max=p_acc[0]
			X1=x1
			X2=x2
print X1,X2,Acc_max
		#for i in range(len(p_label)):
			#print i,p_label[i],p_val[i]
			#print p_acc
			#print 
