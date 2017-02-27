#!/usr/bin/python
# -*- coding: utf-8 -*-
from svmutil import *
from svm import *

y, x = svm_read_problem('/home/awen/Awen/training_set/train.txt')
yt, xt = svm_read_problem('/home/awen/Awen/training_set/test.txt')

prob  = svm_problem(y, x)
param = svm_parameter('-t 0 -c 2 -g 0.03125 -q') #线性核 cost参数(1)   是否估算正确概率(0)
model = svm_train(prob,param)   #训练好的SVM模型
print('test:')
p_label, p_acc, p_val = svm_predict(yt, xt, model)
print p_acc
