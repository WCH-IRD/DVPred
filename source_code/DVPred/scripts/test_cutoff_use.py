#! /usr/bin/env python
# -*-coding: utf-8 -*-
import sys
import math
import numpy as np
from sklearn.metrics import confusion_matrix

#dataset=(train/test)

def opendata(file,dataset,softname,cutoff):
	true_list = []
	pred_list = []
	soft_fu = ['SIFT_score','PROVEAN_score','LRT_score','FATHMM_score','SIFT_RAW','SIFT']

	soft_id = get_soft_id(softname)
	for line in open(file,'r'):
		line = line.strip('\n')
		if "pos" in line:
			continue
		tmp = line.split('\t')

		if tmp[soft_id] == "" :
			continue
		elif tmp[-2] == 'Benign':
			pb = 'B'
		elif tmp[-2] == 'Pathogenic':
			pb = 'P'

		if tmp[soft_id] == "" :
			continue
		elif float(tmp[soft_id]) < cutoff and softname not in soft_fu:
			s_pred = 'B'
		elif float(tmp[soft_id]) >= cutoff and softname not in soft_fu:
			s_pred = 'P'
		elif float(tmp[soft_id]) > cutoff and softname in soft_fu:
			s_pred = 'B'
		elif float(tmp[soft_id]) <= cutoff and softname in soft_fu:
			s_pred = 'P'
		else:
			break
			print ("die\n")

		if dataset in tmp[-3]:
			true_list.append(pb)
			pred_list.append(s_pred)
	return true_list,pred_list

def get_soft_id(name):
	head = ["#chr","pos",".","ref","alt","CADD","ClinPred","MCAP","MetaLR","Polyphen2","PrimateAI","REVEL","SIFT","data","P_B","GBDT"]
	i = 0
	hash={}
	for id in head:
#		print(id)
		hash[id] = i
		i+=1
	return hash[name]

def get_softnum_list(file,dataset,softname):
	soft_id = get_soft_id(softname)
	softnum_list=[]
	for line in open(file,'r'):
		line = line.strip('\n')
		if 'pos' in line:
			continue
		tmp = line.split('\t')
		if tmp[soft_id] == "" :
			next
		else :
			softnum_list.append(float(tmp[soft_id]))
	softnum_list=sorted(list(set(softnum_list)))
	return softnum_list

def get_result(tn,fp,fn,tp):
	if tn==0:
		tn=1
	if fp==0:
		fp=1
	if fn==0:
		fn=1
	if tp==0:
		tp=1
	(tn,fp,fn,tp)=(float(tn),float(fp),float(fn),float(tp))
	tpr=tp/(tp+fn)
	tnr=tn/(fp+tn)
	ppv=tp/(fp+tp)
	npv=tn/(tn+fn)
	gm=(tpr*tnr) ** 0.5
	mk=ppv+npv-1
	bcr=(tpr+tnr)/2
	mcc_fenmu=((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))** 0.5
	mcc_fenzi=(tp*tn)-(fp*fn)
	mcc=mcc_fenzi/mcc_fenmu
	dp0=tpr/(1-tnr)
	dp1=math.log(dp0,10)
	dp_0=tnr/(1-tpr)
	dp2=math.log(dp_0,10)
	dp_chang= (3**0.5)/math.pi
	dp_p=dp_chang*(dp1+dp2)
	f1=(2*tp)/(2*tp+fp+fn)
	acc=(tp+tn)/(tp+tn+fp+fn)
#	print (dp0,dp1,dp_0,dp2,dp_p)
	return(tn,fp,fn,tp,tpr,tnr,ppv,npv,f1,acc,gm,mk,bcr,mcc,dp_p)


(file,dataset,soft,biaodi,cutoff)=(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],float(sys.argv[5]))

soft_id = get_soft_id(soft)
softnum_list = get_softnum_list(file,dataset,soft)


result = open('testing_result.xls','a')

(old_true_list,old_pred_list) = opendata(file,dataset,soft,cutoff)
tn, fp, fn, tp = confusion_matrix(old_true_list,old_pred_list).ravel()
(tn,fp,fn,tp,tpr,tnr,ppv,npv,f1,acc,gm,mk,bcr,mcc,dp_p)=get_result(tn, fp, fn, tp)
print ('%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(soft,biaodi,cutoff,tn,fp,fn,tp,tpr,tnr,ppv,npv,f1,acc,gm,mk,bcr,mcc,dp_p))
result.write('%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(soft,biaodi,cutoff,tn,fp,fn,tp,tpr,tnr,ppv,npv,f1,acc,gm,mk,bcr,mcc,dp_p))
