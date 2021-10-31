#! /usr/bin/env python
# -*-coding: utf-8 -*-
import sys
import math
import numpy as np
from sklearn import metrics
from sklearn.metrics import confusion_matrix
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def opendata(file,dataset,softname):
	true_list = []
	pred_list = []
	ii = 0
	soft_fu = ['SIFT_score','PROVEAN_score','LRT_score','FATHMM_score','SIFT_RAW','SIFT']

	soft_id = get_soft_id(softname)
	for line in open(file,'r'):
		line = line.strip('\n')
		if "pos" in line:
			continue
		tmp = line.split('\t')

		if tmp[soft_id] == "" or tmp[soft_id] == "." :
			continue
		else:
			pb = tmp[-2]

		if tmp[soft_id] == "" or tmp[soft_id] == ".":
			continue
		elif softname == "SIFT":
			s_pred = 0-float(tmp[soft_id])
		else:
			s_pred = float(tmp[soft_id])

		if dataset in tmp[-3]:
			true_list.append(pb)
			pred_list.append(s_pred)
			ii +=1
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

(file,dataset,i)=(sys.argv[1],sys.argv[2],sys.argv[3])
auc=[]
name={}

for softname in ["CADD","ClinPred","MCAP","MetaLR","Polyphen2","PrimateAI","REVEL","SIFT","GBDT"]:
	true_list,pred_list = opendata(file,dataset,softname)
	fpr, tpr, thresholds = metrics.roc_curve(true_list,pred_list,pos_label='Pathogenic')  # M-hM-?M-^YM-iM-^GM-^LM-gM-^ZM-^DM-cM-^@M-^@pos_labelM-cM-^@M-^@M-fM-^LM-^GM-gM-^ZM-^DM-fM-^XM-/M-oM-<M-^LM-fM-^HM-^QM-dM-;M-,M-fM-^JM-^J0M-dM-=M-^\M-dM-8M-:M-fM--M-#M-gM-1M-;M-oM-<M-^LM-dM-9M-^_M-eM-0M-1M-fM-^XM-/Pathogenic
	roc_auc = metrics.auc(fpr, tpr)
	name[roc_auc]=softname
	auc.append(roc_auc)


for roc_auc in sorted(list(set(auc)),reverse=True): 
	softname = name[roc_auc]
	print (softname)
	true_list,pred_list = opendata(file,dataset,softname)
	fpr, tpr, thresholds = metrics.roc_curve(true_list,pred_list,pos_label='Pathogenic')
	plt.plot(fpr, tpr, lw=1,label='%s=%0.4f' % (softname, roc_auc))

plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
#plt.title('Receiver operating characteristic')
plt.legend(loc="lower right")
plt.savefig(sys.argv[2] + sys.argv[3] + '.png')
plt.show()
