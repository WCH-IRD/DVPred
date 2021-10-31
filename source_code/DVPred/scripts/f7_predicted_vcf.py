#-*- coding: utf-8 -*-
from sklearn.model_selection import KFold
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.preprocessing import minmax_scale
from sklearn import metrics
from scipy import interp  
from sklearn.feature_selection import SelectKBest
# from sklearn.feature_selection import SelectPercentile
# from sklearn.feature_selection import f_classif
from sklearn.feature_selection import f_regression
from sklearn import decomposition
import numpy as np
import random
import sys
import os
from sklearn.externals import joblib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import getopt
import warnings
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import accuracy_score
warnings.filterwarnings("ignore")


def load_data(all_file,label_file,train_or_test):

	if train_or_test =='train':
		shuffle_file='/'.join(all_file.split('/')[0:-1]) + "./ML/tmp/shuffle.data"
		print shuffle_file
		shuffle_f=open(shuffle_file,'w')
		lines = open(all_file).readlines()
		random.shuffle(lines)
		for item in lines:
			shuffle_f.write(item)
		shuffle_f.close()
		p=os.popen("awk -F '\t' 'END{ print NF}' " + shuffle_file).readlines()
		NF = int(p[0].strip('\n'))
		os.system('cat '+shuffle_file+' |cut -f 1-7 > prefix.txt')
		os.system('cat '+shuffle_file+' |cut -f 8-'+str(NF-1)+' > feature.txt')
		os.system('cat '+label_file+' |awk \'{if($NF=="Pathogenic"){print 0}else{print 1}}\'> label.txt')
	elif train_or_test =='test':
		shuffle_file = all_file
		p=os.popen("awk -F '\t' 'END{ print NF}' " + shuffle_file).readlines()
		NF = int(p[0].strip('\n'))
		os.system('cat '+shuffle_file+' |cut -f 1-7 > prefix.txt')

        	os.system('cat  ' + shuffle_file + '|cut -f 8-' + str(NF-2) + ' > feature_all_in_scaled_tab.txt')


		if os.path.exists(label_file):
			os.system('cat '+label_file+' |awk \'{if($1=="Pathogenic"){print 0}else{print 1}}\'> label.txt')
		else:
			os.system('awk \'{print "1"}\' '+all_file+' >  label.txt')
        	NF_1 = int(NF-1)
	
	prefix = np.loadtxt('prefix.txt',dtype=np.string_) 
	feature = np.loadtxt('feature_all_in_scaled_tab.txt',dtype=np.float32)  
	label = np.loadtxt('label.txt',dtype=np.int8)  
	return prefix,feature,label,NF_1

def PCA():
    fea = np.loadtxt('feature1.txt', dtype=np.float32)
    fea_scaled = fea
    pca = joblib.load("./models/PCA.m")
    fea_1 = pca.transform(fea_scaled)
    np.savetxt('feature1_scaled.txt', fea_scaled)
    np.savetxt('feature1_1.txt', fea_1)

def LR_train(X_train_scaled,y_train,model_name):
	clf = LogisticRegression(penalty='l2',tol=0.01, C=1.0, solver='liblinear',\
			max_iter=100, multi_class='ovr') 
	clf.fit(X_train_scaled,y_train)
	joblib.dump(clf, model_name)

def SVM_train(c,g,X_train_scaled,y_train,model_name):
	clf = svm.SVC(C=c, gamma=g,probability=True,kernel='rbf')
	clf.fit(X_train_scaled,y_train)
	joblib.dump(clf, model_name)
	
def RandomForest_train(train_x,train_y,n,model_name):
	clf = RandomForestClassifier(n_estimators=n)
	clf.fit(train_x, train_y)
	joblib.dump(clf, model_name)


def GBDT_train(train_x,train_y,n,model_name):
	clf = GradientBoostingClassifier(n_estimators=n)
	clf.fit(train_x, train_y)
	joblib.dump(clf, model_name)

def NN_train(train_x, train_y, layers, model_name):
    clf = MLPClassifier(hidden_layer_sizes=(layers,), activation='relu',
                        shuffle=True, early_stopping=True)
    clf.fit(train_x, train_y)
    joblib.dump(clf, model_name)



# 为特征打分
def selectKfeatures(X_train_scaled,y_train,NF,expand_feature_name):
	NF = NF-8
	select =  SelectKBest(score_func=f_regression(X_train_scaled,y_train), k=NF)
	fea_wei={}
	li_score = select.get_params()['score_func'][0]
	for index,item in enumerate(li_score):
		item=float(item)
		fea_wei[index] = abs(item)
	# print fea_wei
	dict_score= sorted(fea_wei.items(), key=lambda d:d[1], reverse = True)
	a = [ tu[0] for tu in dict_score]
	b = [ tu[1] for tu in dict_score]
	name = expand_feature_name
	print "------------------------------------------------------"
	print 'The scores of all features are listed as the following:'
	for item in dict_score:
		print name[int(item[0])] + '\t' + str(item[1])
	print "------------------------------------------------------"

# 统计
def statistic(f1,f2) :
	os.system('cat '+ f1 +' |sort -u > tmp1.txt')
	single ={line.strip('\n'):0 for line in open("tmp1.txt")}
	staf=open('tmp2.txt','w')
	for line in open(f1):
		line=line.strip('\n')
		single[line] = single[line]+1
	for  (k,v) in  single.items():
		staf.write(str(v)+'\t'+k+"\n")
	staf.close()
	os.system('cat tmp2.txt | sort -k 1 -n -r > '+f2)
	os.system('rm tmp1.txt tmp2.txt')

def calculate_result(actual,pred):  

	target_names = ['Pathogenic', 'Benign']
	report = metrics.classification_report(actual, pred, target_names=target_names,digits=3)
	return report
	#print(metrics.classification_report(actual, pred, target_names=target_names,digits=3))

# 训练过程 见scripts/train_ml.py 脚本
# def ML_train(train_file,model_type,model_name,selectKfea,fea_name_file):
def ML_train(filled_data,label_file,model_type,model_path,Flag_selectKfeatures,expand_feature_name):
	os.system('paste '+ filled_data+' '+label_file+' > ./ML/tmp/all.data')
	prefix,feature,label,NF = load_data('./ML/tmp/all.data',"train")
	print 'Training '+ train_file +' -----------------------'
	X_train_scaled = minmax_scale(feature)
	if Flag_selectKfeatures:
		selectKfeatures(X_train_scaled,label,NF,expand_feature_name)
	model_name = model_type+'.m'
	if model_type=='LR':
		LR_train(X_train_scaled,label,model_name)
	elif model_type=='SVM':
		SVM_train(32,0.1,X_train_scaled,label,model_name)
	elif model_type=='RF':
		RandomForest_train(X_train_scaled,label,200,model_name)
	elif model_type=='GBDT':
		GBDT_train(X_train_scaled,label,200,model_name)
	elif model_type=='NN':
		NN_train(X_train_scaled,label,200,model_name)


def ROC(y_test,y_pred_pro,model_type,out_path):
	fpr, tpr, thresholds = metrics.roc_curve(y_test,y_pred_pro[:, 1]) 
	roc_auc = metrics.auc(fpr, tpr)
	print ("AUC = " + str(roc_auc))
	plt.plot(fpr, tpr, lw=1, label='model:%s (area = %0.4f)' % (model_type,roc_auc))  
	plt.xlim([-0.05, 1.05])  
	plt.ylim([-0.05, 1.05])  
	plt.xlabel('False Positive Rate')  
	plt.ylabel('True Positive Rate')  
	plt.title('Receiver operating characteristic')  
	plt.legend(loc="lower right")  
	plt.savefig(out_path+model_type+'.jpg')
	# plt.show()  

def ML_test(test_file,label_file,model_type,model_path,true_y_test,record_FP,draw_ROC,out_path):

	if true_y_test and record_FP:
		f1=open(out_path + model_type + "test_judge_p_as_b.txt",'a')
		f2=open(out_path + model_type + "test_judge_b_as_p.txt",'a')
	fr=open("./result.txt",'a')
	prof = open(out_path+model_type+'_pro.txt','w')
	prof_all = open(out_path+model_type+'_pro_all.txt','w')
	prefix,X_test,y_test,NF = load_data(test_file,label_file,"test")
	#X_test_scaled = minmax_scale(X_test)
	X_test_scaled = X_test
	clf = joblib.load(model_path)
	y_pred = clf.predict(X_test_scaled)
	y_pred_pro = clf.predict_proba(X_test_scaled) 
	for i in y_pred_pro:
		prof_all.write(str(i) + '\n')
		prof.write(str(i[0]) + '\n')
	prof.close()
	prof_all.close()
	if true_y_test:
		PR_precision, PR_recall, PR_thresholds = precision_recall_curve(y_test,y_pred)
                PR_auc = round(metrics.auc(PR_recall, PR_precision),4)
		print PR_auc
		accuracy = round(accuracy_score(y_test, y_pred),4)
		print accuracy
		print model_type,
		#fr.write("Name--------------  "+model_type+"  -----------------Total")
		fr.write(model_type,)
		report = calculate_result(y_test, y_pred)
		tmp = report.split('\n')
		t2 = tmp[2].split()
		fr.write('\t' + t2[0] + '\t' + t2[1] + '\t' + t2[2] + '\t' + t2[3],)
		report = calculate_result(y_test,y_pred)
		print report
		t3 = tmp[3].split()
		if int(t3[4]) == 0:
			fr.write("\t" + t2[4] + "\n")
		else:
			fpr, tpr, thresholds = metrics.roc_curve(y_test,y_pred_pro[:,0],pos_label=0)
			roc_auc = metrics.auc(fpr, tpr)
			print roc_auc
			roc_auc = round(roc_auc,4)
			fr.write("\t" + str(roc_auc) + '\t' + str(PR_auc) + '\t' + str(accuracy) + '\n')
			p_b_num = t2[4] + '/' + t3[4]
			fr.write('\t' + t3[0] + '\t' + t3[1] + '\t' + t3[2] + '\t' + t3[3] + '\t' + p_b_num + '\n')	
	
	if true_y_test and record_FP:
		for i in range(len(y_test)):
			if y_pred[i]==1 and y_test[i]==0:
				f1.write('\t'.join(list(prefix[i])) + '\n')
			elif y_pred[i]==0 and y_test[i]==1:
				f2.write('\t'.join(list(prefix[i])) + '\n')
		f1.close()
		f2.close()
		statistic(out_path + model_type + "test_judge_p_as_b.txt", out_path + model_type + 'test_p_b.txt')
		statistic(out_path + model_type + "test_judge_b_as_p.txt", out_path + model_type + 'test_b_p.txt')
	if true_y_test and draw_ROC:
		ROC(y_test,y_pred_pro,model_type,out_path)

	return out_path+model_type+'_pro.txt'


def predict_pathogenicity(filled_data,label_file,model_type,model_path,flag_label,flag_FN_FP,flag_fea_score,flag_roc,expand_feature_name,out_path):  
	# ML_train(sys.argv[1],'SVM')
	true_y_test = int(flag_label)
	record_FP = int(flag_FN_FP)
	draw_ROC = int(flag_roc)
	selectKfea = int(flag_fea_score)
	fea_name_file = expand_feature_name
	
	res=[]
	for i in range(len(model_type)):
		# 训练
		# ML_train(filled_data,label_file,model_type[i],model_path[i],Flag_selectKfeatures,expand_feature_name)
		res_pro_file = ML_test(filled_data,label_file,model_type[i],model_path[i],true_y_test,record_FP,draw_ROC,out_path)
		res.append(res_pro_file)
	return res

