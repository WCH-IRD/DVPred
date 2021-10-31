#! /usr/local/bin/python
# -*- coding: utf-8 -*-

from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold
# from sklearn.cross_validation import cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import minmax_scale
from sklearn.preprocessing import maxabs_scale
from sklearn import metrics
from sklearn import svm
from scipy import interp
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import SelectPercentile
from sklearn.feature_selection import f_classif
from sklearn.feature_selection import f_regression
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import recall_score
from sklearn.metrics import precision_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import accuracy_score
import numpy as np
import random
import decimal
import math
import sys
import os
import re
from sklearn.externals import joblib
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
# import matplotlib.pyplot as plt 
import warnings

from mlxtend.classifier import StackingClassifier
from sklearn.neighbors import KNeighborsClassifier # K 近邻 KNN
from sklearn.naive_bayes import MultinomialNB
from sklearn import decomposition
from sklearn.feature_selection import SelectFromModel
from sklearn.linear_model import LassoCV
import xgboost as xgb

warnings.filterwarnings("ignore")


# 特征值和标签分两个文件放，每行一一对应, 标签分配成0,1,2
def load_data(all_file, shuffle_file):
    shuffle_f = open(shuffle_file, 'w')
    lines = open(all_file).readlines()
    random.shuffle(lines)
    for item in lines:
        shuffle_f.write(item)
    shuffle_f.close()

    p = os.popen("awk -F '\t' 'END{ print NF}' " + shuffle_file).readlines()
    NF = int(p[0].strip('\n'))
    os.system('cat ' + shuffle_file + ' |cut -f 1-7 > prefix.txt')
    os.system('cat  ' + shuffle_file + '|cut -f 8-' + str(NF-2) + ' > feature_all_in_scaled_tab.txt')

    os.system('cat ' + shuffle_file + ' |cut -f ' + str(
        NF-1) + ' | awk \'{if($1=="Pathogenic"){print 0}else{print 1}}\' > label.txt')

    print "Input data: all.txt.training"
    prefix = np.loadtxt('prefix.txt', dtype=np.string_)
    feature = np.loadtxt('feature_all_in_scaled_tab.txt', dtype=np.float32)
    label = np.loadtxt('label.txt', dtype=np.int8)
    NF_1 = int(NF-1)
    return prefix, feature, label, NF_1

# 将于gene相关的Consequence、Symbol、Inheritance、Exon、Domain共186列,使用SparsePCA降维为10列
def PCA():
    fea = np.loadtxt('feature1.txt',dtype=np.float32)
    #fea_scaled = minmax_scale(fea)
    fea_scaled = fea
    pca = decomposition.SparsePCA(n_components=10, alpha=0.8)
    fea_1 = pca.fit_transform(fea_scaled)
    joblib.dump(pca, "pca.m")
    np.savetxt('feature1_scaled.txt', fea_scaled)
    np.savetxt('feature1_0.txt', fea_1)
    os.system("sed 's/ /\t/g' feature1_0.txt > feature1_1.txt")

def calculate_result(actual, pred, model_name):
    target_names = ['Pathogenic', 'Benign']
    report = metrics.classification_report(actual, pred,
                                           target_names=target_names, digits=3)
    tmp = report.split('\n')
    print model_name + '      ' + tmp[0]
    print tmp[2]
    print tmp[3]
    p_precision = tmp[2].split()[1]
    p_recall = tmp[2].split()[2]
    p_f1 = tmp[2].split()[3]
    b_precision = tmp[3].split()[1]
    b_recall = tmp[3].split()[2]
    b_f1 = tmp[3].split()[3]
    #return p_precision, p_recall, b_precision, b_recall
    return p_precision, p_recall, p_f1, b_precision, b_recall, b_f1

def statistic(f1, f2):
    os.system('cat ' + f1 + ' |sort -n > f11.txt')
    single = {line.strip('\n'): 0 for line in open("f11.txt")}
    staf = open('f12.txt', 'w')
    for line in open(f1):
        line = line.strip('\n')
        single[line] = single[line] + 1
    for (k, v) in single.items():
        staf.write(str(v) + '\t' + k + "\n")
    staf.close()
    os.system('cat f12.txt | sort -k 8 -n -r > ' + f2)
#    os.system('rm f11.txt f12.txt')


def LR(X_train_scaled, y_train, X_test_scaled):
    clf = LogisticRegression(penalty='l2', tol=0.01, C=1.0, solver='liblinear', \
                             max_iter=100, multi_class='ovr')
    clf.fit(X_train_scaled, y_train)
    joblib.dump(clf, "lr.m")
    y_pred = clf.predict(X_test_scaled)
    y_pred_poss = clf.predict_proba(X_test_scaled)
    return y_pred, y_pred_poss


def SVM(c, g, X_train_scaled, y_train, X_test_scaled):
    clf = svm.SVC(C=c, gamma=g, probability=True, kernel='rbf')
    clf.fit(X_train_scaled, y_train)
    joblib.dump(clf, "svm.m")
    y_pred = clf.predict(X_test_scaled)
    y_pred_pro = clf.predict_proba(X_test_scaled)
    return y_pred, y_pred_pro


def RandomForest(train_x, train_y, test_x, test_y, para_str):
    print para_str
    para_list = re.split('\[|\]|,',para_str)[1:-1]
    #print para_list[0]
    print para_list
    clf = RandomForestClassifier(n_estimators=int(para_list[0]),
				max_depth=float(para_list[1]),
				min_samples_leaf=int(para_list[2]),
				max_leaf_nodes=int(para_list[3]),
				min_weight_fraction_leaf=float(para_list[4]),
				n_jobs=16)
    clf.fit(train_x, train_y)
    joblib.dump(clf, "rf.m")
    y_pred = clf.predict(test_x)
    y_pred_pro = clf.predict_proba(test_x)
    return y_pred, y_pred_pro


def GBDT(train_x, train_y, test_x, test_y, para_str):
    para_list = re.split('\[|\]|,',para_str)[1:-1]
    clf = GradientBoostingClassifier(n_estimators=int(para_list[0]),
				    max_features=int(para_list[1]),
				    max_depth=int(para_list[2]),
				    learning_rate=float(para_list[3]),
				    min_samples_split=int(para_list[4]))
    clf.fit(train_x, train_y)
    joblib.dump(clf, "gbdt.m")
    y_pred = clf.predict(test_x)
    y_pred_pro = clf.predict_proba(test_x)
    return y_pred, y_pred_pro


def NN(train_x, train_y, test_x, test_y, layers):
    clf = MLPClassifier(hidden_layer_sizes=(layers,), activation='relu',
                        shuffle=True, early_stopping=True) #多层感知器分类器模型/MLP neural_network
    clf.fit(train_x, train_y)
    joblib.dump(clf, "nn.m")
    y_pred = clf.predict(test_x)
    y_pred_pro = clf.predict_proba(test_x)
    return y_pred, y_pred_pro


def KNN(train_x, train_y, test_x, test_y, n):
    clf = KNeighborsClassifier(n_neighbors=n)
    clf.fit(train_x, train_y)
    joblib.dump(clf, "knn.m")
    y_pred = clf.predict(test_x)
    y_pred_pro = clf.predict_proba(test_x)
    return y_pred, y_pred_pro

def NB(train_x, train_y, test_x, test_y ,a):
    clf = MultinomialNB(alpha=a)
    #clf = MultinomialNB(fit_prior=True, class_prior=None)
    clf.fit(train_x, train_y)
    joblib.dump(clf, "mnb.m")
    y_pred = clf.predict(test_x)
    y_pred_pro = clf.predict_proba(test_x)
    return y_pred, y_pred_pro


def Stacking(train_x, train_y, test_x, test_y):
    clf1 = KNeighborsClassifier(n_neighbors=4)
    clf2 = RandomForestClassifier(random_state=1)
    clf3 = MultinomialNB(alpha=0.01)
    clf4 = svm.SVC(C=32, gamma=0.1, probability=True, kernel='rbf')
    clf5 = GradientBoostingClassifier(n_estimators=200)
    clf6 = MLPClassifier(hidden_layer_sizes=(20,), activation='relu',
                        shuffle=True, early_stopping=True)
    lr = LogisticRegression()
    clf = StackingClassifier(classifiers=[clf1, clf2, clf3, clf4, clf5, clf6],
                              use_probas=True,
                              average_probas=False,
                              meta_classifier=lr)
    clf.fit(train_x,train_y)
    joblib.dump(clf, "stacking.m")
    y_pred = clf.predict(test_x)
    y_pred_pro = clf.predict_proba(test_x)
    return y_pred, y_pred_pro

def XGBoost(train_x, train_y, test_x, test_yi, para_str):
    para_list = re.split('\[|\]|,',para_str)[1:-1]
    clf = xgb.XGBClassifier(n_estimators=int(para_list[0]),
                            max_depth=int(para_list[1]),
                            min_child_weight=int(para_list[2]),
                            gamma=float(para_list[3]),
			    subsample=float(para_list[4]),
			    colsample_bytree=float(para_list[5]),
			    reg_alpha=float(para_list[6]),
                            learning_rate=float(para_list[7]))
   # clf=xgb.XGBClassifier(max_depth=7,n_estimators=4,learning_rate=0.01,min_child_weight=3,gamma=0.1,subsample=0.9,colsample_bytree=0.8,reg_alpha=0.01)
    clf.fit(train_x, train_y)
    joblib.dump(clf,"xgboost.m")
    y_pred = clf.predict(test_x)
    y_pred_pro = clf.predict_proba(test_x)
    return y_pred, y_pred_pro


def selectKfeatures(X_train_scaled, y_train, NF):
    # ------------------------------- feature scores---------------------------------
    NF = NF - 8
    print NF
    select = GradientBoostingRegressor(n_estimators=700,max_features=20,max_depth=3,learning_rate=0.1,min_samples_split=2).fit(X_train_scaled,y_train)
    li_score = select.feature_importances_
#    select = SelectKBest(score_func=f_regression(X_train_scaled, y_train), k=NF)
    fea_wei = {}
    #li_score = select.get_params()['score_func'][0]
    for index, item in enumerate(li_score):
        item = float(item)
        fea_wei[index] = abs(item)
    dict = sorted(fea_wei.items(), key=lambda d: d[1], reverse=True)
    a = [tu[0] for tu in dict]
    b = [tu[1] for tu in dict]
    name = [item.strip('\n') for item in open(sys.argv[3])]
    for i, j in zip(a, b):
        print(name[int(i)] + '\t' + str(j))

    print "----------------------------------------------------------------------"


def lasso(X_train_scaled, y_train, NF):
    # select from model Lasso
    #clf = Lasso(alpha=0.1)
    clf = LassoCV()
    li_score = clf.fit(X_train_scaled, y_train).coef_
    #print li_score
    fea_wei = {}
    for index, item in enumerate(li_score):
        item = float(item)
        fea_wei[index] = abs(item)
    dict = sorted(fea_wei.items(), key=lambda d: d[1], reverse=True)
    a = [tu[0] for tu in dict]
    b = [tu[1] for tu in dict]
    name = [item.strip('\n') for item in open(sys.argv[3])]
    for i, j in zip(a, b):
        print(name[int(i)] + '\t' + str(j))

    print "----------------------------------------------------------------------"

def main():
    prefix, feature, label, NF = load_data(sys.argv[1], "shuffle.data")
    kf = KFold(n_splits=5, shuffle=False)
    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    f1 = open(sys.argv[2] + "_judge_p_as_b.txt", "w")
    f2 = open(sys.argv[2] + "_judge_b_as_p.txt", "w")
    para_list = sys.argv[6]
    # 分别记录把致病判断为非致病和把非致病判断为致病的样本，第一列是误判的频次，由高到低降序排列。

    # 如果想要对SVM 进行调参，可以从下面挑选c和g(<1)的值。
    # c=[60,80,100,120,140,160,180,200,220,240,260]
    # c=[8,16,32,64,128,256,512]
    # g=[0.0001,0.001,0.01,0.1,0.5,1]
    c = [32]
    g = [0.1]
    folda = 0
    print_score = sys.argv[4]  # 如果输入是1，表示打印特征分值

    for ci in c:
        for gi in g:
            # print("c/g = " + str(ci)+"/"+str(gi) + "\n")
            sum_p_precision = 0
            sum_p_recall = 0
	    sum_p_f1 = 0
            sum_b_precision = 0
            sum_b_recall = 0
            sum_b_f1 = 0
            sum_roc_auc = 0
            for train_index, test_index in kf.split(feature):
                folda = folda + 1

                train_info = prefix[train_index]
                test_info = prefix[test_index]
                X_train, X_test = feature[train_index], feature[test_index]
                y_train, y_test = label[train_index], label[test_index]
                #X_train_scaled = minmax_scale(X_train) #将属性缩放到一个指定的最大和最小值（通常是1-0）之间
                #X_test_scaled = minmax_scale(X_test)
                X_train_scaled = X_train
                X_test_scaled = X_test
                if print_score == "1":
                    #lasso(X_train_scaled, y_train, NF)
                    selectKfeatures(X_train_scaled, y_train, NF)
		    # print 'begin'
                    print_score = "0"
                if sys.argv[2] == 'LR':
                    y_pred, y_pred_pro = LR(X_train_scaled, y_train,
                                            X_test_scaled)
                    p_precision, p_recall, p_f1, b_precision, b_recall, b_f1 = calculate_result(y_test, y_pred, 'LR')
                if sys.argv[2] == 'SVM':
                    y_pred, y_pred_pro = SVM(ci, gi, X_train_scaled, y_train,
                                             X_test_scaled)
                    p_precision, p_recall, p_f1,b_precision, b_recall, b_f1 = calculate_result(y_test, y_pred, 'SVM')
                if sys.argv[2] == 'RandomForest':
                    y_pred, y_pred_pro = RandomForest(X_train_scaled, y_train,
                                                      X_test_scaled, y_test,
                                                      para_list)  # 100 ok
                    p_precision, p_recall, p_f1,b_precision, b_recall, b_f1 = calculate_result(y_test, y_pred, 'RandomForest')
                if sys.argv[2] == 'GBDT':
                    y_pred, y_pred_pro = GBDT(X_train_scaled, y_train,
                                              X_test_scaled, y_test,
                                              para_list)  # > 100
                    p_precision, p_recall, p_f1, b_precision, b_recall, b_f1 = calculate_result(y_test, y_pred, 'GBDT')
		    PR_precision, PR_recall, PR_thresholds = precision_recall_curve(y_test,y_pred)
		    PR_auc = round(metrics.auc(PR_recall, PR_precision),4)
		    accuracy = round(accuracy_score(y_test, y_pred),4)
                if sys.argv[2] == 'NN':
                    y_pred, y_pred_pro = NN(X_train_scaled, y_train,
                                            X_test_scaled, y_test, 20)
                    p_precision, p_recall, p_f1,b_precision, b_recall, b_f1 = calculate_result(y_test, y_pred, 'NN')
                if sys.argv[2] == 'KNN':
                    #n_range = range(1,26)
                    n_range = [4]
                    for n in n_range:
                        print n
                        y_pred, y_pred_pro = KNN(X_train_scaled, y_train,
                                            X_test_scaled, y_test, n)
                        p_precision, p_recall, p_f1,b_precision, b_recall, b_f1 = calculate_result(y_test, y_pred, 'KNN')
                        test_accuracy = []
                        test_accuracy.append(metrics.accuracy_score(y_test, y_pred))
                if sys.argv[2] == 'NB':
                    alphas = np.logspace(-2, 5, num=200)
                    alphas = [0.01]
                    for alpha in alphas:
                        #print alpha
                        y_pred, y_pred_pro = NB(X_train_scaled, y_train,
                                            X_test_scaled, y_test, alpha)
                        p_precision, p_recall, p_f1,b_precision, b_recall, b_f1 = calculate_result(y_test, y_pred, 'NB')
                        test_accuracy = []
                        test_accuracy.append(metrics.accuracy_score(y_test, y_pred))
                if sys.argv[2] == 'Stacking':
                    y_pred,y_pred_pro = Stacking(X_train_scaled, y_train,
                                            X_test_scaled, y_test)
                    p_precision, p_recall, p_f1,b_precision, b_recall, b_f1 = calculate_result(y_test, y_pred, 'Stacking')
		if sys.argv[2] == 'XGBoost':
                    y_pred, y_pred_pro = XGBoost(X_train_scaled, y_train, X_test_scaled, y_test, para_list)
                    p_precision, p_recall, p_f1, b_precision, b_recall, b_f1 = calculate_result(y_test, y_pred, 'XGBoost')
		    PR_precision, PR_recall, PR_thresholds = precision_recall_curve(y_test,y_pred)
                    PR_auc = round(metrics.auc(PR_recall, PR_precision),4)
		    accuracy = round(accuracy_score(y_test, y_pred),4)
                # statistic of sum of p/b precision and recall
                sum_p_precision += float(p_precision)
                sum_p_recall += float(p_recall)
		sum_p_f1 += float(p_f1)
                sum_b_precision += float(b_precision)
                sum_b_recall += float(b_recall)
		sum_b_f1 += float(b_f1)
                # output y pred pro
                fy = open(sys.argv[2] + "_y_pred_pro_" + str(folda) + ".txt",
                          "w")

                # statistic of FN and FP
                for i in range(len(test_index)):
                    fy.write('\t'.join(list(test_info[i])) + '\t' + str(
                        y_test[i]) + '\t' + str(
                        y_pred[i]) + '\t' + str(y_pred_pro[i][0]) + '\n')
                    if y_pred[i] == 1 and y_test[i] == 0:
                        f1.write('\t'.join(list(test_info[i])) + '\n')
                    elif y_pred[i] == 0 and y_test[i] == 1:
                        f2.write('\t'.join(list(test_info[i])) + '\n')
		fy.close()

                fpr, tpr, thresholds = metrics.roc_curve(y_test,
                                                         y_pred_pro[:, 0],
                                                         pos_label=0)  # 这里的　pos_label　指的是，我们把0作为正类，也就是Pathogenic
                print fpr
                print tpr
                print thresholds
                mean_tpr += interp(mean_fpr, fpr,
                                   tpr)  # 对mean_tpr在mean_fpr处进行插值，通过scipy包调用interp()函数
                mean_tpr[0] = 0.0
                roc_auc = metrics.auc(fpr, tpr)
                sum_roc_auc += roc_auc
                plt.plot(fpr, tpr, lw=1,
                         label='ROC fold %d (area = %0.4f)' % (folda, roc_auc))

            ave_p_precision = sum_p_precision/folda
            ave_p_recall = sum_p_recall/folda
	    ave_p_f1 = sum_p_f1/folda
            ave_b_precision = sum_b_precision/folda
            ave_b_recall = sum_b_recall/folda
            ave_b_f1 = sum_b_f1/folda
            ave_roc_auc = round(sum_roc_auc/folda,4)
            #print "\n\nML\tbinary-class\tprecision\trecall\tAUC"
            #print sys.argv[2]
            print sys.argv[2] + "\tPathogenic\t" + str(ave_p_precision) + "\t" + str(ave_p_recall) + "\t" + str(ave_p_f1 ) + "\t" + str(ave_roc_auc) + "\t" + str(PR_auc) + "\t" + str(accuracy)\
                +  "\n \tBenign\t" + str(ave_b_precision) + "\t" + str(ave_b_recall) + "\t" + str(ave_b_f1 )

    ###  add other software roc
    #print sys.argv[5]
    n_list = sys.argv[5].split('/')
    #print n_list
    for n in n_list:
        label_file = open( n + '_label.txt', 'w')
        feature_file = open( n + '_feature.txt', 'w')
        for line in open('./ML/tmp/' + n + '.data', 'r'):
            tmp = line.strip('\n').split('\t')
            if tmp[0]:
                feature_file.write(tmp[0] + '\n')
                if tmp[1] == "Pathogenic":
                    label_file.write(str(0) + '\n')
                else:
                    label_file.write(str(1) + '\n')
        label_file.close()
        feature_file.close()

        y_test_n = np.loadtxt( n + '_label.txt', dtype=np.int8)
        y_pred_n = np.loadtxt( n + '_feature.txt', dtype=np.float32)
        fpr, tpr, thresholds = metrics.roc_curve(y_test_n, y_pred_n, pos_label=0)
        roc_auc = metrics.auc(fpr, tpr)
        plt.plot(fpr, tpr, lw=1, label='ROC ' + n + ' (area = %0.4f)' % roc_auc)


    f1.close()
    f2.close()
    statistic(sys.argv[2] + "_judge_p_as_b.txt", sys.argv[2] + '_p_b.txt')
    statistic(sys.argv[2] + "_judge_b_as_p.txt", sys.argv[2] + '_b_p.txt')
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic')
    plt.legend(loc="lower right")
    plt.savefig(sys.argv[2] + '_5fold.png')
    plt.show()

if __name__ == '__main__':
    # 训练数据在  ~/ML_pred/ML/tmp下面，
    # 注意，被训练数据，最后一列，应该是标签，Pathogenic', 'Benign'其中一个。
    # 训练的命令是：  python train_ml.py ./ML/tmp/filled.data LR ./ML/list/expand_feature.name 1
    # 当前，能跑的data,也就是西南医院给的训练数据，填充完的data，是~/wangyumei/ml/0824/old_all_data/all1.data
    # python train_ml.py ../old_all_data/all1.data LR ./ML/list/expand_feature.name 1
    # 新生成的文件有6个： 模型（*.m） , 错判累计文件（judge_b_as_p.txt和judge_p_as_b.txt），错判文件（ p_b.txt和b_p.txt），ROC曲线（*_5fold.jpg）
    # 生成的模型是 lr.m 或 svm.m  或 rf.m 或 gbdt.m
    # 如果想把新的模型，放到ML_pred这个工具里，就↓ 当然，为了以防万一，可以把旧模型先保存起来……
    # mv lr.m ./models/LR.m
    # mv svm.m ./models/SVM.m
    # mv rf.m ./models/RandomForest.m
    # mv gbdt.m ./models/GBDT.m
    # mv nn.m ./models/NN.m
    #python  scripts/train_ml.py ./ML/tmp/all00.data KNN ./ML/list/expand_feature.name 0 ./ML/tmp/mcap.data

    main()
# os.system('rm prefix.txt feature.txt label.txt shuffle.data')
