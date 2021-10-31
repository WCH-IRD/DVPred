# -*- coding: utf-8 -*-
from sklearn.preprocessing import minmax_scale
import numpy as np
import random
import os
import warnings
from sklearn import decomposition
from sklearn.externals import joblib

warnings.filterwarnings("ignore")


def load_data(all_file, shuffle_file, scaled_file):

    shuffle_file = all_file

    p = os.popen("awk -F '\t' 'END{ print NF}' " + shuffle_file).readlines()
    NF = int(p[0].strip('\n'))
    wc = os.popen("wc -l ML/list/expand_feature.name.gene|awk -F ' ' '{print $1}' -").readlines()
    wcl = int(wc[0].strip()) + 7
    os.system('cat ' + shuffle_file + ' |cut -f 1-7 > prefix.txt')
    os.system('cat ' + shuffle_file + ' |cut -f 8-'+ str(wcl) + ' > feature_gene.txt')
    PCA()
    os.system(
        'cat  ' + shuffle_file + '|cut -f'+ str(wcl+1) +'-' + str(NF - 2) + ' > feature_predict.txt')
    os.system('cat ' + shuffle_file + ' |cut -f ' + str(NF - 1) + ' > label.txt')
    os.system('cat  ' + shuffle_file + '|cut -f ' + str(NF) + ' > name.txt')
    os.system('cat ' + shuffle_file + ' |cut -f ' + str(
        NF-1) + ' | awk \'{if($1=="Pathogenic"){print 0}else{print 1}}\' > label_n.txt')
    os.system(
        "paste prefix.txt feature_pca.txt feature_predict.txt label.txt name.txt > " + scaled_file)
    return scaled_file

def PCA():
    fea = np.loadtxt('feature_gene.txt',dtype=np.float32)
    fea_scaled = minmax_scale(fea)
    #fea_scaled = fea
    pca = decomposition.SparsePCA(n_components=10, alpha=0.8)
    fea_1 = pca.fit_transform(fea_scaled)
    joblib.dump(pca, "pca.m")
    np.savetxt('feature_gene_scaled_0.txt', fea_scaled)
    os.system("sed 's/ /\t/g' feature_gene_scaled_0.txt > feature_gene_scaled.txt")
    np.savetxt('feature_pca_0.txt', fea_1)
    os.system("sed 's/ /\t/g' feature_pca_0.txt > feature_pca.txt")


def grep_data(scaled_in_file, name):
    outfile = './ML/tmp/' + name + '_scaled.data'
    os.system("grep " + name + " " + scaled_in_file + " > " + outfile)
    return outfile


# normal_all_data,training_norm_data,10k_norm_data,morl_norm_data = norm_data(all_data,input_base,"10k","morl")
def norm_data(all_data_in, training_name, test1_name, test2_name):
    scaled_file = load_data(all_data_in, "shuffle.data", "scaled.data")
    tr_file = grep_data(scaled_file, training_name)
    te1_file = grep_data(scaled_file, test1_name)
    te2_file = grep_data(scaled_file, test2_name)
    return scaled_file, tr_file, te1_file, te2_file

