# -*- coding:utf-8 -*- 

import numpy as np
import math
import itertools
import sys
import os
import commands
import re

def consequence_value(consequence_all,Consequence):
	con_list= Consequence.split('&')
	res = ['0' for i in range(len(consequence_all))]
	for co in con_list:
		res[consequence_all.index(co)] = '1'
	return res

def pfam_value(pfam_all,pfam):
	res = ['0' for i in range(len(pfam_all))]
	if pfam in pfam_all:
		res[pfam_all.index(pfam)] = '1'
	return res


def symbol(symbol_all,SYMBOL):
	symbol_list=SYMBOL.split('&')
	#ACTG1&693 SYMBOL + '&' + cDNA_pos
	res = ['0' for i in range(len(symbol_all))]
	if symbol_list[0] in symbol_all:
		ind = symbol_all.index(symbol_list[0])
		res[ind] = symbol_list[1]
	return res,len(res)

def inher_pattern(Inheritance_Pattern,Inheritance):
	inh_p = []
	res = ['0' for i in range(len(Inheritance_Pattern))]
	for item in Inheritance.split('/'):
		res[Inheritance_Pattern.index(item)]='1'
	return res


def exon(exon_list):
	res = ['0' for i in range(3)]
	if int(exon_list[0]) == 1:
		res[0] = '1'
#	if int(exon_list[1]) and int(exon_list[1])-int(exon_list[0]) == 0:
	if int(exon_list[1]) and int(exon_list[1]) - int(exon_list[0]) == 1:
		res[-1] = '1'
#	if not(int(exon_list[0]) == 1) and int(exon_list[1]) and not(int(exon_list[1])-int(exon_list[0]) == 0):
	if not (int(exon_list[0]) == 1) and  not (int(exon_list[1]) - int(exon_list[0]) == 1):
		res[1] = '1'
	return res

def sift(feature):
	return feature.split('(')[-1].split(')')[0]


def domain_value(file,index):
	db = []
	for line in open(file):
		line=line.strip('\n').split('\t')[index]
		if not line == '':
			tmp = line.split('&')
			for item in tmp:
				if 'Pfam_domain' in item:
					ind = item.index(':')
					source = item[ind+1:len(item)]
					if source not in db:
						db.append(source)
	return db

def domain(domain_values,domain_list):
	res = ['0'for i in range(len(domain_values)+1)]
	for item in domain_list:
		if 'Pfam_domain' in item:
			ind = item.index(':')
			source = item[ind+1:len(item)]
			res[domain_values.index(source)] = '1'
			# print 'found Pfam_domain.'
	return res


def missing_val(tmp, fea, ini_fea_dic_new):
	tmp_fea= tmp[ini_fea_dic_new[fea] + 7]
	#print fea
	#print tmp_fea
	code = ['0' for i in range(2)]
	code[0] = tmp_fea
	if not tmp_fea == '':
		code[1] = "1"
	#print code
	tmp[ini_fea_dic_new[fea] + 7] = '\t'.join(code)


# expand_feature  (initial_feature_data,list_feature,list_symbol,list_consequence,'expand_feature_data','expand_feature_name') symbol=gene
#def expand_feature(initial_feature_data,ini_fea,symbol_all,consequence_all,expand_feature_data,expand_feature_name):
def expand_feature(initial_feature_data,ini_fea,symbol_all,consequence_all,domain_all,expand_feature_data,expand_feature_name):
	gene_fea=0
	if (not ini_fea==[]) and os.path.exists(initial_feature_data):
		ini_fea_dic = { value:key for key,value in enumerate(ini_fea)}
		Inheritance_Pattern=[]
		EXON=[]
		sift_v=[]
		expand_feature_name = './ML/list/' + expand_feature_name
		expand_feature_data = './ML/tmp/'+ expand_feature_data

		expand_feature = open(expand_feature_name,'w')
		newf = open(expand_feature_data,'w')
		for line in open(initial_feature_data):
			line = line.strip('\n')
			tmp = line.split('\t')
			#  Consequence---------------------------------------------------------------
			Consequence = tmp[ini_fea_dic['Consequence']+7]
			code = ['0'for i in range(len(consequence_all))]
			if not Consequence=='':
				code = consequence_value(consequence_all,Consequence)
			tmp[ini_fea_dic['Consequence']+7] =  '\t'.join(code)

			# SYMBOL---------------------------------------------------------------------
			SYMBOL = tmp[ini_fea_dic['SYMBOL']+7]
			code =[ '0' for i in range(len(symbol_all))]
			if not SYMBOL=='':
				code,symbol_count = symbol(symbol_all,SYMBOL)
			tmp[ini_fea_dic['SYMBOL']+7] = '\t'.join(code)

			# Inheritance ---------------------------------------------------------------
#			Inheritance_Pattern = ['AD','AR','X-LINK']
			Inheritance_Pattern = ['AD','AR','XL','XLD','XLR','MT']
			Inheritance = tmp[ini_fea_dic['Inheritance']+7]
			code =[ '0' for i in range(len(Inheritance_Pattern))]
			if not Inheritance=='':
				code = inher_pattern(Inheritance_Pattern,Inheritance)
			tmp[ini_fea_dic['Inheritance']+7] = '\t'.join(code)

			#  EXON----------------------------------------------------------------------
			EXON = tmp[ini_fea_dic['EXON']+7]
			code =['0' for i in range(3)]
			if not EXON=='':
				exon_list = EXON.split('/')
				code = exon(exon_list)
			tmp[ini_fea_dic['EXON']+7] = '\t'.join(code)

			# Pfam-----------------------------------------------------------------------
			Pfam = tmp[ini_fea_dic['DOMAINS']+7]
			code = ['0' for i in range(len(domain_all))]
			if not Pfam=='':
				code = pfam_value(domain_all,Pfam)
			tmp[ini_fea_dic['DOMAINS']+7] = '\t'.join(code)


			#  SIFT----------------------------------------------------------------------
#			sift_v = ['tolerated', 'deleterious', 'tolerated_low_confidence', 'deleterious_low_confidence']
#			sift_c = [0.0 for i in range(4)]
#			for i in range(4):
#				if tmp[ini_fea_dic['SIFT']+7].split('(')[0] == sift_v[i]:
#					sift_c[i] = sift_c[i] + float(sift(tmp[ini_fea_dic['SIFT']+7]))
#			sift_str = [str(i) for i in sift_c]
#			tmp[ini_fea_dic['SIFT']+7]='\t'.join(sift_str)


			#  PolyPhen----------------------------------------------------------------------
			if 'Polyphen2_HDIV_score' in ini_fea_dic:
                                #tmp[ini_fea_dic['Polyphen2_HDIV_score']+7] = sift(tmp[ini_fea_dic['Polyphen2_HDIV_score']+7])
                                tmp[ini_fea_dic['Polyphen2_HDIV_score']+7] = (tmp[ini_fea_dic['Polyphen2_HDIV_score']+7]).split('&')[0]
                                if tmp[ini_fea_dic['Polyphen2_HDIV_score']+7] == '.':
                                        tmp[ini_fea_dic['Polyphen2_HDIV_score']+7] = ''


                        if 'Polyphen2_HVAR_score' in ini_fea_dic:
                                tmp[ini_fea_dic['Polyphen2_HVAR_score']+7] = (tmp[ini_fea_dic['Polyphen2_HVAR_score']+7]).split('&')[0]
                                if tmp[ini_fea_dic['Polyphen2_HVAR_score']+7] == '.':
                                        tmp[ini_fea_dic['Polyphen2_HVAR_score']+7] = ''

			# LoFtool------------------------------------------------------------------------
			if 'LoFtool' in ini_fea_dic:
				tmp[ini_fea_dic['LoFtool']+7] = (tmp[ini_fea_dic['LoFtool']+7]).split('(')[0]			
				if tmp[ini_fea_dic['LoFtool']+7] == '-':
					tmp[ini_fea_dic['LoFtool']+7] = ''

			# fitscons_HL--------------------------------------------------------------------
			if 'fitscons' in ini_fea_dic:
				tmp[ini_fea_dic['fitscons']+7] = (tmp[ini_fea_dic['fitscons']+7]).split('&')[0] 

			# phastCons.100way_HL-------------------------------------------------------------
			if 'phastCons.100way' in ini_fea_dic:
				tmp[ini_fea_dic['phastCons.100way']+7] = (tmp[ini_fea_dic['phastCons.100way']+7]).split('&')[0] 

			#Gene_score_Petrovski_HL----------------------------------------------------------
			if 'Gene_score_Petrovski_All' in ini_fea_dic:
				tmp[ini_fea_dic['Gene_score_Petrovski_All']+7] = (tmp[ini_fea_dic['Gene_score_Petrovski_All']+7]).split('(')[0]			

			# IMPACT-------------------------------------------------------------------------
			if 'IMPACT' in ini_fea_dic:
				if tmp[ini_fea_dic['IMPACT']+7] == 'HIGH':
					tmp[ini_fea_dic['IMPACT']+7] = '4'
				elif tmp[ini_fea_dic['IMPACT']+7] == 'MODERATE':
					tmp[ini_fea_dic['IMPACT']+7] = '3'
				elif tmp[ini_fea_dic['IMPACT']+7] == 'MODIFIER':
					tmp[ini_fea_dic['IMPACT']+7] = '2'
				elif tmp[ini_fea_dic['IMPACT']+7] == 'LOW':
					tmp[ini_fea_dic['IMPACT']+7] = '1'

			# FATHMM_score ----------------------------------------------------------
			if 'FATHMM_score' in ini_fea_dic:
				tmp[ini_fea_dic['FATHMM_score']+7] = (tmp[ini_fea_dic['FATHMM_score']+7]).split('&')[0]			
				if tmp[ini_fea_dic['FATHMM_score']+7] == '.':
					tmp[ini_fea_dic['FATHMM_score']+7] = ''


			# FATHMM_score ----------------------------------------------------------
			if 'PROVEAN_score' in ini_fea_dic:
				tmp[ini_fea_dic['PROVEAN_score']+7] = (tmp[ini_fea_dic['PROVEAN_score']+7]).split('&')[0]			
				if tmp[ini_fea_dic['PROVEAN_score']+7] == '.':
					tmp[ini_fea_dic['PROVEAN_score']+7] = ''

	                # all & value--------------------------------------------------------------------
            		for i in ini_fea:
                		if re.findall("&", tmp[ini_fea_dic[i] + 7]):
                    			tmp[ini_fea_dic[i] + 7] = (tmp[ini_fea_dic[i] + 7]).split('&')[0]
				if tmp[ini_fea_dic[i] + 7] == '-':
					tmp[ini_fea_dic[i] + 7] = ''
			# AF AA_AF EA_AF ExAC_AF---------------------------------------------------------
			if 'AF' in ini_fea_dic:
				if '&' in tmp[ini_fea_dic['AF']+7]:
					tmp[ini_fea_dic['AF']+7]=tmp[ini_fea_dic['AF']+7].split('&')[0]
			if 'AA_AF' in ini_fea_dic:
				if '&' in tmp[ini_fea_dic['AA_AF']+7]:
					tmp[ini_fea_dic['AA_AF']+7]=tmp[ini_fea_dic['AA_AF']+7].split('&')[0]
			if 'EA_AF' in ini_fea_dic:
				if '&' in tmp[ini_fea_dic['EA_AF']+7]:
					tmp[ini_fea_dic['EA_AF']+7]=tmp[ini_fea_dic['EA_AF']+7].split('&')[0]
			if 'ExAC_AF' in ini_fea_dic:
				if '&' in tmp[ini_fea_dic['ExAC_AF']+7]:
					tmp[ini_fea_dic['ExAC_AF']+7]=tmp[ini_fea_dic['ExAC_AF']+7].split('&')[0]


			#  PUBMED -----------------------------------------------------------------------
			'''if not tmp[ini_fea_dic['PUBMED']+7]=='0':
				tmp[ini_fea_dic['PUBMED']+7] = '1'
			#  MutationTaster_score ----------------------------------------------------------
			if '&' in tmp[ini_fea_dic['MutationTaster_score']+7]: 
				it = tmp[ini_fea_dic['MutationTaster_score']+7].split('&')
				its = [float(i) for i in it]
				tmp[ini_fea_dic['MutationTaster_score']+7] = str(max(its))
			'''#  gnomad --------------------------------------------------------------------------
			# GERP-----------------------------------------------------------------------------
			if ':' in tmp[ini_fea_dic['GERP']+7]:
				tmp[ini_fea_dic['GERP']+7]=""
			newf.write('\t'.join(tmp) + '\n')
		newf.close()

			
		# ----------------all expanded features-------------------------------------------------
		for item in ini_fea:
			if item == 'Consequence':
				for i in consequence_all:
					gene_fea+=1
					expand_feature.write('Consequence:'+i+'\n')
			elif item == 'SYMBOL':
				for i in symbol_all:
					gene_fea+=1
					expand_feature.write("SYMBOL:"+i+'\n') 
			elif item == 'Inheritance':
				for i in Inheritance_Pattern:
					gene_fea+=1
					expand_feature.write("Inheritance:"+i+'\n') 
			elif item == 'EXON':
				ind = [i for i in range(3)]
				for i in ind:
					expand_feature.write('EXON'+str(i)+'\n')
			elif item == 'DOMAINS':
				for i in domain_all:
					gene_fea+=1
					expand_feature.write('DOMAIN:'+i+'\n')

			else:
				expand_feature.write(item+'\n')

		expand_feature.close()

		# ----------------add val feature-------------------------------------------------
		expand_feature_data = "./ML/tmp/expand_feature.data"
		os.system('mv ./ML/list/expand_feature.name ./ML/list/expand_feature.name.0')
		expand_feature = "./ML/list/expand_feature.name.0"
		p=os.popen("awk -F '\t' 'END{ print NF}' " + expand_feature_data).readlines()
		NF = p[0].strip('\n')
		# print col
		os.environ['col']=NF
		blank_col= os.popen("awk -F '\t' -v c=$col 'BEGIN{for(j=0;j<c;j++){ct[j]=0}}{for(i=1;i<=c;i++){if(ct[i-1]==0){ if($i==\"\"){ct[i-1]=1}}}} END{ for(t in ct){if(ct[t]==1){print t;} }}' " + expand_feature_data).readlines()
		blank=[int(i.strip('\n')) for i in blank_col]
		print blank
		fea_list = []
		for i in blank:
			i -= 6
			fea0 = commands.getoutput('cat ' + expand_feature + '|awk -v c=' + str(i) + ' \'{if(NR==c)print $0}\' ')
			fea_list.append(fea0)
		print fea_list
		#fea_list = ['GERP', 'SIFT_RAW', 'PolyPhen', 'CADD_RAW', 'Eigen-phred', 'MetaLR_score', 'MetaSVM_score', 'MutPred_rankscore', 'MutationTaster_score', 'REVEL_rankscore']
		#fea_list = ['SIFT_RAW', 'Polyphen2_HDIV_score', 'Polyphen2_HVAR_score', 'esp6500siv2_all', 'Gene_score_Aggarwala_All', 'Gene_score_Petrovski_All', 'LoFtool', 'GenicIntolerance_ALL_0.01%', 'GenicIntolerance_%ALL_0.01%', 'GenicIntolerance_%AA_1%', 'GenicIntolerance_ALL_0.1%', 'GenicIntolerance_%ExAC_0.1%popn', 'GenicIntolerance_%ALL_0.1%', 'GenicIntolerance_%ExAC_0.05%popn', 'GenicIntolerance_ALL_1%', 'GenicIntolerance_%ExAC_0.01%', 'Eigen_EIGEN-PC-RAW', 'MetaLR_score', 'MetaSVM_score',  'GERP', 'MutationTaster_score', 'FATHMM_score', 'REVEL_RAW', 'PROVEAN_score', 'VEST4_VEST4', 'fathmm-MKL_coding_score', 'GenicIntolerance_%ALL_1%', 'GenoCanyon_score', 'GenicIntolerance_PP2_ALL_0.1%', 'fitscons', 'GenicIntolerance_%PP2_ALL_0.1%', 'phyloP100way', 'GenicIntolerance_EA_0.1%', 'GenicIntolerance_%EA_0.1%', 'GenicIntolerance_EA_1%', 'GenicIntolerance_%EA_1%', 'GenicIntolerance_AA_0.1%', 'GenicIntolerance_%AA_0.1%', 'GenicIntolerance_AA_1%']
		#fea_list = ['SIFT_RAW', 'Polyphen2_HDIV_score', 'Polyphen2_HVAR_score', 'esp6500siv2_all', 'Gene_score_Aggarwala_All', 'Gene_score_Petrovski_All', 'LoFtool', 'GenicIntolerance_ALL_0.01%', 'GenicIntolerance_%ALL_0.01%', 'GenicIntolerance_%AA_1%', 'GenicIntolerance_ALL_0.1%', 'GenicIntolerance_%ExAC_0.1%popn', 'GenicIntolerance_%ALL_0.1%', 'GenicIntolerance_%ExAC_0.05%popn', 'GenicIntolerance_ALL_1%', 'GenicIntolerance_%ExAC_0.01%', 'Eigen_EIGEN-PC-RAW', 'MetaLR_score', 'MetaSVM_score', 'MutPred_rankscore', 'GERP', 'MutationTaster_score', 'FATHMM_score', 'REVEL_RAW', 'PROVEAN_score', 'VEST4_VEST4', 'fathmm-MKL_score', 'GenicIntolerance_%ALL_1%', 'GenoCanyon', 'GenicIntolerance_PP2_ALL_0.1%', 'fitscons', 'GenicIntolerance_%PP2_ALL_0.1%', 'phyloP100way', 'GenicIntolerance_EA_0.1%', 'GenicIntolerance_%EA_0.1%', 'GenicIntolerance_EA_1%', 'GenicIntolerance_%EA_1%', 'GenicIntolerance_AA_0.1%', 'GenicIntolerance_%AA_0.1%', 'GenicIntolerance_AA_1%','CDTS','M-CAP_score']
	#	fea_list = ['SIFT_RAW', 'Polyphen2_HDIV_score', 'Polyphen2_HVAR_score', 'esp6500siv2_all', 'Gene_score_Aggarwala_All', 'Gene_score_Petrovski_All', 'LoFtool', 'GenicIntolerance_ALL_0.01%', 'GenicIntolerance_%ALL_0.01%', 'GenicIntolerance_%AA_1%', 'GenicIntolerance_ALL_0.1%', 'GenicIntolerance_%ExAC_0.1%popn', 'GenicIntolerance_%ALL_0.1%', 'GenicIntolerance_%ExAC_0.05%popn', 'GenicIntolerance_ALL_1%', 'GenicIntolerance_%ExAC_0.01%', 'Eigen_EIGEN-PC-RAW', 'MetaLR_score', 'MetaSVM_score', 'MutPred_rankscore', 'GERP', 'MutationTaster_score', 'FATHMM_score', 'REVEL_RAW', 'PROVEAN_score', 'VEST4_VEST4', 'fathmm-MKL_score', 'GenicIntolerance_%ALL_1%', 'GenoCanyon', 'GenicIntolerance_PP2_ALL_0.1%', 'fitscons', 'GenicIntolerance_%PP2_ALL_0.1%', 'phyloP100way', 'GenicIntolerance_EA_0.1%', 'GenicIntolerance_%EA_0.1%', 'GenicIntolerance_EA_1%', 'GenicIntolerance_%EA_1%', 'GenicIntolerance_AA_0.1%', 'GenicIntolerance_%AA_0.1%', 'GenicIntolerance_AA_1%','CDTS']
		with open("./ML/list/expand_feature.name", 'w+') as f:
			for line in open("./ML/list/expand_feature.name.0", 'r+'):
				line = line.strip('\n')
				#print line
				if line in fea_list:
					f.write(line + '\n' + line + '_val\n')
				else:
					f.write(line + '\n')
		f.close()

		os.system('mv ./ML/tmp/expand_feature.data ./ML/tmp/expand_feature.data.0')
		expand_feature = [line.strip('\n') for line in open("./ML/list/expand_feature.name.0")]
		ini_fea_dic_new = {value: key for key, value in enumerate(expand_feature)}
		#print ini_fea_dic_new
		with open("./ML/tmp/expand_feature.data", "w+") as newf:
			for line in open("./ML/tmp/expand_feature.data.0", 'r'):
				line = line.strip('\n')
				tmp = line.split('\t')
				for f in fea_list:
					missing_val(tmp, f, ini_fea_dic_new)
				newf.write('\t'.join(tmp) + '\n')
		newf.close()
		gene_fea+=3
		print(gene_fea)
		os.system('head -n '+str(gene_fea)+ ' ./ML/list/expand_feature.name > ./ML/list/expand_feature.name.gene' )
		os.system("sed '1," + str(gene_fea) +"d' ./ML/list/expand_feature.name >./ML/list/expand_feature.name.predict")
		os.system('cp ./ML/list/expand_feature.name ./ML/list/expand_feature.name.all')


		return expand_feature_name,expand_feature_data
	elif ini_fea==[]:
		print ('ERROR: Empty feature list!')
		return '0','Empty_F_list'
		sys.exit()
	elif not os.path.exists(initial_feature_data):
		print ('ERROR: Failed to open initial-feature_data.')
		return '0','NO_ini_fea_data'
		sys.exit()










