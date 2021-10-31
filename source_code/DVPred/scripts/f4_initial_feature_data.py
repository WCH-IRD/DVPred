# -*- coding:utf-8 -*- 
import sys
import os
import re
# from 1.select_csq.py import split_fea

def usage():
	print '  Usage: fetch the specified initial features from single_transcript.vcf  '
	print "Example: python 2.split_csq_to_feature.vcf single_transcript.vcf all_feature selected_feature initial_feature_data.txt"
	#               python 2.split_csq_to_feature.py single_transcript1.vcf $feature initial_feature_name.txt initial_feature_data1.txt inheritance_Pattern.name
	print "Output: initial_feature_data.txt (split by tab)"
def split_fea(string):
	tmp = string.split('|')
	e = enumerate(tmp)
	dict_new = {value:key for key,value in dict(e).items()}
	return dict_new

def find_inher_pattern(symbol,symbol_inhe):
	flag = ''
	for item in symbol_inhe:
		item = item.split('\t')
		if symbol == item[0]:
			flag = item[1]
			break
	return flag



def initial_feature(single_transcript_vcf,csq_field,ini_fea,list_inheritance,initial_feature_data,roc_files):
	if not ini_fea ==[] and  os.path.exists(single_transcript_vcf):
		
		dic_fea = split_fea(csq_field)
		symbol_inhe =list_inheritance

		#roc_files = ["MCAP_RAW", "REVEL_rankscore", "CADD_RAW", "rf_score", "MetaLR_score", "SIFT_RAW", "PolyPhen"]
		for n in roc_files:
			#file_name = n + '_file'
			file_path = './ML/tmp/' + n + '.txt'
			n = open(file_path,'w')
        	for fea in ini_fea:
            		if (fea not in dic_fea) and (fea !=  "Inheritance"):
                		print ('ERROR: feature:"' + fea + '" not in input vcf file !')
                		sys.exit()


		initial_feature_data = './ML/tmp/'+initial_feature_data
		new_file = open(initial_feature_data,'w')
		# error =open('error_split_all.txt','w')
		for line in open(single_transcript_vcf):
			line=line.strip('\n')
			tmp=line.split('\t')
			CSQ = tmp[7].split('CSQ=')[-1]
			# CSQ = CSQ[4:len(CSQ)]   # CSQ=A|blabla
			list_csq = CSQ.split(",")

			if len(list_csq)>1 :
				print 'ERROR！！ more than one transcript left!'
				sys.exit()
			item = CSQ.split('|')
			tt=0
			if  not len(item) == len(dic_fea): # 左边item 是CSQ的区域个数右边。 dic_fea是文件头所有区域的个数，在老的wym-all.vcf中是93个 。
				tt=tt+1
				print line
				print  len(item)

            		#if 'SIFT_RAW' in ini_fea:
                		#print item[dic_fea['SIFT_RAW']]
			# 上面写到检测到区域间隔不一致就break了
			if 'SYMBOL' in ini_fea:
				cDNA_pos='0'
				if item[dic_fea['cDNA_position']]=="":
					cDNA_pos='0'
				elif '?' not in item[dic_fea['cDNA_position']]:
					cDNA_pos = item[dic_fea['cDNA_position']].split('-')[0]
				elif item[dic_fea['cDNA_position']].split('-')[0]=="?" :
					cDNA_pos=item[dic_fea['cDNA_position']].split('-')[1]
				elif item[dic_fea['cDNA_position']].split('-')[1]=="?" :
					cDNA_pos=item[dic_fea['cDNA_position']].split('-')[0]
				Inheritance_Pattern = find_inher_pattern(item[dic_fea['SYMBOL']],symbol_inhe)
				item[dic_fea['SYMBOL']] = item[dic_fea['SYMBOL']] + '&' + cDNA_pos

			# Domains 中提取出Pfam_domain值
			if 'DOMAINS' in ini_fea:
				#print item[dic_fea['DOMAINS']]
				item[dic_fea['DOMAINS']] = re.findall(r'Pfam_domain:(\w+\d+)',item[dic_fea['DOMAINS']])
				if(item[dic_fea['DOMAINS']]):
					item[dic_fea['DOMAINS']] = item[dic_fea['DOMAINS']][0]
				else:
					item[dic_fea['DOMAINS']] = ''
				#print item[dic_fea['DOMAINS']]

			# 一些频率特征的缺失值可以直接补零
			if 'PUBMED' in ini_fea and 'AF' in ini_fea and 'AA_AF' in ini_fea \
				and 'EA_AF' in ini_fea and 'ExAC_AF'in ini_fea:

				if item[dic_fea['PUBMED']]=='':
					item[dic_fea['PUBMED']]='0'
				else: 
					item[dic_fea['PUBMED']]='1'
				if item[dic_fea['AF']] == '':
					item[dic_fea['AF']]='0.000001'
				if item[dic_fea['AA_AF']]=='':
					item[dic_fea['AA_AF']]='0.000001'
				if item[dic_fea['EA_AF']]=='':
					item[dic_fea['EA_AF']]='0.000001'
				if item[dic_fea['ExAC_AF']]=='':
					item[dic_fea['ExAC_AF']]='0.000001'
			if 'rf_score' in ini_fea:
				if item[dic_fea['rf_score']] =='':
					item[dic_fea['rf_score']]='0'
			if 'dpsi_max_tissue' in ini_fea:
				if item[dic_fea['dpsi_max_tissue']] =='':
					item[dic_fea['dpsi_max_tissue']]='0'
				
			if 'gnomad_AF' in dic_fea:
				if item[dic_fea['gnomad_AF']] =='' or item[dic_fea['gnomad_AF']] =='.':
					item[dic_fea['gnomad_AF']]='0.000001'
			if 'gnomad_AF_AFR' in dic_fea:
				if item[dic_fea['gnomad_AF_AFR']] =='' or item[dic_fea['gnomad_AF_AFR']] =='.':
					item[dic_fea['gnomad_AF_AFR']]='0.000001'
			if 'gnomad_AF_AMR' in dic_fea:
				if item[dic_fea['gnomad_AF_AMR']] =='' or item[dic_fea['gnomad_AF_AMR']] =='.':
					item[dic_fea['gnomad_AF_AMR']]='0.000001'
			if 'gnomad_AF_ASJ' in dic_fea:
				if item[dic_fea['gnomad_AF_ASJ']] =='' or item[dic_fea['gnomad_AF_ASJ']] =='.':
					item[dic_fea['gnomad_AF_ASJ']]='0.000001'
			if 'gnomad_AF_EAS' in dic_fea:
				if item[dic_fea['gnomad_AF_EAS']] =='' or item[dic_fea['gnomad_AF_EAS']] =='.':
					item[dic_fea['gnomad_AF_EAS']]='0.000001'
			if 'gnomad_AF_FIN' in dic_fea:
				if item[dic_fea['gnomad_AF_FIN']] =='' or item[dic_fea['gnomad_AF_FIN']] =='.':
					item[dic_fea['gnomad_AF_FIN']]='0.000001'
			if 'gnomad_AF_NFE' in dic_fea:
				if item[dic_fea['gnomad_AF_NFE']] =='' or item[dic_fea['gnomad_AF_NFE']] =='.':
					item[dic_fea['gnomad_AF_NFE']]='0.000001'
			if 'gnomad_AF_OTH' in dic_fea:
				if item[dic_fea['gnomad_AF_OTH']] =='' or item[dic_fea['gnomad_AF_OTH']] =='.':
					item[dic_fea['gnomad_AF_OTH']]='0.000001'
			if 'gnomad_AF_SAS' in dic_fea:
				if item[dic_fea['gnomad_AF_SAS']] =='' or item[dic_fea['gnomad_AF_SAS']] =='.':
					item[dic_fea['gnomad_AF_SAS']]='0.000001'
			if 'esp6500siv2_all' in dic_fea:
				if item[dic_fea['esp6500siv2_all']] == '' or item[dic_fea['esp6500siv2_all']] == '.':
					item[dic_fea['esp6500siv2_all']]='0.000001'
				
			new_csq_list = []

			for fea in ini_fea:
				#print fea
				if not (len(item) == len(dic_fea)):
					print len(item)
					print len(dic_fea)
					print line
					print item
					print '------------'
				else:
					if fea == 'rmsk_simpleRepeat':
						non_repeat_region='0'
						if not item[dic_fea['rmsk_simpleRepeat']] == '':
							non_repeat_region='1'
						new_csq_list.append(non_repeat_region)

					elif fea ==  'Inheritance':
						new_csq_list.append(Inheritance_Pattern)
					elif fea in dic_fea.keys():
						new_csq_list.append(item[int(dic_fea[fea])])

			for n in roc_files:
				if n in dic_fea:
					file_name = n + '_file'
					file_path = './ML/tmp/' + n + '.txt'
					f_r = open(file_path, 'a')
                    			if re.findall("&", item[dic_fea[n]]):
                        			item[dic_fea[n]] = (item[dic_fea[n]]).split('&')[0]
					f_r.write(item[dic_fea[n]] + '\n')
				#print item[dic_fea['MCAP_RAW']]
				#mcap_file.write(item[dic_fea['MCAP_RAW']] + '\n')


			new_csq = '\t'.join(new_csq_list)
			new_rec = '\t'.join(tmp[0:7]) + '\t' +new_csq +'\n'
			new_file.write(new_rec)
		#mcap_file.close()
		#file_name.close()
		new_file.close()

		return initial_feature_data
	elif ini_fea ==[]:
		print ('ERROR: Empty feature list!')
		return 'Empty_F_list'
		sys.exit()
	elif not os.path.exists(single_transcript_vcf):
		print ('ERROR: Failed to create single-transcript vcf.')
		return 'No_single_transcript_vcf'
		sys.exit()

