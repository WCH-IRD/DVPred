# -*- coding:utf-8 -*- 
import sys
import os
def usage():
	print "Usage: select a best one CSQ transcript on the basis of SYMBOL, canonical-(YES), Consequence-(priority) and strand-(+1)"
	print "Example: python select_transcript no_head.vcf $feature single_transcript.vcf"
	print "Output: single_transcript.vcf"
# smaller number corresponds to higher priority

def priority_so(sequence,consequence_all):
	tmp = sequence.split('&')
	if len(tmp) == 1:
		return consequence_all.index(sequence)

	elif len(tmp) > 1:
		tmp_list = []
		for i in tmp:
			tmp_list.append(consequence_all.index(i))
		return (min(tmp_list)-0.5)

def split_fea(string):
	tmp = string.split('|')
	e = enumerate(tmp)
	dict_new = {value:key for key,value in dict(e).items()}
	return dict_new

def select_transcript(no_head_path,csq_field,symbol_all,consequence_all,single_transcript_vcf):
	if (not csq_field =='') and  os.path.exists(no_head_path):
		count = 0
		dic_fea = split_fea(csq_field)
		#dic_fea = {'CDS_position':13,'Allele':0,blabla}
		symbol_count = [0 for i in range(len(symbol_all))]
		#symbol_count = [0,0,0,...0]
		special = 0

		single_transcript_vcf = './ML/tmp/'+single_transcript_vcf
		new_file = open(single_transcript_vcf,'w')
		for line in open(no_head_path):
			line=line.strip('\n')
			line = line.replace('DFNA5','GSDME')
			tmp=line.split('\t')
			#print tmp
			CSQ = tmp[7].split('CSQ=')[-1]
			# CSQ = CSQ[4:len(CSQ)]   # CSQ=A|blabla|blabla...
			CSQ = CSQ.replace("DFNA5","GSDME")
			CSQ = CSQ.replace("DFNB31","WHRN")
			CSQ = CSQ.replace("DFNB59","PJVK")
			CSQ = CSQ.replace("GPR98","ADGRV1")
			list_csq = CSQ.split(",")
			list_pri=[]
			
			for i in range(len(list_csq)):
				csq = list_csq[i]
				item = csq.split('|')
				symbol = item[dic_fea['SYMBOL']]
				canonical = item[dic_fea['CANONICAL']]
				consequence = item[ dic_fea['Consequence']]  # maybe xx&&xx
				strand = item[dic_fea['STRAND']]           # +1 or -1
				item_score = len(item)
				pri = 0
				# 1, symbol
				if symbol in symbol_all:
					symbol_count[symbol_all.index(symbol)] = symbol_count[symbol_all.index(symbol)] +1
					pri = pri + 10000
				# 2, consequence
				pri = pri + (35 - priority_so(consequence,consequence_all))*100
				# 3, canonical = 'YES '
				if canonical == 'YES': # all YES records can be find
					pri = pri + 1
				# 4, strand
				if strand == '1':
					pri = pri + 0.3
				# last, item count
				if item_score == len(csq_field.split('|')) :
					pri =pri + 0.1
				list_pri.append(pri)
		
			if  list_pri.count(max(list_pri)) >=1 :
				top = list_pri.index( max(list_pri) )
				if max(list_pri)>10000:
					count = count+1
				new_csq = "CSQ="
				new_csq = new_csq + list_csq[top]
				new_rec = '\t'.join(tmp[0:7]) + '\t' +new_csq +'\t' + '\t'.join(tmp[8:len(tmp)]) +'\n'
				new_file.write(new_rec)
			else :
				print " ERROR !!!! No acceptable csq for this record!"
				print line
				break
		new_file.close()
		# print count
		return single_transcript_vcf
	elif csq_field=='':
		print ('ERROR: Empty csq list!')
		return 'Empty_csq_list'
		sys.exit()
	elif not os.path.exists(no_head_path):
		print ('ERROR: Failed to create no-header vcf! Please make sure you have right to build dictionary.')
		return 'No_vcf_body'
		sys.exit() 	
