
# -*- coding:utf-8 -*-  
import random
import sys
import os
import numpy as np
import math

def reduce_none(dic):
	newd = dic.copy()
	low = 25.0
	high = -25.0
	count =0
	for key in dic:
		if dic[key]=='':
			count = count+1
			del newd[key]
		else:
			#print dic[key]
			if float(dic[key]) < low:
				low=float(dic[key])
			if float(dic[key]) > high:
				high =float(dic[key])
	if count == len(dic):
		pass
#		print "all 0 !!!!!!!!!!!!!!!!________________-------------------++++++++++++++++++"
	return newd,low,high

def div_to_sets(dic,layers):
	# layer_list中放的元素：每个小集合里的元素个数

	layer_list = [len(dic)/layers for i in range(layers)]
	for i in range(len(dic)%layers):
		layer_list[i] = layer_list[i] + 1
	# layers = layers -1
	# 打乱 全部数据的索引
	ind_list = range(len(dic))
	random.shuffle(ind_list)
	# sets：  [{一个小集合},{又一个小集合},...{很多小集合}]
	sets = [ {} for t in range(layers) ]

	ind = 0
	for i in range(layers): # 数据 分到 几个不同的小集合(称为层)中
		index = ind_list[ind]
		flag = 0
		while len(sets[i]) < layer_list[i]:
			# print 'into while--------'
			index = ind_list[ind]
			sets[i][index] = dic[index]
			ind = ind +1
			if not dic[index]=='':
				flag=1
	return sets,layers,layer_list

def has_non(dic):
	flag = 0
	for it in dic:
		if dic[it]=='':
			flag =1
	return flag

def first_no_0_bit(dec):
	t = 0
	while dec < 1.0:
		dec = dec*10
		t = t+1
	return t
   
def pad(dic):
	newd,low,high = reduce_none(dic)
	# print newd
	#print low
    	if not newd:
        	f = 0
		#print "all 0 !!!!!!!!!!!!!!!!"
        	return f
	newl=[]
	if max(abs(low),abs(high)) >= 10: # （-5～20）（0,～25）（-12.5~7）（-85~15）
		step = 1.0
	elif (abs(low) < 10 and abs(low) >= 0) and ( abs(high) < 10 and abs(high) >= 0 ): #（0~1）（-1.5~1.5）
		step = 0.1
	for i in newd:
		tp = int(float(newd[i])/step)
		newl.append(tp)

	tmp = newd[random.choice( newd.keys() )]
	floa = float(tmp) - step*int(float(tmp)/step)
	inte = float(random.choice( newl ))
	f =inte*step+floa 
	if f > high:
		f=high
	return f

def hot_deck(dic,layers):
	sets,layers,layer_list = div_to_sets(dic,layers)

	rand = [0.0]
	# print layers
	for layer in range(layers):
		# 本层的数据是 sets[layer]
		# 1.随机抽取 n-1 个 0~1之间的小数，排序
		while len(rand) < layer_list[layer]:
			tmp = random.random()
			if tmp not in rand:
				rand.append(tmp)
		rand.append(1.0)
		rand.sort()

		# 2. 分别以 （a1-a0） (a2-a1) ... (an-an-1) 的概率，从本层中，迭代找出所有的缺失值。
		#   每个缺失 估算 缺失值。概率按照当前分布来。
		p =[]
		for ind in range(len(rand)):
			if ind > 0:
				p.append( rand[ind]-rand[ind-1])

		while has_non(sets[layer]):
			tt = [ key for key in sets[layer]]
			# print "len(tt)="
			# print len(tt)
			pkeys = np.random.choice(tt, 1, p )	# 以概率抽样这个层,抽出1个key
			pkey =pkeys[0]
			# 如果有缺失值，就根据概率补上一个值
			if sets[layer][pkey] == '':
				value = pad(sets[layer])
				# print '替换行和替换值分别 为 ：' 
				# print pkey,value
				sets[layer][pkey] = str(value)
		# print len(sets[layer])
		for key in sets[layer]:
			if sets[layer][key]=='':
				print "has none!!!!!!!!!!"
		
	# 结束后合并层，然后 根据 行号 排序
	all_data={}
	for layer in range(layers):
		all_data=mergedic(all_data,sets[layer])

	return all_data


def ave(dic):
	su = 0.0
	count = 0
	no_blank_num = 0
	for i in range(len(dic)):
		#print dic[i],
		if dic[i]:
			no_blank_num += 1
			#print dic[i]
			su += float(dic[i])
		count = count + 1
	if no_blank_num:
		average_value = su*1.0/no_blank_num
	else:
		print "ALL 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		average_value = 0
	average_value = float('%.2f' % average_value)
	print average_value
	for i in range(count):
		if dic[i] == '':
			dic[i] = average_value
	return dic


def mergedic(dic1,dic2):
	dictMerged=dic1.copy()
	dictMerged.update(dic2)
	return dictMerged

def fill_data(expand_feature_data,filled_data, imputation):
	if os.path.exists(expand_feature_data):
		filled_data = './ML/tmp/'+filled_data
		
		p=os.popen("awk -F '\t' 'END{ print NF}' " + expand_feature_data).readlines()
		NF = p[0].strip('\n')
		print NF
		# print col
		os.environ['col']=NF
		blank_col= os.popen("awk -F '\t' -v c=$col 'BEGIN{for(j=0;j<c;j++){ct[j]=0}}{for(i=1;i<=c;i++){if(ct[i-1]==0){ if($i==\"\"){ct[i-1]=1}}}} END{ for(t in ct){if(ct[t]==1){print t;} }}' " + expand_feature_data).readlines()
		blank=[int(i.strip('\n')) for i in blank_col]
		print blank
		#[114, 120, 122, 123, 124, 125, 126, 127, 139]

		for co in blank:
			dic_p={}
			row = 0
			blank_row = 0
			for rec in open(expand_feature_data):
				list_rec = rec.strip('\n').split('\t')
				dic_p[row]=list_rec[co]
				row = row +1
				if list_rec[0] == '':
					blank_row = blank_row+1
			
			layer_ini = int(math.ceil(float(row)/float(200)))
			# print layer_ini
			print "column == "
			print co
			print "--------------------------------------------"
			newf = open(str(co + 1) + '.txt', 'w')
			if imputation == "hot_deck":
				all_data = hot_deck(dic_p,layer_ini)
				for i in all_data:
					newf.write(all_data[i] + '\n')
			elif imputation == "average":
				all_data = ave(dic_p)
				for i in range(len(all_data)):
					newf.write(str(all_data[i]) + '\n')
			newf.close()
		sshf = open('tmp.sh','w')
		ssh = 'paste '
		NF=int(NF)
		for it in range(NF):
			if it not in blank:
				it = str(it)
				sshf.write("blank="+it+"\n")
				sshf.write("awk -F \'\\t\' -v i=\"$blank\" \'{print $(i+1)}\' "+expand_feature_data+" > "+  str(int(it)+1)+'.txt' +"\n")
		blank = sorted(blank)
		for co in range(NF):
			ssh = ssh + str(co+1)+'.txt'+' '
		# print ssh
		sshf.write(ssh+" > "+filled_data +'\n')
		sshrm = 'rm '
		for co in range(NF):
			sshrm = sshrm + str(co+1)+'.txt'+' '
		sshf.write(sshrm +'\n')
		sshf.close()
		os.system('sh tmp.sh')
		os.system('rm tmp.sh')
		return filled_data
		
	else:
		print ('ERROR: Wrong expanded file path!')
		return 'No_exp_fea_data'
		sys.exit()
# def draw_distribution(col):
#	 sift=[]
#	 for c in col:
#		 for line in open(str(c)+'.txt'):
#			 line = line.strip('\n')
#			 tmp = line.split('\t')
#	 #-----------------------------------------------------------------------
#			 fea = tmp[0]
#			 t = float(fea)
#			 sift.append( float(t))
#		 n, bins, patches = plt.hist(np.array(sift), bins=50, normed=1, facecolor='blue', alpha=0.5)
#		 plt.subplots_adjust(left=0.15)  
#		 plt.savefig(str(c)+'.jpg')
#		 plt.show()
