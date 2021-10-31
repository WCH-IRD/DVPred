# encoding:utf-8  

# import pyExcelerator as exl
import csv
import os
import re
# def to_excel(fpath,list_feature,ML_model):  
#     # 读取参数路径文件  
#     lines=open(fpath,'r').readlines()
#     for i in list_feature:
#         if i =='SYMBOL':
#             list_feature[list_feature.index(i)] = 'SYMBOL&cDNA_pos'
#     # 创建workbook  
#     w = exl.Workbook()  
#     # 增加一个sheet页'Sheet1'  
#     ws = w.add_sheet('Sheet1')  
#     # 以'*'分割，获取每行数据  
#     arr_line = [item.strip('\n') for item in lines]
#     rows = len(arr_line)
#     cols = len(arr_line[0].split('\t'))

#     # 头部信息　概率　CHROM POS　ID　REF　ALT　QUAL FILTER
#     head_list = []
#     first_7_col_vcf = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER']
#     csqs = ['single transcript','multiple transcipts']
#     # 家系信息　['FORMAT','PEDIGREE',...]
#     head_list.extend(ML_model)
#     head_list.extend(first_7_col_vcf)
#     head_list.extend(list_feature)
#     head_list.extend(csqs)
#     for j in range(len(head_list)):
#         ws.write(0, j, head_list[j])  
#     # body信息
#     for i in range(1,rows):  
#         # 对行数据进行遍历，获取行数据元素元组
#         cell = arr_line[i].split('\t')
#         for j in range(len(cell)):  
#             # 写入数据  
#             ws.write(i, j, cell[j])  
#     fpath_excel=fpath.replace('vcf','xls')  
#     w.save(fpath_excel)  
#     return fpath_excel

def name_parse(name_file, head=0, scaled=False):
    name_list = []
    n = 0
    for name in open(name_file, 'r'):
        name = name.strip('\n')
        n += 1
        if head > 0:
	    if scaled:
		name = name + '_scaled'
            if n <= head:
                name_list.append(name)
        else:
            name_list.append(name)
    return name_list



#to_csv(input_vcf,output_vcf,csq_field,ML_model)
def to_csv(vcffile,fpath,csq_field,ML_model,name,old_label):
    fpath_csv=fpath.replace('vcf','csv')    
    csvf = open(fpath_csv,'w')
    csv_writer = csv.writer(csvf,dialect='excel')
    # 头部信息　概率　CHROM POS　ID　REF　ALT　QUAL FILTER
    head_list = []
    first_7_col_vcf = ['CHROM:POS','ID','REF','ALT','QUAL','FILTER']
    csqs = csq_field.split('|')
    csqs.extend(['INFO','FORMAT'])
    #csqs.extend(['INFO','multiple transcipts','FORMAT'])
    sample_id = os.popen('grep "#CHROM" ' + vcffile).readlines()[0].split('\t')[9:]
    if sample_id:
        sample_id[-1] = re.sub('\n$', '', sample_id[-1])
    csqs.extend(sample_id)
    # 家系信息　['FORMAT','PEDIGREE',...]
    head_list.extend(ML_model)
    head_list.extend(first_7_col_vcf)
    head_list.extend(csqs)
    csv_writer.writerow(head_list)

    # all csv
    csv_all_f = open(name,'w')
    csv_all_writer = csv.writer(csv_all_f, dialect='excel')
    info1 = name_parse('./ML/list/expand_feature.name.training')  # expand_feature_name
    info2 = name_parse('./ML/list/expand_feature.name.0', 299)  # ./ML/tmp/expand_feature.data.0 1~186
    info3 = name_parse('./ML/list/expand_feature.name.0', 299, scaled=True)  # ./ML/tmp/expand_feature.data.0 1~186
    info4 = name_parse('./lists/initial_features.list')
    info4.append('data_old')
    head_list_all = []
    head_list_all.extend(ML_model)
    head_list_all.extend(first_7_col_vcf)
    csqs = csq_field.split('|')
    head_list_all.extend(csqs)
    head_list_all.append('P/B')
    head_list_all.extend(info1)
    head_list_all.extend(info2)
    head_list_all.extend(info3)
    head_list_all.extend(info4)
    csv_all_writer.writerow(head_list_all)

    # body信息
    for line in open(fpath):
        line_list = line.strip('\n').split('\t')
        csv_writer.writerow(line_list)
    csvf.close()

    os.system('cut -f 8-  ./ML/tmp/initial_feature.data > ./ML/tmp/initial_feature.data.info')
    os.system('paste ' + fpath + ' all.txt.test.merge ./ML/tmp/initial_feature.data.info ' + old_label + '> all.vcf')  # fpath:CDGC_training_anno_HL_pred.vcf
    for tmp in open("all.vcf"):
	#tmp.strip('\n')
        #tmp = re.sub(' ','\t',tmp)
        tmp_list = tmp.strip('\n').split('\t')
        csv_all_writer.writerow(tmp_list)
    csv_all_f.close()

    return fpath_csv
