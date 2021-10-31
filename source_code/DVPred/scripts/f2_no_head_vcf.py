import os
# def no_head(vcffile,non_repeat_region_vcf):
def no_head(vcffile):
	# res1 = os.system("grep -v '#' "+vcffile+"> ./ML/tmp/no_head_1.vcf")
	# res2 = os.system("grep -v '#' "+non_repeat_region_vcf+" | awk -F '\t' '{print $NF}' > ./ML/tmp/no_head_2.vcf")
	# res3 = os.system('paste ./ML/tmp/no_head_1.vcf ./ML/tmp/no_head_2.vcf > ./ML/tmp/no_head.vcf')
	# os.system('rm ./ML/tmp/no_head_1.vcf ./ML/tmp/no_head_2.vcf')
	no_head_path=''
	res1 = os.system("grep -v '#' "+vcffile+"> ./ML/tmp/no_head.vcf")
	# if res1==0 and res2==0 and res3==0:
	if res1==0 :
		no_head_path = './ML/tmp/no_head.vcf'
	return no_head_path