# -*- coding: utf-8 -*-
'''
2018/8/20 Update:
   改变"fw2.write(new_id + '\t' + old_id + ';\n')" 为 "fw2.write(new_id + '\t' + old_id + '\n')"


Usage:
python modify_genblastg_out4muscle.py inputfile outfile
python /share/home/baishenglong/programs/my_script/gene_ann_majianchao/modify_genblastg_out4muscle.py \
    trans_euNOG_hmm.pep.all.fa.gff test.new.gff test.new.gff.id.list

'''
import sys,re

inputfile = sys.argv[1]
outfile = sys.argv[2]
out_idlist = sys.argv[3]

with open(outfile,'w') as fw, open(out_idlist,'w') as fw2:
    pat=re.compile('ID=.*?;')
    for x in open(inputfile):
        #跳过注释行
        if x.startswith('#'):
            pass
        elif x.split()[2] == 'transcript':
            x_replace = x.replace('transcript','mRNA')
            fw.write(x_replace)
            new_id = x.split('=')[1].split(';')[0]
            old_id = x.strip().split('Name=')[1]
            fw2.write(new_id + '\t' + old_id + '\n')
            #fw2.write(new_id + '\t' + old_id + ';\n')
        elif x.split()[2] == 'coding_exon':
            x_replace1 = x.replace('coding_exon', 'CDS')
            x_replace2=pat.sub('',x_replace1)
            fw.write(x_replace2.strip()+';\n')
        else:
            print(x)




