# Usage:
#     sh /share/home/baishenglong/programs/my_script/gene_ann_majianchao/github_genepredict/homolog/little_gene_homolog_annotation.sh <workdir> <genome_file> <homolog_file>

#genome_file="/share/home/baishenglong/database/liuyang_mosses_genome_BGI/moss_genome/10_moss_assembly/Anomodon_attenuatus.fasta"
#workdir="/share/home/baishenglong/huangjinling/01.HET_search/09.ten_mosses_genome"
#homolog_file="$workdir/final_het.fasta"

workdir=$1
genome_file=$2
homolog_file=$3

genome_name=$(echo $genome_file|rev |cut -d '/' -f 1|rev|sed 's/.fasta//g')

mkdir -p $workdir/$genome_name
outdir=$workdir/$genome_name
ref=$genome_file

mkdir $outdir/bin $outdir/input $outdir/result;

ln -s $ref $outdir/input/ref.fa
cp ~/programs/my_script/gene_ann_majianchao/github_genepredict/homolog/genblastg_cut_query.pl $outdir/bin/genblastg_cut_query.pl

ln -s $homolog_file $outdir/input/homology.fa

cd $outdir/result
ln -s /share/home/baishenglong/programs/genblast/genBlast_v138_linux_x86_64/blastall ./
ln -s /share/home/baishenglong/programs/genblast/genBlast_v138_linux_x86_64/formatdb ./
ln -s /share/home/baishenglong/programs/genblast/genBlast_v138_linux_x86_64/alignscore.txt ./

# Declare environment variables of genblastg to avoid errors caused by inconsistent blast versions
export PATH=/share/home/baishenglong/programs/genblast/genBlast_v138_linux_x86_64/:$PATH
# Format database
$outdir/result/formatdb -i $outdir/input/ref.fa -p F -o T

# Run genblastg
## cut query file into small subfiles
perl /share/home/baishenglong/programs/my_script/gene_ann_majianchao/github_genepredict/homolog/common_bin/fastaDeal.pl -cutf 10 ../input/homology.fa -outdir ./;

for x in `ls homology.fa.cut`;do
echo "/share/home/baishenglong/programs/genblast/genBlast_v138_linux_x86_64/genblast_v138_linux_x86_64 -f F -e 1e-2 -p genblastg -g T -h 0 -j 0 -v 1 -s 0 -r 5 -norepair -gff -q ./homology.fa.cut/$x -t ../input/ref.fa  -o ./homology.fa.cut/$x.genblastg.out"
done > homology.fa.genblastg.shell

sh homology.fa.genblastg.shell

## merge genblastg result
#rm core* # remove core file
cat homology.fa.cut/*gff > homology.fa.genblastg.gff

## Output of Homologous Annotation Results
grep transcript homology.fa.genblastg.gff |wc -l
 
## Modify the format of genblastg output GFF file to adapt to the subsequent muscle comparison script and extract PEP sequence
python /share/home/baishenglong/programs/my_script/gene_ann_majianchao/modify_genblastg_out4muscle.py homology.fa.genblastg.gff homology.fa.genblastg.new.gff homology.fa.genblastg.new.gff.id.list
/share/home/baishenglong/programs/gffread-0.9.12/gffread/gffread homology.fa.genblastg.new.gff -g $outdir/input/ref.fa -x homology.fa.genblastg.new.gff.cds -y homology.fa.genblastg.new.gff.pep

## Run muscle

perl /share/home/baishenglong/programs/my_script/gene_ann_majianchao/github_genepredict/homolog/genBlastA/run_muscle.pl homology.fa.genblastg.new.gff.id.list homology.fa.genblastg.new.gff.pep ../input/homology.fa
sh homology.fa.genblastg.new.gff.id.list.muscle.sh
#perl ../bin/genblastg_cut_query.pl --verbose --pep_cut 1000 --run bsub --queue normal --step 6 ../input/homology.fa $outdir/input/ref.fa 1>out.log 2>err.log

## Screening genes based on muscle results
for i in homology.fa.genblastg.new.gff.id.list.cut*/*; do 
    for j in $i/*;do 
        for k in $j/*.muscle;do 
            perl /share/home/baishenglong/programs/my_script/gene_ann_majianchao/github_genepredict/homolog/muscle_identity.pl $k;
        done;
    done;
done > homology.fa.genblastg.new.gff.ident.list;
awk '$3 >= 30 && $6 >= 50 && $7 >= 10 && $8 >= 25' homology.fa.genblastg.new.gff.ident.list > homology.fa.genblastg.new.gff.ident.list.filter;
perl /share/home/baishenglong/programs/my_script/gene_ann_majianchao/github_genepredict/homolog/algnRat2identity_gff.pl homology.fa.genblastg.new.gff.ident.list.filter homology.fa.genblastg.new.gff > homology.fa.genblastg.new.gff.ident

## Removal of redundancy caused by repeated prediction of the same section due to query similarity
python /share/home/baishenglong/programs/my_script/gene_ann_majianchao/filter_genewise_gff_redundancy.py homology.fa.genblastg.new.gff.ident 0 homology.fa.genblastg.new.gff.ident.filter

# Extraction of filtered gene sets
/share/home/baishenglong/programs/gffread-0.9.12/gffread/gffread homology.fa.genblastg.new.gff.ident.filter -g $outdir/input/ref.fa -x filter.cds -y filter.pep
# Modify Sequence Name
sed "s/>/>${genome_name}_/g" filter.pep > filter.rename.pep
