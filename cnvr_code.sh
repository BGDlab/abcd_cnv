#! /bin/bash

#################
# CNVR analysis #
#################
# https://github.com/CAG-CNV/ParseCNV2
# https://parsecnv.sourceforge.net/
module load perl
module load R/4.0
module load plink/1.9-20210416

pheno=(anxdep withdep somatic social thought attention rulebreak aggressive totprob totalcomp fluid flanker list cardsort reading pattern picture picvocab)

for i in ${pheno[@]}
do

mkdir /project/bbl-bgd/zhisha/abcd_cnv/cnvr_asso_results/${i}

cp /project/bbl-bgd/zhisha/abcd_cnv/pheno/${i}.txt /project/bbl-bgd/zhisha/abcd_cnv/cnvr_asso_results/${i}/${i}.txt

cd /project/bbl-bgd/zhisha/abcd_cnv/cnvr_asso_results/${i}/
perl /project/bbl-bgd/software/ParseCNV2/ParseCNV2_RevisedABCD.pl \
-i /project/bbl-bgd/zhisha/abcd_cnv/cnv_data/cnv_input.rawcnv \
-q ${i}.txt \
-b hg19 \
-batch /project/bbl-bgd/zhisha/abcd_cnv/cnv_data/batch_id.txt \
-stat linear -no_freq \
-m 1 \
-o ${i}

done



