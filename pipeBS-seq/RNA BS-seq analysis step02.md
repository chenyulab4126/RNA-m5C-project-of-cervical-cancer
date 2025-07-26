### Differential m5C analysis of the two groups of samples

```bash
paste sample1.geneid.txt sample2.geneid.txt sample3.geneid.txt sample4.geneid.txt | cut -f1,2,4,6,8 > Counts.txt

perl /home/data/software/estimateSizeFactors.pl -cf Counts.txt -cn sample1,sample2

#!/bin/bash
source /home/data/miniconda3/etc/profile.d/conda.sh
conda activate m5c

DIR=$(cd $(dirname $0); pwd)
echo $DIR

meRanCompare \
-fa $DIR/03calldata/sample1.result_FDR_0.01.txt, sample2.result_FDR_0.01.txt \
-fb $DIR/03calldata/sample3.result_FDR_0.01.txt, sample4.result_FDR_0.01.txt \
-na normal \
-nb cancer \
-sfa 0.9389,0.9773 \
-sfb 1.4595,0.8816 \
-sig 0.01 \
-fdr 0.05 
conda deactivate

#!/bin/bash
source /home/data/miniconda3/etc/profile.d/conda.sh
conda activate m5c

DIR=$(cd $(dirname $0); pwd)
echo $DIR

gtf=/home/data/database/genome/Ensembl/Homo_sapiens/GRCh38_release102/Homo_sapiens.GRCh38.102.chr.gtf

meRanAnnotate -p 8 -t $DIR/05comparedata/intersect_normal_cancer_FDR_0.05.txt -f 'gene|transcript' -g $gtf -gtf -o $DIR/06comparedata_anno/intersect_normal_cancer_FDR_0.05_Annotated.txt

meRanAnnotate -p 8 -t $DIR/05comparedata/intersect_normal_cancer.txt -f 'gene|transcript' -g $gtf -gtf -o $DIR/06comparedata_anno/intersect_normal_cancer_Annotated.txt

meRanAnnotate -p 8 -t $DIR/05comparedata/uniq_normal.txt -f 'gene|transcript' -g $gtf -gtf -o $DIR/06comparedata_anno/uniq_normal_Annotated.txt

meRanAnnotate -p 8 -t $DIR/05comparedata/uniq_cancer.txt -f 'gene|transcript' -g $gtf -gtf -o $DIR/06comparedata_anno/uniq_cancer_Annotated.txt
conda deactivate
```

