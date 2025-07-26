### **Preparation before analysis**

```bash
ls *.fq.gz >id.txt
cat id.txt | while read id; do echo $id | cut -d "." -f 1; done >id2.txt
cat id2.txt | uniq > id3.txt

#!/bin/bash
mkdir 00rawdata
mkdir 01cleandata
mkdir 02mapdata
mkdir 03calldata
mkdir 04calldata_anno
mkdir 05comparedata
mkdir 06comparedata_anno
mkdir 07methrate
```

### Aligned to human reference genome

```bash
#!/bin/bash
source /home/data/miniconda3/etc/profile.d/conda.sh
conda activate m5c

DIR=$(cd $(dirname $0); pwd)
echo $DIR
softwaredir=/home/data/liujiejie/software

mkdir $DIR/UID/fastqc
mkdir $DIR/UID/multiqc
ID=bis_seq_cc
sample=$DIR/UID/id.txt
index=/home/data/database/genome/Ensembl/Homo_sapiens/GRCh38_release102/BSGsgenomeIDX
gtf=/home/data/database/genome/Ensembl/Homo_sapiens/GRCh38_release102/Homo_sapiens.GRCh38.102.chr.gtf
genome=/home/data/database/genome/Ensembl/Homo_sapiens/GRCh38_release102/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa

fastqc $DIR/UID/*.fastq.gz -o $DIR/UID/fastqc -t 10 
multiqc $DIR/UID/fastqc/*.zip -o $DIR/UID/multiqc -n clean

cat $sample | while read id
do
mkdir $DIR/02mapdata/${id}
meRanGs align -o $DIR/02mapdata/${id} -f $DIR/UID/${id}.dedup.R1.fastq.gz -r $DIR/UID/${id}.dedup.R2.fastq.gz -t 20  -S RNA-BSseq.sam -deleteSAM -id $index -star_outFilterMultimapNmax 10  -star_outFilterMatchNminOverLread 0 -star_outFilterScoreMinOverLread 0 -star_outFilterMatchNmin 0 -star_outFilterMismatchNmax 2

htseq-count -f bam $DIR/02mapdata/${id}/RNA-BSseq_sorted.bam -i gene_id $gtf > $DIR/05comparedata/${id}.geneid.txt

meRanCall -p 32 -o $DIR/03calldata/${id}.result.txt -bam $DIR/02mapdata/${id}/RNA-BSseq_sorted.bam -f $genome -mBQ 0 -cr 0.99 -fs5 6 -rl 75 -sc 10 -ei 0.1 -fdr 0.01 -tref -bed63 -np -md 20 -mcov 20 -C_cutoff 3 -mr 0 -snr 0.9 >/dev/null 2>&1

meRanAnnotate -p 8 -t $DIR/03calldata/${id}.result.txt -f 'gene|transcript' -g $gtf -gtf -o $DIR/04calldata_anno/${id}_Annotated.txt

done
conda deactivate
```

### Calculation of methylation rate

```bash
meRanGh mkbsidx \
-t 4 \
-fa luci.fa \
-id ./BSGhlucigenomeIDX
```

```bash
#!/bin/bash
source /home/bio/biosoft/miniconda3/etc/profile.d/conda.sh
conda activate m5c

DIR=$(cd $(dirname $0); pwd)
echo $DIR
softwaredir=/home/data/software

ID=bis_seq_cc
sample=$DIR/UID/id.txt
index=/home/data/database/genome/Ensembl/Homo_sapiens/GRCh38_release102/luci_BS
genome=/home/data/database/genome/Ensembl/Homo_sapiens/GRCh38_release102/Luci.fa

cat $sample | while read id
do
mkdir $DIR/07methrate/${id}
meRanGh align -o $DIR/07methrate/${id} -f $DIR/UID/${id}.dedup.R1.fastq.gz -r $DIR/UID/${id}.dedup.R2.fastq.gz -t 12 -S RNA-BSseq.sam -id $index -bg -hisat2_ignore-quals

luc_bam=$DIR/07methrate/${id}/RNA-BSseq_sorted.bam
meRanCall -p 32 -fs5 6 -fs3 0 -rs5 0 -rs3 0 -s $luc_bam -f $genome -rl 150 -ccr -mBQ 0 -tref -cSeqID Luciferase > $DIR/07methrate/${id}.txt
done
conda deactivate
```



