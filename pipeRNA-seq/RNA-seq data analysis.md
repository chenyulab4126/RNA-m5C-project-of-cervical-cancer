```bash
#!/bin/bash
source /home/data/miniconda3/etc/profile.d/conda.sh
conda activate rna

DIR=$(cd $(dirname $0); pwd)
echo $DIR
softwaredir=/home/data/software

ID=rna_seq
sample=$DIR/01cleandata/id3.txt
index=/home/data/database/genome/Ensembl/Homo_sapiens/GRCh38_release102/star_index_2.7.11b
gtf=/home/data/database/genome/Ensembl/Homo_sapiens/GRCh38_release102/Homo_sapiens.GRCh38.102.chr.gtf


cat $sample | while read id
do
STAR --runThreadN 18 \
    --quantMode TranscriptomeSAM \
    --genomeDir $index \
    --sjdbGTFfile $gtf \
    --sjdbScore 1 \
    --readFilesIn $DIR/01cleandata/${id}.dedup.R1.fastq.gz $DIR/01cleandata/${id}.dedup.R2.fastq.gz \
    --outSAMunmapped Within \
    --outFilterType BySJout \
    --outSAMattributes NH HI AS NM MD \
    --outFilterMultimapNmax 10 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outReadsUnmapped Fastx \
    --outFileNamePrefix $DIR/02map/${id}_\
    --readFilesCommand gunzip -c

samtools sort  -o $DIR/02map/${id}.sort.bam $DIR/02map/${id}_Aligned.out.sam -@ 2 -m 2G
samtools index $DIR/02map/${id}.sort.bam
rm $DIR/02map/${id}_Aligned.out.sam
rm $DIR/02map/${id}_Aligned.toTranscriptome.out.bam
gzip -c $DIR/02map/${id}_Unmapped.out.mate1 > $DIR/02map/${id}.r1.unmap.fastq.gz
gzip -c $DIR/02map/${id}_Unmapped.out.mate2 > $DIR/02map/${id}.r2.unmap.fastq.gz
rm $DIR/02map/${id}_Unmapped.out.mate1
rm $DIR/02map/${id}_Unmapped.out.mate2
rm -rf $DIR/02map/${id}__STARgenome
done

#featureCount
featureCounts -T 5 -p -t exon -g gene_name -a $gtf -o $DIR/03count/all.id.txt $DIR/02map/*.sort.bam

conda deactivate
```

