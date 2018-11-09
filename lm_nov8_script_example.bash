#!/bin/bash
#$ -pe smp 2
#$ -l h_vmem=30G
#$ -e "error.txt"
#$ -o "output.txt"
#$ -cwd

set -e 

fastq=$1

name_exp=$2
echo $name_exp

full_fastq="`readlink -f $fastq`"


#TRIM ADAPTER

zcat $full_fastq | /netapp/home/lcalviel/miniconda2/bin/cutadapt -m 18 -a AGATCGGAAGAGC -o trim_cut.fastq - 1> cutadapt.log

# MAP TO CONTAMINANTS/SMALL RNAs

less trim_cut.fastq  | /netapp/home/lcalviel/software/bowtie2-2.3.4.1-linux-x86_64/bowtie2 -q --local -p 2 -x /netapp/home/lcalviel/Annotation/bowtie2_contam/rRNA --al trim_cut_rRNA.fastq --un trim_cut_norRNA.fastq -U - | /netapp/home/lcalviel/software/bin/samtools view -bS - | /netapp/home/lcalviel/software/bin/samtools sort  -m 20G  - -o $name_exp"_rRNA_bowtie.bam"

/netapp/home/lcalviel/scripts/bef_sf_work/make_5p_bw.R $name_exp"_rRNA_bowtie.bam"
/netapp/home/lcalviel/software/bowtie2-2.3.4.1-linux-x86_64/bowtie2 -q  --local -p 2 -x /netapp/home/lcalviel/Annotation/bowtie2_contam/sno_miRNA --al trim_cut_sno_miRNA.fastq --un trim_cut_norsnomiRNA.fastq -U trim_cut_norRNA.fastq | /netapp/home/lcalviel/software/bin/samtools view -bS - | /netapp/home/lcalviel/software/bin/samtools sort  -m 20G  - -o $name_exp"_snomiRNA_bowtie.bam"

/netapp/home/lcalviel/software/bowtie2-2.3.4.1-linux-x86_64/bowtie2 -q  --local -p 2 -x /netapp/home/lcalviel/Annotation/bowtie2_contam/tRNA --al trim_cut_tRNA.fastq --un "$name_exp"_clean.fastq -U trim_cut_norsnomiRNA.fastq | /netapp/home/lcalviel/software/bin/samtools view -bS - | /netapp/home/lcalviel/software/bin/samtools sort  -m 20G  - -o $name_exp"_tRNA_bowtie.bam"
full_path_norRNA="`readlink -f "$name_exp"_clean.fastq`"

rm -f trim_cut_sno_miRNA.fastq trim_cut_norsnomiRNA.fastq trim_cut_norRNA.fastq trim_cut_tRNA.fastq trim_cut_rRNA.fastq

ls *.fastq | xargs gzip

full_path_norRNA="`readlink -f "$name_exp"_clean.fastq.gz`"


mkdir "starmapp_"$name_exp/

cd "starmapp_"$name_exp/

#MAP READS TO GENOME AND SPLICE JUNCTIONS

/netapp/home/lcalviel/software/STAR-2.6.0a/bin/Linux_x86_64/STAR --genomeDir /netapp/home/lcalviel/Annotation/Star_Index29/ --readFilesIn $full_path_norRNA --runThreadN 2 --alignEndsType EndToEnd --readFilesCommand zcat --outFilterMismatchNmax 2 --outFilterMultimapNmax 20 --chimScoreSeparation 10 --chimScoreMin 20 --chimSegmentMin 15 --outSAMattributes NH HI AS nM NM MD XS --outFilterIntronMotifs RemoveNoncanonicalUnannotated --alignSJoverhangMin 500 --outFileNamePrefix $name_exp"_" --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outWigType bedGraph --outWigNorm RPM --outSAMstrandField intronMotif --outSAMmultNmax 1 --outMultimapperOrder Random --limitBAMsortRAM 25000000000 --outFilterMismatchNoverLmax 0.1

#MAP LONG READS WITH NEW SPLICE DISCOVERY

/netapp/home/lcalviel/software/STAR-2.6.0a/bin/Linux_x86_64/STAR --genomeDir /netapp/home/lcalviel/Annotation/Star_Index99/ --readFilesIn $full_fastq_1 $full_fastq_2 --readFilesCommand zcat --alignSJoverhangMin 8 --runThreadN 2 --outFilterMultimapNmax 20 --peOverlapNbasesMin 10 --chimJunctionOverhangMin 20 --chimSegmentMin 20 --outSAMattributes NH HI AS nM NM MD XS --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFileNamePrefix $name_exp"_" --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --outWigNorm RPM --outSAMstrandField intronMotif --outSAMmultNmax 1 --outMultimapperOrder Random --limitBAMsortRAM 40000000000

bamfull="`readlink -f $name_exp"_"Aligned.sortedByCoord.out.bam`"

/netapp/home/lcalviel/software/bin/samtools index $bamfull


#COLLECT COVERAGE PROFILES AROUND DIFFERENT REGIONS

/netapp/home/lcalviel/scripts/bef_sf_work/SaTAnn_Ribo-seQC/analysis_qc_mod_jun18_2018_latest_modwynt.R BSgenome.Homo.sapiens.gencode25 /netapp/home/lcalviel/Annotation/genc25 $bamfull

#COUNT READS OVER REGIONS

/netapp/home/lcalviel/scripts/count_summary.R BSgenome.Homo.sapiens.gencode25 /netapp/home/lcalviel/Annotation/genc25 $bamfull single 1

#IDENTIFY READ MUTATIONS 

/netapp/home/lcalviel/scripts/mutations_BAM_RL.R $bamfull single 1

