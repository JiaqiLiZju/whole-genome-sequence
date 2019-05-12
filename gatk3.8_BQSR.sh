#!/bin/bash
###set parameter
fq1=$1
fq2=$2
speice_name=$3

###set env
out_dir=out
##define software path/jar file
trimmomatic=~/jiaqiLi/fqQC/software/Trimmomatic-0.38/trimmomatic-0.38.jar
bwa=/home/ggj/jiaqiLi/whole_genome_seq/software/bwa-0.7.15/bwa
samtools=/home/ggj/jiaqiLi/whole_genome_seq/software/samtools-1.3.1/samtools
picard=/home/ggj/jiaqiLi/whole_genome_seq/software/GenomeAnalysisTK-3.8-1/picard.jar
gatk=/home/ggj/jiaqiLi/whole_genome_seq/software/GenomeAnalysisTK-3.8-1/GenomeAnalysisTK.jar
snpeff=/home/ggj/jiaqiLi/whole_genome_seq/software/snpEff/snpEff.jar
#ref_file
reference=/home/ggj/jiaqiLi/whole_genome_seq/zebrafish/reference_byLI/Danio_rerio.GRCz10.dna.toplevel.fa
knownSites=/home/ggj/jiaqiLi/whole_genome_seq/zebrafish/reference_byLI/danio_rerio.vcf
chromosome=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 )
#sequence information
RGID=L006
#thread_num
thread_num=8

###prepare ref index
time $bwa index $reference && echo -e "\n\nbwa index done\n\n"
time $samtools faidx $reference && echo -e "\n\nsamtools index done\n\n"
time $samtools dict \
	-R $reference -o ${refrence}.dict \
	&& echo -e "\n\ngatk dict done\n\n"

#prepare for human variation
#$gatk IndexFeatureFile --Feature-file $GATK_bundle/hapmap_3.3.hg38.vcf
#$gatk IndexFeatureFile --Feature-file $GATK_bundle/1000G_omni2.5.hg38.vcf
#$gatk IndexFeatureFile --Feature-file $GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf
#$gatk IndexFeatureFile --Feature-file $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf
#$gatk IndexFeatureFile --Feature-file $GATK_bundle/dbsnp_146.hg38.vcf

###work flow
##mkdir
if [ ! -d "$out_dir"  ]; then
    mkdir $out_dir
fi

##fast quality control
echo "quality control begin"
if [ ! -d "$out_dir/cleanfq"  ]; then
    mkdir $out_dir/cleanfq
fi
#fastqc of original fq
time fastqc $fq1 $fq2 -o $out_dir

#trimmomatic set as illumina
#put the adapter file in the work dir&&change the shell
time java -jar ${trimmomatic} PE -phred33 -trimlog log_trim \
	$fq1 $fq2 \
	$out_dir/cleanfq/$speice_name.pair.1.fq $out_dir/cleanfq/$speice_name.trimUnpair.1.fq \
	$out_dir/cleanfq/$speice_name.pair.2.fq $out_dir/cleanfq/$speice_name.trimUnpair.2.fq \
	ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
	SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:50 \
	&& echo -e "\n\ntrimmomatic done\n\n"

#fastqc of trimmed fq
time fastqc $out_dir/cleanfq/$speice_name.pair.1.fq $out_dir/cleanfq/$speice_name.pair.2.fq -o $out_dir/cleanfq \
	&& echo -e "\n\nquality control done\n\n"

##bwa align with reference and sort to bam file
if [ ! -d "$out_dir/bwa_v11"  ]; then
    mkdir $out_dir/bwa_v11
fi

#align with ref_path genome
time $bwa mem -t $thread_num -M -R "@RG\tID:${RGID}\tPL:ILLUMINA\tSM:${speice_name}" \
	$reference $out_dir/cleanfq/$speice_name.pair.1.fq $out_dir/cleanfq/$speice_name.pair.2.fq \
	| $samtools view -b -o $out_dir/bwa_v11/$speice_name.bam \
	&& echo -e "\n\nbwa mem done\n\n"
time $samtools sort -@ $thread_num -m 8G -O bam -o $out_dir/bwa_v11/$speice_name.sorted.bam $out_dir/bwa_v11/$speice_name.bam \
	&& echo -e "\n\nbwa sort done\n\n"
	
#mark duplicate
time java -jar ${picard} MarkDuplicates \
	I=$out_dir/bwa_v11/$speice_name.sorted.bam O=$out_dir/bwa_v11/$speice_name.sorted.markdup.bam \
	M=$out_dir/bwa_v11/$speice_name.sorted.markdup_metrics.txt \
	CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT \
	&& echo -e "\n\n markdup done\n\n"
	
time $samtools index $out_dir/bwa_v11/$speice_name.sorted.markdup.bam && echo -e "\n\nsamtools index done\n\n"

##apply realiagnment
time java -jar $gatk -T RealignerTargetCreator \
	-R $reference \
	-known $knownSites \
	-I $out_dir/bwa_v11/$speice_name.sorted.markdup.bam \
	-o $out_dir/bwa_v11/realign_interval.list && echo -e "\n\nrealignment done\n\n"

time java -jar $gatk -T IndelRealigner \
	-R $reference \
	-known $knownSites \
	-targetIntervals $out_dir/bwa_v11/realign_interval.list \
	-I $out_dir/bwa_v11/$speice_name.sorted.markdup.bam \
	-o $out_dir/bwa_v11/$speice_name.sorted.markdup.realign.bam \
	&& echo -e "\n\nrealignment bam done\n\n"

##BQSR
time java -jar $gatk -T BaseRecalibrator \
	-R $reference \
	--knownSites $knownSites \
	-I $out_dir/bwa_v11/$speice_name.sorted.markdup.realign.bam -o $out_dir/bwa_v11/recal_data.table \
	&& echo -e "\n\nrecal done\n\n"

time java -jar $gatk -T PrintReads \
	-R $reference \
	-BQSR $out_dir/bwa_v11/recal_data.table \
	-I $out_dir/bwa_v11/$speice_name.sorted.markdup.realign.bam -o $out_dir/bwa_v11/$speice_name.sorted.markdup.realign.BQSR.bam \
	&& echo -e "\n\nBQSR done\n\n"