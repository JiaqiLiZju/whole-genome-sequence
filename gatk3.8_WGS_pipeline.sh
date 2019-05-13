#!/bin/bash
###set parameter
fq1=/share/home/guoguoji/RAWDATA/AML-genome-seq/Merged-lane/H_R1.fastq
fq2=/share/home/guoguoji/RAWDATA/AML-genome-seq/Merged-lane/H_R2.fastq
speice_name=AML

###set env
out_dir=out
##define software path/jar file
trimmomatic=/share/home/guoguoji/tools/Biology-tools-package/trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar
bwa=/share/home/guoguoji/tools/Biology-tools-package/bwa/bwa-0.7.15/bwa
samtools=/share/apps/samtools/1.8/bin/samtools
picard=/share/home/guoguoji/tools/picard-tools-1.119/picard-1.119.jar
gatk=/share/home/guoguoji/tools/Biology-tools-package/gatk/GenomeAnalysisTK-3.8-1/GenomeAnalysisTK.jar
snpeff=/share/home/guoguoji/Desktop/callsnptools/snpEff/snpEff.jar
#ref_file
reference=/share/home/guoguoji/tools/BWA_Reference_Human/Homo_sapiens.GRCh38.fa
knownSites1=/share/home/guoguoji/Desktop/callsnptools/human_reference/Mills_and_1000G_gold_standard.indels.hg38.vcf
knownSites2=/share/home/guoguoji/Desktop/callsnptools/human_reference/dbsnp_138.hg38.vcf.gz
#chromosome=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 )
#sequence information
RGID=L005
#thread_num
thread_num=8

###prepare ref index
#time $bwa index $reference && echo -e "\n\nbwa index done\n\n"
#time $samtools faidx $reference && echo -e "\n\nsamtools index done\n\n"
#time $samtools dict \
#	$reference -o ${refrence}.dict \
#	&& echo -e "\n\ngatk dict done\n\n"

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
#time fastqc $fq1 $fq2 -o $out_dir

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
if [ ! -d "$out_dir/bwa"  ]; then
    mkdir $out_dir/bwa
fi

#align with ref_path genome
time $bwa mem -t $thread_num -M -R "@RG\tID:${RGID}\tPL:ILLUMINA\tSM:${speice_name}" \
	$reference $out_dir/cleanfq/$speice_name.pair.1.fq $out_dir/cleanfq/$speice_name.pair.2.fq \
	| $samtools view -b -o $out_dir/bwa/$speice_name.bam \
	&& echo -e "\n\nbwa mem done\n\n"
time $samtools sort -@ $thread_num -m 8G -O bam -o $out_dir/bwa/$speice_name.sorted.bam $out_dir/bwa/$speice_name.bam \
	&& echo -e "\n\nbwa sort done\n\n"
	
#mark duplicate
time java -jar ${picard} MarkDuplicates \
	I=$out_dir/bwa/$speice_name.sorted.bam O=$out_dir/bwa/$speice_name.sorted.markdup.bam \
	M=$out_dir/bwa/$speice_name.sorted.markdup_metrics.txt \
	CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT \
	&& echo -e "\n\n markdup done\n\n"
	
time $samtools index $out_dir/bwa/$speice_name.sorted.markdup.bam && echo -e "\n\nsamtools index done\n\n"

##apply realiagnment
time java -jar $gatk -T RealignerTargetCreator \
	-R $reference \
	-known $knownSites1 \
	-known $knownSites2 \
	-I $out_dir/bwa/$speice_name.sorted.markdup.bam \
	-o $out_dir/bwa/realign_interval.list && echo -e "\n\nrealignment done\n\n"

time java -jar $gatk -T IndelRealigner \
	-R $reference \
	-known $knownSites1 \
	-known $knownSites2 \
	-targetIntervals $out_dir/bwa/realign_interval.list \
	-I $out_dir/bwa/$speice_name.sorted.markdup.bam \
	-o $out_dir/bwa/$speice_name.sorted.markdup.realign.bam \
	&& echo -e "\n\nrealignment bam done\n\n"

##BQSR
time java -jar $gatk -T BaseRecalibrator \
	-R $reference \
	-known $knownSites1 \
	-known $knownSites2 \
	-I $out_dir/bwa/$speice_name.sorted.markdup.realign.bam -o $out_dir/bwa/recal_data.table \
	&& echo -e "\n\nrecal done\n\n"

time java -jar $gatk -T PrintReads \
	-R $reference \
	-BQSR $out_dir/bwa/recal_data.table \
	-I $out_dir/bwa/$speice_name.sorted.markdup.realign.bam -o $out_dir/bwa/$speice_name.sorted.markdup.realign.BQSR.bam \
	&& echo -e "\n\nBQSR done\n\n"

##gatk single hard sort
if [! -d "$out_dir/gatk_single_hard" ]; then
	mkdir $out_dir/gatk_single_hard
fi

time java -jar $gatk -T HaplotypeCaller \
	-R $reference \
	--dbsnp $knownSites2 \
	-I dedup_split.bam \
	-dontUseSoftClippedBases \
	-stand_call_conf 20.0 \
	-o dedup_realign_BQSR.vcf
time java -jar $gatk \
	-T VariantFiltration \
	-R $reference \
	-V dedup_realign_BQSR.vcf \
	-window 35 \
	-cluster 3 \
	-filterName FS -filter "FS > 30.0" \
	-filterName QD -filter "QD < 2.0" \
	-o dedup_realign_BQSR_filtered.vcf

#annotate snv
time java -jar $snpeff GRCh38.92 -i vcf \
	$out_dir/gatk_single_hard/${specie}.filter.vcf.gz > $out_dir/gatk_single_hard/${specie}.filter.anpeff.vcf.gz