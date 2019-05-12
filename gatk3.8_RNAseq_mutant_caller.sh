#/bin/bash
##software path
picard=/home/ggj/jiaqiLi/whole_genome_seq/software/GenomeAnalysisTK-3.8-1/picard.jar
gatk=/home/ggj/jiaqiLi/whole_genome_seq/software/GenomeAnalysisTK-3.8-1/GenomeAnalysisTK.jar
snpeff=/home/ggj/jiaqiLi/whole_genome_seq/software/snpEff/snpEff.jar
##reference
reference=/home/ggj/jiaqiLi/bulk-RNA-seq/tea_proj/workspace/reference/tea_genome.fa
ref_vcf=/home/ggj/jiaqiLi/bulk-RNA-seq/tea_proj/workspace/reference/tea_variation.vcf

###picard Add read groups, sort, mark duplicates, and create index
java -jar $picard AddOrReplaceReadGroups \
	I=Aligned.sortedByCoord.out.bam O=Aligned.sortedByCoord.RG_added_sorted.bam \
	SO=coordinate RGID=T1 RGLB=rna RGPL=illumina RGPU=hiseq RGSM=T1

java -jar $picard MarkDuplicates \
	I=Aligned.sortedByCoord.RG_added_sorted.bam O=dedup.bam \
	CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=dedup.metrics

###gatk
##mapq
java -jar $gatk -T SplitNCigarReads \
	-R $reference \
	-I dedup.bam -o dedup_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
##IndelRealign
#java -jar $gatk -T RealignerTargetCreator \
#	-R $reference \
#	-I dedup_split.bam -o dedup_realign_interval.list 

#java -jar $gatk -T IndelRealigner \
#	-R $reference \
#	-I dedup_split.bam -o dedup_realign.bam -targetIntervals dedup_realign_interval.list
##BQSR
#java -jar $gatk -T BaseRecalibrator \
#	-R $reference \
#	--knownSites $ref_vcf \
#	-I dedup_realign.bam \
#	-o recal_data.table

#java -jar $gatk -T PrintReads \
#	-R $reference \
#	-I dedup_realign.bam -BQSR recal_data.table -o dedup_realign_BQSR.bam

#mutant caller
java -jar $gatk -T HaplotypeCaller \
	-R $reference \
	-I dedup_split.bam \
	-dontUseSoftClippedBases \
	-stand_call_conf 20.0 \
	-o dedup_realign_BQSR.vcf
java -jar $gatk \
	-T VariantFiltration \
	-R $reference \
	-V dedup_realign_BQSR.vcf \
	-window 35 \
	-cluster 3 \
	-filterName FS -filter "FS > 30.0" \
	-filterName QD -filter "QD < 2.0" \
	-o dedup_realign_BQSR_filtered.vcf

java -Xmx4G -jar $snpeff T001 -i vcf \
	$out_dir/${speice_name}_${i}_mutect2.vcf.gz > $out_dir/${speice_name}_${i}_mutect2.snpeff.vcf \