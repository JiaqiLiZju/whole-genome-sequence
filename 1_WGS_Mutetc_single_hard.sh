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
gatk=/home/ggj/jiaqiLi/whole_genome_seq/software/gatk-4.0.5.1/gatk
#VEP=/your_path_to/ensembl-vep/vep
#ref_file
reference=/home/ggj/jiaqiLi/whole_genome_seq/zebrafish/reference/zebrafish_genomic.fna.gz
chromosome=( 	CM002885.2, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, \
				chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, \
				chr20, chr21, chr22, chr23, chr24, chr25 )
#export GATK_bundle=GATK_bundle_PATH
#sequence information
RGID=L005
#thread_num
thread_num=16

###prepare ref index
time $bwa index $reference && echo "bwa index done"
time $samtools faidx $reference && echo "samtools index done"
time $gatk CreateSequenceDictionary \
	-R $reference -O ${refrence}.dict \
	&& echo "gatk dict done"

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
	&& echo "trimmomatic done"

#fastqc of trimmed fq
time fastqc $out_dir/cleanfq/$speice_name.pair.1.fq $out_dir/cleanfq/$speice_name.pair.2.fq -o $out_dir/cleanfq \
	&& echo "quality control done"

##bwa align with reference and sort to bam file
if [ ! -d "$out_dir/bwa"  ]; then
    mkdir $out_dir/bwa
fi
#align with ref_path genome
time $bwa mem -t $thread_num -M -R "@RG\tID:$RGID\tPL:ILLUMINA\tSM:$speice_name" \
	$reference $out_dir/cleanfq/$speice_name.pair.1.fq $out_dir/cleanfq/$speice_name.pair.2.fq \
	| $samtools view -b -o $out_dir/bwa/$speice_name.bam \
	&& echo "bwa mem done"
time $samtools sort -@ $thread_num -m 4G -O bam -o $out_dir/bwa/$speice_name.sorted.bam $out_dir/bwa/$speice_name.bam \
	&& echo "bwa sort done"
	
##gatk
#mark duplicate
time $gatk MarkDuplicates \
	-I $out_dir/bwa/$speice_name.sorted.bam \
	-M $out_dir/bwa/$speice_name.sorted.markdup_metrics.txt \
	-O $out_dir/bwa/$speice_name.sorted.markdup.bam \
	&& echo "** markdup done **"

##apply realiagnment

##gatk BQSR
time $gatk BaseRecalibrator \
    -R $ref_path/$ref_file \
    -I $out_dir/bwa/$speice_name.sorted.markdup.bam \
    --known-sites $GATK_bundle/1000G_phase1.indels.hg38.vcf \
	--known-sites $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf \
	--known-sites $GATK_bundle/dbsnp_146.indels.hg38.vcf \
	-O $out_dir/bwa/$speice_name.sorted.markdup.recal.table \
	&& echo "$speice_name.sorted.markdup.recal.table done"
	
time $gatk ApplyBQSR \
	--bqsr-recal-file $out_dir/bwa/$speice_name.sorted.markdup.recal.table \
	-R $ref_path/$ref_file \
	-I $out_dir/bwa/$speice_name.sorted.markdup.bam \
	-O $out_dir/bwa/$speice_name.sorted.markdup.BQSR.bam \
	&& echo "BQSR.bam done"
	
time $samtools sort -@ $thread_num -m 4G -O bam -o $out_dir/bwa/$speice_name.sorted.markdup.BQSR.bam \
	&& echo "bwa sort done"
	
time $samtools index $out_dir/bwa/$speice_name.sorted.markdup.BQSR.bam \
	&& echo "sorted bam index done"
	
rm $out_dir/bwa/$speice_name.bam
rm $out_dir/bwa/$speice_name.sorted.bam	

##gatk single hard sort
if [! -d "$out_dir/gatk" ]; then
	mkdir $out_dir/gatk_single_hard
fi

#use LINUX FIFO to Concurrent
#chromosome=( 	chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, \
				chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, \
				chr20, chr21, chr22, chr23, chr24, chr25 )
TMPFIFO=/tmp/$$.fifo			#声明管道名称，'$$'表示脚本当前运行的进程PID
mkfifo $TMPFIFO                #创建管道
exec 5<>${TMPFIFO}             #创建文件标示符“5”，调用进程内部执行一个可执行文件
rm -rf ${TMPFIFO}              #清除创建的管道文件
#为并发线程创建同样个数的占位
for((i=1;i<=$thread_num;i++))
do
   echo ;               
   #借用read命令一次读取一行的特性，使用一个echo默认输出一个换行符，来确保每一行只有一个线程占位；
done >&5                
#将占位信息写入管道
for i in ${chromosome[@]}; do
	read -u5
	#从文件描述符管道中，获取一个管道的线程占位然后开始执行操作；read中 -u 后面跟fd，表示从文件描述符中读入，该文件描述符可以是exec新开启的
	{
	time $gatk HaplotypeCaller \
		-R $reference \
		--emit-ref-confidence GVCF \
		-I $out_dir/bwa/${specie}.sorted.markdup.bam \
		-L $i \
		-O $out_dir/gatk_single_hard/${specie}.${i}.gvcf \
		&& echo "gvcf done" && \
	#call variety
	time $gatk GenotypeGVCFs \
		-R $reference \
		-V $out_dir/gatk_single_hard/${specie}.${i}.gvcf \
		-O $out_dir/gatk_single_hard/${specie}.${i}.vcf \
		&& echo "vcf done"
	sleep 1
	echo "" >&5
	#任务执行完后在fd5中写入一个占位符，以保证这个线程执行完后，线程继续保持占位
	} &
done
wait                #等待父进程的子进程都执行结束后再结束父进程          
exec 5>&-           #关闭fd5的管道
#merge chromosome file
merge_vcfs=""
for i in ${chromosome[@]}; do
	merge_vcfs=${merge_vcfs}" -I $out_dir/gatk_single_hard/${specie}.${i}.vcf \\"\n
done 
time $gatk merge_vcfs ${merge_vcfs} -O $out_dir/gatk_single_hard/${specie}.vcf \
	&& echo "merge vcfs done"

time bgzip -f $out_dir/gatk_single_hard/${specie}.vcf
time tabix -p vcf $out_dir/gatk_single_hard/${specie}.vcf.gz

##SNP INDEL select
# 使用SelectVariants，选出SNP
time $gatk SelectVariants \
    -select-type SNP \
    -V $out_dir/gatk_single_hard/${specie}.vcf.gz \
    -O $out_dir/gatk_single_hard/${specie}.snp.vcf.gz

# 为SNP作硬过滤
time $gatk VariantFiltration \
    -V $out_dir/gatk_single_hard/${specie}.snp.vcf.gz \
    --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "Filter" \
    -O $out_dir/gatk_single_hard/${specie}.snp.fiter.vcf.gz

# 使用SelectVariants，选出Indel
time $gatk SelectVariants \
    -select-type INDEL \
    -V $out_dir/gatk_single_hard/${specie}.vcf.gz \
    -O $out_dir/gatk_single_hard/${specie}.indel.vcf.gz

# 为Indel作过滤
time $gatk VariantFiltration \
    -V $out_dir/gatk_single_hard/${specie}.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "Filter" \
    -O $out_dir/gatk_single_hard/${specie}.fiter.indel.vcf.gz

# 重新合并过滤后的SNP和Indel
time $gatk MergeVcfs \
    -I $out_dir/gatk_single_hard/${specie}.snp.fiter.vcf.gz \
    -I $out_dir/gatk_single_hard/${specie}.fiter.indel.vcf.gz \
    -O $out_dir/gatk_single_hard/${specie}.filter.vcf.gz

# 删除无用中间文件
#rm -f 

##variety annotation
## VEP annotation
time $VEP --fasta $reference \
 --vcf --merged --fork 10 --hgvs --force_overwrite --everything \
   --offline --dir_cache /your_path_to/ensembl-vep/cachedir \
   -i $out_dir/population/$out_name.snp.VQSR.vcf.gz \
   -o $out_dir/population/$out_name.snp.VQSR.VEP.vcf.gz

time $VEP --fasta $reference \
 --vcf --merged --fork 10 --hgvs --force_overwrite --everything \
   --offline --dir_cache /your_path_to/ensembl-vep/cachedir \
   -i $out_dir/population/$out_name.indel.VQSR.vcf.gz \
   -o $out_dir/population/$out_name.indel.VQSR.VEP.vcf.gz