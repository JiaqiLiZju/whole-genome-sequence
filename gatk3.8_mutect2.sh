#!/bin/bash
###set parameter
speice_name=$1

###set env
out_dir=mutect
##define software path/jar file
gatk=/home/ggj/jiaqiLi/whole_genome_seq/software/GenomeAnalysisTK-3.8-1/GenomeAnalysisTK.jar
snpeff=/home/ggj/jiaqiLi/whole_genome_seq/software/snpEff/snpEff.jar
#ref_file
reference=/home/ggj/jiaqiLi/whole_genome_seq/zebrafish/reference_byLI/Danio_rerio.GRCz10.dna.toplevel.fa
knownSites=/home/ggj/jiaqiLi/whole_genome_seq/zebrafish/reference_byLI/danio_rerio.vcf
chromosome=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25)
#thread_num
thread_num=8

###work flow
##mkdir
if [ ! -d "$out_dir"  ]; then
    mkdir $out_dir
fi

##Mutect2
#use LINUX FIFO to Concurrent
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
    time java -jar $gatk -T MuTect2 \
	    -R $reference \
		--dbsnp $knownSites \
	    --output_mode EMIT_VARIANTS_ONLY \
		-I:normal tail8/out/bwa_v11/zebrafish_tail8.sorted.markdup.realign.BQSR.bam \
		-I:tumor tumor8/out/bwa_v11/zebrafish_tumor8.sorted.markdup.realign.BQSR.bam \
		-o $out_dir/${speice_name}_${i}_mutect2.vcf.gz \
		-L $i
	time java -Xmx4G -jar $snpeff GRCz10.86 -i vcf \
		$out_dir/${speice_name}_${i}_mutect2.vcf.gz > $out_dir/${speice_name}_${i}_mutect2.snpeff.vcf \
	&& mv snpEff_genes.txt snpEff_genes_${i}.txt && mv snpEff_summary.html snpEff_summary_${i}.html
	sleep 1
	echo "" >&5
	#任务执行完后在fd5中写入一个占位符，以保证这个线程执行完后，线程继续保持占位
	} &
done
wait                #等待父进程的子进程都执行结束后再结束父进程          
exec 5>&-           #关闭fd5的管道

