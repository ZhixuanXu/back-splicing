
#BWA + CIRI2 to identify circRNA
bwa mem –T 19 ref.fa read1.fq read2.fq 1> aln-pe.sam 2> aln-pe.log
perl CIRI2.pl -I test.sam -O outfile -F ref.fa -A ref.gtf
perl CIRI2.pl -I test.sam -O outfile2 -F ref.fa

#STAR + CIRCexplorer2 to identify circRNA
STAR --chimSegmentMin 10 --runThreadN 27 --genomeDir STAR_index --readFilesIn read_1.fastq read_2.fastq
CIRCexplorer2 parse -t STAR Chimeric.out.junction > CIRCexplorer2_parse.log
CIRCexplorer2 annotate -r ref.txt -g ref.fa -b back_spliced_junction.bed -o circularRNA_known.txt > CIRCexplorer2_annotate.log

#Tophat2 or STAR to identify linear splicing
tophat2 -p 100 -G annotation.gtf -o dir/ bowtie2_index read_1.fastq read_2.fastq
STAR --chimSegmentMin 10 --runThreadN 100 --genomeDir STAR_index --readFilesIn read_1.fastq read_2.fastq
