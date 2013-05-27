#!/bin/bash
#$ -S /bin/bash
#$ -cwd
bwa=/home/usuyama/genomon/bin/bwa-0.5.10/bwa
ref=/home/usuyama/genomon/ref/hg19_bwa-0.5.10/hg19.fasta
samtools view -Sb $1.sam > $1.bam
samtools sort $1.bam $1.sort
samtools index $1.sort.bam
range=`samtools view $1.sort.bam | head -n1 | ruby -ane 'puts "#{$F[2]}:#{$F[3].to_i}-#{$F[3].to_i + 800}"'`
wd=`pwd`
base=`pwd | ruby -ane 'puts $F[0].split("/")[-3..-1].join("_")'`
cat << EOT > igv.batch
snapshotDirectory tmp/igvss
new
load $wd/tumor.sort.bam
load $wd/normal.sort.bam
genome hg19
goto $range
collapse
snapshot ${base}_${range}.png
EOT
