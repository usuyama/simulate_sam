#!/bin/bash
#$ -S /bin/bash
#$ -cwd
map=/home/usuyama/codes/hapmuc_sim2/to_bam.sh
ref=/home/usuyama/codes/hapmuc_sim2/simulated.20130527.fasta
sim=/home/usuyama/codes/hapmuc_sim2/sim_sam.rb
mutation_call=/home/usuyama/codes/hapmuc_sc/mutation_call_for_sim.sh
test_type=$1
tp=$2
depth=$3
count=$4
out_f=$5
a=$6
b=${7}
k=${8}

echo SGE_TASK_ID:$SGE_TASK_ID
echo SGE_TASK_FIRST:$SGE_TASK_FIRST
echo SGE_TASK_LAST:$SGE_TASK_LAST
echo SGE_TASK_STEPSIZE:$SGE_TASK_STEPSIZE

id=$SGE_TASK_ID # 1から
if test "$id" = "" ;then
  id=1
fi
let start_num=(id-1)*count+1
let end_num=id*count

gout=$out_f/${a}_${b}_${k}/$test_type/tp${tp}
mkdir -p $gout
for i in `seq $start_num $end_num`
do
  out=$gout/$i
  mkdir $out
  cd $out
  ruby $sim $test_type $depth $tp ./ 20 $a $b $k
  sh $map normal
  sh $map tumor
  sh $mutation_call $out/tumor.sort.bam $out/normal.sort.bam $out $ref
  rm normal.sam tumor.sam normal.bam tumor.bam normal.bam.bai tumor.bam.bai mc.glf.txt mc.lower_bounds.txt cand_somatic_mutations cigar.variants.txt cigar.libraries.txt
  rm igv.batch normal.pileup tumor.pileup normal.sort.bam.bai tumor.sort.bam.bai windows.igv
  rm -r pics
  rm mc.haplotypes.txt
  rm mc.hap2_haplotypes.txt  mc.mm.txt  mc.non-mm.txt
done
