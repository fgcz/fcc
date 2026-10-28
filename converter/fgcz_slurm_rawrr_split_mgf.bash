#!/bin/bash
# Christian Panse <cp@fgcz.ethz.ch>
#
#SBATCH --mem-per-cpu=4096M
#SBATCH --cpus-per-task=16
#SBATCH -N 1 # count of nodes                   
#SBATCH -n 1 # cpu count / number of tasks
#SBATCH --partition=prx
#SBATCH --nodelist=fgcz-c-072
#SBATCH --job-name=splitMGF
#SBATCH -o /tmp/slurm.%x.%j.out
#SBATCH -e /tmp/slurm.%x.%j.err

set -e
set -x

input=$1
output=$2


SCRATCH="/scratch/fcc/converterID-88-rawrr.$SLURM_JOB_ID/"

mkdir -p $(dirname $output) \
&& mkdir -p $SCRATCH \
&& cd $SCRATCH/ \
&& ln -svf $input \
&& echo "source('/home/cpanse/src/gitlab.bfabric.org/proteomics/qc/inst/scripts/rawrr_mgf_extract_top5_most_intense_ions.R'); mgf_split('$(basename $input)')" | nice -19 R -q --no-save \
&& cp *.mgf $(dirname $output)
