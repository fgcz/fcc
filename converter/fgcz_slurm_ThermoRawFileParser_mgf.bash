#!/bin/bash
#
# Christian Panse <cp@fgcz.ethz.ch> 2025-10-30
#SBATCH --mem-per-cpu=1024M
#SBATCH --cpus-per-task=1
#SBATCH -N 1 # count of nodes                   
#SBATCH -n 1 # cpu count / number of tasks
#SBATCH --partition=zip
#SBATCH --nodelist=fgcz-r-035
#SBATCH --job-name=TRFPmgf
#SBATCH -o /tmp/slurm.%x.%j.out
#SBATCH -e /tmp/slurm.%x.%j.err


input=$1
output=$2

#input=/home/cpanse/.cache/R/ExperimentHub/18e7d5bfb4a3b_4590.raw
#output=/tmp/oo.mgf

mkdir -p $(dirname $output) \
&& /home/biobeamer/fcc/converter/ThermoRawFileParser-2.0/ThermoRawFileParser \
--input=$input \
--format=0 \
--mgfPrecursor \
--stdout > $output
