#!/bin/bash
# Christian Panse <cp@fgcz.ethz.ch>


#SBATCH --mem-per-cpu=1024M
#SBATCH --cpus-per-task=1
#SBATCH -N 1 # count of nodes                   
#SBATCH -n 1 # cpu count / number of tasks
#SBATCH --partition=zip
#SBATCH --nodelist=fgcz-r-035
#SBATCH --job-name="timsTOF-zipper"
#SBATCH -o /tmp/slurm.%x.%j.out
#SBATCH -e /tmp/slurm.%x.%j.err

## TODO(cp):
# srun -p prx --nodelist=fgcz-r-038 /src/gitlab.bfabric.org/proteomics/qc/inst/reports/tims_qc_multi/fgcz_slurm_timsTOF_qc_report.bash $1


declare -r ZIP_OPTIONS="-0 -@ -m"

test -x ${ZIP} || { echo "zip not available"; exit 1; }

timstof_zip(){
  local source_dir="$1"
  local output_zip_file="$2"

  if [[ -f ${output_zip_file} ]];
  then
	  echo "output zip already existing"
	  ls -l ${output_zip_file}
	  return
  fi

  if [ -d ${source_dir} ];
  then
  	mkdir -p `dirname ${output_zip_file}` || { echo "ERROR: could not make `dirname $output_zip_file`"; exit 1; }
  
  	test `/bin/ls ${source_dir} | /usr/bin/wc -l` -gt 0  || { echo "Warning; no files exit 0."; exit 0; }

    set -x
  	cd ${source_dir}/.. \
    && pwd \
    && echo `basename $source_dir` \
    && find `basename $source_dir` -type f \
      | zip ${ZIP_OPTIONS} ${output_zip_file} \
  	&& cd - \
  	&& unzip -l ${output_zip_file} || { echo "ERROR: $0 failed in '$source_dir'."; exit 1; }
  fi
}

echo "sbatch /home/cpanse/src/gitlab.bfabric.org/proteomics/qc/inst/reports/tims_qc_multi/fgcz_slurm_timsTOF_qc_report.bash $2" >> /tmp/qcqcqcqc.log
sbatch /home/cpanse/src/gitlab.bfabric.org/proteomics/qc/inst/reports/tims_qc_multi/fgcz_slurm_timsTOF_qc_report.bash $2

timstof_zip $1 $2


exit 0

