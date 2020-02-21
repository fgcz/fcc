#!/bin/bash
# Christian Panse <cp@fgcz.ethz.ch>


declare -r ZIP_OPTIONS="-0 -@ -m"

test -x ${ZIP} || { echo "zip not available"; exit 1; }

timstof_zip(){
  local source_dir="$1"
  local output_zip_file="$2"

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

timstof_zip $1 $2
exit 0

