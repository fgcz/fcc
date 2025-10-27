#!/bin/bash
# Witold E. Wolski <wew@fgcz.ethz.ch> 2015-07-08
# Christian Panse <cp@fgcz.ethz.ch> 2025-10-27
#SBATCH --mem-per-cpu=4096M
#SBATCH --cpus-per-task=4
#SBATCH -N 1 # count of nodes                   
#SBATCH -n 1 # cpu count / number of tasks
#SBATCH --partition=zip
#SBATCH --nodelist=fgcz-r-035
#SBATCH --job-name=rawrr
#SBATCH -o /tmp/slurm.%x.%j.out
#SBATCH -e /tmp/slurm.%x.%j.err


my_source=$1
my_target=$2
ZIP=/usr/bin/zip
test -x $ZIP || { echo "zip not available"; exit 1; }




if [ -d $my_source ];
then
	[[ $(/bin/ls ${my_source} | /usr/bin/wc -l) -eq 0 ]] && { echo "Warning: empty directory ${my_source} done."; exit 0; }

	echo "waters zip is sleeping for 2 seconds ..."
	sleep 2
	echo "waters zip continues!"

	mkdir -p `dirname $my_target` || { echo "ERROR: could not make `dirname $my_target`"; exit 1; }


	cd $my_source/.. \
	&& ls -l \
	&& pwd \
	&& $ZIP -rm $my_target `basename $my_source` \
	&& cd - \
	&& unzip -l $my_target || { echo "ERROR: $0 failed in '$my_source'."; exit 1; }
fi

exit 0

