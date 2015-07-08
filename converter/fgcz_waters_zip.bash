#!/bin/bash
# wew@fgcz.ethz.ch 2015-07-08


my_source=$1
my_target=$2

if [ -d $my_source ];
then
	mkdir -p `dirname $my_target` || { echo "ERROR: could not make `dirname $my_target`"; exit 1; }

	cd $my_source/.. \
	&& pwd \
	&& zip -r $my_target `basename $my_source` \
	&& cd - \
	&& unzip -l $my_target || { echo "ERROR: $0 failed."; exit 1; }
fi

exit 0

