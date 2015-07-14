#!/bin/bash
# Witold E. Wolski <wew@fgcz.ethz.ch> 2015-07-08

my_source=$1
my_target=$2
ZIP=/usr/bin/zip
test -x $ZIP || { echo "zip not available"; exit 1; }


if [ -d $my_source ];
then
	mkdir -p `dirname $my_target` || { echo "ERROR: could not make `dirname $my_target`"; exit 1; }

	test `/bin/ls $my_source | /usr/bin/wc -l` -gt 0  || { echo "Warning; no files exit 0."; exit 0; }

	cd $my_source/.. \
	&& ls -l \
	&& pwd \
	&& $ZIP -rm $my_target `basename $my_source` \
	&& cd - \
	&& unzip -l $my_target || { echo "ERROR: $0 failed in '$my_source'."; exit 1; }
fi

exit 0

