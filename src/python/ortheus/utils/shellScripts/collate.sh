#!/bin/sh


fragmentNo=$1
prefix=$2
suffix=$3
lineNo=$4
ortheusPath=$5
numberOfReps=$6

rep=0

#echo $fragmentNo
#echo $prefix
#echo $suffix 
#echo $lineNo 
#echo $ortheusPath
#echo $numberOfReps

while [ ${rep} -lt ${numberOfReps} ]
do
	python ${ortheusPath}/src/python/ortheus/utils/shellScripts/collate.py $fragmentNo $prefix ${suffix}.${rep} $lineNo >> collate.sh.out
	#echo ${rep}
	rep=`expr ${rep} + 1`
done
echo >> collate.sh.out
cat collate.sh.out
rm collate.sh.out
