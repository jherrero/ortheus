#!/bin/sh

ortheusPath=$1
errorPath=$2
suffix=$3
orthArgs1=$4
orthArgs2=$5
numberOfReps=$6
passToOrth=$7
extra=$8

rep=0

echo $ortheusPath 
echo $errorPath 
echo $suffix 
echo $orthArgs1 
echo $orthArgs2 
echo $numberOfReps 
echo $passToOrth

while [ ${rep} -lt ${numberOfReps} ]
do
	for i in chunk*.tar
	do
		bsub -P birneygroup -q long -R"select[mem>7000] rusage[mem=7000]" -M 7000000 -e ${errorPath} ${ortheusPath}/src/python/ortheus/utils/shellScripts/script${extra}.sh ${i} ${ortheusPath} ${suffix}.${rep} "${orthArgs1}" "${orthArgs2}" $passToOrth 
		sleep 1
	done
	echo ${rep}
	rep=`expr ${rep} + 1`
done
