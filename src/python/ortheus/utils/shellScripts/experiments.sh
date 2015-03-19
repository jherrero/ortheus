#!/bin/sh

#ortheus/src/python/ortheus/utils/shellScripts/experiments.sh /lustre/work1/ensembl/bpaten/ortheus /lustre/work1/ensembl/bpaten/ortheusTest/error 5 /lustre/work1/ensembl/bpaten/ortheusTest

ortheusPath=$1
errorPath=$2
numberOfReps=$3
experimentPath=$4

masterScriptPath=${ortheusPath}/src/python/ortheus/utils/shellScripts/masterScript.sh

cd ${ortheusPath}
rm ortheusExt.tar
cd ..
cp ortheusExt.tar ${ortheusPath}
cd ${experimentPath}

for var in 100
do
	orthArgs='#-i '${var}
	${masterScriptPath} ${ortheusPath} ${errorPath} cross_Valid "${orthArgs}" "${orthArgs}" ${numberOfReps} 1 2
done

for var in 500 
#200 100 50 10 5 2 1
do
	orthArgs='#-i '${var}
	${masterScriptPath} ${ortheusPath} ${errorPath} samples.${var} "${orthArgs}" "${orthArgs}" ${numberOfReps} 1 
done
 
for var in 3000 1500 750 500
do
	orthArgs='-n '${var}
	#${masterScriptPath} ${ortheusPath} ${errorPath} memoryGaps.${var} "${orthArgs}" "${orthArgs}" ${numberOfReps} 0 
done

for var in 1000 500 200 100 50 25 10 0
do
	orthArgs1='-o '2000' -n '5000
	orthArgs2='-o '${var}' -n '`expr 3000 + ${var}` 
	echo ${orthArgs1} ${orthArgs2}
	#${masterScriptPath} ${ortheusPath} ${errorPath} traceGap.${var} "${orthArgs1}" "${orthArgs2}" ${numberOfReps} 0 
done

#
for var in 50 
#0 5 10 20 30 50
do
	orthArgs='#-j '${var}
	${masterScriptPath} ${ortheusPath} ${errorPath} relax.${var} "${orthArgs}" "${orthArgs}" ${numberOfReps} 1 
done
#
var1=0.056
var2=0.084
orthArgs1='#-o '${var1}' -q '${var1}
orthArgs2='#-o '${var2}' -q '${var2}
#${masterScriptPath} ${ortheusPath} ${errorPath} params.${var1}.${var2} "${orthArgs1}" "${orthArgs2}" ${numberOfReps} 1 
#
var1=0.056
var2=0.0373333
orthArgs1='#-o '${var1}' -q '${var1}
orthArgs2='#-o '${var2}' -q '${var2}
#${masterScriptPath} ${ortheusPath} ${errorPath} params.${var1}.${var2} "${orthArgs1}" "${orthArgs2}" ${numberOfReps} 1 
#