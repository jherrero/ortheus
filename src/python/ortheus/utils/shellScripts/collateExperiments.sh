#!/bin/sh


ortheusPath=$1
lineNo=$2
numberOfReps=$3
fragmentNo=$4

masterScriptPath=${ortheusPath}/src/python/ortheus/utils/shellScripts/collate.sh

suffix='tar.out.cross_Valid'
rm chunk.${suffix}.out
for var in 100
do
	echo $1
	echo ${var} >> chunk.${suffix}.out
	${masterScriptPath} ${fragmentNo} chunk ${suffix} ${lineNo} ${ortheusPath} ${numberOfReps} >> chunk.${suffix}.out
done

suffix='tar.out.samples'
#rm chunk.${suffix}.out
for var in 500 200 100 50 10 5 2 1
do
	echo $1
	#echo ${var} >> chunk.${suffix}.out
	#${masterScriptPath} ${fragmentNo} chunk ${suffix}.${var} ${lineNo} ${ortheusPath} ${numberOfReps} >> chunk.${suffix}.out
done

suffix=tar.out.memoryGaps
#rm chunk.${suffix}.out
for var in 3000 1500 750 500
do
	echo $1
	#echo ${var} >> chunk.${suffix}.out
	#${masterScriptPath} ${fragmentNo} chunk ${suffix}.${var} ${lineNo} ${ortheusPath} ${numberOfReps} >> chunk.${suffix}.out
done

suffix=tar.out.traceGap
#rm chunk.${suffix}.out
for var in 1000 500 200 100 50 25 10 0
do
	echo $1
	#echo ${var} >> chunk.${suffix}.out
	#${masterScriptPath} ${fragmentNo} chunk ${suffix}.${var} ${lineNo} ${ortheusPath} ${numberOfReps} >> chunk.${suffix}.out
done

suffix=tar.out.relax
#rm chunk.${suffix}.out
for var in 0 5 10 20 30 50
do
	echo $1
	#echo ${var} >> chunk.${suffix}.out
	#${masterScriptPath} ${fragmentNo} chunk ${suffix}.${var} ${lineNo} ${ortheusPath} ${numberOfReps} >> chunk.${suffix}.out
done

var1=0.056
var2=0.084
suffix=tar.out.params.${var1}.${var2}
#rm chunk.${suffix}.out
#echo ${var1} ${var2} >> chunk.${suffix}.out
#${masterScriptPath} $fragmentNo chunk ${suffix} ${lineNo} ${ortheusPath} ${numberOfReps} >> chunk.${suffix}.out

var1=0.056
var2=0.0373333
suffix=tar.out.params.${var1}.${var2}
#rm chunk.${suffix}.out
#echo ${var1} ${var2} >> chunk.${suffix}.out
#${masterScriptPath} $fragmentNo chunk ${suffix} ${lineNo} ${ortheusPath} ${numberOfReps} >> chunk.${suffix}.out
