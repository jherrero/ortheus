#!/bin/sh

chunkName=$1
ortheusPath=$2
suffix=$3
orthArgs1=$4
orthArgs2=$5
passToOrth=$6

echo $chunkName >&2
echo $ortheusPath >&2
echo $suffix >&2
echo $orthArgs1 >&2
echo $orthArgs2 >&2
echo $passToOrth >&2

tree="((((((((human:0.006969,chimp:0.009727):0.025291,(colobus_monkey:0.014610,(baboon:0.008968,macaque:0.014961):0.011019):0.024581):0.023649,(dusky_titi:0.025441,(owl_monkey:0.012081,marmoset:0.031665):0.008905):0.029440):0.066673,(mouse_lemur:0.058520,galago:0.123757):0.030984):0.018405,((rat:0.081244,mouse:0.072818):0.238435,rabbit:0.204316):0.021892):0.023260,(((cow:0.164728,(dog:0.107805):0.049576):0.004663,rfbat:0.142667):0.010883,((hedgehog:0.013056):0.216601,shrew:0.267086):0.057014):0.033242):0.028346,armadillo:0.161811):0.016015,(elephant:0.103584,tenrec:0.251557):0.035);"

chunkingDir=/tmp/chunkName_${chunkName}.${suffix}

back=`pwd`
outFile=${chunkName}.out.${suffix}
mkdir $chunkingDir
cp ${chunkName} $chunkingDir
cp ${ortheusPath}/ortheusExt.tar $chunkingDir
cd $chunkingDir
tar -xf ${chunkName}     
tar -xf ortheusExt.tar
cd ortheus
make all
cd ..
cd chunks
ortheusPath='../ortheus'
PYTHONPATH="${chunkingDir}/ortheus/src/python/"

if [ $passToOrth -eq 1 ]
then
	python ${ortheusPath}/Ortheus.py -b -l '#-w ' -e *.fa -d ${tree} -f out1.mfa -l "$orthArgs1" -h out1.score -m "java -server -Xmx1700m " > one
	python ${ortheusPath}/Ortheus.py -b -l '#-w ' -e *.fa -d ${tree} -f out2.mfa -l "$orthArgs2" -h out2.score -m "java -server -Xmx1700m " >> one
fi

if [ $passToOrth -eq 0 ]
then
	python ${ortheusPath}/Ortheus.py -b -l '#-w ' -e *.fa -d ${tree} -f out1.mfa $orthArgs1 -h out1.score -m "java -server -Xmx1700m " > one
	python ${ortheusPath}/Ortheus.py -b -l '#-w ' -e *.fa -d ${tree} -f out2.mfa $orthArgs2 -h out2.score -m "java -server -Xmx1700m " >> one
fi

python ${ortheusPath}/src/python/ortheus/utils/CompareAncestors.py -b -d out1.mfa out1.score -e out2.mfa out2.score -f ${tree} -g *.ranges > ${outFile}
grep "total_time" one >> ${outFile}

rcp $outFile $back
cd $back
rm -r $chunkingDir
