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
ortheusPath="${chunkingDir}/ortheus"
PYTHONPATH="${chunkingDir}/ortheus/src/python"

echo $PYTHONPATH >&2

tree1="((((((((human:0.006969,chimp:0.009727):0.025291,((baboon:0.008968,macaque:0.014961):0.011019):0.024581):0.023649,(dusky_titi:0.025441,(owl_monkey:0.012081,marmoset:0.031665):0.008905):0.029440):0.066673,(mouse_lemur:0.058520,galago:0.123757):0.030984):0.018405,((rat:0.081244,mouse:0.072818):0.238435,rabbit:0.204316):0.021892):0.023260,(((cow:0.164728,(dog:0.107805):0.049576):0.004663,rfbat:0.142667):0.010883,((hedgehog:0.013056):0.216601,shrew:0.267086):0.057014):0.033242):0.028346,armadillo:0.161811):0.016015,(elephant:0.103584,tenrec:0.251557):0.035);"

tree2="((((((((human:0.006969,chimp:0.009727):0.025291,(colobus_monkey:0.014610,(baboon:0.008968,macaque:0.014961):0.011019):0.024581):0.023649,((owl_monkey:0.012081,marmoset:0.031665):0.008905):0.029440):0.066673,(mouse_lemur:0.058520,galago:0.123757):0.030984):0.018405,((rat:0.081244,mouse:0.072818):0.238435,rabbit:0.204316):0.021892):0.023260,(((cow:0.164728,(dog:0.107805):0.049576):0.004663,rfbat:0.142667):0.010883,((hedgehog:0.013056):0.216601,shrew:0.267086):0.057014):0.033242):0.028346,armadillo:0.161811):0.016015,(elephant:0.103584,tenrec:0.251557):0.035);"

tree3="((((((((human:0.006969,chimp:0.009727):0.025291,(colobus_monkey:0.014610,(baboon:0.008968,macaque:0.014961):0.011019):0.024581):0.023649,(dusky_titi:0.025441,(owl_monkey:0.012081,marmoset:0.031665):0.008905):0.029440):0.066673,(mouse_lemur:0.058520,galago:0.123757):0.030984):0.018405,((rat:0.081244,mouse:0.072818):0.238435):0.021892):0.023260,(((cow:0.164728,(dog:0.107805):0.049576):0.004663,rfbat:0.142667):0.010883,((hedgehog:0.013056):0.216601,shrew:0.267086):0.057014):0.033242):0.028346,armadillo:0.161811):0.016015,(elephant:0.103584,tenrec:0.251557):0.035);"

tree4="((((((((human:0.006969,chimp:0.009727):0.025291,(colobus_monkey:0.014610,(baboon:0.008968,macaque:0.014961):0.011019):0.024581):0.023649,(dusky_titi:0.025441,(owl_monkey:0.012081,marmoset:0.031665):0.008905):0.029440):0.066673,(mouse_lemur:0.058520,galago:0.123757):0.030984):0.018405,((rat:0.081244,mouse:0.072818):0.238435,rabbit:0.204316):0.021892):0.023260,(((cow:0.164728,(dog:0.107805):0.049576):0.004663):0.010883,((hedgehog:0.013056):0.216601,shrew:0.267086):0.057014):0.033242):0.028346,armadillo:0.161811):0.016015,(elephant:0.103584,tenrec:0.251557):0.035);"

tree5="((((((((human:0.006969,chimp:0.009727):0.025291,(colobus_monkey:0.014610,(baboon:0.008968,macaque:0.014961):0.011019):0.024581):0.023649,(dusky_titi:0.025441,(owl_monkey:0.012081,marmoset:0.031665):0.008905):0.029440):0.066673,(mouse_lemur:0.058520,galago:0.123757):0.030984):0.018405,((rat:0.081244,mouse:0.072818):0.238435,rabbit:0.204316):0.021892):0.023260,(((cow:0.164728,(dog:0.107805):0.049576):0.004663,rfbat:0.142667):0.010883,((hedgehog:0.013056):0.216601,shrew:0.267086):0.057014):0.033242):0.028346):0.016015,(elephant:0.103584,tenrec:0.251557):0.035);"

string1=colobus_baboon

string2=dusky_owl_monkey

string3=rabbit_mouse_lemur

string4=rfbat_dog

string5=armadillo_elephant

seq1=002

seq2=005

seq3=012

seq4=015

seq5=018

#8
anc1=12 

#14
anc2=24

#20
anc3=40

#28
anc4=56

#38
anc5=72

clo1=10   
#003

clo2=22	  
#006

clo3=34	  
#008

clo4=58   
#014

clo5=74   
#019

echo `ls -l` >&2 

${ortheusPath}/src/python/ortheus/utils/shellScripts/script3.sh ${anc1} ${clo1} ${seq1}.*.fa ${string1} "${tree1}" "${orthArgs1}" ${ortheusPath} ${passToOrtheus} ${outFile}
${ortheusPath}/src/python/ortheus/utils/shellScripts/script3.sh ${anc2} ${clo2} ${seq2}.*.fa ${string2} "${tree2}" "${orthArgs1}" ${ortheusPath} ${passToOrtheus} ${outFile}
${ortheusPath}/src/python/ortheus/utils/shellScripts/script3.sh ${anc3} ${clo3} ${seq3}.*.fa ${string3} "${tree3}" "${orthArgs1}" ${ortheusPath} ${passToOrtheus} ${outFile}
${ortheusPath}/src/python/ortheus/utils/shellScripts/script3.sh ${anc4} ${clo4} ${seq4}.*.fa ${string4} "${tree4}" "${orthArgs1}" ${ortheusPath} ${passToOrtheus} ${outFile}
${ortheusPath}/src/python/ortheus/utils/shellScripts/script3.sh ${anc5} ${clo5} ${seq5}.*.fa ${string5} "${tree5}" "${orthArgs1}" ${ortheusPath} ${passToOrtheus} ${outFile}

rcp $outFile $back
cd $back
rm -r $chunkingDir
