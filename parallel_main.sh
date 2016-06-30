#!/bin/bash

#Script, for calling  PhySca (Main.py) parallel on different cores
#aufzurufen
#
#
#make a directory for each call
if [ $1 = "-h" ] || [ $1 = "--help" ]
	then
		echo "parallel_main.sh [-h](-nhx<nhx_tree> | -nf <newick_tree>)
    	(-a <adjacencies>/ -m <markers>)
    	[-out <output>]
    	[-pN <number of processes>]
    	[-alpha <alpha>] [-s <Z>] [-x <x>] [-kT <kT>]
    	[-phySca <phySca>]
    	
    	
		-h,--help	show this help message and exit
		-nhx <nhx_tree>	tree file in newick (NHX) format
		-nf <newick_tree>	path to the file with NEWICK-tree
		-a <adjacencies>	adjacencies in extant genomes
		-m <markers>	marker order of extant genomes
		-output <output>	path to directory for preprocessing output,
							default=./testlauf
		-pN <processnumber>	number of processes used for sampling. 
							Max: [number of cpu], default=1
		-alpha <alpha>	alpha parameter in objective function, [0,1], default=0.0
		-x <x>	assign potential adjacencies by weight threshold,[0,1], default=0
		-kT <kT> deClone constant, default=0.1
		-s <Z>	sample Z solutions for given set of parameters, default=0
		-phySca <phySca>	path to directory with main program (PhySca), default=./"
    	
    	exit
fi


out="${PWD}/testlauf"
for param in "$@"
	do
		parameterArray+=(${param})
	done

for((p=0;p<$#;p+=1))
	do
		echo ${parameterArray[$p]}
		case ${parameterArray[$p]} in
		-alpha) alpha=${parameterArray[$((p+1))]}
		;;
		-pN) processNumber=${parameterArray[$((p+1))]}
		;;
		-x) x=${parameterArray[$((p+1))]}
		;;
		-s) sampleNumber=${parameterArray[$((p+1))]}
		;;
		-nhx) nhx=${parameterArray[$((p+1))]}
		;;
		-nf) nf=${parameterArray[$((p+1))]}
		;;
		-a) adjacencies=${parameterArray[$((p+1))]}
		;;
		-m) markers=${parameterArray[$((p+1))]}
		;;
		-kT) kT=${parameterArray[$((p+1))]}
		;;
		-out) out=${parameterArray[$((p+1))]}
		;;
		-phySca) phySca=${parameterArray[$((p+1))]}
		;;
		*)
		;;
		esac
	done
 
if [ -z $nhx ] && [ -z $nf ]
	then
		echo "Error: No tree for preprocessing given"
		exit
	elif [ -z $nhx  ]
		then
			treeparam="-nf"
			tree="${nf}"

	else
		treeparam="-nhx"
		tree="${nhx}"
fi

if [ -z $markers ] && [ -z $adjacencies ]
	then
		echo "Error:No markers or adjacencies given"
		exit
	elif [ -z $markers ]
		then
			MAparam="-a"
			MA="${adjacencies}"
			
	else
		MAparam="-m"
		MA="${markers}"
		
fi

if [ -z $processNumber ]
	then
		echo "Using only one process"
		processNumber="1"
fi

if [ -z $x ]
	then
		echo "default: x=0 "
		x="0"
fi

if [ -z $kT ]
	then
		echo "default: kT=0.1 "
		kT="0.1"
fi

if [ -z $alpha ]
	then
		echo "default: alpha=0.0"
		alpha="0.0"
fi

if [ -z $sampleNumber ]
	then
		echo "No Samples"
		sampleNumber="0"
fi

if [ -z $out ]
    then
        echo "No directory for preprocessing output given."
        exit
    elif [ -d $out ]
	     then
		    echo "directory ${out} is already there."
	else
	    echo "create directory ${out}"
	    mkdir ${out}
				
fi

if [ -z $phySca ]
    then
        echo "Assuming PhySca scripts in current directory"
        phySca=${PWD}
    elif [ -d $phySca ]
	     then
		    echo "PhySca scripts are located in ${phySca}"
	else
	    echo "directory for PhySca scripts ${phySca} doesnt' exist."
	    exit

fi
#preprocessing: weighting With DeClone
python ${phySca}/weightingWithDeClone.py ${MAparam} ${MA} ${treeparam} ${tree} -kT ${kT} -out ${out}

extant="${out}/extant_adjacencies"
internal="${out}/weighted_internal_adjacencies"
if [ -z $nhx ]
    then
        nhx_tree="${out}/nhx_tree"
    else
        nhx_tree=$tree
fi
#create directories for each sample run
for ((i=0;i<$processNumber;i+=1))
	do 
		if [ -d ${PWD}/test_${i} ]
			then
				echo "directory ${PWD}/test_${i} is already there."
			else
				echo "create directory ${PWD}/test_${i}"
				mkdir ${PWD}/test_${i}
				
		fi
		directories[i]="${PWD}/test_${i}"
	done


START=$(date +%s.%N)

for ((i=0;i<$processNumber;i+=1))
	do
		cd ${directories[i]}
		if [ $i -eq 0 ]
			then
			#PhySca call with 0th sample
				python ${phySca}/Main.py -tree $nhx_tree -alpha $alpha -s $sampleNumber -x $x -extant $extant -internal $internal -out ${PWD} -sc $i > ${PWD}/log.txt&
		else
		    #PhySca call without 0th sample
			python ${phySca}/Main.py -tree $nhx_tree -alpha $alpha -s $sampleNumber -x $x -extant $extant -internal $internal -out ${PWD} -sk -sc $(($i*$sampleNumber)) > ${PWD}/log.txt&
		
		fi
		cd ..
	done
	wait

END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo ${DIFF}

#Now the different outputs need to be put together
#Common output folder
if [ -d ${PWD}/output ]
			then
				echo "directory ${PWD}/output is already there."
			else
				echo "create directory ${PWD}/output"
				mkdir ${PWD}/output
fi
outputDir=${PWD}/output

printf "Combined log file of ${i} runs, each with ${sampleNumber} samples" > $outputDir/log.txt
printf "Combined SCJ distances of ${i} runs, each with ${sampleNumber} samples\n" > $outputDir/SCJ_distances

#First of all: files doubled_scaffolds_*, undoubled_scaffolds_*, reconstructed_adjacencies_*
#then the files, which needs to be combined: SCJ_distances, log.txt, statistic_allSampled_ReconstructedAdjacencies
#and conflict (only for the first process)
for ((i=0;i<$processNumber;i+=1))
	do
	    if [ $i -eq 0 ]
			then
		        mv ${directories[i]}/conflicts $outputDir
	    fi
        mv ${directories[i]}/doubled_scaffolds* $outputDir
        mv ${directories[i]}/undoubled_scaffolds* $outputDir
        mv ${directories[i]}/reconstructed_adjacencies* $outputDir
        cat ${directories[i]}/log.txt >> $outputDir/log.txt
        rm ${directories[i]}/log.txt
        cat ${directories[i]}/SCJ_distances >> $outputDir/temp_distances
        rm ${directories[i]}/SCJ_distances
        fileList[i]=${directories[i]}/statistic_allSampled_ReconstructedAdjacencies

    done
flist=$(printf ",%s" "${fileList[@]}")
python ${phySca}/putStatisticFilesTogether.py ${flist:1} ${outputDir}/statistic_allSampled_ReconstructedAdjacencies

#"  " in the sort statement is an tab
cat ${outputDir}/temp_distances | sort -t " " -nk1 >> ${outputDir}/SCJ_distances
rm ${outputDir}/temp_distances
for ((i=0;i<$processNumber;i+=1))
    do
       rm  ${directories[i]}/statistic_allSampled_ReconstructedAdjacencies
       rm -r ${directories[i]}
    done

echo "Done"
