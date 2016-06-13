from ete2 import Tree
import getAdjacencies
import globalAdjacencyGraph
import SR
import scaffolding
import argparse
import time

from calculate_SCJ import calculate_SCJ
import multiprocessing

#run method, does one sampling
def runSample(params):
        #retrieving the given parameter
        ccs=params[0]
        tree=params[1]
        extantAdjacencies=params[2]
        adjacencyProbs=params[3]
        alpha=params[4]
        i=params[5]
        extantAdjacencies_species_adj=params[6]
        outputDirectory=params[7]
        allSampleReconstructionStatistic={}
        dict_SCJ={}
        lock = multiprocessing.Lock()

        #start sampling method like in the Main.py
        jointLabels, first = SR.enumJointLabelings(ccs)
        validLabels, validAtNode = SR.validLabels(jointLabels, first)

        lock.acquire()
        topDown = SR.sampleLabelings(tree, ccs, validAtNode, extantAdjacencies, adjacencyProbs, alpha)
        lock.release()
        reconstructedAdj = SR.reconstructedAdjacencies(topDown)
        SR.outputReconstructedAdjacencies(reconstructedAdj, outputDirectory+"reconstructed_adjacencies_" + str(i))

        for node in reconstructedAdj:
            # count for each adjaency on each internal node, how often this adjacencies over all samples occurs there
            for adjacency in reconstructedAdj[node]:
                lock.acquire()
                if (node,adjacency) in allSampleReconstructionStatistic:
                    allSampleReconstructionStatistic[(node,adjacency)] += 1
                else:
                    allSampleReconstructionStatistic.update({(node,adjacency):1})
                lock.release()

        scaffolds = scaffolding.scaffoldAdjacencies(reconstructedAdj)
        undoubled = scaffolding.undoubleScaffolds(scaffolds)
        scaffolding.outputUndoubledScaffolds(undoubled, outputDirectory+"undoubled_scaffolds_" + str(i))
        scaffolding.outputScaffolds(scaffolds, outputDirectory+"doubled_scaffolds_" + str(i))
        scaffolding.sanityCheckScaffolding(undoubled)

        lock.acquire()
        scj = calculate_SCJ(tree, reconstructedAdj, extantAdjacencies_species_adj)
        dict_SCJ.update({'Sample_' + str(i): scj})
        lock.release()
        return (allSampleReconstructionStatistic,dict_SCJ)

