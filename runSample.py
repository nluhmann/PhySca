from ete2 import Tree
import getAdjacencies
import globalAdjacencyGraph
import SR
import scaffolding
import argparse
import time

from calculate_SCJ import calculate_SCJ
import multiprocessing

def runSample(params):
        ccs=params[0]
        tree=params[1]
        extantAdjacencies=params[2]
        adjacencyProbs=params[3]
        alpha=params[4]
        i=params[5]
        extantAdjacencies_species_adj=params[6]
        allSampleReconstructionStatistic={}
        dict_SCJ={}
        lock = multiprocessing.Lock()

        jointLabels, first = SR.enumJointLabelings(ccs)
        validLabels, validAtNode = SR.validLabels(jointLabels, first)

        lock.acquire()
        topDown = SR.sampleLabelings(tree, ccs, validAtNode, extantAdjacencies, adjacencyProbs, alpha)
        lock.release()
        reconstructedAdj = SR.reconstructedAdjacencies(topDown)
        SR.outputReconstructedAdjacencies(reconstructedAdj, "reconstructed_adjacencies_" + str(i))

        for node in reconstructedAdj:
            # count for each adjaency on each internal node, how often this adjacencies over all samples occurs there
            for adjacency in reconstructedAdj[node]:
                #if node in allSampleReconstructionStatistic:
                #    # print allSampleReconstructionStatistic[node]
                #    if adjacency in allSampleReconstructionStatistic[node]:
                #        allSampleReconstructionStatistic[node][adjacency] += 1
                #    else:
                #        dict_adj = {adjacency: 1}
                #        allSampleReconstructionStatistic[node].update(dict_adj)
                #        # print allSampleReconstructionStatistic[node][adjacency]
                #else:
                #    allSampleReconstructionStatistic.update({node: {adjacency: 1}})
                lock.acquire()
                if (node,adjacency) in allSampleReconstructionStatistic:
                    allSampleReconstructionStatistic[(node,adjacency)] += 1
                else:
                    allSampleReconstructionStatistic.update({(node,adjacency):1})
                lock.release()
        scaffolds = scaffolding.scaffoldAdjacencies(reconstructedAdj)
        undoubled = scaffolding.undoubleScaffolds(scaffolds)
        scaffolding.outputUndoubledScaffolds(undoubled, "undoubled_scaffolds_" + str(i))
        scaffolding.outputScaffolds(scaffolds, "doubled_scaffolds_" + str(i))
        scaffolding.sanityCheckScaffolding(undoubled)

        lock.acquire()
        scj = calculate_SCJ(tree, reconstructedAdj, extantAdjacencies_species_adj)
        dict_SCJ.update({'Sample_' + str(i): scj})
        lock.release()
        return (allSampleReconstructionStatistic,dict_SCJ)

def sample(ccs, tree, extantAdjacencies, adjacencyProbs,alpha,sampling, samplesize,extantAdjacencies_species_adj):
    controlTree=tree
    pool = multiprocessing.Pool(processes=samplesize)
    tasks=(( ccs, tree,extantAdjacencies, adjacencyProbs, alpha, i,
              extantAdjacencies_species_adj) for i in range(0,sampling))
    results=pool.map_async(runSample, tasks)
    pool.close()
    pool.join()
    output=results.get()


    allSampleReconstructionStatistic={}
    dict_SCJ={}
    if tree==controlTree:
        print "Tree untouched!"
    else:
        print "------Tree Changed-------"
    for tuple in output:
        tempRS = tuple[0]
        tempSCJ = tuple[1]
        for key in tempRS:
            if key in allSampleReconstructionStatistic:
                allSampleReconstructionStatistic[key] +=tempRS[key]
            else:
                allSampleReconstructionStatistic.update({key:tempRS[key]})
        for key in tempSCJ:
            dict_SCJ[key]=tempSCJ[key]

    return allSampleReconstructionStatistic, dict_SCJ