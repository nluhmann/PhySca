from ete2 import Tree
import getAdjacencies
import globalAdjacencyGraph
import SR
import scaffolding
import argparse
import time

from calculate_SCJ import calculate_SCJ
import threading, Queue


#thread class
class  samplingThread(threading.Thread):

    lock = threading.Lock()
    allSampleReconstructionStatistic={}
    dict_SCJ={}
    input_q = Queue.Queue()
    output_q = Queue.Queue()

    #constructor of samplingThread
    def __init__(self):
        threading.Thread.__init__(self)
        self.stoprequest = threading.Event()
    #overridden join method
    def join(self, timeout=None):
        self.stoprequest.set()
        super(samplingThread, self).join(timeout)
    #run method of samplingThread
    def run(self):
            #input=(ccs, tree, extantAdjacencies, adjacencyProbs, alpha,i,  extantAdjacencies_species_adj)
            while not self.stoprequest.isSet():
                try:

                    #get input
                    params = samplingThread.input_q.get(True, 0.05)
                    ccs=params[0]
                    tree=params[1]
                    extantAdjacencies=params[2]
                    adjacencyProbs=params[3]
                    alpha=params[4]
                    i=params[5]
                    extantAdjacencies_species_adj=params[6]

                    jointLabels, first = SR.enumJointLabelings(ccs)
                    validLabels, validAtNode = SR.validLabels(jointLabels, first)

                    samplingThread.lock.acquire()
                    topDown = SR.sampleLabelings(tree, ccs, validAtNode, extantAdjacencies, adjacencyProbs, alpha)
                    samplingThread.lock.release()
                    reconstructedAdj = SR.reconstructedAdjacencies(topDown)
                    SR.outputReconstructedAdjacencies(reconstructedAdj, "reconstructed_adjacencies_" + str(i))

                    for node in reconstructedAdj:
                        # count for each adjaency on each internal node, how often this adjacencies over all samples occurs there
                        for adjacency in reconstructedAdj[node]:
                            samplingThread.lock.acquire()
                            if (node,adjacency) in self.allSampleReconstructionStatistic:
                                self.allSampleReconstructionStatistic[(node,adjacency)] += 1
                            else:
                                self.allSampleReconstructionStatistic.update({(node,adjacency):1})
                            samplingThread.lock.release()
                    scaffolds = scaffolding.scaffoldAdjacencies(reconstructedAdj)
                    undoubled = scaffolding.undoubleScaffolds(scaffolds)
                    scaffolding.outputUndoubledScaffolds(undoubled, "undoubled_scaffolds_" + str(i))
                    scaffolding.outputScaffolds(scaffolds, "doubled_scaffolds_" + str(i))
                    scaffolding.sanityCheckScaffolding(undoubled)

                    samplingThread.lock.acquire()
                    scj = calculate_SCJ(tree, reconstructedAdj, extantAdjacencies_species_adj)
                    self.dict_SCJ.update({'Sample_' + str(i): scj})
                    samplingThread.lock.release()

                    samplingThread.output_q.put((self.allSampleReconstructionStatistic, self.dict_SCJ))

                    samplingThread.input_q.task_done()
                except Queue.Empty:
                    continue
#create samplesize worker threads which working on sampling jobs
def sample(ccs, tree, extantAdjacencies, adjacencyProbs, alpha, sampling, samplesize,  extantAdjacencies_species_adj):
    allSampleReconstructionStatistic = {}
    dict_SCJ = {}

    #create a pool of  workers (number=samplesize)
    pool = [samplingThread() for i in range(0,samplesize)]
    for thread in pool:
        thread.start()

    # create jobs (number=sampling) shared by the workers
    for i in range(0, sampling):
           samplingThread.input_q.put((ccs, tree, extantAdjacencies, adjacencyProbs, alpha, i, extantAdjacencies_species_adj))

    #receive output of jobs
    for i in range(0,sampling):

        output=samplingThread.output_q.get()
        threadStatistic=output[0]
        threadSCJ=output[1]

        for key in threadStatistic:
            allSampleReconstructionStatistic[key]=threadStatistic[key]

        for key in threadSCJ:
            dict_SCJ[key] = threadSCJ[key]

    #join all threads in pool
    for thread in pool:
        thread.join()

    return allSampleReconstructionStatistic,dict_SCJ
