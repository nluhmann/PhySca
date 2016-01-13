__author__ = 'Nina'


#scaffold reconstructed adjacencies at each internal node




def scaffoldAdjacencies(reconstructedAdj):
    print "Scaffolding..."
    scaffoldsPerNode = {}
    for node in reconstructedAdj:
        adjacencies = reconstructedAdj[node]
        #key: (leftmost marker in scaffold, rightmost marker in scaffold), value: [marker]
        scaffold = {}
        for adj in adjacencies:
            inserted = False
            #check if adjacency can extend an existing scaffold, extend only one scaffold if there are two possibilities
            for scaff in scaffold.keys():
                if adj[0] == scaff[0]:
                    markerList = scaffold.pop(scaff)
                    other = getOtherExtremity(adj[1])
                    markerList.insert(0,adj[1])
                    markerList.insert(0,other)
                    scaffold[(other,scaff[1])] = markerList
                    inserted = True
                    break
                elif adj[0] == scaff[1]:
                    #remove and return current scaffold
                    markerList = scaffold.pop(scaff)
                    #extend markerList with other extremity of the current adjacency
                    right = getOtherExtremity(adj[1])
                    markerList.append(adj[1])
                    markerList.append(right)
                    #save new scaffold with updated key
                    scaffold[(scaff[0],right)] = markerList
                    #adjacency has been inserted
                    inserted = True
                    #get out as soon as one scaffold has been extended
                    break
                elif adj[1] == scaff[0]:
                    markerList = scaffold.pop(scaff)
                    other = getOtherExtremity(adj[0])
                    markerList.insert(0, adj[0])
                    markerList.insert(0,other)
                    #save new scaffold with updated key
                    scaffold[(other,scaff[1])] = markerList
                    inserted = True
                    break
                elif adj[1] == scaff[1]:
                    markerList = scaffold.pop(scaff)
                    other = getOtherExtremity(adj[0])
                    markerList.append(adj[0])
                    markerList.append(other)
                    #save new scaffold with updated key
                    scaffold[(scaff[0],other)] = markerList
                    inserted = True
                    break
            if not inserted:
                #create new scaffold
                left = getOtherExtremity(adj[0])
                right = getOtherExtremity(adj[1])
                scaffold[(left,right)] = [left,adj[0],adj[1],right]

        scaffold2 = mergeScaffolds(scaffold)
        scaffoldsPerNode[node] = scaffold2
    return scaffoldsPerNode


def mergeScaffolds(scaffold):
    #then merge scaffolds that have the same border
        merg = True
        while merg:
            merg = False
            #case 1: two scaffolds have the same left border (find by sorting keylist!)
            keys = scaffold.keys()

            keys.sort(key=lambda tup: tup[0])
            for i in range(0,len(keys)-1):
                if keys[i][0] == keys[i+1][0]-1 and keys[i][0] % 2 == 1:
                    #merge scaffolds
                    firstList = scaffold.pop(keys[i])
                    secondList = scaffold.pop(keys[i+1])
                    #reverse first list, then chop of two last elements, then merge lists
                    rev = firstList[::-1]
                    revChopped = rev[:-2]
                    revChopped += secondList
                    scaffold[(revChopped[0],revChopped[-1])] = revChopped
                    merg = True

            #case 2: two scaffolds have the same right border (find by sorting keylist!)
            keys = scaffold.keys()
            keys.sort(key=lambda tup: tup[1])
            for i in range(0,len(keys)-1):
                if keys[i][1] == keys[i+1][1]-1 and keys[i][1] % 2 == 1:
                    #merge scaffolds
                    firstList = scaffold.pop(keys[i])
                    secondList = scaffold.pop(keys[i+1])

                    #reverse first list, then chop of two last elements, then merge lists
                    rev2 = firstList[::-1]
                    revChopped2 = rev2[2:]
                    secondList += revChopped2
                    scaffold[(secondList[0],secondList[-1])] = secondList
                    merg = True


        #case 3,4: the higher,lower extremity is at the right border, the lower,higher at the left (we have to compare all of them here)
        merged = True
        while merged:
            merged = False
            keys = scaffold.keys()
            for i in range(0,len(keys)-1):
                changed = False
                for j in range(i+1,len(keys)):
                    if keys[i][1] == keys[j][0]+1 and keys[i][1] % 2 == 0:
                        firstList = scaffold.pop(keys[i])
                        secondList = scaffold.pop(keys[j])
                        secondChop = secondList[2:]
                        firstList += secondChop
                        scaffold[(firstList[0],firstList[-1])] = firstList
                        changed = True
                        break
                    elif keys[i][0] == keys[j][1]-1 and keys[i][0] % 2 == 1:
                        firstList = scaffold.pop(keys[i])
                        secondList = scaffold.pop(keys[j])
                        firstChop = firstList[2:]
                        secondList += firstChop
                        scaffold[(secondList[0],secondList[-1])] = secondList
                        changed = True
                        break
                    elif keys[i][1] == keys[j][0]-1 and keys[i][1] % 2 == 1:
                        firstList = scaffold.pop(keys[i])
                        secondList = scaffold.pop(keys[j])
                        secondChop = secondList[2:]
                        firstList += secondChop
                        scaffold[(firstList[0],firstList[-1])] = firstList
                        changed = True
                        break
                    elif keys[i][0] == keys[j][1]+1 and keys[i][0] % 2 == 0:
                        firstList = scaffold.pop(keys[i])
                        secondList = scaffold.pop(keys[j])
                        firstChop = firstList[2:]
                        secondList += firstChop
                        scaffold[(secondList[0],secondList[-1])] = secondList
                        changed = True
                        break
                if changed:
                    merged = True
                    break
        return scaffold


def getOtherExtremity(extremity):
    if int(extremity) % 2 == 0:
        other = int(extremity) - 1
    elif int(extremity) % 2 == 1:
        other = int(extremity) + 1
    else:
        print "Houston there is a problem!"
    return other


def undoubleScaffolds(scaffoldsPerNode):
    undoubled = {}
    for node in scaffoldsPerNode:
        undoubled[node] = []
        scaffolds = scaffoldsPerNode[node]
        for sc in scaffolds:
            scaffold = []
            scaff = scaffolds[sc]
            for i in range(0,len(scaff)-1,2):
                first = int(scaff[i])
                second = int(scaff[i+1])

                if first - second == -1:
                    mark = str(second/2)

                elif first - second == 1:
                    marker = first/2
                    mark = "-"+str(marker)
                else:
                    print "Houston there is a problem!"
                scaffold.append(mark)

            undoubled[node].append(scaffold)
    return undoubled

def outputUndoubledScaffolds(scaffolds,output):
    file = open(output,"w")
    for node in scaffolds:
        file.write(">"+node+" "+str(len(scaffolds[node]))+"\n")
        counter = 1
        for scaff in scaffolds[node]:
            file.write("#"+str(counter)+"\n")
            file.write(" ".join(str(x) for x in scaff)+" $\n")
            counter = counter + 1
    file.close()

def outputScaffolds(scaffolds,output):
    file = open(output,"w")
    for node in scaffolds:
        file.write(">"+node+" "+str(len(scaffolds[node]))+"\n")
        counter = 1
        for scaff in scaffolds[node]:
            file.write("#"+str(counter)+"\n")
            file.write(" ".join(str(x) for x in scaffolds[node][scaff])+" $\n")
            counter = counter + 1
    file.close()


def sanityCheckScaffolding(undoubledScaffolds):
    for node in undoubledScaffolds:
        scaffolds = undoubledScaffolds[node]
        testHash = {}
        for scaff in scaffolds:
            for elem in scaff:
                if "-" in str(elem):
                    elem = elem[1:]
                if elem in testHash:
                    print node
                    print elem
                else:
                    testHash[elem] = 1

def compareAdjacencyLists(list1, list2):
    set1 = set(list1)
    set2 = set(list1)
    diff1 = set1.difference(set2)
    diff2 = set2.difference(set1)