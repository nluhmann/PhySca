import sys,os,subprocess,tempfile
import argparse
import random

##########################
### PARSING PARAMETERS ###
##########################

parser = argparse.ArgumentParser(description='Weights given tree in nhx-format with DeClone. Also converts tree in NEWICK-format into tree in nhx-format')
groupFormat = parser.add_mutually_exclusive_group(required=True)
groupFormat.add_argument("-nhx", "--nhx_Tree", type=str, help="path to the file with nhx-tree")
groupFormat.add_argument("-nf",'--Newick', type=str,help='path to the file with NEWICK-tree')
parser.add_argument("-i","--ignore_weights",action='store_const', const=True, help="boolean, for either ignore or consider edge length/weights, when parsing Newick Tree into nhx Tree")
parser.add_argument("-sm","--set_minimum", type=float, help="minimal value for any edge length, when parsing Newick Tree into nhx Tree", default=0.0)
groupAM = parser.add_mutually_exclusive_group(required=True)
groupAM.add_argument("-a","--adjacencies",type=str, help="path to adjacency-file")
groupAM.add_argument("-m","--markers",type=str,help="path to marker-file")
groupAM.add_argument("-f","--families",type=str,help="path to family file")
parser.add_argument("-kT",type=float,help="deClone constant", default=0.1)
parser.add_argument("-jP","--just_Parse",action='store_const', const=True, help="boolean, for either just parse the Newick-file or run DeClone after it.")
parser.add_argument("-out","--output",type=str, help="output directory, current directory as default", default=".")

args = parser.parse_args()


# create output directory if not existing
if not os.path.isdir(args.output):
    os.makedirs(args.output)


#if not args.adjacencies and not args.markers:
#    parser.error('Error: wrong parameter number or usage.')
#if not args.nhx_Tree and not args.Newick:
#    parser.error('Error: wrong parameter number or usage.')


#global variables
nhxFileOut = args.output+'/nhx_tree'
listOfExtWeightOut =  args.output+'/extant_adjacencies'
listOfIntWeightOut =  args.output+'/weighted_internal_adjacencies'
singleLeafAdjOut =  args.output+'/single_leaf_adjacencies'
listOfPotentialExtant = args.output+'/potential_extant_adjacencies'








#################################################
### TREE MANIPULATION AND WEIGHTING FUNCTIONS ###
#################################################
  
#convert NEWICK tree notation with weights into nhx notation
#parameter: file - the path to the file with the tree in NEWICK-notation
#           ignore_weights - boolean, if the parsed nhx-tree shall include weights or not
#			set_minimum - str for minimal edge length (default: 0.0) 
def parse_NEWICK_to_nhx(file,ignore_weights,minimal_branch_length=0.0):
    print "Parsing Newick into nhx format"
    nhx_tree = ''  # empty tree, will be filled according to nhx-format
    charpos = 0  # current charpos in line
    numberUnamedNodes = 0  # current unamed internal nodes
    numberOfSpecies = -1  # numberOfSpecies is equal to the number of colons
    nameStartPos = 0  # startposition of the last species name
    weightStartPos=0 #startposition of the last weight
    isLeaf = False  # current species is a leaf or not
    newick_signs = [';', ',', '(', ')', ':']  # signs which have a distinct function in the NEWICK-Format
    isWeight=False
    includesNoWeight=False
    f = open(file, "r")  # file contains only one line
    for line in f:
		# if a : is in the Newick-Tree the tree is very likely weighted
        if ':' not in line:
            includesNoWeight=True
        #for each character
        for char in line:
            if char == '(':
                nameStartPos = charpos #the current species name begins
                isLeaf = True  # a species directly after a ( is a leaf
                isWeight=False # after a (, there's never a weight
            elif char == ')':
                weightEndPos=charpos
                if not includesNoWeight:
                    nhx_tree=check_edge_length(minimal_branch_length,nhx_tree,line,weightStartPos,weightEndPos)
                nameEndPos = charpos  # the current species name ends
                if line[charpos - 1] not in newick_signs:  # if there was a species id directly before this ) ...
                    # ...the nhx-block needs to be added
                    nhx_tree += '[&&NHX:D=?:Ev='
                    if isLeaf:  # was the last species a leaf or an internal node?
                        nhx_tree += 'Extant'
                    else:
                        nhx_tree += 'Spec'
                    numberOfSpecies += 1  # the number of species is increased
                    lastSpeciesName=construct_species_name(line,nameStartPos,nameEndPos,numberUnamedNodes)
                    nhx_tree += ':S=' + str(numberOfSpecies) + ':ND=' + str(lastSpeciesName) + ']'
                isLeaf = False  # the next species directly behind a ) isn't a leaf
                nameStartPos=charpos
                isWeight=False # after a ), there's never directly a weight
            elif char == ',':
                weightEndPos=charpos
                if not includesNoWeight:
                    nhx_tree=check_edge_length(minimal_branch_length,nhx_tree,line,weightStartPos,weightEndPos)
                nameEndPos  = charpos  # the current species name ended
                if line[charpos - 1] not in newick_signs:  # if there was a species id directly before this ( ...
                    # ...the nhx-block needs to be added
                    nhx_tree += '[&&NHX:D=?:Ev='
                    if isLeaf:  # was the last species a leaf or an internal node
                        nhx_tree += 'Extant'
                    else:
                        nhx_tree += 'Spec'
                    numberOfSpecies += 1  # the number of species is increased
                    lastSpeciesName=construct_species_name(line,nameStartPos,nameEndPos,numberUnamedNodes)
                    nhx_tree += ':S=' + str(numberOfSpecies) + ':ND=' + str(lastSpeciesName) + ']'
                    isLeaf = True
                    nameStartPos  = charpos  # the current species name begins
                    isWeight=False
            elif char == ':':
                weightStartPos=charpos
                #setting temporarily a species name to check if it's empty
                tempSpeciesName=line[nameStartPos+1:charpos]
                #print 'temp: '+'start: '+str(nameStartPos)+' end: '+str(charpos)+' name: '+str(tempSpeciesName)
                if len(tempSpeciesName) == 0:
                    numberUnamedNodes += 1
                    tempSpeciesName = "node" + str(numberUnamedNodes)
                    nhx_tree+=tempSpeciesName
                if ignore_weights:
                    isWeight=True
            elif char == ';':  # the end of the Newick-tree is reached and a root node needs to be added
                if line[charpos-1] != ')': #if a identifier for root node is given
                    nameEndPos=charpos
                    lastSpeciesName=construct_species_name(line,nameStartPos,nameEndPos,numberUnamedNodes)
                    nhx_tree+="[&&NHX:D=?:Ev=Spec:S=" + str(numberOfSpecies + 1) + ":ND="+ str(lastSpeciesName) +"]"

                else:# if no identifier for root node is given
                    nhx_tree += "root[&&NHX:D=?:Ev=Spec:S=" + str(numberOfSpecies + 1) + ":ND=root]"
                isWeight=False
            #if the given Newick -tree doesn't include :, the check for unnamed internal nodes is here
            elif includesNoWeight and line[charpos-1]==')':
                #checking if the current internal node is unamed:
                if line[charpos]==',' or line[charpos]==')':
                    numberUnamedNodes += 1
                    tempSpeciesName = "node" + str(numberUnamedNodes)
                    nhx_tree+=tempSpeciesName

            #if the nhx-Tree shall not include weights, skip adding the current char to the tree, if it belongs to a weight    
            if not isWeight:
                nhx_tree += char  # append the current char to the nhx_tree
            charpos += 1
        nhx_tree.strip()
        print nhx_tree
    f.close()
    #write the nhx-tree into a file
    f = open(nhxFileOut,"w")
    f.write(nhx_tree)
    f.close()

def check_edge_length(minimal_branch_length,nhx_tree,line,weightStartPos,weightEndPos):
    #if current edge length is smaller than demanded minimal edge length
    #replace it with minimal edge length
    weight_str=line[weightStartPos+1:weightEndPos]
    weight=float(weight_str)
    if weight<minimal_branch_length:
       last_chars=len(weight_str)
       nhx_tree=nhx_tree[:-last_chars] #removing last characters from nhx_tree
       nhx_tree+=str(minimal_branch_length) # and replacing it with the minimal_branch_length
    return nhx_tree

#constructs the name of the current species and gives it back.
def construct_species_name(line,nameStartPos,nameEndPos,numberUnamedNodes):
    lastSpeciesName=line[nameStartPos+1:nameEndPos]
    if ':' in lastSpeciesName:
        lastSpeciesName=lastSpeciesName.split(':')[0]
        # if the current species is an unnamed internal node, name it node_number
        if len(lastSpeciesName) == 0:
            lastSpeciesName = "node" + str(numberUnamedNodes)
    return lastSpeciesName

#retrieve extant adjacencies
def readAdjacencyFile(file):
    print "Collect extant adjacencies from provided adjacency file"
    # keys: (left marker, right marker), value: [(species,chromosome),...]
    adjacencies = {}
    lookUpHash = {}
    speciesList = set()
    f = open(file,"r")
    for line in f:
        fields = line.rstrip("\n").split("#")
        adjas = fields[0].split(":")
        left = adjas[1].split(" ")[0]
        right = adjas[1].split(" ")[1]
        specieslist = []
        for extant in fields[1:]:
            species = extant.split(":")[0].split(".")[0]
            speciesList.add(species)
            specieslist.append((species,"mockchrom"))
        specieslist = set(specieslist)
        adjacencies[(left,right)] = specieslist

        if left in lookUpHash:
            lookUpHash[left].add(right)
        else:
            lookUpHash[left] = set()
            lookUpHash[left].add(right)
        if right in lookUpHash:
            lookUpHash[right].add(left)
        else:
            lookUpHash[right] = set()
            lookUpHash[right].add(left)

    return adjacencies, lookUpHash, list(speciesList)

def getAllAdjacenciesFromFamilyFile(file):
    print "Collect extant adjacencies from provided families file"
    # keys: (left marker, right marker), value: [(species,chromosome),...]
    adjacencies = {}
    species_hash = {}

    lookUpHash = {}
    f = open(file, "r")
    for line in f:
        if line.startswith(">"):
            fam_id = line.rstrip("\n")[1:]
        elif not line == "" and not line == "\n":
            species = line.split(":")[0].split(".")[0]
            x = line.rstrip("\n").split(" ")
            ori = x[1]
            coords = x[0].split(":")[1]
            s = coords.split("-")
            start = int(s[0])
            stop = int(s[1])
            tup = (fam_id, ori, start, stop)
            if species in species_hash:
                species_hash[species].append(tup)
            else:
                species_hash[species] = []
                species_hash[species].append(tup)
    f.close()

    speciesList = species_hash.keys()
    #sort all families for all species according to start position in genome
    for spec in species_hash:
        sorted_by_start = sorted(species_hash[spec], key=lambda tup: tup[2])
        #then adjacencies are contiguous entries in the sorted hash, but be careful with the orientation!
        for i in range(0, len(sorted_by_start)):
            #assume circular chromosomes
            #TODO: this should be a parameter here
            if i == len(sorted_by_start) -1:
                first = sorted_by_start[-1]
                second = sorted_by_start[0]
            else:
                first = sorted_by_start[i]
                second = sorted_by_start[i + 1]
            doubleFirst = doubleMarker(first[0])
            doubleSecond = doubleMarker(second[0])
            if first[1] == "+":
                left = doubleFirst[1]
            else:
                left = doubleFirst[0]
            if second[1] == "+":
                right = doubleSecond[0]
            else:
                right = doubleSecond[1]

            if left < right:
                adj = (str(left),str(right))
                rev = (str(right), str(left))
            else:
                rev = (str(left), str(right))
                adj = (str(right), str(left))

            if adj in adjacencies:
                adjacencies[adj].append((spec,"mockchrom"))
            elif rev in adjacencies:
                adjacencies[rev].append((spec, "mockchrom"))
            else:
                adjacencies[adj] = []
                adjacencies[adj].append((spec, "mockchrom"))

            if left in lookUpHash:
                lookUpHash[left].add(right)
            else:
                lookUpHash[left] = set()
                lookUpHash[left].add(right)

            if right in lookUpHash:
                lookUpHash[right].add(left)
            else:
                lookUpHash[right] = set()
                lookUpHash[right].add(left)


    return adjacencies, lookUpHash, speciesList

# double marker id to account for orientation
def doubleMarker(marker):
    if "-" in marker:
        return int(marker[1:]) * 2, (int(marker[1:]) * 2) - 1
    else:
        return (int(marker) * 2) - 1, int(marker) * 2

# double marker and find adjacencies in chromosomes
def findAdjacencies(speciesHash):
    print "Collect extant adjacencies from marker file..."
    # keys: (left marker, right marker), value: [(species,chromosome),...]
    adjacencies = {}
    lookUpHash = {}
    speciesList = []
    for species in speciesHash:
        speciesList.append(species)
        chromosomes = speciesHash[species]
        for chrom in chromosomes:
            orderList = chromosomes[chrom]
            for i in range(0,len(orderList)-1):
                first = orderList[i]
                second = orderList[i+1]
                doubleFirst = doubleMarker(first)
                doubleSecond = doubleMarker(second)

                #take right extrem of first and left extrem of second for adj tuple
                f = doubleFirst[1]
                s = doubleSecond[0]
                adj = (f,s)
                rev = (s,f)

                if adj in adjacencies:
                    adjacencies[adj].append((species,chrom))
                elif rev in adjacencies:
                    adjacencies[rev].append((species,chrom))
                else:
                    adjacencies[adj] = [(species,chrom)]

                if f in lookUpHash:
                    lookUpHash[f].add(s)
                else:
                    lookUpHash[f] = set()
                    lookUpHash[f].add(s)

                if s in lookUpHash:
                    lookUpHash[s].add(f)
                else:
                    lookUpHash[s] = set()
                    lookUpHash[s].add(f)
    return adjacencies, lookUpHash, speciesList

#for each extant adjacency, use declone to compute probability that the
#adjacency is present in an adjacency forest sampled randomly from a
#Boltzmann distribution. Then assign adjacency if it is above threshold to internal node of the tree.
def deCloneProbabilities(extantAdjacencies, kT, listOfInternalNodes, treefile):
    print "Compute probabilities with DeClone..."
    locateDeClone = os.path.dirname(sys.argv[0]) 
    adjacenciesPerNode = {}
    singleLeafAdj={}
    extantWeightedAdjacencies={}
    for adj in extantAdjacencies:
        for species in extantAdjacencies[adj]:
            node=species[0]
            if node in extantWeightedAdjacencies:
                extantWeightedAdjacencies[node].add(adj)
            else:
                adjset = set()
                adjset.add(adj)
                extantWeightedAdjacencies[node] = adjset
    for adjacency in extantAdjacencies:
        #produce list of extant adjacencies for declone
        path = os.path.dirname(os.path.realpath(__file__))
        tmpfile=tempfile.NamedTemporaryFile(delete=True) #create a temporary named file (appears with a arbitrary name in the directory)
        species = extantAdjacencies[adjacency]

        for spec in species:
            tmpfile.write(spec[0]+" "+spec[0]+"\n")
        tmpfile.seek(0) #go to the beginning of the tmpfile
        command = locateDeClone+'/DeClone -t1 '+treefile+' -t2 '+treefile+' -a '+tmpfile.name+' -i -kT '+str(kT)
        #use declone to compute probs
        output = subprocess.check_output(command, shell=True)
        #output is just matrix with probabilities
        #each line of the output should contain max one number greater 0, save for internal nodes
        lines = output.split("\n")
        for line in lines:
            if not line == "" and not line.startswith("\t") and not line.startswith(">"):
                node = line.split("\t")[0]  #for each internal node,find the prob for the current adj to be at this node
                probs = line.split("\t")[1]
                probs = probs.rstrip(" ")
                probArray = probs.split(" ")
                probArrayFl = [float(x) for x in probArray]
                probability = max(probArrayFl, key=float)
                #print probability
                if node in listOfInternalNodes: #if node is an internal one
                    if len(species)>1: #if an adjacency occurs just in one external leaf,
                                       # it's ignored for internal nodes (just evolved at this leaf)
                        if node in adjacenciesPerNode:
                            adjacenciesPerNode[node].add((adjacency,probability))
                        else:
                            adjset = set()
                            adjset.add((adjacency,probability))
                            adjacenciesPerNode[node] = adjset
                    else:
                            #ignored adjacencies with only one leaf occuring in
                            singleLeafAdj.update({adjacency:(species[0],probability)})


        tmpfile.close() #tmpfile is closed and immediately deleted

    #ignored adjacencies are written into a special file
    f=open(singleLeafAdjOut,'w')
    print('Removed '+str(len(singleLeafAdj))+' Adjacencies occurring just in one external node/leaf')
    for adj in singleLeafAdj:
        f.write('('+str(adj[0])+','+str(adj[1])+')'+'\t'+str(singleLeafAdj[adj][0])+'\t'+str(singleLeafAdj[adj][1])+'\n')
    f.close()


    print "Generating output..."
    # output format: >internal node  adjacency  weight
    #                >external node  adjacency

    #internal nodes and root
    file=open(listOfIntWeightOut, 'w')
    for node in adjacenciesPerNode:
        for adj_weight in adjacenciesPerNode[node]: #for each adjacency tuple with weight
            file.write('>'+str(node)+'\t') #write the node's name
            adj_L=str(adj_weight[0][0])
            adj_R=str(adj_weight[0][1])
            weight=str(adj_weight[1])
            file.write('('+adj_L+','+adj_R+')\t') #write the adjacency tuple
            file.write(weight+'\n') #write the adjacency weight

    file.close()
    #external leaves
    file=open(listOfExtWeightOut, 'w')
    for leaf in extantWeightedAdjacencies:
        for adj in extantWeightedAdjacencies[leaf]:
            file.write('>' + str(leaf) + '\t' + str(adj) +  '\n')
    file.close()

#get the internal nodes form the given nhx-treefile
#   -seperatingChar: either ':', if tree includes weights or '[', if not
#   -treefile: the path to the treefile
def get_internal_nodes_from_treefile(seperatingChar,treefile):
    print "Get internal nodes from treefile"
    listOfInternalNodes=[]
    #read the nhx-tree
    file=open(treefile,'r')
    line=file.readline() #treefile consists only of one line
    file.close()
    #split the tree at each ')', because each internal node's name starts right after one
    line_splitted_at_closed_brackets=line.strip().split(')')
    for element in line_splitted_at_closed_brackets[1:]:
        internal_node=element.strip().split(seperatingChar)[0]
        #in case of the root node the '[&&NHX...'-block needs to be removed
        listOfInternalNodes.append(internal_node.split('[')[0])
    return listOfInternalNodes

def read_Marker_file(marker):
 #compute adjacencies, compute weights with DeClone

    print "Number of undoubled marker: "
    # read undoubled marker, for each species
    species_marker_order = {}
    file = open(marker, "r")
    species = ""
    for line in file:
        # new species
        if line.startswith(">"):
            if not species == "":
                species_marker_order[species] = chromosomes
                print species
                print markerCount
            species = line.split(" ")[0][1:].rstrip("\n")
            chromosomes = {}
            markerCount = 0
            # new chromosome
        elif line.startswith("#"):
            chrom = line.rstrip("\n")[2:]
        elif not line == "\n":
            order = line.rstrip("\n")[:-2].split(" ")
            markerCount += len(order)
            chromosomes[chrom] = order
    species_marker_order[species] = chromosomes
    print species
    print markerCount
    file.close()
    return species_marker_order

# SE: here it is easy to find all potential adjacencies by computing the global set of adjacencies,
# then compare for each set of extant adjacencies and check if a potential adjacency induces a conflict
def identify_potential_extant_adjacencies(extantAdjacencies, lookUpHash, speciesList):
    potential_adjacencies = {}
    print "Identify potential extant adjacencies"
    # extant Adjacencies: keys: (left marker, right marker), value: [(species,chromosome),...]
    # collect all different adjacencies present in any extant adjacency
    #species_present = {}
    #for elem in extantAdjacencies:
    #    species_present = list(set(species_present).union(extantAdjacencies[elem]))
    #print "Find adjacencies"
    # check for each adjacency, if it can be a potential adjacency for all species not present
    # for this, check all other adjacencies that include one of the extremities and check if they are present in the species

    for adj in extantAdjacencies:
        if len(extantAdjacencies[adj]) > 1:
            species_missing = [x for x in speciesList if x not in extantAdjacencies[adj]]
            for spec in species_missing:
                first_extrem = adj[0]
                second_extrem = adj[1]
                good = True
                first_adj = lookUpHash[first_extrem]
                for elem in first_adj:
                    if (elem,first_extrem) in extantAdjacencies:
                        x = extantAdjacencies[(elem,first_extrem)]
                    else:
                        x = extantAdjacencies[(first_extrem, elem)]
                    for a in x:
                        if spec in a:
                            good = False
                            break
                    if not good:
                        break
                if good:
                    second_adj = lookUpHash[second_extrem]
                    for elem in second_adj:
                        if (elem,second_extrem) in extantAdjacencies:
                            x = extantAdjacencies[(elem,second_extrem)]
                        else:
                            x = extantAdjacencies[(second_extrem,elem)]
                        for a in x:
                            if spec in a:
                                good = False
                                break
                    if not good:
                        break
                # for el in extantAdjacencies:
                #     if not el == adj:
                #         if first_extrem in el or second_extrem in el:
                #             x = extantAdjacencies[el]
                #             for a in x:
                #                 if spec[0] in a:
                #                     good = False
                #                     break
                if good:
                    if adj in potential_adjacencies:
                        potential_adjacencies[adj].add(spec)
                    else:
                        potential_adjacencies[adj] = set()
                        potential_adjacencies[adj].add(spec)
    # write output
    file = open(listOfPotentialExtant, 'w')
    for adj in potential_adjacencies:
        for spec in potential_adjacencies[adj]:
            file.write('>' + str(spec) + '\t' + str(adj) + '\n')
    file.close()

    return potential_adjacencies

#TESTMETHOD
def simFragmentedExtantGenomes(extantAdjacencies,numFrags):
    #choose a couple of adjacencies and remove a random element in their species list
    for i in range (1,numFrags):
        k = random.choice(extantAdjacencies.keys())
        x = extantAdjacencies[k]
        if len(x) > 1:
            x.pop(random.randrange(len(x)))

    return extantAdjacencies







#################
### COMPUTING ###
#################

if args.nhx_Tree:
    # get List of internal nodes from treefile
    listOfInternalNodes=get_internal_nodes_from_treefile(':',args.nhx_Tree)
    if args.adjacencies:
        extantAdjacencies, lookUpHash = readAdjacencyFile(args.adjacencies)
    elif args.markers:
        extantAdjacencies, lookUpHash, speciesList =findAdjacencies(read_Marker_file(args.markers))
    elif args.families:
        extantAdjacencies=getAllAdjacenciesFromFamilyFile(args.families)
    else:
        parser.error('Error: wrong parameter number or usage.')
    identify_potential_extant_adjacencies(extantAdjacencies, lookUpHash, speciesList)
    deCloneProbabilities(extantAdjacencies, args.kT, listOfInternalNodes, args.nhx_Tree)

if args.Newick:

    if args.ignore_weights:
        parse_NEWICK_to_nhx(args.Newick,True,args.set_minimum)
        listOfInternalNodes=get_internal_nodes_from_treefile('[',nhxFileOut) #when no weights given, after a node's name comes a [
    else:
        parse_NEWICK_to_nhx(args.Newick,False,args.set_minimum)
        listOfInternalNodes=get_internal_nodes_from_treefile(':',nhxFileOut) #when weights given, after a node's name comes a :
    if args.just_Parse:
        print 'Parsed Newick-File to nhx-File.'
    else:
        if args.adjacencies:
            extantAdjacencies, lookUpHash = readAdjacencyFile(args.adjacencies)
        elif args.markers:
            extantAdjacencies, lookUpHash, speciesList = findAdjacencies(read_Marker_file(args.markers))
        elif args.families:
            extantAdjacencies = getAllAdjacenciesFromFamilyFile(args.families)

            ###!!! sim method, should be commented when you are using this script!!!
            #ToDo separate script for simulating extant fragmented genomes
            #extantAdjacencies = simFragmentedExtantGenomes(extantAdjacencies,1000)
        else:
            parser.error('Error: wrong parameter number or usage.')
        #TODO add parameter so this is not always done
        identify_potential_extant_adjacencies(extantAdjacencies, lookUpHash, speciesList)
        deCloneProbabilities(extantAdjacencies, args.kT, listOfInternalNodes, nhxFileOut)
