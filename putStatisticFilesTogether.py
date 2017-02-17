import sys



def putFilesTogether(dictionaries):
    statistic = {}
    for dic in dictionaries:
        for spec in dic:
            if spec in statistic:
                for adj in dic[spec]:
                    if adj in statistic[spec]:

                        temp=statistic[spec][adj]
                        temp+=dic[spec][adj]
                        statistic[spec][adj]=temp
                    else:
                        #adj ins't in statistic[spec] by now
                        statistic[spec].update({adj: dic[spec][adj]})
            else:
                #species isn't in statistic by now

                for adj in dic[spec]:
                    if spec in statistic:
                        statistic[spec].update({adj:dic[spec][adj]})
                    else:
                        statistic[spec]={}
                        statistic[spec].update({adj: dic[spec][adj]})
    return statistic

def readFile(filePath):
    dic={}
    f=open(filePath,"r")
    line=f.readline()
    while line:
        split=line.split("\t")
        species=split[0]
        if (species != ">Header:"):
            gene=split[1]
            number=int(split[2].strip())
            if species in dic:
                dic[species].update({gene:number})
            else:
                dic[species]={gene:number}
        line=f.readline()
    return dic

def writeOutput(statistic,output):
    f = open(output, "w")
    for species in statistic:
        for gene in statistic[species]:
            number=statistic[species][gene]
            f.write(str(species)+"\t"+gene+"\t"+str(number)+"\n")
    f.close()

if len(sys.argv) == 3:
    file_list=sys.argv[1].split(',')
    stat_list=[]
    for file in file_list:
        stat_list.append(readFile(file))
    writeOutput(putFilesTogether(stat_list),sys.argv[2])
else:
    print('Error: Wrong Parameter number!')