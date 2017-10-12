import numpy as np
import pandas as pd
from string import ascii_uppercase
import cPickle as pkl
import time
import argparse
import os.path
from collections import OrderedDict

#Import the csv
def importPCR_instructions(fname):
    rawdf = pd.read_csv(fname)
    return rawdf

# Generate a new source plate as dictionaries
def newsourceplate():
    coord = list()
    for char in ascii_uppercase[0:16]:
        for col in range(24):
            coord.append(char+str(col+1))
            sourcepl = OrderedDict.fromkeys(coord)
    return sourcepl

def getoligos(rawdf):
    oligos = OrderedDict()
    count = 0

    while pd.isnull(rawdf.iloc[count,0]) == False:
        oligos[rawdf.iloc[count,0]] = rawdf.iloc[count,3]
        count += 1

    return oligos

def checkOligoSourcePlate(oligos, sname, platetype):
    date = time.strftime("%Y%m%d")
    coord = list()
    for char in ascii_uppercase[0:16]:
        for col in range(24):
            coord.append(char+str(col+1))

    if platetype is None:
        DV = 4 #Assume using a low dead volume plate
        MV = 12
    elif platetype == 'LDVplate':
        DV = 4 #dead volume for a regular dead volume plate
        MV = 12
    elif platetype == 'rDVplate':
        DV = 15 #dead volume for a regular dead volume plate
        MV = 60

    if sname is None: #The user does not have a source plate saved yet, so create one
        sourcepl = newsourceplate()

        userinstr = OrderedDict()
        for idx,oligo in enumerate(oligos):
            sourcepl[coord[idx]] = {oligos[oligo]: MV}
            userinstr[oligos[oligo]] = 'Add '+str(MV)+ ' ul to well ' +str(coord[idx])

        allsourcepl = OrderedDict({'Source[1]': sourcepl})

    else: #import the users source plate dictionary
        allsourcepl = pkl.load( open(sname, "rb" ) )
        userinstr = OrderedDict()

        #Check if plasmids already exist in the dictionary
        for plate in allsourcepl:
            for oligo in oligos:
                for key, value in allsourcepl[plate].iteritems():
                    if value.keys()[0] == oligo:
                        #Check the volume of the plasmid in the well
                        if value.values() > DV: #Enough volume, return all clear on that plasmid
                            userinstr[oligos[oligo]] = [key,'Oligo already exists with sufficient volume']
                        else:
                            userinstr[oligos[oligo]] = [key,'Add '+str(MV)+ ' ul to well '+key]
                            allsourcepl[plate][key][allsourcepl[plate][key].keys()[0]] += MV
            #Check to see if the plate is full. Save the empty plate
            if not all(allsourcepl[plate].values()):
                emptyplate = plate

        for oligo, instr in userinstr.iteritems():
            if instr == None: #Find all the plasmids that do not yet exist in the counts, then add them to the most recent plate
                counter = 0
                while allsourcepl[emptyplate][coord[counter]] != None:
                    counter += 1

                well = coord[counter]
                allsourcepl[emptyplate][well] = {oligos[oligo]: MV}
                userinstr[oligos[oligo]] = [well,'Add '+str(MV)+ ' ul to well '+well]

        if not all(userinstr.values()):
            sourcepl = newsourceplate()
            for oligo, instr in userinstr.iteritems():
                if instr == None: #Find all the plasmids that do not yet exist in the counts, then add them to the most recent plate
                    counter = 0
                    while sourcepl[coord[counter]] != None:
                        counter += 1

                    well = coord[counter]
                    sourcepl[well] = {oligos[oligo]: MV}
                    userinstr[oligos[oligo]] = [well,'Add '+str(MV)+ ' ul to well '+well]
            allsourcepl['Source[' + str(len(allsourcepl)+1) + ']'] = sourcepl

    return allsourcepl, userinstr

def makeOligoPicklists(rawdf, oligos, allsourcepl, dname):
    #Define 96 well plate coordinates
    coord_ninetysix = list()
    for char in ascii_uppercase[0:8]:
        for col in range(12):
            coord_ninetysix.append(char+str(col+1))

    coord_threeeightyfour = list()
    for char in ascii_uppercase[0:16]:
        for col in range(24):
            coord_threeeightyfour.append(char+str(col+1))

    #Define an empty picklist
    oligopicklist = pd.DataFrame(columns = ['Source Plate Name', 'Source Well', 'Destination Plate Name', 'Destination Well', 'Transfer Volume'])

    destinationloc = OrderedDict()

    row = 0
    while rawdf.iloc[row,0] != 'PCR':
        row +=1

    j=0
    for idx, i in enumerate(rawdf.iloc[row+1:,:].index):
        #Associate the oligo with a oAN number and locate that oligo in the source plate
        o1 = oligos[rawdf.iloc[i,1]]
        for plate in allsourcepl:
            for well in allsourcepl[plate]:
                if allsourcepl[plate][well] is None:
                    continue
                if allsourcepl[plate][well].keys()[0] == o1:
                    oligopicklist.loc[j] = [plate,well,'Destination[1]',coord_ninetysix[idx],100] #hardcode in to do less than 96 GGs at once for now
                    allsourcepl[plate][well][o1] -= .1
                    j+=1

        o2 = oligos[rawdf.iloc[i,2]]
        for plate in allsourcepl:
            for well in allsourcepl[plate]:
                if allsourcepl[plate][well] is None:
                    continue
                if allsourcepl[plate][well].keys()[0] == o2:
                    oligopicklist.loc[j] = [plate,well,'Destination[1]',coord_ninetysix[idx],100] #hardcode in to do less than 96 GGs at once for now
                    allsourcepl[plate][well][o2] -= .1
                    j+=1

        destinationloc[rawdf.iloc[i,0].rstrip()] = coord_ninetysix[idx]

    date = time.strftime("%Y%m%d")

    opdest = date + 'oligoPicklist'
    counter = 1
    while os.path.isfile(opdest+'.csv'):
        opdest = opdest + '_' + str(counter)
        counter += 1

    oligopicklist.to_csv(opdest+'.csv',index=False)
    #Save the source plates to a file
    dumpdest = date + dname
    counter = 1
    while os.path.isfile(dumpdest+'.p'):
        dumpdest = dumpdest + '_' + str(counter)
        counter += 1

    pkl.dump(allsourcepl, open(dumpdest +'.p', "wb"),-1)

    return destinationloc, oligopicklist

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', action='store', dest = 'fname', type=str,
        help = 'CSV file of Dueber lab output')
    parser.add_argument('-s', action='store', dest = 'sname', type=str, default = None,
        help = 'File name of source plate')
    parser.add_argument('-p', action='store', dest = 'platetype', type=str, default = None,
        help = 'Plate Type')
    parser.add_argument('-d', action='store', dest = 'dname', type=str, default = 'oligoSourcePlates.p',
        help = 'Destination for new source plate file')

    args = vars(parser.parse_args())
    rawdf = importPCR_instructions(args['fname'])
    oligos = getoligos(rawdf)
    allsourcepl, userinstr = checkOligoSourcePlate(oligos,args['sname'],args['platetype'])
    destinationloc, oligopicklist = makeOligoPicklists(rawdf,oligos,allsourcepl,args['dname'])

    print userinstr
    print destinationloc

    print "{:<8} {:<15}".format('Plasmid','Instructions')
    for k, v in userinstr.iteritems():
        print "{:<8} {:<15}".format(k, v)

    print "{:<8} {:<15}".format('Plasmid','Location')
    for k, v in destinationloc.iteritems():
        print "{:<8} {:<15}".format(k, v)
