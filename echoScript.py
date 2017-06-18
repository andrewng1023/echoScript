# Test script for generating Echo picklist

# Import an csv file that is the output of the dueber lab website in "Table"
# mode. This csv file should contain all the information needed to generate a
# protocol for the user as well as a picklist for the Echo

# # Master Mix Recipe
#
# Perform GG reactions in a 2 ul total volume:
#
# 0.1 ul BsaI/BsmBI (Esp3I)
#
# 0.1 T4 DNA Ligase
#
# 0.2 10x Ligase Buffer + PEG
#
# 1.2 ul H20 (assume ~4 parts per rxn, can modify this if there is larger variance between parts)

# Do rxns in a 2 ul volume (Prepare a master mix based on the number of reactions,
# also accounting for the minimum dead volume of the plate

import numpy as np
import pandas as pd
from string import ascii_uppercase
import cPickle as pkl
import time
import argparse
import os.path

#Import the csv
def import_instructions(fname):
    rawdf = pd.read_csv(fname,names=['GoldenGate','Marker','Pieces','Enzyme','Sequence'])
    return rawdf

# First extract the plasmids that will be required for the assembly
def getparts(rawdf):
    allparts = list()

    for i in rawdf.index:
        for part in rawdf.loc[i].Pieces.split(','):
            allparts.append(part.strip())

    counts = [[x,allparts.count(x)] for x in set(allparts)]
    return counts

# Generate a new source plate as dictionaries
def newsourceplate():
    coord = list()
    for char in ascii_uppercase[0:16]:
        for col in range(24):
            coord.append(char+str(col+1))
            sourcepl = dict.fromkeys(coord)
    return sourcepl

# Now import existing source plate, or if none then generate a new source plate

def checksourceplate(counts, sname, platetype):
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

        userinstr = {count[0]:None for count in counts}
        for idx,count in enumerate(counts):
            sourcepl[coord[idx]] = {count[0]: 25}
            userinstr[count[0]] = 'Add '+str(MV)+ ' ul to well ' +str(coord[idx])

        allsourcepl = {'Source[1]': sourcepl}

    else: #import the users source plate dictionary
        allsourcepl = pkl.load( open(sname, "rb" ) )
        userinstr = {count[0]:None for count in counts}

        #Check if plasmids already exist in the dictionary
        for plate in allsourcepl:
            for count in counts:
                for key, value in allsourcepl[plate].iteritems():
                    if value == None:
                        continue
                    elif value.keys()[0] == count[0]:
                        #Check the volume of the plasmid in the well
                        if value.values() > DV + .1*count[1]: #Enough volume, return all clear on that plasmid
                            userinstr[count[0]] = [key,'Plasmid already exists with sufficient volume']
                        else:
                            userinstr[count[0]] = [key,'Add '+str(MV)+ ' ul to well '+key]
                            allsourcepl[plate][key][allsourcepl[plate][key].keys()[0]] += MV
            #Check to see if the plate is full. Save the empty plate
            if not all(allsourcepl[plate].values()):
                emptyplate = plate

        for plasmid, instr in userinstr.iteritems():
            if instr == None: #Find all the plasmids that do not yet exist in the counts, then add them to the most recent plate
                counter = 0
                while allsourcepl[emptyplate][coord[counter]] != None:
                    counter += 1

                well = coord[counter]
                allsourcepl[emptyplate][well] = {plasmid: MV}
                userinstr[plasmid] = [well,'Add '+str(MV)+ ' ul to well '+well]

        if not all(userinstr.values()):
            sourcepl = newsourceplate()
            for plasmid, instr in userinstr.iteritems():
                if instr == None: #Find all the plasmids that do not yet exist in the counts, then add them to the most recent plate
                    counter = 0
                    while sourcepl[coord[counter]] != None:
                        counter += 1

                    well = coord[counter]
                    sourcepl[well] = {plasmid: MV}
                    userinstr[plasmid] = [well,'Add '+str(MV)+ ' ul to well '+well]
            allsourcepl['Source[' + str(len(allsourcepl)+1) + ']'] = sourcepl

    return allsourcepl, userinstr

def makepicklists(rawdf, allsourcepl,dname):
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
    plasmidpicklist = pd.DataFrame(columns = ['Source Plate Name', 'Source Well', 'Destination Plate Name', 'Destination Well', 'Transfer Volume'])
    mmpicklist = pd.DataFrame(columns = ['Source Plate Name', 'Source Well', 'Destination Plate Name', 'Destination Well', 'Transfer Volume'])

    destinationloc = dict()

    j=0
    for i in rawdf.index:
        for part in rawdf.loc[i].Pieces.split(','):
            part = part.strip()
            #Find the plasmid in the source plate
            for plate in allsourcepl:
                for well in allsourcepl[plate]:
                    if allsourcepl[plate][well] is None:
                        continue
                    if allsourcepl[plate][well].keys()[0] == part:
                        plasmidpicklist.loc[j] = [plate,well,'Destination[1]',coord_ninetysix[i],100] #hardcode in to do less than 96 GGs at once for now
                        allsourcepl[plate][well][part] -= .1
                        j+=1
        destinationloc[rawdf.loc[i].GoldenGate.rstrip()] = coord_ninetysix[i]

        #Assume a regular dead volume plate
        if rawdf.loc[i].Enzyme == 'BsaI':
            mmpicklist.loc[i] = ['Source[1]',coord_threeeightyfour[i/29],'Destination[1]',coord_ninetysix[i],1800]
        else:
            mmpicklist.loc[i] = ['Source[1]',coord_threeeightyfour[i/29 + 24],'Destination[1]',coord_ninetysix[i],1800]

    date = time.strftime("%Y%m%d")

    ppdest = date + 'plasmidpicklist'
    counter = 1
    while os.path.isfile(ppdest+'.csv'):
        ppdest = ppdest + '_' + str(counter)
        counter += 1

    mmdest = date + 'mmpicklist'
    counter = 1
    while os.path.isfile(mmdest+'.csv'):
        mmdest = mmdest + '_' + str(counter)
        counter += 1

    plasmidpicklist.to_csv(ppdest+'.csv',index=False)
    mmpicklist.to_csv(mmdest+'.csv',index=False)
    #Save the source plates to a file
    dumpdest = date + dname
    counter = 1
    while os.path.isfile(dumpdest+'.p'):
        dumpdest = dumpdest + '_' + str(counter)
        counter += 1

    pkl.dump(allsourcepl, open(dumpdest +'.p', "wb"),-1)

    return destinationloc, plasmidpicklist, mmpicklist

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', action='store', dest = 'fname', type=str,
        help = 'CSV file of Dueber lab output')
    parser.add_argument('-s', action='store', dest = 'sname', type=str, default = None,
        help = 'File name of source plate')
    parser.add_argument('-p', action='store', dest = 'platetype', type=str, default = None,
        help = 'Plate Type')
    parser.add_argument('-d', action='store', dest = 'dname', type=str, default = 'sourceplates',
        help = 'Destination for new source plate file (no need for .p ext)')

    args = vars(parser.parse_args())
    rawdf = import_instructions(args['fname'])
    counts = getparts(rawdf)
    allsourcepl, userinstr = checksourceplate(counts,args['sname'],args['platetype'])
    destinationloc, plasmidpicklist, mmpicklist = makepicklists(rawdf,allsourcepl,args['dname'])

    print "{:<8} {:<15}".format('Plasmid','Instructions')
    for k, v in userinstr.iteritems():
        print "{:<8} {:<15}".format(k, v)

    print "{:<8} {:<15}".format('Plasmid','Location')
    for k, v in destinationloc.iteritems():
        print "{:<8} {:<15}".format(k, v)
