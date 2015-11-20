#!/usr/bin/env python


# Copyright 2015 Ryan Wick

# This file is part of GetSPAdesContigGraph.

# GetSPAdesContigGraph is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free 
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.

# GetSPAdesContigGraph is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.

# You should have received a copy of the GNU General Public License along with
# GetSPAdesContigGraph.  If not, see <http:# www.gnu.org/licenses/>.


from __future__ import division
from __future__ import print_function
import sys
import subprocess
import os
import argparse
import datetime
import shutil
import itertools


def main():
    startTime = datetime.datetime.now()
    args = getArguments()


    links = loadGraph(args.graph)
    contigs = loadContigs(args.contigs)
    paths = loadPaths(args.paths)

    addGraphSegmentsToContigs(contigs, paths)

    for contig in contigs:
        print(contig.fullname)
        print(contig.startingSegment)
        print(contig.endingSegment)
        print("")

    quit()





    # Create the output graph
    outputLines = []


    for contig in contigs:
        outputLines.append(contig.getHeaderWithEdges())
        outputLines += contig.sequenceLines

    if args.output != "":
        outputFile = open(args.output, 'w')
        for line in outputLines:
            outputFile.write(line + '\n')
    else:
        for line in outputLines:
            print(line)






def getArguments():
    parser = argparse.ArgumentParser(description='GetSPAdesContigGraph: a tool for creating a contig graph using a SPAdes assembly graph')
    parser.add_argument('graph', help='The assembly_graph.fastg file made by SPAdes')
    parser.add_argument('contigs', help='A contigs or scaffolds fasta file made by SPAdes')
    parser.add_argument('paths', help='The paths file which corresponds to the contigs or scaffolds file')
    parser.add_argument('-o', '--output', action='store', help='Save the output graph to this file (default: write graph to stdout)', default="")

    return parser.parse_args()




# This function takes a contig filename and returns a list of Contig objects.
def loadContigs(contigFilename):

    contigs = []
    contigFile = open(contigFilename, 'r')

    name = ''
    sequenceLines = []
    for line in contigFile:

        strippedLine = line.strip()

        # Skip empty lines.
        if len(strippedLine) == 0:
            continue

        # Header lines indicate the start of a new contig.
        if strippedLine[0] == '>':

            # If a contig is currently stored, save it now.
            if len(name) > 0:
                contig = Contig(name, sequenceLines)
                contigs.append(contig)
                name = ''
                sequenceLines = []

            name = strippedLine[1:]

        # If not a header line, we assume this is a sequence line.
        else:
            sequenceLines.append(strippedLine)

    # Save the last contig.
    if len(name) > 0:
        contig = Contig(name, sequenceLines)
        contigs.append(contig)

    return contigs





# This function takes a graph filename and returns a dictionary of the graph links.
# The dictionary key is the starting graph segment.
# The dictionary value is a list of the ending graph segments.
def loadGraph(graphFilename):

    links = {}
    graphFile = open(graphFilename, 'r')

    for line in graphFile:
        line = line.strip()

        # Skip empty lines.
        if len(line) == 0:
            continue

        # Skip lines that aren't headers
        if line[0] != '>':
            continue

        # Remove the last semicolon
        if line[-1] == ';':
            line = line[:-1]

        lineParts = line.split(":")

        startingSegment = lineParts[0]
        start = getNumberFromFullSequenceName(startingSegment)
        if start not in links:
            links[start] = []
        startRevComp = getOppositeSequenceNumber(start)
        if startRevComp not in links:
            links[start] = []

        # If there are no outgoing links from this sequence, there's nothing else to do.
        if len(lineParts) < 2:
            continue

        endingSegments = lineParts[1].split(",")
        ends = []
        for endingSegment in endingSegments:
            ends.append(getNumberFromFullSequenceName(endingSegment))

        # Add the links to the dictionary in the forward direction.
        for end in ends:
            if end not in links[start]:
                links[start].append(end)

        # Add the links to the dictionary in the reverse direction:
        for end in ends:
            endRevComp = getOppositeSequenceNumber(end)
            if endRevComp not in links:
                links[endRevComp] = []
            if startRevComp not in links[endRevComp]:
                links[endRevComp].append(startRevComp)

    return links





# This function takes a path filename and returns a dictionary.
# The dictionary key is the contig name.
# The dictionary value is a tuple:
#   the first graph segment in the path
#   the last graph segment in the path
def loadPaths(pathFilename):

    paths = {}

    pathFile = open(pathFilename, 'r')

    contigName = ''
    pathString = ''
    for line in pathFile:

        line = line.strip()

        # Skip empty lines.
        if len(line) == 0:
            continue

        # Lines starting with 'NODE' are the start of a new path
        if len(line) > 3 and line[0:4] == 'NODE':

            # If a contig is currently stored, save it now.
            if len(contigName) > 0:
                paths[contigName] = getFirstAndLastSegmentFromPathString(pathString)
                contigName = ''
                pathString = ''

            contigName = line

        # If not a node name line, we assume this is a path line.
        else:
            pathLine = line
            if pathLine[-1] == ';':
                pathLine = pathLine[0:-1]
            if len(pathString) == 0:
                pathString = pathLine
            else:
                pathString += ',' + pathLine

    # Save the last contig.
    if len(contigName) > 0:
        paths[contigName] = getFirstAndLastSegmentFromPathString(pathString)

    return paths




def getFirstAndLastSegmentFromPathString(pathString):
    pathStringParts = pathString.split(',')
    firstPart = pathStringParts[0]
    lastPart = pathStringParts[-1]
    return (firstPart, lastPart)





def getNumberFromFullSequenceName(fullSequenceName):
    number = fullSequenceName.split("_")[1]
    if fullSequenceName[-1] == "'":
        number += "-"
    else:
        number += "+"
    return number




def getOppositeSequenceNumber(number):
    if number[-1] == "+":
        return number[0:-1] + "-"
    else:
        return number[0:-1] + "+"




# This function uses the contents of the paths dictionary to add
# first/last graph segment information to the contigs.
def addGraphSegmentsToContigs(contigs, paths):
    for contig in contigs:
        firstSegment, lastSegment = paths[contig.fullname]
        contig.addGraphSegments(firstSegment, lastSegment)





# This class holds a contig: its name, sequence and length.
class Contig:

    def __init__(self, name, sequenceLines):
        self.fullname = name
        nameParts = name.split("_")
        self.number = int(nameParts[1])
        self.sequenceLines = sequenceLines

    def __str__(self):
        return self.fullname

    def __repr__(self):
        return self.fullname

    def addGraphSegments(self, startingSegment, endingSegment):
        self.startingSegment = startingSegment
        self.endingSegment = endingSegment

    def getHeaderWithEdges(self):
        headerWithEdges = '>' + self.fullname

        # ADD THE EDGES HERE

        headerWithEdges += ";"
        return headerWithEdges




# This class holds a contig: its name, sequence and length.
class Path:

    def __init__(self, name, pathLines):
        self.firstSegment = "DSSDF" # TO DO: get the actual first segment
        self.lastSegment = "DFGDFG" # TO DO: get the actual last segment








# Standard boilerplate to call the main() function to begin the program.
if __name__ == '__main__':
    main()