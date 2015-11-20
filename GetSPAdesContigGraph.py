#!/usr/bin/env python


# Copyright 2015 Ryan Wick

# This file is part of GetSPAdesContigGraph.

# GetSPAdesContigGraph is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.

# GetSPAdesContigGraph is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details.

# You should have received a copy of the GNU General Public License along with
# GetSPAdesContigGraph.  If not, see <http:# www.gnu.org/licenses/>.


from __future__ import division
from __future__ import print_function
import sys
import os
import argparse


def main():
    args = getArguments()

    links = loadGraph(args.graph)
    contigs = loadContigs(args.contigs)
    paths = loadPaths(args.paths)

    addLinksToContigs(contigs, paths, links)

    output = []
    for contig in contigs:
        output.append(contig.getHeaderWithLinks())
        output.append(contig.getSequenceWithLineBreaks())

    if args.output != "":
        outputFile = open(args.output, 'w')
        for line in output:
            outputFile.write(line)
    else:
        for line in output:
            print(line, end='')






def getArguments():
    parser = argparse.ArgumentParser(description='GetSPAdesContigGraph: a tool for creating a contig graph using a SPAdes assembly graph')
    parser.add_argument('graph', help='The assembly_graph.fastg file made by SPAdes')
    parser.add_argument('contigs', help='A contigs or scaffolds fasta file made by SPAdes')
    parser.add_argument('paths', help='The paths file which corresponds to the contigs or scaffolds file')
    parser.add_argument('-o', '--output', action='store', help='Save the output graph to this file (default: write graph to stdout)', default="")

    return parser.parse_args()




# This function takes a contig filename and returns a list of Contig objects.
# It expects a file of just foward contigs, but it creates both forward and
# reverse complement Contig objects.
def loadContigs(contigFilename):

    contigs = []

    contigFile = open(contigFilename, 'r')

    name = ''
    sequence = ''
    for line in contigFile:

        strippedLine = line.strip()

        # Skip empty lines.
        if len(strippedLine) == 0:
            continue

        # Header lines indicate the start of a new contig.
        if strippedLine[0] == '>':

            # If a contig is currently stored, save it now.
            if len(name) > 0:
                contig = Contig(name, sequence)
                contigs.append(contig)
                contigs.append(getReverseComplementContig(contig))
                name = ''
                sequence = ''

            name = strippedLine[1:]

        # If not a header line, we assume this is a sequence line.
        else:
            sequence += strippedLine

    # Save the last contig.
    if len(name) > 0:
        contig = Contig(name, sequence)
        contigs.append(contig)
        contigs.append(getReverseComplementContig(contig))

    return contigs





# This function takes a graph filename and returns a dictionary of the graph
# links.
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

        # If there are no outgoing links from this sequence, there's nothing
        # else to do.
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




# This function uses the contents of the paths dictionary to add link
# information to the contigs.
def addLinksToContigs(contigs, paths, links):

    # First store the first and last graph segments for each contig.
    for contig in contigs:
        firstSegment, lastSegment = paths[contig.fullname]
        contig.addGraphSegments(firstSegment, lastSegment)

    # Now for each contig, we take the last graph segment, find the segments
    # that it leads to, and then find the contigs which start with that next
    # segment.  These make up the links to the current contig.
    for contig1 in contigs:
        endingSegment = contig1.endingSegment
        followingSegments = links[endingSegment]
        linkedContigs = []
        for followingSegment in followingSegments:
            for contig2 in contigs:
                startingSegment = contig2.startingSegment
                if followingSegment == startingSegment:
                    linkedContigs.append(contig2)
        contig1.addLinkedContigs(linkedContigs)






# This function takes a contig and returns its reverse complement contig.
def getReverseComplementContig(contig):

    revCompContigFullname = contig.fullname

    # Add or remove the ' as necessary
    if revCompContigFullname[-1] == "'":
        revCompContigFullname = revCompContigFullname[0:-1]
    else:
        revCompContigFullname = revCompContigFullname + "'"

    revCompSequence = getReverseComplement(contig.sequence)

    return Contig(revCompContigFullname, revCompSequence)





def getReverseComplement(forwardSequence):

    reverseComplement = ""
    for i in reversed(range(len(forwardSequence))):
        base = forwardSequence[i]

        if base == 'A': reverseComplement += 'T'
        elif base == 'T': reverseComplement += 'A'
        elif base == 'G': reverseComplement += 'C'
        elif base == 'C': reverseComplement += 'G'
        elif base == 'a': reverseComplement += 't'
        elif base == 't': reverseComplement += 'a'
        elif base == 'g': reverseComplement += 'c'
        elif base == 'c': reverseComplement += 'g'
        elif base == 'R': reverseComplement += 'Y'
        elif base == 'Y': reverseComplement += 'R'
        elif base == 'S': reverseComplement += 'S'
        elif base == 'W': reverseComplement += 'W'
        elif base == 'K': reverseComplement += 'M'
        elif base == 'M': reverseComplement += 'K'
        elif base == 'r': reverseComplement += 'y'
        elif base == 'y': reverseComplement += 'r'
        elif base == 's': reverseComplement += 's'
        elif base == 'w': reverseComplement += 'w'
        elif base == 'k': reverseComplement += 'm'
        elif base == 'm': reverseComplement += 'k'
        elif base == 'B': reverseComplement += 'V'
        elif base == 'D': reverseComplement += 'H'
        elif base == 'H': reverseComplement += 'D'
        elif base == 'V': reverseComplement += 'B'
        elif base == 'b': reverseComplement += 'v'
        elif base == 'd': reverseComplement += 'h'
        elif base == 'h': reverseComplement += 'd'
        elif base == 'v': reverseComplement += 'b'
        elif base == 'N': reverseComplement += 'N'
        elif base == 'n': reverseComplement += 'n'
        elif base == '.': reverseComplement += '.'
        elif base == '-': reverseComplement += '-'
        elif base == '?': reverseComplement += '?'
        else: reverseComplement += 'N'

    return reverseComplement



# This class holds a contig: its name, sequence and length.
class Contig:

    def __init__(self, name, sequence):
        self.fullname = name
        nameParts = name.split("_")
        self.number = int(nameParts[1])
        self.sequence = sequence
        self.linkedContigs = []

    def __str__(self):
        return self.fullname

    def __repr__(self):
        return self.fullname

    def addGraphSegments(self, startingSegment, endingSegment):
        self.startingSegment = startingSegment
        self.endingSegment = endingSegment

    def addLinkedContigs(self, linkedContigs):
        self.linkedContigs = linkedContigs

    # This function produces a FASTG header for the contig, with links and a
    # line break at the end.
    def getHeaderWithLinks(self):
        headerWithEdges = '>' + self.fullname

        if len(self.linkedContigs) > 0:
            headerWithEdges += ":"
            for linkedContig in self.linkedContigs:
                headerWithEdges += linkedContig.fullname + ','
            headerWithEdges = headerWithEdges[0:-1]

        headerWithEdges += ";\n"
        return headerWithEdges

    def getSequenceWithLineBreaks(self):
        sequenceRemaining = self.sequence
        returnSequence = ""
        while len(sequenceRemaining) > 60:
            returnSequence += sequenceRemaining[0:60] + '\n'
            sequenceRemaining = sequenceRemaining[60:]
        returnSequence += sequenceRemaining + '\n'
        return returnSequence








# Standard boilerplate to call the main() function to begin the program.
if __name__ == '__main__':
    main()