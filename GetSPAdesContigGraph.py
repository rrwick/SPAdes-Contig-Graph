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
import shutil
import subprocess
from distutils import spawn


def main():
    args = getArguments()

    # Load in the user-specified files.
    print("Loading graph...     ", end="")
    sys.stdout.flush()
    links = loadGraphLinks(args.graph)
    print("done\nLoading contigs...   ", end="")
    sys.stdout.flush()
    contigs = loadContigs(args.contigs)
    print("done\nLoading paths...     ", end="")
    sys.stdout.flush()
    paths = loadPaths(args.paths)
    print("done")
    sys.stdout.flush()

    # Add the paths to each contig object, so each contig knows its graph path,
    # and add the links to each contig object, turning the contigs into a
    # graph.
    print("Building graph...    ", end="")
    sys.stdout.flush()
    addPathsToContigs(contigs, paths)
    addLinksToContigs(contigs, links)
    print("done")

    # If the user chose to prioritise connections, then a number of graph
    # modifications are carried out.
    if (args.connection_priority):

        print("Removing duplicates... ", end="")
        sys.stdout.flush()
        contigs = removeDuplicateContigs(contigs)
        print("done")

        print("Splitting contigs... ", end="")
        sys.stdout.flush()

        if not isBlastInstalled():
            print("Error: could not find BLAST program", file=sys.stderr)
            quit()

        segmentSequences = loadGraphSequences(args.graph)
        graphOverlap = getGraphOverlap(links, segmentSequences)
        contigs = splitContigs(contigs, links, segmentSequences, graphOverlap)

        # Since the contigs have changed, we need to redo the linking.
        addLinksToContigs(contigs, links)
        print("done")

    # Output the graph to file
    print("Saving graph...      ", end="")
    sys.stdout.flush()
    outputFile = open(args.output, 'w')
    for contig in contigs:
        outputFile.write(contig.getHeaderWithLinks())
        outputFile.write(contig.getSequenceWithLineBreaks())
    print("done")

    # If the user asked for a paths file, save that to file too.
    if args.paths_out != '':
        print("Saving paths...      ", end="")
        sys.stdout.flush()
        outputPathsFile = open(args.paths_out, 'w')
        for contig in contigs:
            outputPathsFile.write(contig.fullname + '\n')
            outputPathsFile.write(contig.path.getPathsWithLineBreaks())
        print("done")
        



def getArguments():
    parser = argparse.ArgumentParser(description='GetSPAdesContigGraph: a tool for creating a FASTG contig graph from a SPAdes assembly')

    parser.add_argument('graph', help='The assembly_graph.fastg file made by SPAdes')
    parser.add_argument('contigs', help='A contigs or scaffolds fasta file made by SPAdes')
    parser.add_argument('paths', help='The paths file which corresponds to the contigs or scaffolds file')
    parser.add_argument('output', help='The graph file made by this program')
    parser.add_argument('-p', '--paths_out', action='store', help='Save the paths to this file (default: do not save paths)', default='')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-c', '--connection_priority', action='store_true', help='Prioritise graph connections over segment length')
    group.add_argument('-l', '--length_priority', action='store_true', help='Prioritise segment length over graph connections')

    return parser.parse_args()




def isBlastInstalled():
    makeblastdbPath = spawn.find_executable('makeblastdb')
    blastnPath = spawn.find_executable('blastn')
    return makeblastdbPath != None and blastnPath != None


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
def loadGraphLinks(graphFilename):

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

        lineParts = line.split(':')

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

        endingSegments = lineParts[1].split(',')
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



# This function takes a graph filename and returns a dictionary of the graph
# sequences.
# The dictionary key is the graph segment names.
# The dictionary value is the graph segment sequence string.
def loadGraphSequences(graphFilename):

    sequences = {}

    graphFile = open(graphFilename, 'r')

    name = ''
    sequence = ''
    for line in graphFile:

        strippedLine = line.strip()

        # Skip empty lines.
        if len(strippedLine) == 0:
            continue

        # Header lines indicate the start of a new contig.
        if strippedLine[0] == '>':

            # If a sequence is currently stored, save it now.
            if len(name) > 0:
                sequences[name] = sequence
                name = ''
                sequence = ''

            if strippedLine[-1] == ';':
                strippedLine = strippedLine[:-1]

            lineParts = strippedLine.split(':')
            segmentFullName = lineParts[0]
            name = getNumberFromFullSequenceName(segmentFullName)

        # If not a header line, we assume this is a sequence line.
        else:
            sequence += strippedLine

    # Save the last contig.
    if len(name) > 0:
        sequences[name] = sequence

    return sequences





# This function takes a path filename and returns a dictionary.
# The dictionary key is the contig name.
# The dictionary value is a Path object.
def loadPaths(pathFilename):

    paths = {}

    pathFile = open(pathFilename, 'r')

    contigName = ''
    pathSegments = []
    for line in pathFile:

        line = line.strip()

        # Skip empty lines.
        if len(line) == 0:
            continue

        # Lines starting with 'NODE' are the start of a new path
        if len(line) > 3 and line[0:4] == 'NODE':

            # If a path is currently stored, save it now.
            if len(contigName) > 0:
                paths[contigName] = Path(pathSegments)
                contigName = ''
                pathSegments = []

            contigName = line

        # If not a node name line, we assume this is a path line.
        else:
            pathLine = line

            # Replace a semicolon at the end of a line with a placeholder
            # segment called 'gap'
            if pathLine[-1] == ';':
                pathLine = pathLine[0:-1]
                pathLine += ',gap'

            # Add the path to the path list
            if len(pathLine) > 0:
                pathSegments.extend(pathLine.split(','))

    # Save the last contig.
    if len(contigName) > 0:
        paths[contigName] = Path(pathSegments)

    return paths





def getNumberFromFullSequenceName(fullSequenceName):
    number = fullSequenceName.split('_')[1]
    if fullSequenceName[-1] == "'":
        number += '-'
    else:
        number += '+'
    return number




def getOppositeSequenceNumber(number):
    if number[-1] == '+':
        return number[0:-1] + '-'
    else:
        return number[0:-1] + '+'



# This function adds the path information to the contig objects, so each contig
# knows its graph path.
def addPathsToContigs(contigs, paths):
    for contig in contigs:
        contig.addPath(paths[contig.fullname])




# This function uses the contents of the contig paths to add link
# information to the contigs.
def addLinksToContigs(contigs, links):

    # For each contig, we take the last graph segment, find the segments that
    # it leads to, and then find the contigs which start with that next
    # segment.  These make up the links to the current contig.
    for contig1 in contigs:
        endingSegment = contig1.getEndingSegment()
        followingSegments = links[endingSegment]
        linkedContigs = []
        for followingSegment in followingSegments:
            for contig2 in contigs:
                startingSegment = contig2.getStartingSegment()
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
    revCompContig = Contig(revCompContigFullname, revCompSequence)

    # Get the path's reverse complement.
    revCompPath = getReverseComplementPath(contig.path)
    revCompContig.path = revCompPath

    return revCompContig


def getReverseComplementPath(path):
    revSegmentList = []
    for segment in reversed(path.segmentList):
        if segment == 'gap':
            revSegmentList.append(segment)
        else:
            revSegmentList.append(getOppositeSequenceNumber(segment))
    return Path(revSegmentList)




def getReverseComplement(forwardSequence):

    reverseComplement = ''
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



# This function splits contigs as necessary to maintain all graph connections.
# Specifically, it looks for graph segments which are connected to the end of
# one contig and occur in the middle of a second contig.  In such cases, the
# second contig is split to allow for the connection.
def splitContigs(contigs, links, segmentSequences, graphOverlap):

    # # Create a reverse links dictionary.
    # reverseLinks = {}
    # for start, ends in links.iteritems():
    #     for end in ends:
    #         if end not in reverseLinks:
    #             reverseLinks[end] = []
    #         reverseLinks[end].append(start)

    # # Compile lists of all segments which reside on contigs dead ends.
    # deadEndEndSegments = []
    # deadEndStartSegments = []
    # for contig in contigs:
    #     if not contig.linkedContigs:
    #         deadEndEnd = contig.getEndingSegment()
    #         deadEndEndSegments.append(deadEndEnd)
    #         deadEndStartSegments.append(getOppositeSequenceNumber(deadEndEnd))

    # # Find all graph segments which are connected to these dead end segments.
    # # These will need to be on contig ends, to allow for these connections.
    # segmentsWhichMustBeOnContigEnds = []
    # for segment in deadEndStartSegments:
    #     if segment in reverseLinks:
    #         segmentsWhichMustBeOnContigEnds.extend(reverseLinks[segment])
    # segmentsWhichMustBeOnContigStarts = []
    # for segment in deadEndEndSegments:
    #     if segment in links:
    #         segmentsWhichMustBeOnContigStarts.extend(links[segment])

    # Find all missing links.  A missing link is defined as a link in the
    # assembly graph which is not represented somewhere in the contig graph.
    # 'Represented somewhere' includes both within a contig and between two
    # linked contigs.
    missingLinks = []
    for start, ends in links.iteritems():
        for end in ends:
            if (not linkIsInContigs(start, end, contigs)) and (not linkIsBetweenContigs(start, end, contigs)):
                missingLinks.append((start,end))

    # In order for these links to be present in the graph, we need to split
    # contigs such that the start segments of missing links are on the ends
    # of contigs and the end segments of missing links are on the starts of
    # contigs.
    segmentsWhichMustBeOnContigEnds = []
    segmentsWhichMustBeOnContigStarts = []
    for missingLink in missingLinks:
        segmentsWhichMustBeOnContigEnds.append(missingLink[0])
        segmentsWhichMustBeOnContigStarts.append(missingLink[1])

    # Now split contigs, as necessary.
    newPositiveContigs = []
    nextContigNumber = 1
    for contig in contigs:

        # We only split positive contigs and will make the negative complement
        # after we are done.
        if not contig.isPositive():
            continue

        # This will contain the locations at which the contig must be split.
        # It is a list of integers which are indices for path segments that
        # must be at the start of the contig.
        splitPoints = []
        for segment in segmentsWhichMustBeOnContigStarts:
            splitPoints.extend(contig.findSegmentLocations(segment))
        for segment in segmentsWhichMustBeOnContigEnds:
            splitPoints.extend(contig.findSegmentLocationsPlusOne(segment))

        # Remove duplicates and sort split points.
        splitPoints = sorted(list(set(splitPoints)))

        # If the first split point is zero, then remove it, as there is no need
        # to split a contig at its start.
        if splitPoints and splitPoints[0] == 0:
            splitPoints = splitPoints[1:]

        # If the last split point is one past the end, then remove it.
        if splitPoints and splitPoints[-1] == contig.getSegmentCount():
            splitPoints = splitPoints[:-1]

        # If there are splits to be done, then we make the new contigs!
        if splitPoints:
            contig.determineAllSegmentLocations(segmentSequences, graphOverlap)

            for splitPoint in reversed(splitPoints):
                contigPart1, contigPart2 = splitContig(contig, splitPoint, nextContigNumber)
                newPositiveContigs.append(contigPart2)
                contig = contigPart1
                nextContigNumber += 1

            newPositiveContigs.append(contig)

        # If there weren't any split points, then we don't have to split the
        # contig.
        else:
            newPositiveContigs.append(contig)

        # At this point, the last contig added will have a number of 0, so we
        # need to renumber it.
        newPositiveContigs[-1].renumber(nextContigNumber)
        nextContigNumber += 1

    # Now we make the reverse complements for all of our new contigs.
    newContigs = []
    for contig in newPositiveContigs:
        newContigs.append(contig)
        newContigs.append(getReverseComplementContig(contig))

    return newContigs



# This function takes a contig and returns two contigs, split at the split
# point.  The split point is an index for the segment for the segment in the
# contig's path.
def splitContig(contig, splitPoint, nextContigNumber):

    # The indices are bit confusing here, as the contigCoordinates are 1-based
    # with an inclusive end (because that's how BLAST does it).  To get to
    # 0-based and exclusive end (for Python), we subtract one from the start.
    newContig1Path = Path(contig.path.segmentList[:splitPoint])
    newContig1PathContigCoordinates = contig.path.contigCoordinates[:splitPoint]
    newContig1SeqStart = newContig1PathContigCoordinates[0][0] - 1
    newContig1SeqEnd = newContig1PathContigCoordinates[-1][1]
    newContig1Sequence = contig.sequence[newContig1SeqStart:newContig1SeqEnd]

    newContig2Path = Path(contig.path.segmentList[splitPoint:])
    newContig2PathContigCoordinates = contig.path.contigCoordinates[splitPoint:]
    newContig2SeqStart = newContig2PathContigCoordinates[0][0] - 1
    newContig2SeqEnd = newContig2PathContigCoordinates[-1][1]
    newContig2Sequence = contig.sequence[newContig2SeqStart:newContig2SeqEnd]

    # Give the next contig number to the second piece, as the first one may
    # have to be split further.  It will be renumber, if necessary, later.
    newContig1Name = 'NODE_0_length_' + str(len(newContig1Sequence)) + '_cov_' + str(contig.cov)
    newContig2Name = 'NODE_' + str(nextContigNumber) + '_length_' + str(len(newContig2Sequence)) + '_cov_' + str(contig.cov)

    # The first contig is the one that will potentially be split again, so it
    # still needs to have contig coordinates in its path.
    newContig1 = Contig(newContig1Name, newContig1Sequence)
    newContig1.addPath(newContig1Path)
    newContig1.path.contigCoordinates = newContig1PathContigCoordinates

    newContig2 = Contig(newContig2Name, newContig2Sequence)
    newContig2.addPath(newContig2Path)

    return newContig1, newContig2




# This function tests a single overlap between two sequences.
def doesOverlapWork(s1, s2, overlap):
    endOfS1 = s1[-overlap:]
    startOfS2 = s2[:overlap]
    return endOfS1 == startOfS2


# This function tries the possibleOverlaps in the given list and returns
# those that work.
def getPossibleOverlaps(s1, s2, possibleOverlaps):

    newPossibleOverlaps = []

    for possibleOverlap in possibleOverlaps:
        if doesOverlapWork(s1, s2, possibleOverlap):
            newPossibleOverlaps.append(possibleOverlap)

    return newPossibleOverlaps


# This function figures out what the graph overlap size is.
def getGraphOverlap(links, segmentSequences):

    if len(links) == 0:
        return 0

    # Determine the shortest segment in the graph, as this will be the maximum
    # possible overlap.
    segmentLengths = []
    for key, value in segmentSequences.items():
        segmentLengths.append(len(value))
    shortestSegmentSequence = min(segmentLengths)

    # Now we loop through each overlap looking at the segment pairs.
    possibleOverlaps = range(1, shortestSegmentSequence + 1)
    for start, ends in links.iteritems():
        s1 = segmentSequences[start]
        for endingSegment in ends:
            s2 = segmentSequences[endingSegment]
            possibleOverlaps = getPossibleOverlaps(s1, s2, possibleOverlaps)

            # If no overlaps work, then we return 0.
            # This shouldn't happen, as every SPAdes graph should have overlaps.
            if len(possibleOverlaps) == 0:
                return 0

            # If only one overlap works, then that's our answer!
            if len(possibleOverlaps) == 1:
                return possibleOverlaps[0]

    # If the code gets here, that means we have tried every segment pair and
    # there are still multiple possible overlaps.  This shouldn't happen in
    # anything but tiny graphs or seriously pathological cases.
    print("Error: failed to correctly determine graph overlap", file=sys.stderr)
    return 0



def saveSequenceToFastaFile(sequence, sequenceName, filename):
    fastaFile = open(filename, 'w')
    fastaFile.write('>' + sequenceName)
    fastaFile.write('\n')
    fastaFile.write(sequence)
    fastaFile.write('\n')

# This function looks through all contigs to see if the link is present in
# their paths.
def linkIsInContigs(start, end, contigs):
    for contig in contigs:
        if contig.path.containsLink(start, end):
            return True
    return False

# This function looks through all pairs of contigs to see if the link exists
# between their ends.
def linkIsBetweenContigs(start, end, contigs):
    for contig1 in contigs:
        for contig2 in contig1.linkedContigs:
            if contig1.getEndingSegment() == start and contig2.getStartingSegment() == end:
                return True
    return False



# This function returns a reduced list of contigs which has no cases where one
# contig is contained entirely within another.
def removeDuplicateContigs(contigs):

    # Sort the contigs from big to small.  This ensures that containing contigs
    # will be encountered before contained contigs.
    sortedContigs = sorted(contigs, key=lambda contig: len(contig.sequence), reverse=True)

    contigsNoDuplicates = []
    for contig1 in contigs:
        for contig2 in contigsNoDuplicates:
            if contig2.contains(contig1):
                break
        else:
            contigsNoDuplicates.append(contig1)

    return contigsNoDuplicates












# This class holds a contig: its name, sequence and length.
class Contig:

    def __init__(self, name, sequence):
        self.fullname = name
        nameParts = name.split('_')
        self.number = int(nameParts[1])
        covString = nameParts[5]
        if covString[-1] == "'":
            covString = covString[:-1]
        self.cov = float(covString)
        self.sequence = sequence
        self.linkedContigs = [] # filled by addLinkedContigs
        self.path = Path()

    def __str__(self):
        return self.fullname

    def __repr__(self):
        return self.fullname

    def renumber(self, newNumber):
        positive = self.isPositive()
        self.number = newNumber
        self.fullname = 'NODE_' + str(newNumber) + '_length_' + str(len(self.sequence)) + '_cov_' + str(self.cov)
        if not positive:
            self.fullname += "'"

    def addPath(self, path):
        self.path = path

    def getStartingSegment(self):
        return self.path.getFirstSegment()

    def getEndingSegment(self):
        return self.path.getLastSegment()

    def addLinkedContigs(self, linkedContigs):
        self.linkedContigs = linkedContigs

    def isPositive(self):
        return self.fullname[-1] != "'"

    # This function produces a FASTG header for the contig, with links and a
    # line break at the end.
    def getHeaderWithLinks(self):
        headerWithEdges = '>' + self.fullname

        if len(self.linkedContigs) > 0:
            headerWithEdges += ':'
            for linkedContig in self.linkedContigs:
                headerWithEdges += linkedContig.fullname + ','
            headerWithEdges = headerWithEdges[0:-1]

        headerWithEdges += ';\n'
        return headerWithEdges

    def getSequenceWithLineBreaks(self):
        sequenceRemaining = self.sequence
        returnSequence = ''
        while len(sequenceRemaining) > 60:
            returnSequence += sequenceRemaining[0:60] + '\n'
            sequenceRemaining = sequenceRemaining[60:]
        returnSequence += sequenceRemaining + '\n'
        return returnSequence

    def getSegmentCount(self):
        return self.path.getSegmentCount()

    def findSegmentLocations(self, s):
        return self.path.findSegmentLocations(s)

    def findSegmentLocationsPlusOne(self, s):
        segmentLocations = self.path.findSegmentLocations(s)
        return [x+1 for x in segmentLocations]

    # This function determines the start and end coordinates of each of the
    # contig's segments, in contig sequence coordinates.  This information is
    # necessary before we can split a contig.
    def determineAllSegmentLocations(self, segmentSequences, graphOverlap):

        # Create a temporary directory for doing BLAST searches.
        if not os.path.exists('GetSPAdesContigGraph-temp'):
            os.makedirs('GetSPAdesContigGraph-temp')

        saveSequenceToFastaFile(self.sequence, self.fullname, 'GetSPAdesContigGraph-temp/contig.fasta')

        # Create a BLAST database for the contig.
        makeblastdbCommand = ['makeblastdb', '-dbtype', 'nucl', '-in', 'GetSPAdesContigGraph-temp/contig.fasta']
        p = subprocess.Popen(makeblastdbCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()

        # TO DO: CHECK THAT makesblastdb ran okay

        for i in range(len(self.path.segmentList)):
            segment = self.path.segmentList[i]

            # Don't deal with assembly gaps just yet - we'll give them contig
            # start/end coordinates after we've finished with the real
            # segments.
            if segment == "gap":
                continue

            segmentSequence = segmentSequences[segment]
            segmentLength = len(segmentSequence)

            saveSequenceToFastaFile(segmentSequence, segment, 'GetSPAdesContigGraph-temp/segment.fasta')

            # BLAST for the segment in the contig
            sys.stdout.flush()
            blastnCommand = ['blastn', '-task', 'blastn', '-db', 'GetSPAdesContigGraph-temp/contig.fasta', '-query', 'GetSPAdesContigGraph-temp/segment.fasta', '-outfmt', '6 length pident sstart send qstart qend']
            p = subprocess.Popen(blastnCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()

            # TO DO: CHECK THAT blastn ran okay

            # Save the alignments in Python objects.
            alignmentStrings = out.splitlines()
            blastAlignments = []
            for alignmentString in alignmentStrings:
                alignment = BlastAlignment(alignmentString, segmentLength)
                blastAlignments.append(alignment)

            # Only keep alignments that are very high identity and cover the
            # vast majority of the query.  The matches will usually be perfect,
            # but if SPAdes was run with MismatchCorrector, then there can be
            # some differences.
            filteredAlignments = []
            for blastAlignment in blastAlignments:
                fractionQueryPresent = blastAlignment.getAlignmentQueryLength() / segmentLength
                if fractionQueryPresent >= 0.999 and blastAlignment.percentIdentity >= 99.9:
                    filteredAlignments.append(blastAlignment)
            filteredAlignments.sort()

            segmentOccurrencesInPath = self.path.segmentList.count(segment)
            # TO DO: CHECK HERE TO MAKE SURE THE FILTERED BLAST HITS AND THE SEGMENT OCCURRENCES ARE THE SAME!  IF NOT, I NEED TO INTELLIGENTLY DEAL WITH THE DISCREPANCY!

            # Determine which occurrence of this segment in this path we are
            # currently on, and use that to identify the correct BLAST
            # alignment.
            segmentOccurrencesInPathBefore = self.path.segmentList[:i].count(segment)
            correctBlastAlignment = filteredAlignments[segmentOccurrencesInPathBefore]

            # Use the BLAST alignment to determine the query's start and end
            # coordinates in the contig.
            segmentStartInContig = correctBlastAlignment.getQueryStartInReference()
            segmentEndInContig = correctBlastAlignment.getQueryEndInReference()
            self.path.contigCoordinates[i] = (segmentStartInContig, segmentEndInContig)

            os.remove('GetSPAdesContigGraph-temp/segment.fasta')

        shutil.rmtree('GetSPAdesContigGraph-temp')

        # Now we have to go back and assign contig start/end positions for any
        # gap segments.
        segmentCount = len(self.path.segmentList)
        for i in range(segmentCount):
            segment = self.path.segmentList[i]
            if segment != 'gap':
                continue

            gapStartInContig = 1
            gapEndInContig = len(self.sequence)

            if i > 0:
                gapStartInContig = self.path.contigCoordinates[i - 1][1] - graphOverlap + 1
            if i < segmentCount - 1:
                gapEndInContig = self.path.contigCoordinates[i + 1][0] + graphOverlap - 1

            if gapStartInContig < 1:
                gapStartInContig = 1
            if gapStartInContig > len(self.sequence):
                gapStartInContig = len(self.sequence)

            self.path.contigCoordinates[i] = (gapStartInContig, gapEndInContig)

    # This function returns true if this contig contains the entirety of the
    # other contig.  It will also return true if the two contigs are identical.
    # http://stackoverflow.com/questions/3847386/testing-if-a-list-contains-another-list-with-python
    def contains(self, otherContig):
        small = otherContig.path.segmentList
        big = self.path.segmentList

        for i in xrange(len(big)-len(small)+1):
            for j in xrange(len(small)):
                if big[i+j] != small[j]:
                    break
            else:
                return True

        return False



# This class holds a path: the lists of graph segments making up a contig.
class Path:
    def __init__(self, segmentList = []):
        self.segmentList = segmentList
        self.contigCoordinates = [(0,0) for i in range(len(segmentList))]

    def getFirstSegment(self):
        return self.segmentList[0]

    def getLastSegment(self):
        return self.segmentList[-1]

    def __str__(self):
        return str(self.segmentList) + ', ' + str(self.contigCoordinates) 

    def __repr__(self):
        return str(self.segmentList) + ', ' + str(self.contigCoordinates) 

    def getSegmentCount(self):
        return len(self.segmentList)

    def findSegmentLocations(self, s):
        locations = []
        for i in range(len(self.segmentList)):
            if s == self.segmentList[i]:
                locations.append(i)
        return locations

    def getPathsWithLineBreaks(self):
        output = ''
        for segment in self.segmentList:
            if segment == 'gap':
                output += ';\n'
            else:
                output += segment + ','
        return output[:-1] + '\n'

    # This function looks for a particular link within the path.  It returns
    # true if the link is present, false if not.
    def containsLink(self, start, end):

        # A path must have at least 2 segments to possibly contain the link.
        if len(self.segmentList) < 2:
            return False

        for i in range(len(self.segmentList) - 1):
            s1 = self.segmentList[i]
            s2 = self.segmentList[i + 1]
            if s1 == start and s2 == end:
                return True
        return False






class BlastAlignment:

    def __init__(self, blastString, queryLength):
        blastStringParts = blastString.split("\t")

        self.percentIdentity = float(blastStringParts[1])

        self.referenceStart = int(blastStringParts[2])
        self.referenceEnd = int(blastStringParts[3])

        self.queryStart = int(blastStringParts[4])
        self.queryEnd = int(blastStringParts[5])

        self.queryLength = queryLength

    def getAlignmentQueryLength(self):
        return self.queryEnd - self.queryStart + 1

    def getQueryStartInReference(self):
        queryMissingBasesAtStart = self.queryStart - 1
        return self.referenceStart - queryMissingBasesAtStart

    def getQueryEndInReference(self):
        queryMissingBasesAtEnd = self.queryLength - self.queryEnd
        return self.referenceEnd + queryMissingBasesAtEnd


    def __lt__(self, other):
        return self.getQueryStartInReference() < other.getQueryStartInReference()

    def __str__(self):
        return "reference: " + str(self.referenceStart) + " to " + str(self.referenceEnd) + \
               ", query: " + str(self.queryStart) + " to " + str(self.queryEnd) + \
               ", identity: " + str(self.percentIdentity) + "%"

    def __repr__(self):
        return "reference: " + str(self.referenceStart) + " to " + str(self.referenceEnd) + \
               ", query: " + str(self.queryStart) + " to " + str(self.queryEnd) + \
               ", identity: " + str(self.percentIdentity) + "%"






# Standard boilerplate to call the main() function to begin the program.
if __name__ == '__main__':
    main()