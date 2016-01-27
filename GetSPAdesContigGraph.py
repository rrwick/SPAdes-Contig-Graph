#!/usr/bin/env python
"""
GetSPAdesContigGraph

This is a tool to combine the assembly graph and contigs made by the SPAdes
assembler.  For more information, go to:
https://github.com/rrwick/GetSPAdesContigGraph

Author: Ryan Wick
email: rrwick@gmail.com

"""
from __future__ import division
from __future__ import print_function
import sys
import os
import argparse
import shutil
import subprocess
from distutils import spawn

def main():
    args = get_arguments()

    links = load_graph_links(args.graph)
    contigs = load_contigs(args.contigs)
    paths = load_paths(args.paths, links)
    build_graph(contigs, paths, links)

    if args.connection_priority:
        check_for_blast()
        segment_sequences, segment_depths = load_graph_sequences(args.graph)
        graph_overlap = get_graph_overlap(links, segment_sequences)
        contigs = split_contigs(contigs, links, segment_sequences, graph_overlap)
        contigs = merge_contigs(contigs, graph_overlap)
        recalculate_contig_depths(contigs, segment_sequences, segment_depths, graph_overlap)
        renumber_contigs(contigs)

    save_graph_to_file(contigs, args.output)
    if args.paths_out:
        save_paths_to_file(contigs, args.paths_out)

def get_arguments():
    parser = argparse.ArgumentParser(description='GetSPAdesContigGraph: a tool for creating FASTG contig graphs from SPAdes assemblies')
    parser.add_argument('graph', help='The assembly_graph.fastg file made by SPAdes')
    parser.add_argument('contigs', help='A contigs or scaffolds fasta file made by SPAdes')
    parser.add_argument('paths', help='The paths file which corresponds to the contigs or scaffolds file')
    parser.add_argument('output', help='The graph file made by this program')
    parser.add_argument('-p', '--paths_out', action='store', help='Save the paths to this file (default: do not save paths)', default='')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-c', '--connection_priority', action='store_true', help='Prioritise graph connections over segment length')
    group.add_argument('-l', '--length_priority', action='store_true', help='Prioritise segment length over graph connections')
    return parser.parse_args()

def save_graph_to_file(contigs, graph_filename):
    print('Saving graph.......... ', end='')
    sys.stdout.flush()
    output_file = open(graph_filename, 'w')
    for contig in contigs:
        output_file.write(contig.getHeaderWithLinks())
        output_file.write(contig.getSequenceWithLineBreaks())
    print('done')

def save_paths_to_file(contigs, paths_filename):
    print('Saving paths.......... ', end='')
    sys.stdout.flush()
    output_paths_file = open(paths_filename, 'w')
    for contig in contigs:
        output_paths_file.write(contig.fullname + '\n')
        output_paths_file.write(contig.path.getPathsWithLineBreaks())
    print('done')

def check_for_blast():
    makeblastdb_path = spawn.find_executable('makeblastdb')
    blastn_path = spawn.find_executable('blastn')
    blast_installed = (makeblastdb_path != None and blastn_path != None)
    if not blast_installed:
        print('Error: could not find BLAST program', file=sys.stderr)
        quit()

def load_contigs(contig_filename):
    print('Loading contigs....... ', end='')
    sys.stdout.flush()
    try:
        contigs = load_contigs_2(contig_filename)
    except Exception:
        print('\nError: could not load ' + contig_filename, file=sys.stderr)
        quit()
    print('done')
    return contigs




# This function takes a contig filename and returns a list of Contig objects.
# It expects a file of just foward contigs, but it creates both forward and
# reverse complement Contig objects.
def load_contigs_2(contig_filename):

    contigs = []

    contigFile = open(contig_filename, 'r')

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





def load_graph_links(graph_filename):
    print('Loading graph......... ', end='')
    sys.stdout.flush()
    try:
        links = load_graph_links_2(graph_filename)
    except Exception:
        print('\nError: could not load ' + graph_filename, file=sys.stderr)
        quit()
    print('done')
    return links


# This function takes a graph filename and returns a dictionary of the graph
# links.
# The dictionary key is the starting graph segment.
# The dictionary value is a list of the ending graph segments.
def load_graph_links_2(graph_filename):

    links = {}
    graphFile = open(graph_filename, 'r')

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



def build_graph(contigs, paths, links):
    # Add the paths to each contig object, so each contig knows its graph path,
    # and add the links to each contig object, turning the contigs into a
    # graph.
    print('Building graph........ ', end='')
    sys.stdout.flush()
    addPathsToContigs(contigs, paths)
    addLinksToContigs(contigs, links, True, False)
    print('done')




def load_graph_sequences(graph_filename):
    try:
        segmentSequences, segmentDepths = load_graph_sequences_2(graph_filename)
    except Exception:
        print('\nError: could not determine graph sequences and depths', file=sys.stderr)
        quit()
    return segmentSequences, segmentDepths



# This function takes a graph filename and returns two dictionaries:
#  Graph sequences:
#     The dictionary key is the graph segment names.
#     The dictionary value is the graph segment sequence string.
#  Graph depth:
#     The dictionary key is the graph segment names.
#     The dictionary value is the graph segment read depth (cov).
def load_graph_sequences_2(graph_filename):

    sequences = {}
    depths = {}

    graphFile = open(graph_filename, 'r')

    name = ''
    sequence = ''
    depth = 1.0

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
                depths[name] = depth

                name = ''
                sequence = ''
                depth = 1.0

            if strippedLine[-1] == ';':
                strippedLine = strippedLine[:-1]

            lineParts = strippedLine.split(':')
            segmentFullName = lineParts[0]
            name, depth = getNumberAndDepthFromFullSequenceName(segmentFullName)

        # If not a header line, we assume this is a sequence line.
        else:
            sequence += strippedLine

    # Save the last contig.
    if len(name) > 0:
        sequences[name] = sequence
        depths[name] = depth

    return sequences, depths




def load_paths(pathFilename, links):
    print('Loading paths......... ', end='')
    sys.stdout.flush()
    try:
        paths = load_paths_2(pathFilename, links)
    except Exception:
        print('\nError: could not load ' + pathFilename, file=sys.stderr)
        quit()
    print('done')
    return paths




# This function takes a path filename and returns a dictionary.
# The dictionary key is the contig name.
# The dictionary value is a Path object.
def load_paths_2(pathFilename, links):

    paths = {}

    pathFile = open(pathFilename, 'r')

    contigNumber = ''
    pathSegments = []
    for line in pathFile:

        line = line.strip()

        # Skip empty lines.
        if len(line) == 0:
            continue

        # Lines starting with 'NODE' are the start of a new path
        if len(line) > 3 and line[0:4] == 'NODE':

            # If a path is currently stored, save it now.
            if len(contigNumber) > 0:
                paths[contigNumber] = Path(pathSegments)
                contigNumber = ''
                pathSegments = []

            positive = line[-1] != "'"
            lineParts = line.split('_')
            contigNumber = lineParts[1]
            if positive:
                contigNumber += '+'
            else:
                contigNumber += '-'

        # If not a node name line, we assume this is a path line.
        else:
            pathLine = line

            # Replace a semicolon at the end of a line with a placeholder
            # segment called 'gap_SEG' where SEG will be the name of the
            # preceding segment.
            if pathLine[-1] == ';':
                pathLine = pathLine[0:-1]
                if contigNumber.endswith('+'):
                    pathLine += ',gap_POS'
                else:
                    pathLine += ',gap_NEG'

            # Add the path to the path list
            if len(pathLine) > 0:
                pathSegments.extend(pathLine.split(','))

    # Save the last contig.
    if len(contigNumber) > 0:
        paths[contigNumber] = Path(pathSegments)

    # Now we go through the paths and rename our gap segments.  Anything named
    # gap_POS will get the name of the preceding path segment.  Anything named
    # gap_NEG will get the name of the following path segment.
    for path in paths.itervalues():
        for i in range(len(path.segmentList)):
            segment = path.segmentList[i]
            if segment == 'gap_POS':
                path.segmentList[i] = 'gap_' + path.segmentList[i-1]
            elif segment == 'gap_NEG':
                path.segmentList[i] = 'gap_' + path.segmentList[i+1]

    # Now we must go through all of the paths we just made and add any links
    # for new gap segments.
    for path in paths.itervalues():
        for i in range(len(path.segmentList) - 1):
            j = i + 1
            s1 = path.segmentList[i]
            s2 = path.segmentList[j]
            if s1.startswith('gap') or s2.startswith('gap'):
                if s1 in links:
                    links[s1].append(s2)
                else:
                    links[s1] = [s2]

    return paths




def getNumberFromFullSequenceName(fullSequenceName):
    nameParts = fullSequenceName.split('_')
    number = nameParts[1]
    if fullSequenceName[-1] == "'":
        number += '-'
    else:
        number += '+'

    return number


def getNumberAndDepthFromFullSequenceName(fullSequenceName):
    nameParts = fullSequenceName.split('_')
    number = nameParts[1]
    if fullSequenceName[-1] == "'":
        number += '-'
    else:
        number += '+'

    depthString = nameParts[5]
    if depthString[-1] == "'":
        depthString = depthString[:-1]
    depth = float(depthString)

    return number, depth




def getOppositeSequenceNumber(number):
    if number[-1] == '+':
        return number[0:-1] + '-'
    else:
        return number[0:-1] + '+'



# This function adds the path information to the contig objects, so each contig
# knows its graph path.
def addPathsToContigs(contigs, paths):
    for contig in contigs:
        contig.addPath(paths[contig.getNumberWithSign()])




# This function uses the contents of the contig paths to add link
# information to the contigs.
def addLinksToContigs(contigs, links, clear, deadEndsOnly):

    if clear:
        for contig in contigs:
            contig.outgoingLinkedContigs = []
            contig.incomingLinkedContigs = []

    # If we are only adding links for dead ends, then we make sets to easily
    # tell if a contig has no connections.
    if deadEndsOnly:
        noOutgoingLinks = set()
        noIncomingLinks = set()
        for contig in contigs:
            if not contig.outgoingLinkedContigs:
                noOutgoingLinks.add(contig)
            if not contig.incomingLinkedContigs:
                noIncomingLinks.add(contig)

    # For each contig, we take the last graph segment, find the segments that
    # it leads to, and then find the contigs which start with that next
    # segment.  These make up the links to the current contig.
    for contig1 in contigs:
        endingSegment = contig1.getEndingSegment()
        followingSegments = links[endingSegment]

        for followingSegment in followingSegments:
            for contig2 in contigs:
                startingSegment = contig2.getStartingSegment()
                if followingSegment == startingSegment:
                    if deadEndsOnly and contig1 not in noOutgoingLinks and contig2 not in noIncomingLinks:
                        continue
                    contig1.outgoingLinkedContigs.append(contig2)
                    contig2.incomingLinkedContigs.append(contig1)

    # This process can create duplicate linked contigs, so we now go through
    # them to remove the duplicates.
    for contig in contigs:
        contig.outgoingLinkedContigs = list(set(contig.outgoingLinkedContigs))
        contig.incomingLinkedContigs = list(set(contig.incomingLinkedContigs))






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




def split_contigs(contigs, links, segmentSequences, graph_overlap):
    print('Splitting contigs..... ', end='')
    sys.stdout.flush()
    contigs = splitContigs2(contigs, links, segmentSequences, graph_overlap)
    addLinksToContigs(contigs, links, False, True)
    print('done')
    return contigs



# This function splits contigs as necessary to maintain all graph connections.
# Specifically, it looks for graph segments which are connected to the end of
# one contig and occur in the middle of a second contig.  In such cases, the
# second contig is split to allow for the connection.
def splitContigs2(contigs, links, segmentSequences, graph_overlap):

    # Find all missing links.  A missing link is defined as a link in the
    # assembly graph which is not represented somewhere in the contig graph.
    # 'Represented somewhere' includes both within a contig and between two
    # linked contigs.
    linksInContigs = set()
    for contig in contigs:
        linksInContig = contig.getLinksInThisContigAndToOtherContigs()
        for link in linksInContig:
            linksInContigs.add(link)
    allGraphLinks = set()
    for start, ends in links.iteritems():
        for end in ends:
            allGraphLinks.add((start, end))
    missingLinks = []
    for graphLink in allGraphLinks:
        if graphLink not in linksInContigs:
            missingLinks.append(graphLink)

    # In order for these links to be present in the graph, we need to split
    # contigs such that the start segments of missing links are on the ends
    # of contigs and the end segments of missing links are on the starts of
    # contigs.
    segmentsWhichMustBeOnContigEnds = []
    segmentsWhichMustBeOnContigStarts = []
    for missingLink in missingLinks:
        segmentsWhichMustBeOnContigEnds.append(missingLink[0])
        segmentsWhichMustBeOnContigStarts.append(missingLink[1])

    # Create a reverse links dictionary.
    reverseLinks = {}
    for start, ends in links.iteritems():
        for end in ends:
            if end not in reverseLinks:
                reverseLinks[end] = []
            reverseLinks[end].append(start)

    # Compile lists of all segments which reside on contigs dead ends.
    deadEndEndSegments = []
    deadEndStartSegments = []
    for contig in contigs:
        if not contig.outgoingLinkedContigs:
            deadEndEnd = contig.getEndingSegment()
            deadEndEndSegments.append(deadEndEnd)
            deadEndStartSegments.append(getOppositeSequenceNumber(deadEndEnd))

    # Find all graph segments which are connected to these dead end segments.
    # These will need to be on contig ends, to allow for these connections.
    segmentsWhichMustBeOnContigEnds = []
    for segment in deadEndStartSegments:
        if segment in reverseLinks:
            segmentsWhichMustBeOnContigEnds.extend(reverseLinks[segment])
    segmentsWhichMustBeOnContigStarts = []
    for segment in deadEndEndSegments:
        if segment in links:
            segmentsWhichMustBeOnContigStarts.extend(links[segment])

    # Remove any duplicates from the segments lists just made.
    segmentsWhichMustBeOnContigStarts = list(set(segmentsWhichMustBeOnContigStarts))
    segmentsWhichMustBeOnContigEnds = list(set(segmentsWhichMustBeOnContigEnds))

    # Before we split the contigs, we need to remember all of the links present
    # so they can be remade.
    linksBeforeSplits = {}
    for contig in contigs:
        start = contig.getNumberWithSign()
        ends = []
        for outgoingLinkedContig in contig.outgoingLinkedContigs:
            ends.append(outgoingLinkedContig.getNumberWithSign())
        linksBeforeSplits[start] = ends

    # Now split contigs, as necessary.
    newPositiveContigs = []
    nextContigNumber = 1
    oldNumbersToNewNumbers = {}
    for contig in contigs:

        # We only split positive contigs and will make the negative complement
        # after we are done.
        if not contig.isPositive():
            continue

        # We need to keep track of the mapping from old contig numbers to new
        # numbers, as this will be needed for reconnecting contigs.
        contigNumber = contig.number
        oldNumbersToNewNumbers[contigNumber] = []

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
            contig.determineAllSegmentLocations(segmentSequences, graph_overlap)

            for splitPoint in reversed(splitPoints):
                contigPart1, contigPart2 = splitContig(contig, splitPoint, nextContigNumber)

                # If the split was successful, then both contigPart1 and
                # contigPart2 are new Contig objects.  If unsuccessful, then
                # they are None.
                if contigPart1 is not None:
                    newPositiveContigs.append(contigPart2)
                    oldNumbersToNewNumbers[contigNumber] = [contigPart2.number] + oldNumbersToNewNumbers[contigNumber]
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
        oldNumbersToNewNumbers[contigNumber] = [newPositiveContigs[-1].number] + oldNumbersToNewNumbers[contigNumber]

    # Now we make the reverse complements for all of our new contigs.
    newContigs = []
    for contig in newPositiveContigs:
        newContigs.append(contig)
        newContigs.append(getReverseComplementContig(contig))

    # Now we have to put together the links in new contig numbers.  First we
    # Create the internal links between parts of a split contig.
    linksAfterSplits = {}
    for newNum in oldNumbersToNewNumbers.itervalues():
        for i in range(len(newNum) - 1):
            start = str(newNum[i]) + '+'
            end = str(newNum[i+1]) + '+'
            linksAfterSplits[start] = [end]
            start = str(newNum[i+1]) + '-'
            end = str(newNum[i]) + '-'
            linksAfterSplits[start] = [end]

    # Add the external links between contigs.  To do this we need to
    # translate old contig numbers to new contig numbers.
    for start, ends in linksBeforeSplits.iteritems():
        startSign = start[-1]
        startNum = int(start[:-1])
        if startSign == '+':
            newStartNum = oldNumbersToNewNumbers[startNum][-1]
        else:
            newStartNum = oldNumbersToNewNumbers[startNum][0]
        newStart = str(newStartNum) + startSign
        newEnds = []
        for end in ends:
            endSign = end[-1]
            endNum = int(end[:-1])
            if endSign == '+':
                newEndNum = oldNumbersToNewNumbers[endNum][0]
            else:
                newEndNum = oldNumbersToNewNumbers[endNum][-1]
            newEnd = str(newEndNum) + endSign
            newEnds.append(newEnd)
        linksAfterSplits[newStart] = newEnds

    # Also make the links in reverse direction.
    reverseLinksAfterSplits = {}
    for start, ends in linksAfterSplits.iteritems():
        for end in ends:
            if end in reverseLinksAfterSplits:
                reverseLinksAfterSplits[end].append(start)
            else:
                reverseLinksAfterSplits[end] = [start]

    # Now we add the links back to our new contigs.
    newContigDict = {}
    for contig in newContigs:
        newContigDict[contig.getNumberWithSign()] = contig
    for contig in newContigs:
        contig.incomingLinkedContigs = []
        contig.outgoingLinkedContigs = []
        contigNum = contig.getNumberWithSign()
        if contigNum in linksAfterSplits:
            for outgoingNum in linksAfterSplits[contigNum]:
                contig.outgoingLinkedContigs.append(newContigDict[outgoingNum])
        if contigNum in reverseLinksAfterSplits:
            for incomingNum in reverseLinksAfterSplits[contigNum]:
                contig.incomingLinkedContigs.append(newContigDict[incomingNum])

    return newContigs



# This function takes a contig and returns two contigs, split at the split
# point.  The split point is an index for the segment for the segment in the
# contig's path.
def splitContig(contig, splitPoint, nextContigNumber):

    # Determine the new contig paths.
    newContig1Path = Path(contig.path.segmentList[:splitPoint])
    newContig2Path = Path(contig.path.segmentList[splitPoint:])

    # Get the coordinates for the new contig paths.
    newContig1PathContigCoordinates = contig.path.contigCoordinates[:splitPoint]
    newContig2PathContigCoordinates = contig.path.contigCoordinates[splitPoint:]
    newContig1FirstSegmentCoordinates = newContig1PathContigCoordinates[0]
    newContig1LastSegmentCoordinates = newContig1PathContigCoordinates[-1]
    newContig2FirstSegmentCoordinates = newContig2PathContigCoordinates[0]
    newContig2LastSegmentCoordinates = newContig2PathContigCoordinates[-1]

    # Check to see if any of the important coordinates are absent, as will be
    # the case if this program was unable to find the segment in the contig.
    # In this case, the split fails.
    if newContig1FirstSegmentCoordinates[0] is None or \
       newContig1LastSegmentCoordinates[1] is None or \
       newContig2FirstSegmentCoordinates[0] is None or \
       newContig2LastSegmentCoordinates[1] is None:
        return None, None

    # Determine exact coordinate for the new segments.
    # The indices are bit confusing here, as the contigCoordinates are 1-based
    # with an inclusive end (because that's how BLAST does it).  To get to
    # 0-based and exclusive end (for Python), we subtract one from the start.
    newContig1SeqStart = newContig1FirstSegmentCoordinates[0] - 1
    newContig1SeqEnd = newContig1LastSegmentCoordinates[1]
    newContig1Sequence = contig.sequence[newContig1SeqStart:newContig1SeqEnd]
    newContig2SeqStart = newContig2FirstSegmentCoordinates[0] - 1
    newContig2SeqEnd = newContig2LastSegmentCoordinates[1]
    newContig2Sequence = contig.sequence[newContig2SeqStart:newContig2SeqEnd]

    # Give the next contig number to the second piece, as the first one may
    # have to be split further.  It will be renumbered, if necessary, later.
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
def get_graph_overlap(links, segmentSequences):

    if len(links) == 0:
        return 0

    # Determine the shortest segment in the graph, as this will be the maximum
    # possible overlap.
    segmentLengths = []
    for sequence in segmentSequences.itervalues():
        segmentLengths.append(len(sequence))
    shortestSegmentSequence = min(segmentLengths)

    # Now we loop through each overlap looking at the segment pairs.
    possibleOverlaps = range(1, shortestSegmentSequence + 1)
    for start, ends in links.iteritems():
        if start.startswith('gap'):
            continue
        s1 = segmentSequences[start]
        for endingSegment in ends:
            if endingSegment.startswith('gap'):
                continue
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
    print('Error: failed to correctly determine graph overlap', file=sys.stderr)
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
        for contig2 in contig1.outgoingLinkedContigs:
            if contig1.getEndingSegment() == start and contig2.getStartingSegment() == end:
                return True
    return False


def merge_contigs(contigs, graph_overlap):
    print('Merging contigs....... ', end='')
    sys.stdout.flush()
    contigs = mergeIdenticalContigs(contigs)
    contigs = mergeLinearRuns(contigs, graph_overlap)
    print('done')
    return contigs



# This function looks for contigs which are made of the exact same graph
# segments as each other.
def mergeIdenticalContigs(contigs):

    contigGroups = []

    # Group contigs into collections with the exact same segment list.
    for contig in contigs:
        for contigGroup in contigGroups:
            if contig.path.segmentList == contigGroup[0].path.segmentList:
                contigGroup.append(contig)
                break
        else:
            contigGroups.append([contig])

    # For the first contig in each group, give it the linked contigs of its
    # entire group.
    oldNumberToNewContig = {}
    newContigs = []
    for contigGroup in contigGroups:
        allIncomingLinkedContigNumbers = set()
        allOutgoingLinkedContigNumbers = set()
        firstContigInGroup = contigGroup[0]
        for contig in contigGroup:
            for incomingLinkedContig in contig.incomingLinkedContigs:
                allIncomingLinkedContigNumbers.add(incomingLinkedContig.getNumberWithSign())
            for outgoingLinkedContig in contig.outgoingLinkedContigs:
                allOutgoingLinkedContigNumbers.add(outgoingLinkedContig.getNumberWithSign())
            oldNumberToNewContig[contig.getNumberWithSign()] = firstContigInGroup
        firstContigInGroup.incomingLinkedContigs = list(allIncomingLinkedContigNumbers)
        firstContigInGroup.outgoingLinkedContigs = list(allOutgoingLinkedContigNumbers)
        newContigs.append(firstContigInGroup)

    # Now for each of the new contigs, we must convert the linked contig lists
    # from numbers to actual contigs.
    for contig in newContigs:
        newIncomingLinkedContigs = set()
        newOutgoingLinkedContigs = set()
        for incomingLinkedContigOldNumber in contig.incomingLinkedContigs:
            newIncomingLinkedContigs.add(oldNumberToNewContig[incomingLinkedContigOldNumber])
        for outgoingLinkedContigOldNumber in contig.outgoingLinkedContigs:
            newOutgoingLinkedContigs.add(oldNumberToNewContig[outgoingLinkedContigOldNumber])
        contig.incomingLinkedContigs = list(newIncomingLinkedContigs)
        contig.outgoingLinkedContigs = list(newOutgoingLinkedContigs)

    return newContigs


# This function simplifies the graph by merging simple linear runs together.
def mergeLinearRuns(contigs, graph_overlap):

    mergeHappened = True
    while mergeHappened:

        contigDict = {}
        for contig in contigs:
            contigDict[contig.getNumberWithSign()] = contig

        for contig in contigs:

            # Make sure that this contig has one downstream contig and that
            # only has one upstream contig.
            if len(contig.outgoingLinkedContigs) != 1:
                continue
            nextContig = contig.outgoingLinkedContigs[0]
            if len(nextContig.incomingLinkedContigs) != 1 or \
               nextContig.incomingLinkedContigs[0] != contig:
                continue

            # Make sure this contig isn't just looping back to itself.
            if contig == nextContig:
                continue

            # Make sure that the reverse complement contigs are also properly
            # simple and linear.
            revCompContig = contigDict[getOppositeSequenceNumber(contig.getNumberWithSign())]
            revCompNextContig = contigDict[getOppositeSequenceNumber(nextContig.getNumberWithSign())]
            if len(revCompNextContig.outgoingLinkedContigs) != 1 or \
               len(revCompContig.incomingLinkedContigs) != 1 or \
               revCompNextContig.outgoingLinkedContigs[0] != revCompContig or \
               revCompContig.incomingLinkedContigs[0] != revCompNextContig:
                continue

            # Make sure that the contig is not simply looping back onto its own
            # reverse complement.
            if nextContig == revCompContig or contig == revCompNextContig:
                continue

            # If the code got here, then we can merge contig and nextContig
            # (and their reverse complements).
            contigs = mergeTwoContigs(contigs, graph_overlap, contig, nextContig, revCompContig, revCompNextContig)
            break

        else:
            mergeHappened = False

    return contigs



def mergeTwoContigs(allContigs, graph_overlap, contig1, contig2, contig1RevComp, contig2RevComp):

    largestContigNumber = 0
    for contig in allContigs:
        contigNum = contig.number
        if contigNum > largestContigNumber:
            largestContigNumber = contigNum
    mergedContigNum = largestContigNumber + 1

    mergedContigSequence = contig1.sequence[:-graph_overlap] + contig2.sequence
    mergedContigSequenceRevComp = contig2RevComp.sequence + contig1RevComp.sequence[graph_overlap:]

    mergedContigName = "NODE_" + str(mergedContigNum) + "_length_" + str(len(mergedContigSequence)) + "_cov_1.0"
    mergedContigNameRevComp = "NODE_" + str(mergedContigNum) + "_length_" + str(len(mergedContigSequenceRevComp)) + "_cov_1.0'"

    mergedContig = Contig(mergedContigName, mergedContigSequence)
    mergedContigRevComp = Contig(mergedContigNameRevComp, mergedContigSequenceRevComp)

    # Copy the connections over to the new merged contigs.
    mergedContig.incomingLinkedContigs = contig1.incomingLinkedContigs
    mergedContig.outgoingLinkedContigs = contig2.outgoingLinkedContigs
    mergedContigRevComp.incomingLinkedContigs = contig2RevComp.incomingLinkedContigs
    mergedContigRevComp.outgoingLinkedContigs = contig1RevComp.outgoingLinkedContigs

    # Copy the path segments over to the new merged contigs.
    mergedContig.path.segmentList = contig1.path.segmentList + contig2.path.segmentList
    mergedContigRevComp.path.segmentList = contig2RevComp.path.segmentList + contig1RevComp.path.segmentList

    # Build a new list of contigs, and while we're looping through, we can fix
    # up any links to the merged contig.
    newContigs = []
    for contig in allContigs:
        if contig == contig1 or contig == contig2 or contig == contig1RevComp or contig == contig2RevComp:
            continue

        newIncomingLinkedContigs = []
        for incomingLinkedContig in contig.incomingLinkedContigs:
            if incomingLinkedContig == contig2:
                newIncomingLinkedContigs.append(mergedContig)
            elif incomingLinkedContig == contig1RevComp:
                newIncomingLinkedContigs.append(mergedContigRevComp)
            else:
                newIncomingLinkedContigs.append(incomingLinkedContig)
        contig.incomingLinkedContigs = newIncomingLinkedContigs

        newOutgoingLinkedContigs = []
        for outgoingLinkedContig in contig.outgoingLinkedContigs:
            if outgoingLinkedContig == contig1:
                newOutgoingLinkedContigs.append(mergedContig)
            elif outgoingLinkedContig == contig2RevComp:
                newOutgoingLinkedContigs.append(mergedContigRevComp)
            else:
                newOutgoingLinkedContigs.append(outgoingLinkedContig)
        contig.outgoingLinkedContigs = newOutgoingLinkedContigs

        newContigs.append(contig)

    newContigs.append(mergedContig)
    newContigs.append(mergedContigRevComp)

    return newContigs



# This function is run after removeDuplicateContigs.  It goes through all of
# the contigs' links and removes ones that are no longer valid (because their
# destination contig has been removed as a duplicate).
def removeBogusLinks(contigs):

    contigNumSet = set()
    for contig in contigs:
        contigNumSet.add(contig.number)

    for contig in contigs:
        newIncomingLinkedContigs = []
        newOutgoingLinkedContigs = []

        for incomingLinkedContig in contig.incomingLinkedContigs:
            if incomingLinkedContig.number in contigNumSet:
                newIncomingLinkedContigs.append(incomingLinkedContig)
        for outgoingLinkedContig in contig.outgoingLinkedContigs:
            if outgoingLinkedContig.number in contigNumSet:
                newOutgoingLinkedContigs.append(outgoingLinkedContig)

        contig.incomingLinkedContigs = newIncomingLinkedContigs
        contig.outgoingLinkedContigs = newOutgoingLinkedContigs



# This function recalculates contig depths using the depths of the graph
# segments which make up the contigs.
# For this function I tried to copy what SPAdes does: it seems to define a
# contig's depth as the weighted average of the graph segments in the contig.
# The weight for the average is the segment length minus the overlap.
# Notably, SPAdes does not seem to divide up the available depth for a segment
# which appears in multiple places in the contigs, so I don't do that either.
def recalculate_contig_depths(contigs, sequences, depths, graph_overlap):

    print('Calculating depth..... ', end='')
    sys.stdout.flush()

    for contig in contigs:
        totalLength = 0
        totalDepthTimesLength = 0.0

        for segment in contig.path.segmentList:
            if segment.startswith('gap'):
                continue

            # Get the segment depth and length.  In some odd cases, SPAdes does
            # not save both segments in a complementary pair, so we may have to
            # look for the complementary segment.
            if segment in depths:
                depth = depths[segment]
            else:
                depth = depths[getOppositeSequenceNumber(segment)]
            adjustedDepth = depth

            if segment in sequences:
                length = len(sequences[segment])
            else:
                length = len(sequences[getOppositeSequenceNumber(segment)])
            adjustedLength = length - graph_overlap

            totalLength += adjustedLength
            totalDepthTimesLength += adjustedDepth * adjustedLength

        if totalLength > 0:
            finalDepth = totalDepthTimesLength / totalLength
            contig.cov = finalDepth
            contig.rebuildFullName()

    print('done')




def renumber_contigs(contigs):

    print('Renumbering contigs... ', end='')
    sys.stdout.flush()

    # Sort from contigs big to small, so contig 1 is the largest.
    positiveContigs = []
    for contig in contigs:
        if contig.isPositive():
            positiveContigs.append(contig)
    sortedContigs = sorted(positiveContigs, key=lambda contig: len(contig.sequence), reverse=True)

    # Create the new number mapping.
    oldNumToNewNum = {}
    i = 1
    for contig in sortedContigs:
        oldNumToNewNum[contig.number] = i
        i += 1

    # Assign new numbers and create complement contigs.
    for contig in contigs:
        contig.renumber(oldNumToNewNum[contig.number])

    # Sort the contigs using their new numbers.
    contigs.sort(key=lambda contig: (contig.number, contig.getSign()))

    print('done')


# This function takes blast alignments and returns the best one.
# 'Best' is defined first as covering the entirety of the alignment, being high
# identity and having an appropriate start position.
def getBestBlastAlignment(blastAlignments, segmentLength, expectedReferenceStart):

    if not blastAlignments:
        return None

    # Sort by alignment length and identity.
    sortedAlignments = sorted(blastAlignments, key=lambda alignment: (alignment.getQueryLength(), alignment.percentIdentity), reverse=True)

    # If we have an expected reference start to work with, then we find the
    # first alignment in the list which occurs very near the expected reference
    # start.
    if expectedReferenceStart is not None:
        for alignment in sortedAlignments:
            discrepancy = abs(alignment.getQueryStartInReference() - expectedReferenceStart)
            if discrepancy < 5:
                return alignment

    # If we don't have an expected reference start to work with, then we are
    # more limited.  We can only give a result if there is a single full length
    # alignment.
    else:
        fullLengthAlignments = []
        for alignment in sortedAlignments:
            fractionPresent = alignment.getQueryLength() / segmentLength
            if fractionPresent == 1.0:
                fullLengthAlignments.append(alignment)
        if len(fullLengthAlignments) == 1:
            return fullLengthAlignments[0]

    # If the code got here, then no good match was found.
    return None












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
        self.outgoingLinkedContigs = []
        self.incomingLinkedContigs = []
        self.path = Path()

    def __str__(self):
        return self.fullname

    def __repr__(self):
        return self.fullname

    def renumber(self, newNumber):
        self.number = newNumber
        self.rebuildFullName()

    def rebuildFullName(self):
        positive = self.isPositive()
        self.fullname = 'NODE_' + str(self.number) + '_length_' + str(len(self.sequence)) + '_cov_' + str(self.cov)
        if not positive:
            self.fullname += "'"

    def getNumberWithSign(self):
        num = str(self.number)
        if self.isPositive():
            num += '+'
        else:
            num += '-'
        return num

    def getSign(self):
        if self.isPositive():
            return '+'
        else:
            return '-'

    def addPath(self, path):
        self.path = path

    def getStartingSegment(self):
        return self.path.getFirstSegment()

    def getEndingSegment(self):
        return self.path.getLastSegment()

    def isPositive(self):
        return self.fullname[-1] != "'"

    # This function produces a FASTG header for the contig, with links and a
    # line break at the end.
    def getHeaderWithLinks(self):
        headerWithEdges = '>' + self.fullname

        if len(self.outgoingLinkedContigs) > 0:
            headerWithEdges += ':'
            for linkedContig in self.outgoingLinkedContigs:
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
    def determineAllSegmentLocations(self, segmentSequences, graph_overlap):

        # Create a temporary directory for doing BLAST searches.
        if not os.path.exists('GetSPAdesContigGraph-temp'):
            os.makedirs('GetSPAdesContigGraph-temp')

        saveSequenceToFastaFile(self.sequence, self.fullname, 'GetSPAdesContigGraph-temp/contig.fasta')

        # Create a BLAST database for the contig.
        makeblastdbCommand = ['makeblastdb', '-dbtype', 'nucl', '-in', 'GetSPAdesContigGraph-temp/contig.fasta']
        p = subprocess.Popen(makeblastdbCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()

        # Check that makeblastdb ran okay
        if len(err) > 0:
            print('\nmakeblastdb encountered an error:', file=sys.stderr)
            print(err, file=sys.stderr)
            quit()

        # This value tracks where we expect the next segment to start, in
        # contig coordinates.  When it is set to None, that means we don't
        # know.
        expectedSegmentStartInContig = 1

        for i in range(len(self.path.segmentList)):
            segment = self.path.segmentList[i]

            # Don't deal with assembly gaps just yet - we'll give them contig
            # start/end coordinates after we've finished with the real
            # segments.
            if segment.startswith('gap'):
                expectedSegmentStartInContig = None
                continue

            segmentSequence = segmentSequences[segment]
            segmentLength = len(segmentSequence)

            saveSequenceToFastaFile(segmentSequence, segment, 'GetSPAdesContigGraph-temp/segment.fasta')

            # BLAST for the segment in the contig
            sys.stdout.flush()
            blastnCommand = ['blastn', '-task', 'blastn', '-db', 'GetSPAdesContigGraph-temp/contig.fasta', '-query', 'GetSPAdesContigGraph-temp/segment.fasta', '-outfmt', '6 length pident sstart send qstart qend']
            p = subprocess.Popen(blastnCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()

            # Check that blastn ran okay
            if len(err) > 0:
                print('\nblastn encountered an error:', file=sys.stderr)
                print(err, file=sys.stderr)
                quit()

            # Save the alignments in Python objects.
            alignmentStrings = out.splitlines()
            blastAlignments = []
            for alignmentString in alignmentStrings:
                alignment = BlastAlignment(alignmentString, segmentLength)
                blastAlignments.append(alignment)

            # Get the alignment (if any) that best matches the segment sequence
            # and position.
            bestAlignment = getBestBlastAlignment(blastAlignments, segmentLength, expectedSegmentStartInContig)

            # If we found an alignment, use it to determine the query's start
            # and end coordinates in the contig.
            if bestAlignment is not None:
                segmentStartInContig = bestAlignment.getQueryStartInReference()
                segmentEndInContig = bestAlignment.getQueryEndInReference()

            # If we failed to find an alignment, then we don't have segment
            # coordinates and will be unable to split the contig at this point.
            else:
                segmentStartInContig = None
                segmentEndInContig = None

            self.path.contigCoordinates[i] = (segmentStartInContig, segmentEndInContig)

            # Update the expected segment start for the next segment.
            # Hopefully we can use the alignment to do this precisely.
            if bestAlignment is not None:
                expectedSegmentStartInContig = (bestAlignment.referenceEnd - graph_overlap + 1)

            # If there is no alignment with which we can predict the next
            # segment start, we can guess using the segment length.
            elif expectedSegmentStartInContig is not None:
                expectedSegmentStartInContig += (segmentLength - graph_overlap)

            os.remove('GetSPAdesContigGraph-temp/segment.fasta')

        shutil.rmtree('GetSPAdesContigGraph-temp')

        # Now we have to go back and assign contig start/end positions for any
        # gap segments.
        segmentCount = len(self.path.segmentList)
        for i in range(segmentCount):
            segment = self.path.segmentList[i]
            if not segment.startswith('gap'):
                continue

            gapStartInContig = 1
            gapEndInContig = len(self.sequence)

            if i > 0:
                previousSegmentEnd = self.path.contigCoordinates[i - 1][1]
                if previousSegmentEnd is not None:
                    gapStartInContig = previousSegmentEnd - graph_overlap + 1
                    if gapStartInContig < 1:
                        gapStartInContig = 1
                else:
                    gapStartInContig = None

            if i < segmentCount - 1:
                nextSegmentStart = self.path.contigCoordinates[i + 1][0]
                if nextSegmentStart is not None:
                    gapEndInContig = nextSegmentStart + graph_overlap - 1
                    if gapEndInContig > len(self.sequence):
                        gapEndInContig = len(self.sequence)
                else:
                    gapEndInContig = None

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

    # This function returns a list of the links from this contig to any
    # other contig.
    def getLinksToOtherContigs(self):
        linksToOtherContigs = []
        for outgoingLinkedContig in self.outgoingLinkedContigs:
            linksToOtherContigs.append((self.getEndingSegment(), outgoingLinkedContig.getStartingSegment()))
        for incomingLinkedContig in self.incomingLinkedContigs:
            linksToOtherContigs.append((incomingLinkedContig.getEndingSegment(), self.getStartingSegment()))
        return list(set(linksToOtherContigs))

    def getLinksInThisContigAndToOtherContigs(self):
        links = self.path.getAllLinks()
        links.extend(self.getLinksToOtherContigs())
        return list(set(links))






# This class holds a path: the lists of graph segments making up a contig.
class Path:
    def __init__(self, segmentList=[]):
        self.segmentList = segmentList
        self.contigCoordinates = [(0, 0) for i in range(len(segmentList))]

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
            if segment.startswith('gap'):
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

    # This function returns a list of tuple for each link in the path.
    def getAllLinks(self):
        links = []
        for i in range(len(self.segmentList) - 1):
            s1 = self.segmentList[i]
            s2 = self.segmentList[i + 1]
            links.append((s1, s2))
        return links







class BlastAlignment:

    def __init__(self, blastString, queryLength):
        blastStringParts = blastString.split('\t')

        self.percentIdentity = float(blastStringParts[1])

        self.referenceStart = int(blastStringParts[2])
        self.referenceEnd = int(blastStringParts[3])

        self.queryStart = int(blastStringParts[4])
        self.queryEnd = int(blastStringParts[5])

        self.queryLength = queryLength

    def getQueryLength(self):
        return self.queryEnd - self.queryStart + 1

    def getReferenceLength(self):
        return self.referenceEnd - self.referenceStart + 1

    def getQueryStartInReference(self):
        queryMissingBasesAtStart = self.queryStart - 1
        return self.referenceStart - queryMissingBasesAtStart

    def getQueryEndInReference(self):
        queryMissingBasesAtEnd = self.queryLength - self.queryEnd
        return self.referenceEnd + queryMissingBasesAtEnd


    def __lt__(self, other):
        return self.getQueryStartInReference() < other.getQueryStartInReference()

    def __str__(self):
        return 'reference: ' + str(self.referenceStart) + ' to ' + str(self.referenceEnd) + \
               ', query: ' + str(self.queryStart) + ' to ' + str(self.queryEnd) + \
               ', identity: ' + str(self.percentIdentity) + '%'

    def __repr__(self):
        return 'reference: ' + str(self.referenceStart) + ' to ' + str(self.referenceEnd) + \
               ', query: ' + str(self.queryStart) + ' to ' + str(self.queryEnd) + \
               ', identity: ' + str(self.percentIdentity) + '%'






# Standard boilerplate to call the main() function to begin the program.
if __name__ == '__main__':
    main()
