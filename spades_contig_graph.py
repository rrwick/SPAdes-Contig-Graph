#!/usr/bin/env python
"""
SPAdes Contig Graph

This is a tool to combine the assembly graph and contigs made by the SPAdes
assembler.  For more information, go to:
https://github.com/rrwick/spades_contig_graph

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
    parser = argparse.ArgumentParser(description='SPAdes Contig Graph: a tool for creating FASTG contig graphs from SPAdes assemblies')
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
        output_file.write(contig.get_header_with_links())
        output_file.write(contig.get_sequence_with_line_breaks())
    print('done')

def save_paths_to_file(contigs, paths_filename):
    print('Saving paths.......... ', end='')
    sys.stdout.flush()
    output_paths_file = open(paths_filename, 'w')
    for contig in contigs:
        output_paths_file.write(contig.fullname + '\n')
        output_paths_file.write(contig.path.get_paths_with_line_breaks())
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

def load_contigs_2(contig_filename):
    """
    This function takes a contig filename and returns a list of Contig objects.
    It expects a file of just foward contigs, but it creates both forward and
    reverse complement Contig objects.
    """
    contigs = []
    contig_file = open(contig_filename, 'r')
    name = ''
    sequence = ''

    for line in contig_file:
        line = line.strip()
        if not line:
            continue
        if line[0] == '>': # Header line = start of new contig
            if len(name) > 0: # If a contig is currently in memory...
                contig = Contig(name, sequence)
                contigs.append(contig)
                contigs.append(make_reverse_complement_contig(contig))
                name = ''
                sequence = ''
            name = line[1:]
        else: # Not a header line = sequence
            sequence += line
    if len(name) > 0: # Save the last contig
        contig = Contig(name, sequence)
        contigs.append(contig)
        contigs.append(make_reverse_complement_contig(contig))

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

def load_graph_links_2(graph_filename):
    """
    This function takes a graph filename and returns a dictionary of the graph
    links.
    The dictionary key is the starting graph segment.
    The dictionary value is a list of the ending graph segments.
    """
    links = {}
    graph_file = open(graph_filename, 'r')

    for line in graph_file:
        line = line.strip()
        if not line:
            continue
        if line[0] != '>':
            continue
        if line[-1] == ';':
            line = line[:-1]

        line_parts = line.split(':')
        start = get_number_from_sequence_name(line_parts[0])
        if start not in links:
            links[start] = []
        start_rev_comp = get_opposite_sequence_number(start)
        if start_rev_comp not in links:
            links[start] = []

        if len(line_parts) < 2:
            continue

        ending_segments = line_parts[1].split(',')
        ends = []
        for ending_segment in ending_segments:
            ends.append(get_number_from_sequence_name(ending_segment))

        # Add the links to the dictionary in the forward direction.
        for end in ends:
            if end not in links[start]:
                links[start].append(end)

        # Add the links to the dictionary in the reverse direction.
        for end in ends:
            end_rev_comp = get_opposite_sequence_number(end)
            if end_rev_comp not in links:
                links[end_rev_comp] = []
            if start_rev_comp not in links[end_rev_comp]:
                links[end_rev_comp].append(start_rev_comp)

    return links

def build_graph(contigs, paths, links):
    """
    Add the paths to each contig object, so each contig knows its graph path,
    and add the links to each contig object, turning the contigs into a
    graph.
    """
    print('Building graph........ ', end='')
    sys.stdout.flush()
    add_paths_to_contigs(contigs, paths)
    add_links_to_contigs(contigs, links, True, False)
    print('done')

def load_graph_sequences(graph_filename):
    try:
        segment_sequences, segment_depths = load_graph_sequences_2(graph_filename)
    except Exception:
        print('\nError: could not determine graph sequences and depths', file=sys.stderr)
        quit()
    return segment_sequences, segment_depths

def load_graph_sequences_2(graph_filename):
    """
    This function takes a graph filename and returns two dictionaries:
    Graph sequences:
      The dictionary key is the graph segment names.
      The dictionary value is the graph segment sequence string.
    Graph depth:
      The dictionary key is the graph segment names.
      The dictionary value is the graph segment read depth (cov).
    """
    sequences = {}
    depths = {}

    graph_file = open(graph_filename, 'r')

    name = ''
    sequence = ''
    depth = 1.0

    for line in graph_file:
        line = line.strip()
        if not line:
            continue

        # Header lines indicate the start of a new contig.
        if line[0] == '>':

            # If a sequence is currently stored, save it now.
            if len(name) > 0:
                sequences[name] = sequence
                depths[name] = depth
                name = ''
                sequence = ''
                depth = 1.0

            if line[-1] == ';':
                line = line[:-1]

            line_parts = line.split(':')
            name, depth = get_number_and_depth_from_sequence_name(line_parts[0])

        # If not a header line, we assume this is a sequence line.
        else:
            sequence += line

    # Save the last contig.
    if len(name) > 0:
        sequences[name] = sequence
        depths[name] = depth

    return sequences, depths

def load_paths(path_filename, links):
    print('Loading paths......... ', end='')
    sys.stdout.flush()
    try:
        paths = load_paths_2(path_filename, links)
    except Exception:
        print('\nError: could not load ' + path_filename, file=sys.stderr)
        quit()
    print('done')
    return paths

def load_paths_2(path_filename, links):
    """
    This function takes a path filename and returns a dictionary.
    The dictionary key is the contig name.
    The dictionary value is a Path object.
    """
    paths = {}
    path_file = open(path_filename, 'r')
    contig_number = ''
    path_segments = []

    for line in path_file:
        line = line.strip()
        if not line:
            continue

        # Lines starting with 'NODE' are the start of a new path
        if len(line) > 3 and line[0:4] == 'NODE':

            # If a path is currently stored, save it now.
            if len(contig_number) > 0:
                paths[contig_number] = Path(path_segments)
                contig_number = ''
                path_segments = []

            positive = line[-1] != "'"
            line_parts = line.split('_')
            contig_number = line_parts[1]
            if positive:
                contig_number += '+'
            else:
                contig_number += '-'

        # If not a node name line, we assume this is a path line.
        else:
            path_line = line

            # Replace a semicolon at the end of a line with a placeholder
            # segment called 'gap_SEG' where SEG will be the name of the
            # preceding segment.
            if path_line[-1] == ';':
                path_line = path_line[0:-1]
                if contig_number.endswith('+'):
                    path_line += ',gap_POS'
                else:
                    path_line += ',gap_NEG'

            # Add the path to the path list
            if len(path_line) > 0:
                path_segments.extend(path_line.split(','))

    # Save the last contig.
    if len(contig_number) > 0:
        paths[contig_number] = Path(path_segments)

    # Now we go through the paths and rename our gap segments.  Anything named
    # gap_POS will get the name of the preceding path segment.  Anything named
    # gap_NEG will get the name of the following path segment.
    for path in paths.itervalues():
        for i in range(len(path.segment_list)):
            segment = path.segment_list[i]
            if segment == 'gap_POS':
                path.segment_list[i] = 'gap_' + path.segment_list[i-1]
            elif segment == 'gap_NEG':
                path.segment_list[i] = 'gap_' + path.segment_list[i+1]

    # Now we must go through all of the paths we just made and add any links
    # for new gap segments.
    for path in paths.itervalues():
        for i in range(len(path.segment_list) - 1):
            j = i + 1
            s_1 = path.segment_list[i]
            s_2 = path.segment_list[j]
            if s_1.startswith('gap') or s_2.startswith('gap'):
                if s_1 in links:
                    links[s_1].append(s_2)
                else:
                    links[s_1] = [s_2]

    return paths

def get_number_from_sequence_name(sequence_name):
    name_parts = sequence_name.split('_')
    number = name_parts[1]
    if sequence_name[-1] == "'":
        number += '-'
    else:
        number += '+'
    return number

def get_number_and_depth_from_sequence_name(sequence_name):
    name_parts = sequence_name.split('_')
    number = name_parts[1]
    if sequence_name[-1] == "'":
        number += '-'
    else:
        number += '+'
    depth_string = name_parts[5]
    if depth_string[-1] == "'":
        depth_string = depth_string[:-1]
    depth = float(depth_string)
    return number, depth

def get_opposite_sequence_number(number):
    if number[-1] == '+':
        return number[0:-1] + '-'
    else:
        return number[0:-1] + '+'

def add_paths_to_contigs(contigs, paths):
    """
    This function adds the path information to the contig objects, so each contig
    knows its graph path.
    """
    for contig in contigs:
        contig.path = paths[contig.get_number_with_sign()]

def add_links_to_contigs(contigs, links, clear, dead_ends_only):
    """
    This function uses the contents of the contig paths to add link
    information to the contigs.
    """
    if clear:
        for contig in contigs:
            contig.outgoing_linked_contigs = []
            contig.incoming_linked_contigs = []

    # If we are only adding links for dead ends, then we make sets to easily
    # tell if a contig has no connections.
    if dead_ends_only:
        no_outgoing_links = set()
        no_incoming_links = set()
        for contig in contigs:
            if not contig.outgoing_linked_contigs:
                no_outgoing_links.add(contig)
            if not contig.incoming_linked_contigs:
                no_incoming_links.add(contig)

    # For each contig, we take the last graph segment, find the segments that
    # it leads to, and then find the contigs which start with that next
    # segment.  These make up the links to the current contig.
    for contig_1 in contigs:
        ending_segment = contig_1.path.get_last_segment()
        following_segments = links[ending_segment]

        for following_segment in following_segments:
            for contig_2 in contigs:
                starting_segment = contig_2.path.get_first_segment()
                if following_segment == starting_segment:
                    if dead_ends_only and contig_1 not in no_outgoing_links and contig_2 not in no_incoming_links:
                        continue
                    contig_1.outgoing_linked_contigs.append(contig_2)
                    contig_2.incoming_linked_contigs.append(contig_1)

    # This process can create duplicate linked contigs, so we now go through
    # them to remove the duplicates.
    for contig in contigs:
        contig.outgoing_linked_contigs = list(set(contig.outgoing_linked_contigs))
        contig.incoming_linked_contigs = list(set(contig.incoming_linked_contigs))

def make_reverse_complement_contig(contig):
    """
    This function takes a contig and returns its reverse complement contig.
    """
    rev_comp_contig_fullname = contig.fullname

    # Add or remove the ' as necessary
    if rev_comp_contig_fullname[-1] == "'":
        rev_comp_contig_fullname = rev_comp_contig_fullname[0:-1]
    else:
        rev_comp_contig_fullname = rev_comp_contig_fullname + "'"

    rev_comp_sequence = make_reverse_complement_sequence(contig.sequence)
    rev_comp_contig = Contig(rev_comp_contig_fullname, rev_comp_sequence)
    rev_comp_contig.path = make_reverse_complement_path(contig.path)

    return rev_comp_contig

def make_reverse_complement_path(path):
    rev_segment_list = []
    for segment in reversed(path.segment_list):
        rev_segment_list.append(get_opposite_sequence_number(segment))
    return Path(rev_segment_list)

def make_reverse_complement_sequence(forward_sequence):
    rev_comp = ''
    for i in reversed(range(len(forward_sequence))):
        base = forward_sequence[i]
        if base == 'A': rev_comp += 'T'
        elif base == 'T': rev_comp += 'A'
        elif base == 'G': rev_comp += 'C'
        elif base == 'C': rev_comp += 'G'
        elif base == 'a': rev_comp += 't'
        elif base == 't': rev_comp += 'a'
        elif base == 'g': rev_comp += 'c'
        elif base == 'c': rev_comp += 'g'
        elif base == 'R': rev_comp += 'Y'
        elif base == 'Y': rev_comp += 'R'
        elif base == 'S': rev_comp += 'S'
        elif base == 'W': rev_comp += 'W'
        elif base == 'K': rev_comp += 'M'
        elif base == 'M': rev_comp += 'K'
        elif base == 'r': rev_comp += 'y'
        elif base == 'y': rev_comp += 'r'
        elif base == 's': rev_comp += 's'
        elif base == 'w': rev_comp += 'w'
        elif base == 'k': rev_comp += 'm'
        elif base == 'm': rev_comp += 'k'
        elif base == 'B': rev_comp += 'V'
        elif base == 'D': rev_comp += 'H'
        elif base == 'H': rev_comp += 'D'
        elif base == 'V': rev_comp += 'B'
        elif base == 'b': rev_comp += 'v'
        elif base == 'd': rev_comp += 'h'
        elif base == 'h': rev_comp += 'd'
        elif base == 'v': rev_comp += 'b'
        elif base == 'N': rev_comp += 'N'
        elif base == 'n': rev_comp += 'n'
        elif base == '.': rev_comp += '.'
        elif base == '-': rev_comp += '-'
        elif base == '?': rev_comp += '?'
        else: rev_comp += 'N'
    return rev_comp

def split_contigs(contigs, links, segment_sequences, graph_overlap):
    print('Splitting contigs..... ', end='')
    sys.stdout.flush()
    contigs = split_contigs_2(contigs, links, segment_sequences, graph_overlap)
    add_links_to_contigs(contigs, links, False, True)
    print('done')
    return contigs

def split_contigs_2(contigs, links, segment_sequences, graph_overlap):
    """
    This function splits contigs as necessary to maintain all graph connections.
    Specifically, it looks for graph segments which are connected to the end of
    one contig and occur in the middle of a second contig.  In such cases, the
    second contig is split to allow for the connection.
    """

    # Find all missing links.  A missing link is defined as a link in the
    # assembly graph which is not represented somewhere in the contig graph.
    # 'Represented somewhere' includes both within a contig and between two
    # linked contigs.
    links_in_contigs = set()
    for contig in contigs:
        links_in_contig = contig.get_links_in_this_contig_and_to_other_contigs()
        for link in links_in_contig:
            links_in_contigs.add(link)
    all_graph_links = set()
    for start, ends in links.iteritems():
        for end in ends:
            all_graph_links.add((start, end))
    missing_links = []
    for graph_link in all_graph_links:
        if graph_link not in links_in_contigs:
            missing_links.append(graph_link)

    # In order for these links to be present in the graph, we need to split
    # contigs such that the start segments of missing links are on the ends
    # of contigs and the end segments of missing links are on the starts of
    # contigs.
    segments_which_must_be_on_contig_ends = []
    segments_which_must_be_on_contig_starts = []
    for missing_link in missing_links:
        segments_which_must_be_on_contig_ends.append(missing_link[0])
        segments_which_must_be_on_contig_starts.append(missing_link[1])

    # Create a reverse links dictionary.
    reverse_links = {}
    for start, ends in links.iteritems():
        for end in ends:
            if end not in reverse_links:
                reverse_links[end] = []
            reverse_links[end].append(start)

    # Compile lists of all segments which reside on contigs dead ends.
    dead_end_end_segments = []
    dead_end_start_segments = []
    for contig in contigs:
        if not contig.outgoing_linked_contigs:
            dead_end_end = contig.path.get_last_segment()
            dead_end_end_segments.append(dead_end_end)
            dead_end_start_segments.append(get_opposite_sequence_number(dead_end_end))

    # Find all graph segments which are connected to these dead end segments.
    # These will need to be on contig ends, to allow for these connections.
    segments_which_must_be_on_contig_ends = []
    for segment in dead_end_start_segments:
        if segment in reverse_links:
            segments_which_must_be_on_contig_ends.extend(reverse_links[segment])
    segments_which_must_be_on_contig_starts = []
    for segment in dead_end_end_segments:
        if segment in links:
            segments_which_must_be_on_contig_starts.extend(links[segment])

    # Remove any duplicates from the segments lists just made.
    segments_which_must_be_on_contig_starts = list(set(segments_which_must_be_on_contig_starts))
    segments_which_must_be_on_contig_ends = list(set(segments_which_must_be_on_contig_ends))

    # Before we split the contigs, we need to remember all of the links present
    # so they can be remade.
    links_before_splits = {}
    for contig in contigs:
        start = contig.get_number_with_sign()
        ends = []
        for outgoing_linked_contig in contig.outgoing_linked_contigs:
            ends.append(outgoing_linked_contig.get_number_with_sign())
        links_before_splits[start] = ends

    # Now split contigs, as necessary.
    new_positive_contigs = []
    next_contig_number = 1
    old_nums_to_new_nums = {}
    for contig in contigs:

        # We only split positive contigs and will make the negative complement
        # after we are done.
        if not contig.is_positive():
            continue

        # We need to keep track of the mapping from old contig numbers to new
        # numbers, as this will be needed for reconnecting contigs.
        contig_number = contig.number
        old_nums_to_new_nums[contig_number] = []

        # This will contain the locations at which the contig must be split.
        # It is a list of integers which are indices for path segments that
        # must be at the start of the contig.
        split_points = []
        for segment in segments_which_must_be_on_contig_starts:
            split_points.extend(contig.find_segment_locations(segment))
        for segment in segments_which_must_be_on_contig_ends:
            split_points.extend(contig.find_segment_locations_plus_one(segment))

        # Remove duplicates and sort split points.
        split_points = sorted(list(set(split_points)))

        # If the first split point is zero, then remove it, as there is no need
        # to split a contig at its start.
        if split_points and split_points[0] == 0:
            split_points = split_points[1:]

        # If the last split point is one past the end, then remove it.
        if split_points and split_points[-1] == contig.path.get_segment_count():
            split_points = split_points[:-1]

        # If there are splits to be done, then we make the new contigs!
        if split_points:
            contig.determine_all_segment_locations(segment_sequences, graph_overlap)

            for split_point in reversed(split_points):
                contig_part_1, contig_part_2 = split_contig(contig, split_point, next_contig_number)

                # If the split was successful, then both contig_part_1 and
                # contig_part_2 are new Contig objects.  If unsuccessful, then
                # they are None.
                if contig_part_1 is not None:
                    new_positive_contigs.append(contig_part_2)
                    old_nums_to_new_nums[contig_number] = [contig_part_2.number] + old_nums_to_new_nums[contig_number]
                    contig = contig_part_1
                    next_contig_number += 1

            new_positive_contigs.append(contig)

        # If there weren't any split points, then we don't have to split the
        # contig.
        else:
            new_positive_contigs.append(contig)

        # At this point, the last contig added will have a number of 0, so we
        # need to renumber it.
        new_positive_contigs[-1].renumber(next_contig_number)
        next_contig_number += 1
        old_nums_to_new_nums[contig_number] = [new_positive_contigs[-1].number] + old_nums_to_new_nums[contig_number]

    # Now we make the reverse complements for all of our new contigs.
    new_contigs = []
    for contig in new_positive_contigs:
        new_contigs.append(contig)
        new_contigs.append(make_reverse_complement_contig(contig))

    # Now we have to put together the links in new contig numbers.  First we
    # Create the internal links between parts of a split contig.
    links_after_splits = {}
    for new_num in old_nums_to_new_nums.itervalues():
        for i in range(len(new_num) - 1):
            start = str(new_num[i]) + '+'
            end = str(new_num[i+1]) + '+'
            links_after_splits[start] = [end]
            start = str(new_num[i+1]) + '-'
            end = str(new_num[i]) + '-'
            links_after_splits[start] = [end]

    # Add the external links between contigs.  To do this we need to
    # translate old contig numbers to new contig numbers.
    for start, ends in links_before_splits.iteritems():
        start_sign = start[-1]
        start_num = int(start[:-1])
        if start_sign == '+':
            new_start_num = old_nums_to_new_nums[start_num][-1]
        else:
            new_start_num = old_nums_to_new_nums[start_num][0]
        new_start = str(new_start_num) + start_sign
        new_ends = []
        for end in ends:
            end_sign = end[-1]
            end_num = int(end[:-1])
            if end_sign == '+':
                new_end_num = old_nums_to_new_nums[end_num][0]
            else:
                new_end_num = old_nums_to_new_nums[end_num][-1]
            new_end = str(new_end_num) + end_sign
            new_ends.append(new_end)
        links_after_splits[new_start] = new_ends

    # Also make the links in reverse direction.
    reverse_links_after_splits = {}
    for start, ends in links_after_splits.iteritems():
        for end in ends:
            if end in reverse_links_after_splits:
                reverse_links_after_splits[end].append(start)
            else:
                reverse_links_after_splits[end] = [start]

    # Now we add the links back to our new contigs.
    new_contig_dict = {}
    for contig in new_contigs:
        new_contig_dict[contig.get_number_with_sign()] = contig
    for contig in new_contigs:
        contig.incoming_linked_contigs = []
        contig.outgoing_linked_contigs = []
        contig_num = contig.get_number_with_sign()
        if contig_num in links_after_splits:
            for outgoing_num in links_after_splits[contig_num]:
                contig.outgoing_linked_contigs.append(new_contig_dict[outgoing_num])
        if contig_num in reverse_links_after_splits:
            for incoming_num in reverse_links_after_splits[contig_num]:
                contig.incoming_linked_contigs.append(new_contig_dict[incoming_num])

    return new_contigs

def split_contig(contig, split_point, next_contig_number):
    """
    This function takes a contig and returns two contigs, split at the split
    point.  The split point is an index for the segment for the segment in the
    contig's path.
    """
    # Determine the new contig paths.
    new_contig_1_path = Path(contig.path.segment_list[:split_point])
    new_contig_2_path = Path(contig.path.segment_list[split_point:])

    # Get the coordinates for the new contig paths.
    new_contig_1_path_coordinates = contig.path.contig_coordinates[:split_point]
    new_contig_2_path_coordinates = contig.path.contig_coordinates[split_point:]
    new_contig_1_first_segment_coordinates = new_contig_1_path_coordinates[0]
    new_contig_1_last_segment_coordinates = new_contig_1_path_coordinates[-1]
    new_contig_2_first_segment_coordinates = new_contig_2_path_coordinates[0]
    new_contig_2_last_segment_coordinates = new_contig_2_path_coordinates[-1]

    # Check to see if any of the important coordinates are absent, as will be
    # the case if this program was unable to find the segment in the contig.
    # In this case, the split fails.
    if new_contig_1_first_segment_coordinates[0] is None or \
       new_contig_1_last_segment_coordinates[1] is None or \
       new_contig_2_first_segment_coordinates[0] is None or \
       new_contig_2_last_segment_coordinates[1] is None:
        return None, None

    # Determine exact coordinate for the new segments.
    # The indices are bit confusing here, as the contig coordinates are 1-based
    # with an inclusive end (because that's how BLAST does it).  To get to
    # 0-based and exclusive end (for Python), we subtract one from the start.
    new_contig_1_seq_start = new_contig_1_first_segment_coordinates[0] - 1
    new_contig_1_seq_end = new_contig_1_last_segment_coordinates[1]
    new_contig_1_sequence = contig.sequence[new_contig_1_seq_start:new_contig_1_seq_end]
    new_contig_2_seq_start = new_contig_2_first_segment_coordinates[0] - 1
    new_contig_2_seq_end = new_contig_2_last_segment_coordinates[1]
    new_contig_2_sequence = contig.sequence[new_contig_2_seq_start:new_contig_2_seq_end]

    # Give the next contig number to the second piece, as the first one may
    # have to be split further.  It will be renumbered, if necessary, later.
    new_contig_1_name = 'NODE_0_length_' + str(len(new_contig_1_sequence)) + '_cov_' + str(contig.cov)
    new_contig_2_name = 'NODE_' + str(next_contig_number) + '_length_' + str(len(new_contig_2_sequence)) + '_cov_' + str(contig.cov)

    # The first contig is the one that will potentially be split again, so it
    # still needs to have contig coordinates in its path.
    new_contig_1 = Contig(new_contig_1_name, new_contig_1_sequence)
    new_contig_1.path = new_contig_1_path
    new_contig_1.path.contig_coordinates = new_contig_1_path_coordinates

    new_contig_2 = Contig(new_contig_2_name, new_contig_2_sequence)
    new_contig_2.path = new_contig_2_path

    return new_contig_1, new_contig_2

def does_overlap_work(s_1, s_2, overlap):
    """
    Test a single overlap between two sequences.
    """
    return s_1[-overlap:] == s_2[:overlap]

def get_possible_overlaps(s_1, s_2, possible_overlaps):
    """
    Try the possible overlaps in the given list and return those that work.
    """
    new_possible_overlaps = []
    for possible_overlap in possible_overlaps:
        if does_overlap_work(s_1, s_2, possible_overlap):
            new_possible_overlaps.append(possible_overlap)
    return new_possible_overlaps

def get_graph_overlap(links, segment_sequences):
    """
    Figure out the graph overlap size.
    """
    if not links:
        return 0

    # Determine the shortest segment in the graph, as this will be the maximum
    # possible overlap.
    segment_lengths = []
    for sequence in segment_sequences.itervalues():
        segment_lengths.append(len(sequence))
    shortest_segment_sequence = min(segment_lengths)

    # Now we loop through each overlap looking at the segment pairs.
    possible_overlaps = range(1, shortest_segment_sequence + 1)
    for start, ends in links.iteritems():
        if start.startswith('gap'):
            continue
        s_1 = segment_sequences[start]
        for ending_segment in ends:
            if ending_segment.startswith('gap'):
                continue
            s_2 = segment_sequences[ending_segment]
            possible_overlaps = get_possible_overlaps(s_1, s_2, possible_overlaps)

            # If no overlaps work, then we return 0.
            # This shouldn't happen, as every SPAdes graph should have overlaps.
            if len(possible_overlaps) == 0:
                return 0

            # If only one overlap works, then that's our answer!
            if len(possible_overlaps) == 1:
                return possible_overlaps[0]

    # If the code gets here, that means we have tried every segment pair and
    # there are still multiple possible overlaps.  This shouldn't happen in
    # anything but tiny graphs or seriously pathological cases.
    print('Error: failed to correctly determine graph overlap', file=sys.stderr)
    return 0

def save_sequence_to_fasta_file(sequence, sequence_name, filename):
    fasta_file = open(filename, 'w')
    fasta_file.write('>' + sequence_name)
    fasta_file.write('\n')
    fasta_file.write(sequence)
    fasta_file.write('\n')

def merge_contigs(contigs, graph_overlap):
    print('Merging contigs....... ', end='')
    sys.stdout.flush()
    contigs = merge_identical_contigs(contigs)
    contigs = merge_linear_runs(contigs, graph_overlap)
    print('done')
    return contigs

def merge_identical_contigs(contigs):
    """
    Find contigs which are made of the exact same graph segments as each other
    and merge them.
    """
    contig_groups = []

    # Group contigs into collections with the exact same segment list.
    for contig in contigs:
        for contig_group in contig_groups:
            if contig.path.segment_list == contig_group[0].path.segment_list:
                contig_group.append(contig)
                break
        else:
            contig_groups.append([contig])

    # For the first contig in each group, give it the linked contigs of its
    # entire group.
    old_num_to_new_contig = {}
    new_contigs = []
    for contig_group in contig_groups:
        all_incoming_linked_contig_numbers = set()
        all_outgoing_linked_contig_numbers = set()
        first_contig_in_group = contig_group[0]
        for contig in contig_group:
            for incoming_linked_contig in contig.incoming_linked_contigs:
                all_incoming_linked_contig_numbers.add(incoming_linked_contig.get_number_with_sign())
            for outgoing_linked_contig in contig.outgoing_linked_contigs:
                all_outgoing_linked_contig_numbers.add(outgoing_linked_contig.get_number_with_sign())
            old_num_to_new_contig[contig.get_number_with_sign()] = first_contig_in_group
        first_contig_in_group.incoming_linked_contigs = list(all_incoming_linked_contig_numbers)
        first_contig_in_group.outgoing_linked_contigs = list(all_outgoing_linked_contig_numbers)
        new_contigs.append(first_contig_in_group)

    # Now for each of the new contigs, we must convert the linked contig lists
    # from numbers to actual contigs.
    for contig in new_contigs:
        new_incoming_linked_contigs = set()
        new_outgoing_linked_contigs = set()
        for incoming_linked_contig_old_number in contig.incoming_linked_contigs:
            new_incoming_linked_contigs.add(old_num_to_new_contig[incoming_linked_contig_old_number])
        for outgoing_linked_contig_old_number in contig.outgoing_linked_contigs:
            new_outgoing_linked_contigs.add(old_num_to_new_contig[outgoing_linked_contig_old_number])
        contig.incoming_linked_contigs = list(new_incoming_linked_contigs)
        contig.outgoing_linked_contigs = list(new_outgoing_linked_contigs)

    return new_contigs

def merge_linear_runs(contigs, graph_overlap):
    """
    Find and merge simple linear runs of contigs.
    """
    merge_happened = True
    while merge_happened:

        contig_dict = {}
        for contig in contigs:
            contig_dict[contig.get_number_with_sign()] = contig

        for contig in contigs:

            # Make sure that this contig has one downstream contig and that
            # only has one upstream contig.
            if len(contig.outgoing_linked_contigs) != 1:
                continue
            next_contig = contig.outgoing_linked_contigs[0]
            if len(next_contig.incoming_linked_contigs) != 1 or \
               next_contig.incoming_linked_contigs[0] != contig:
                continue

            # Make sure this contig isn't just looping back to itself.
            if contig == next_contig:
                continue

            # Make sure that the reverse complement contigs are also properly
            # simple and linear.
            rev_comp_contig = contig_dict[get_opposite_sequence_number(contig.get_number_with_sign())]
            rev_comp_next_contig = contig_dict[get_opposite_sequence_number(next_contig.get_number_with_sign())]
            if len(rev_comp_next_contig.outgoing_linked_contigs) != 1 or \
               len(rev_comp_contig.incoming_linked_contigs) != 1 or \
               rev_comp_next_contig.outgoing_linked_contigs[0] != rev_comp_contig or \
               rev_comp_contig.incoming_linked_contigs[0] != rev_comp_next_contig:
                continue

            # Make sure that the contig is not simply looping back onto its own
            # reverse complement.
            if next_contig == rev_comp_contig or contig == rev_comp_next_contig:
                continue

            # If the code got here, then we can merge contig and next_contig
            # (and their reverse complements).
            contigs = merge_two_contigs(contigs, graph_overlap, contig, next_contig, rev_comp_contig, rev_comp_next_contig)
            break

        else:
            merge_happened = False

    return contigs

def merge_two_contigs(all_contigs, graph_overlap, contig_1, contig_2, contig_1_rev_comp, contig_2_rev_comp):
    largest_contig_number = 0
    for contig in all_contigs:
        contig_num = contig.number
        if contig_num > largest_contig_number:
            largest_contig_number = contig_num
    merged_contig_num = largest_contig_number + 1

    merged_contig_sequence = contig_1.sequence[:-graph_overlap] + contig_2.sequence
    merged_contig_sequence_rev_comp = contig_2_rev_comp.sequence + contig_1_rev_comp.sequence[graph_overlap:]

    merged_contig_name = "NODE_" + str(merged_contig_num) + "_length_" + str(len(merged_contig_sequence)) + "_cov_1.0"
    merged_contig_name_rev_comp = "NODE_" + str(merged_contig_num) + "_length_" + str(len(merged_contig_sequence_rev_comp)) + "_cov_1.0'"

    merged_contig = Contig(merged_contig_name, merged_contig_sequence)
    merged_contig_rev_comp = Contig(merged_contig_name_rev_comp, merged_contig_sequence_rev_comp)

    # Copy the connections over to the new merged contigs.
    merged_contig.incoming_linked_contigs = contig_1.incoming_linked_contigs
    merged_contig.outgoing_linked_contigs = contig_2.outgoing_linked_contigs
    merged_contig_rev_comp.incoming_linked_contigs = contig_2_rev_comp.incoming_linked_contigs
    merged_contig_rev_comp.outgoing_linked_contigs = contig_1_rev_comp.outgoing_linked_contigs

    # Copy the path segments over to the new merged contigs.
    merged_contig.path.segment_list = contig_1.path.segment_list + contig_2.path.segment_list
    merged_contig_rev_comp.path.segment_list = contig_2_rev_comp.path.segment_list + contig_1_rev_comp.path.segment_list

    # Build a new list of contigs, and while we're looping through, we can fix
    # up any links to the merged contig.
    new_contigs = []
    for contig in all_contigs:
        if contig == contig_1 or contig == contig_2 or contig == contig_1_rev_comp or contig == contig_2_rev_comp:
            continue

        new_incoming_linked_contigs = []
        for incoming_linked_contig in contig.incoming_linked_contigs:
            if incoming_linked_contig == contig_2:
                new_incoming_linked_contigs.append(merged_contig)
            elif incoming_linked_contig == contig_1_rev_comp:
                new_incoming_linked_contigs.append(merged_contig_rev_comp)
            else:
                new_incoming_linked_contigs.append(incoming_linked_contig)
        contig.incoming_linked_contigs = new_incoming_linked_contigs

        new_outgoing_linked_contigs = []
        for outgoing_linked_contig in contig.outgoing_linked_contigs:
            if outgoing_linked_contig == contig_1:
                new_outgoing_linked_contigs.append(merged_contig)
            elif outgoing_linked_contig == contig_2_rev_comp:
                new_outgoing_linked_contigs.append(merged_contig_rev_comp)
            else:
                new_outgoing_linked_contigs.append(outgoing_linked_contig)
        contig.outgoing_linked_contigs = new_outgoing_linked_contigs

        new_contigs.append(contig)

    new_contigs.append(merged_contig)
    new_contigs.append(merged_contig_rev_comp)

    return new_contigs


def recalculate_contig_depths(contigs, sequences, depths, graph_overlap):
    """
    This function recalculates contig depths using the depths of the graph
    segments which make up the contigs.
    For this function I tried to copy what SPAdes does: it seems to define a
    contig's depth as the weighted average of the graph segments in the contig.
    The weight for the average is the segment length minus the overlap.
    Notably, SPAdes does not seem to divide up the available depth for a
    segment which appears in multiple places in the contigs, so I don't do that
    either.
    """
    print('Calculating depth..... ', end='')
    sys.stdout.flush()

    for contig in contigs:
        total_length = 0
        total_depth_times_length = 0.0

        for segment in contig.path.segment_list:
            if segment.startswith('gap'):
                continue

            # Get the segment depth and length.  In some odd cases, SPAdes does
            # not save both segments in a complementary pair, so we may have to
            # look for the complementary segment.
            if segment in depths:
                depth = depths[segment]
            else:
                depth = depths[get_opposite_sequence_number(segment)]
            adjusted_depth = depth

            if segment in sequences:
                length = len(sequences[segment])
            else:
                length = len(sequences[get_opposite_sequence_number(segment)])
            adjusted_length = length - graph_overlap

            total_length += adjusted_length
            total_depth_times_length += adjusted_depth * adjusted_length

        if total_length > 0:
            final_depth = total_depth_times_length / total_length
            contig.cov = final_depth
            contig.rebuild_full_name()

    print('done')

def renumber_contigs(contigs):
    print('Renumbering contigs... ', end='')
    sys.stdout.flush()

    # Sort from contigs big to small, so contig 1 is the largest.
    positive_contigs = []
    for contig in contigs:
        if contig.is_positive():
            positive_contigs.append(contig)
    sorted_contigs = sorted(positive_contigs, key=lambda contig: len(contig.sequence), reverse=True)

    # Create the new number mapping.
    old_nums_to_new_nums = {}
    i = 1
    for contig in sorted_contigs:
        old_nums_to_new_nums[contig.number] = i
        i += 1

    # Assign new numbers and sort using them.
    for contig in contigs:
        contig.renumber(old_nums_to_new_nums[contig.number])
    contigs.sort(key=lambda contig: (contig.number, contig.get_sign()))

    print('done')

def get_best_blast_alignment(blast_alignments, segment_length, expected_reference_start):
    """
    Find and return the best of the given BLAST alignments.  'Best' is defined
    as covering the entirety of the alignment, being high identity and having
    an appropriate start position.
    """
    if not blast_alignments:
        return None

    # Sort by alignment length and identity.
    sorted_alignments = sorted(blast_alignments, key=lambda alignment: (alignment.get_query_length(), alignment.percent_identity), reverse=True)

    # If we have an expected reference start to work with, then we find the
    # first alignment in the list which occurs very near the expected reference
    # start.
    if expected_reference_start is not None:
        for alignment in sorted_alignments:
            discrepancy = abs(alignment.get_query_start_in_reference() - expected_reference_start)
            if discrepancy < 5:
                return alignment

    # If we don't have an expected reference start to work with, then we are
    # more limited.  We can only give a result if there is a single full length
    # alignment.
    else:
        full_length_alignments = []
        for alignment in sorted_alignments:
            fraction_present = alignment.get_query_length() / segment_length
            if fraction_present == 1.0:
                full_length_alignments.append(alignment)
        if len(full_length_alignments) == 1:
            return full_length_alignments[0]

    # If the code got here, then no good match was found.
    return None


# This class holds a contig: its name, sequence and length.
class Contig:

    def __init__(self, name, sequence):
        self.fullname = name
        name_parts = name.split('_')
        self.number = int(name_parts[1])
        cov_string = name_parts[5]
        if cov_string[-1] == "'":
            cov_string = cov_string[:-1]
        self.cov = float(cov_string)
        self.sequence = sequence
        self.outgoing_linked_contigs = []
        self.incoming_linked_contigs = []
        self.path = Path()

    def __str__(self):
        return self.fullname

    def __repr__(self):
        return self.fullname

    def renumber(self, new_number):
        self.number = new_number
        self.rebuild_full_name()

    def rebuild_full_name(self):
        """
        Remake the contig's full name using its current number, sequence and
        depth.
        """
        positive = self.is_positive()
        self.fullname = 'NODE_' + str(self.number) + '_length_' + str(len(self.sequence)) + '_cov_' + str(self.cov)
        if not positive:
            self.fullname += "'"

    def get_number_with_sign(self):
        """
        Return the contig number in string form with a + or - at the end.
        """
        num = str(self.number)
        if self.is_positive():
            num += '+'
        else:
            num += '-'
        return num

    def get_sign(self):
        if self.is_positive():
            return '+'
        else:
            return '-'

    def is_positive(self):
        return self.fullname[-1] != "'"

    def get_header_with_links(self):
        """
        Produce a FASTG header for the contig, with links and a line break at
        the end.
        """
        header_with_edges = '>' + self.fullname

        if len(self.outgoing_linked_contigs) > 0:
            header_with_edges += ':'
            for linked_contig in self.outgoing_linked_contigs:
                header_with_edges += linked_contig.fullname + ','
            header_with_edges = header_with_edges[0:-1]

        header_with_edges += ';\n'
        return header_with_edges

    def get_sequence_with_line_breaks(self):
        sequence_remaining = self.sequence
        sequence_with_line_breaks = ''
        while len(sequence_remaining) > 60:
            sequence_with_line_breaks += sequence_remaining[0:60] + '\n'
            sequence_remaining = sequence_remaining[60:]
        sequence_with_line_breaks += sequence_remaining + '\n'
        return sequence_with_line_breaks

    def find_segment_locations(self, segment):
        return self.path.find_segment_locations(segment)

    def find_segment_locations_plus_one(self, segment):
        segment_locations = self.path.find_segment_locations(segment)
        return [x+1 for x in segment_locations]

    def determine_all_segment_locations(self, segment_sequences, graph_overlap):
        """
        Determine the start and end coordinates of each segment in the contig.
        This information is necessary before the contig can be split.
        """
        # Create a temporary directory for doing BLAST searches.
        if not os.path.exists('spades_contig_graph-temp'):
            os.makedirs('spades_contig_graph-temp')

        save_sequence_to_fasta_file(self.sequence, self.fullname, 'spades_contig_graph-temp/contig.fasta')

        # Create a BLAST database for the contig.
        makeblastdb_command = ['makeblastdb', '-dbtype', 'nucl', '-in', 'spades_contig_graph-temp/contig.fasta']
        process = subprocess.Popen(makeblastdb_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()

        # Check that makeblastdb ran okay
        if len(err) > 0:
            print('\nmakeblastdb encountered an error:', file=sys.stderr)
            print(err, file=sys.stderr)
            quit()

        # This value tracks where we expect the next segment to start, in
        # contig coordinates.  When it is set to None, that means we don't
        # know.
        expected_segment_start_in_contig = 1

        for i in range(len(self.path.segment_list)):
            segment = self.path.segment_list[i]

            # Don't deal with assembly gaps just yet - we'll give them contig
            # start/end coordinates after we've finished with the real
            # segments.
            if segment.startswith('gap'):
                expected_segment_start_in_contig = None
                continue

            segment_sequence = segment_sequences[segment]
            segment_length = len(segment_sequence)

            save_sequence_to_fasta_file(segment_sequence, segment, 'spades_contig_graph-temp/segment.fasta')

            # BLAST for the segment in the contig
            sys.stdout.flush()
            blastn_command = ['blastn', '-task', 'blastn', '-db', 'spades_contig_graph-temp/contig.fasta', '-query', 'spades_contig_graph-temp/segment.fasta', '-outfmt', '6 length pident sstart send qstart qend']
            process = subprocess.Popen(blastn_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = process.communicate()

            # Check that blastn ran okay
            if len(err) > 0:
                print('\nblastn encountered an error:', file=sys.stderr)
                print(err, file=sys.stderr)
                quit()

            # Save the alignments in Python objects.
            alignment_strings = out.splitlines()
            blast_alignments = []
            for alignment_string in alignment_strings:
                alignment = BlastAlignment(alignment_string, segment_length)
                blast_alignments.append(alignment)

            # Get the alignment (if any) that best matches the segment sequence
            # and position.
            best_alignment = get_best_blast_alignment(blast_alignments, segment_length, expected_segment_start_in_contig)

            # If we found an alignment, use it to determine the query's start
            # and end coordinates in the contig.
            if best_alignment is not None:
                segment_start_in_contig = best_alignment.get_query_start_in_reference()
                segment_end_in_contig = best_alignment.get_query_end_in_reference()

            # If we failed to find an alignment, then we don't have segment
            # coordinates and will be unable to split the contig at this point.
            else:
                segment_start_in_contig = None
                segment_end_in_contig = None

            self.path.contig_coordinates[i] = (segment_start_in_contig, segment_end_in_contig)

            # Update the expected segment start for the next segment.
            # Hopefully we can use the alignment to do this precisely.
            if best_alignment is not None:
                expected_segment_start_in_contig = (best_alignment.reference_end - graph_overlap + 1)

            # If there is no alignment with which we can predict the next
            # segment start, we can guess using the segment length.
            elif expected_segment_start_in_contig is not None:
                expected_segment_start_in_contig += (segment_length - graph_overlap)

            os.remove('spades_contig_graph-temp/segment.fasta')

        shutil.rmtree('spades_contig_graph-temp')

        # Now we have to go back and assign contig start/end positions for any
        # gap segments.
        segment_count = len(self.path.segment_list)
        for i in range(segment_count):
            segment = self.path.segment_list[i]
            if not segment.startswith('gap'):
                continue

            gap_start_in_contig = 1
            gap_end_in_contig = len(self.sequence)

            if i > 0:
                previous_segment_end = self.path.contig_coordinates[i - 1][1]
                if previous_segment_end is not None:
                    gap_start_in_contig = previous_segment_end - graph_overlap + 1
                    if gap_start_in_contig < 1:
                        gap_start_in_contig = 1
                else:
                    gap_start_in_contig = None

            if i < segment_count - 1:
                next_segment_start = self.path.contig_coordinates[i + 1][0]
                if next_segment_start is not None:
                    gap_end_in_contig = next_segment_start + graph_overlap - 1
                    if gap_end_in_contig > len(self.sequence):
                        gap_end_in_contig = len(self.sequence)
                else:
                    gap_end_in_contig = None

            self.path.contig_coordinates[i] = (gap_start_in_contig, gap_end_in_contig)

    def get_links_to_other_contigs(self):
        """
        Return a list of the links from this contig to any other contig.
        """
        links_to_other_contigs = []
        for outgoing_linked_contig in self.outgoing_linked_contigs:
            links_to_other_contigs.append((self.path.get_last_segment(), outgoing_linked_contig.path.get_first_segment()))
        for incoming_linked_contig in self.incoming_linked_contigs:
            links_to_other_contigs.append((incoming_linked_contig.path.get_last_segment(), self.path.get_first_segment()))
        return list(set(links_to_other_contigs))

    def get_links_in_this_contig_and_to_other_contigs(self):
        links = self.path.get_all_links()
        links.extend(self.get_links_to_other_contigs())
        return list(set(links))


# This class holds a path: the lists of graph segments making up a contig.
class Path:
    def __init__(self, segment_list=[]):
        self.segment_list = segment_list
        self.contig_coordinates = [(0, 0) for i in range(len(segment_list))]

    def get_first_segment(self):
        return self.segment_list[0]

    def get_last_segment(self):
        return self.segment_list[-1]

    def __str__(self):
        return str(self.segment_list) + ', ' + str(self.contig_coordinates)

    def __repr__(self):
        return str(self.segment_list) + ', ' + str(self.contig_coordinates)

    def get_segment_count(self):
        return len(self.segment_list)

    def find_segment_locations(self, segment):
        locations = []
        for i in range(len(self.segment_list)):
            if segment == self.segment_list[i]:
                locations.append(i)
        return locations

    def get_paths_with_line_breaks(self):
        output = ''
        for segment in self.segment_list:
            if segment.startswith('gap'):
                output += ';\n'
            else:
                output += segment + ','
        return output[:-1] + '\n'

    # This function returns a list of tuple for each link in the path.
    def get_all_links(self):
        links = []
        for i in range(len(self.segment_list) - 1):
            s_1 = self.segment_list[i]
            s_2 = self.segment_list[i + 1]
            links.append((s_1, s_2))
        return links


class BlastAlignment:

    def __init__(self, blast_string, query_length):
        blast_string_parts = blast_string.split('\t')

        self.percent_identity = float(blast_string_parts[1])

        self.reference_start = int(blast_string_parts[2])
        self.reference_end = int(blast_string_parts[3])

        self.query_start = int(blast_string_parts[4])
        self.query_end = int(blast_string_parts[5])

        self.query_length = query_length

    def get_query_length(self):
        return self.query_end - self.query_start + 1

    def get_query_start_in_reference(self):
        query_missing_bases_at_start = self.query_start - 1
        return self.reference_start - query_missing_bases_at_start

    def get_query_end_in_reference(self):
        query_missing_bases_at_end = self.query_length - self.query_end
        return self.reference_end + query_missing_bases_at_end


    def __lt__(self, other):
        return self.get_query_start_in_reference() < other.get_query_start_in_reference()

    def __str__(self):
        return 'reference: ' + str(self.reference_start) + ' to ' + str(self.reference_end) + \
               ', query: ' + str(self.query_start) + ' to ' + str(self.query_end) + \
               ', identity: ' + str(self.percent_identity) + '%'

    def __repr__(self):
        return 'reference: ' + str(self.reference_start) + ' to ' + str(self.reference_end) + \
               ', query: ' + str(self.query_start) + ' to ' + str(self.query_end) + \
               ', identity: ' + str(self.percent_identity) + '%'


# Standard boilerplate to call the main() function to begin the program.
if __name__ == '__main__':
    main()
