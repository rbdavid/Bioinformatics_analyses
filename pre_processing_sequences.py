#!/home/rbdavid/bin/python

import sys
import re
from Bio import SeqIO

original_sequence_file = sys.argv[1]
output_file_name_root = sys.argv[2]
regular_expression = sys.argv[3]     # regular expression formatting
avg_sequence_length = int(sys.argv[4])
sequence_range_value = int(sys.argv[5])
avg_expression_length = int(sys.argv[6])
expression_range_value = int(sys.argv[7])
number_residues_pre_expression = int(sys.argv[5])
prob_XBZ = float(sys.argv[7])

# Code outline:
# Overall goal:
#    Loop over all sequences and check to see if I want to keep it or drop it
#    I want to keep only sequences of the helicase domain of NS3
#    therefore, I need to remove sequences that are too short to be complete sequences of NS3h or too poorly resolved to be of any use; 
#    Additionally, I want to remove segments of sequences that do not correspond to non-NS3h domains; sequences that have abnormal distances between motifs are to be cut; sequences that are originally too long need to be truncated down to a good size.
# Code break down:
# 1) Read in the original sequence file (assumed to be fasta)
# 2) Loop through all sequences, performing a set of tests on each sequence to output only the desired sequences. To do: 
#    ID the 'start' of the domain of interest
#    check to see if sequence is of expected size; if too short, chuck it; if too long, truncate the sequence
#    check to see if this truncated sequence has poor resolution (too many 'X', 'B', 'Z' residues), if so, chuck it
#    check to see if the sequence is identical to other sequences already analyzed,
#        if not, append the truncated sequence to a dictionary where the sequence is the key and accession number is the value
#        if so, append the accession number to the original sequence and chuck the redundant sequence
# 3) Output a new fasta file that holds only the unique, truncated sequences of the NS3 helicase

# objective 1
sequences = list(SeqIO.parse(original_sequence_file,'fasta'))
nSeqs = len(sequences)
minimum_sequence_length = avg_sequence_length - sequence_range_value
maximum_sequence_length = avg_sequence_length + sequence_range_value + 1
minimum_expression_length = avg_expression_length - expression_range_value
maximum_expression_length = avg_expression_length + expression_range_value + 1

# objective 2
sequence_dictionary = {}
with open(output_file_name_root+'.summary.txt','w') as W:
        for i in range(nSeqs):
                sequence = str(sequences[i].seq)
                last_found = -1
                found = []
                while True:
                        try:
                                search_string = re.search(regular_expression,sequence[last_found+1:]).group(0)
                                zeroth_index_search_string = sequence.find(search_string)
                                zeroth_index_sequence = zeroth_index_search_string - number_residues_pre_expression
                                # elements in found: index of first character in search_string, search_string, length of sequence x residues before search_string starts, truncated sequence string
                                if zeroth_index_sequence < 0:   # 2nd element in found will correspond to the length of the whole sequence, since search_string is found early in sequence
                                        found.append([zeroth_index_search_string,search_string,len(sequence[:]),sequence[:]])
                                else:   # 2nd element in found will correspond to the length of the truncated sequence (have only removed residues well before search_string instance)
                                        found.append([zeroth_index_search_string,search_string,len(sequence[zeroth_index_sequence:]),sequence[zeroth_index_sequence:]])
                        except AttributeError:
                                break
                        last_found = sequence.find(search_string,last_found+1)

                if last_found == -1:
                        W.write('%-5d %s Regular Expression input not found, removed from population of sequences.\n'%(i,sequences[i].id))
                        continue

                for j in found:
                        if j[2] < minimum_sequence_length:  # ignore search_strings that would create truncated sequences that are shorter than the sequences of interest are expected to be...
                                W.write('%-5d %s Regular Expression input found but remaining sequence is too short (length = %d) to be the desired sequence.\n'%(i,sequences[i].id,j[2]))
                                continue
                        elif len(j[1]) < minimum_expression_length or len(j[1]) > maximum_expression_length:    # ignore search_strings that have a length outside of an expected range of expression lengths; the number of residues between motifs X and Y are generally conserved...
                                W.write('%-5d %s Regular Expression input found but length of the search string is unexpected (length = %d) (expected length within %d to %d residues).\n'%(i,sequences[i].id,len(j[1]),minimum_expression_length,maximum_expression_length))
                                continue
                        elif j[2] > maximum_sequence_length:    # truncate (post-search_string) segments of the sequence that are beyond the expected length of sequences; only want the sequence segment j[3][:avg_sequence_length+1] 
                                sequence = j[3][:avg_sequence_length+1]
                                W.write('%-5d %s Regular Expression input found but length truncated (pre-search_string) sequence is longer than expected length; truncating (post-search_string) segment by assuming the sequence should be approximately the average length. '%(i,sequences[i].id))
                        else:
                                sequence = j[3]
                                W.write('%-5d %s Regular Expression input found; sequence length is within expected range. '%(i,sequences[i].id))
                        
                        if float(sequence.count('X') + sequence.count('B') + sequence.count('Z'))/len(sequence) <= prob_XBZ:    # set a resolution cutoff for sequences; if sequence resolution is larger than prob_XBZ
                                if sequence not in sequence_dictionary:
                                        sequence_dictionary[sequence] = sequences[i].id
                                        W.write('This sequence string is the first instance of such a sequence; Adding sequence to sequence_dictionary.\n')
                                else:
                                        sequence_dictionary[sequence] += '___' + sequences[i].id
                                        W.write('This sequence string has been seen before; appending sequence info to previous sequence_dictionary value.\n')
                        else:
                                W.write('This sequence string is poor resolution.\n')

                print 'Finished analyzing sequence', i

# objective 3
with open(output_file_name_root + '.fasta','w') as output_file:
        for sequence in sequence_dictionary:
                output_file.write('>'+sequence_dictionary[sequence]+'\n'+sequence+'\n')

