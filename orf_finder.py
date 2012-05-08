#
# Program to find open reading frames.
#
# Prints a list of open reading frames in all six frames for a sequence file
#   in FASTA format.
#
# Aniket Schneider 
# Feb 2012
#
# Updated revcomp to accept lower case characters
#
# Based in part on orf_skel.py by Jeff Parker
#


import string
import sys
import cs190FileUtil

MIN_ORF_LENGTH = 300   # minimum ORF length in bp, can be overridden

def revcomp(seq):
    '''
    Returns the reverse complement of seq.

    Non-DNA nucleotides are not complemented.
    '''
    return seq[::-1].translate(string.maketrans("ACGTacgt","TGCAtgca"))


def find_all_orfs(seq, minlen=0):
    '''
    Returns records of all orfs at least minlen in length in seq.

    Only the longest overlapping ORF in a given reading frame is recorded.
    '''
    orf_list = []
    for frame in xrange(3):
        orf_list.extend(
                find_inframe_orfs(seq, frame, reverse=False, minlen=minlen))
        
    seq = revcomp(seq)
    
    for frame in xrange(3):
        orf_list.extend(
                find_inframe_orfs(seq, frame, reverse=True, minlen=minlen))

    return orf_list

        

def find_inframe_orfs(seq, frame, reverse=False, minlen=0):
    '''
    Prints a list of orfs in seq in the given reading frame.

    To search the reverse direction, pass in the reverse complement of the 
    sequence and set reverse to True.

    frame is an offset between 0 and 2 corresponding to frames +1 to +3.
    reverse specifies whether printed sequence indices should be in reference 
    to the end or the beginning of seq.
    Note: running off the end of the strand is counted as a translation stop.
    '''
    length = len(seq)
    in_orf = False
    start = ["ATG","atg"]
    stops = ["TAA","TAG","TGA","taa","tag","tga"]
    orf_list = []
    
    for i in xrange(frame, length, 3):
        if not in_orf:
            if seq[i:i+3] in start:
                if not reverse:
                    orf_start = i + 1
                else: # if reverse
                    orf_start = length - i
                in_orf = True

        else: # if in_orf == True
            if seq[i:i+3] in stops or i >= (length - 3):
                in_orf = False
                
                if not reverse:
                    orf_end = i + 3
                    orf_len = orf_end - orf_start + 1
                    if orf_len >= minlen:
                        orf_list.append( {'frame': frame+1, 
                            'start': orf_start, 
                            'end': orf_end, 
                            'length': orf_len,
                            'seq': seq[ orf_start-1 : orf_end ]} )
                else: # if reverse
                    orf_end = length - i - 2
                    orf_len = orf_start - orf_end + 1
                    if orf_len >= minlen:
                        # note reversal of start and end for reverse strand
                        orf_list.append( {'frame': -(frame+1), 
                            'start': orf_end, 
                            'end': orf_start, 
                            'length': orf_len,
                            'seq': seq[ length-orf_start : i + 3 ]} )

    return orf_list
    
# Execute when run as a script
if __name__ == '__main__':
    if (len(sys.argv) < 2):
        print "Usage: python", sys.argv[0], "<filename in FASTA format> [<min ORF length>]"
    else:
        fileName = sys.argv[1]
        if (len(sys.argv) > 2):             # 2nd arg should be an integer
            try:
                minlen = int(sys.argv[2])    # Convert string to integer
            except ValueError:              # try-except catches errors
                print "\n\tExpecting an integer to define min ORF length, found",
                print sys.argv[2]
                exit()
        else:
            minlen = MIN_ORF_LENGTH

        print "ORF must be at least", minlen, "Base pairs long"

        text = cs190FileUtil.readFastaFile(fileName)
    
        # Time to start finding ORFs!
        orf_list = find_all_orfs(text, minlen)

        for orf in orf_list:
            print "Frame {frame} Start {start} End {end} Len {length}".format(
                frame=orf['frame'],start=orf['start'],end=orf['end'],length=orf['length'])
