# cs190FileUtil.py 
#
# Jeff Parker        Feb 2009
#
# Utilities to take a file name of a fasta file and return processed string
#
# Usage
#     text = readFastaFile(sys.argv[1])

import string

# Remove \n, change to upper case
def prepare(text):
    text = string.replace(text, " ", "")
    text = string.replace(text, "\n", "")
    text = text.upper()
    return text


def readFastaFile(fileName):
    f = open(fileName, 'r')

    # First line in Fasta format is header
    line = f.readline()
    print "Saw", line

    text = f.read()
    f.close()
    
    # Remove end of lines
    text = prepare(text)
    return text
