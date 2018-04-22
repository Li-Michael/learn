import sys

def overlap(a, b, min_length=3):
    """ return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length' 
        character long. If no such overlap exists, return 0. """
    assert a != '' and b != ''
    start = 0   # start all the way at the left
    for i in range(min_length, len(a)):
        start =  a.find(b[:min_length], start)  # look for b's suffix in a
        if start == -1:   ## no more occurrences to left 
            return 0
        elif b.startswith(a[start:]):  ## found occurrence; check if full prefix/suffix match
            return len(a) -start
        start += 1  ## move just past previous match
        
#from itertools import permutations,product
import itertools
def naive_overlap(reads, k=3):
    """ find all overlaps in reads """
    overlaps = {}
    for a,b in itertools.permutations(reads, 2):
        overlen = overlap(a, b, min_length=k)
        if overlen >0:
            overlaps[(a, b)] = overlen
    return overlaps


def kmer_overlap(reads, k=3):
    """ use k-mer dictionary to find all overlaps in reads. faster """
    kmer_index = {}
    kmer_prefix = {}
    overlaps = {}
    
    for read in reads:
        for i in range(len(read)-k+1):
            if read[i:i+k] not in kmer_index:
                kmer_index[read[i:i+k]] = set()
            if read[:k] not in kmer_prefix:
                kmer_prefix[read[:k]] = set()
            kmer_index[read[i:i+k]].add(read)
            kmer_prefix[read[:k]].add(read)

    for key in kmer_prefix:
        for a,b in itertools.product(kmer_index[key], kmer_prefix[key]):
            if a!= b and (a,b) not in overlaps:
                overlen = overlap(a,b,min_length=k)
                if overlen > 0:
                    overlaps[(a,b)] = overlen

    return overlaps
       
def scs(reads):
    """ Returns shortest common superstring of given 
        strings, which must be the same length. 
        Use brute force to solve the shortest common superstring """
    superstring = None
    for ss in itertools.permutations(reads):
        sstring = ss[0]      # superstring starts as first string
        for i in range(1,len(ss)):
            # overlap adjacent strings A and B in the permutationi
            # and  add non-overlapping portion of B to superstring
            sstring += ss[i][overlap(sstring,ss[i], min_length=1):]
        if superstring is None or len(sstring) < len(superstring):    
            superstring = sstring   # found shorter superstring
    return superstring


def pick_max_overlap(reads, k):
    reada, readb = None, None
    best_olen = 0
    for a,b in itertools.permutations(reads, 2):
        olen = overlap(a, b, min_length=k)
        if olen > best_olen:
            reada, readb = a, b
            best_olen = olen
    return reada, readb, best_olen


def greedy_scs(reads, k=3):
    """ greedy shortest common superstring """
    reada, readb, olen = pick_max_overlap(reads, k)
    while olen > 0:
        reads.remove(reada)
        reads.remove(readb)
        reads.append(reada+readb[olen:])
        reada, readb, olen = pick_max_overlap(reads, k)
    return ''.join(reads)
        



