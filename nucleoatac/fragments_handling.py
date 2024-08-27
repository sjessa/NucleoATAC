"""
Utilities for handling fragments files as input to NucleoATAC.

@author: Selin Jessa, adapted from Alicia Schep
"""

import numpy as np

def getAllFragmentSizesFromFragsFile(fragments, lower, upper):
    """
    Calculate fragment size distribution from a fragments file in BED format.
    Replaces the getAllFragmentSizes function from pyatac.fragments, to use
    a fragments file as input instead of a BAM file.

    Parameters:
    - fragments: Path to the BED file containing fragments
    - lower: Lower bound for fragment sizes to consider
    - upper: Upper bound for fragment sizes to consider
    
    """

    # initialize an array to hold the size counts within the specified range
    sizes = np.zeros(upper - lower, dtype=np.float)

    # open file
    with open(fragments, 'r') as f:
        for line in f:
            
			# split the line into columns
            cols = line.strip().split()
            start = int(cols[1])
            end = int(cols[2])

            # calculate fragment size
            ilen = end - start
            
            # only count fragments within the specified size range
            if lower <= ilen < upper:
                sizes[ilen - lower] += 1

    return sizes
