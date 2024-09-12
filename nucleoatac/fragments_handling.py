"""
Utilities for handling fragments files as input to NucleoATAC.

@author: Selin Jessa, Greenleaf Lab, Stanford University. Adapted from Alicia Schep.

NOTE: in general, for fragments, we will assume they're already shifted +4/-4.
So unlike the BAM-handling counterparts to the functions here, we won't take in an `atac`
argument, and no shifting will be done for fragments.
"""

import numpy as np
import pysam

def getAllFragmentSizesFromFragsFile(fragments, lower, upper):
    """
    Calculate fragment size distribution from a fragments file in BED format.
    Adapted from pyatac.fragments.getAllFragmentSizes, to use
    a fragments file as input instead of a BAM file.

    Args:
    - fragments: path to the compressed, tabix-indexed BED file containing fragments. Tabix index
    expected at `fragments + '.tbi'`.
    - lower: lower bound for fragment sizes to consider
    - upper: upper bound for fragment sizes to consider
    
    Returns:
    - sizes: NumPy array containing the fragment size distribution within the specified range.
    """

    # initialize an array to hold the size counts within the specified range
    sizes = np.zeros(upper - lower, dtype=np.float)

    # open tabix-indexed frags file
    tbx = pysam.TabixFile(fragments)

    # fetch to get an iterator over all fragments
    for row in tbx.fetch(parser=pysam.asBed()):
            
        # extract coords of fragment
        start = int(row.start)
        end = int(row.end)

        # calculate fragment size
        ilen = end - start
        
        # only count fragments within the specified size range
        if lower <= ilen < upper:
            sizes[ilen - lower] += 1

    tbx.close()

    return sizes




def getAllFragmentSizesFromFragsFileFromChunkList(chunks, fragments, lower, upper):
    """
    Calculate fragment size distribution from a tabix-indexed fragments file in BED format,
    using a ChunkList to specify the regions of interest. Adapted from pyatac.fragments.getFragmentSizesFromChunkList.

    Args:
    - chunks: ChunkList object containing regions of interest.
    - fragments: path to the tabix-indexed BED file containing fragments (e.g., .bed.gz).
    - lower: lower bound for fragment sizes to consider.
    - upper: upper bound for fragment sizes to consider.

    Returns:
    - sizes: NumPy array containing the fragment size distribution within the specified range.
    """
    
    # initialize an array to hold the size counts within the specified range
    sizes = np.zeros(upper - lower, dtype=np.float)

    # open the tabix-indexed fragments file using pysam
    tbx = pysam.TabixFile(fragments)

    # iterate over each chunk in the ChunkList
    for chunk in chunks:
        
        # fetch fragments in the specified chunk region
        for row in tbx.fetch(chunk.chrom, max(0, chunk.start - upper), chunk.end + upper, parser=pysam.asBed()):

            # extract coords of fragment
            fragment_start = int(row.start)
            fragment_end = int(row.end)

            # calculate fragment size
            ilen = fragment_end - fragment_start

            # calculate the center of the fragment
            center = fragment_start + (ilen - 1) / 2

            # check if fragment size is within the specified bounds and if the center is within the chunk
            if lower <= ilen < upper and chunk.start <= center < chunk.end:
                sizes[ilen - lower] += 1

    tbx.close()

    return sizes




def makeFragmentMatFromFragments(fragments, chrom, start, end, lower, upper):
    """
    Generate a fragment matrix from a fragments file in BED format. This function
    produces a 2D matrix (fragment matrix) where each row represents a different
    fragment size (or length), and each column represents a position within a
    specified genomic region. The matrix records the number of fragments of a
    specific size that start or overlap at each position. Adapted from
    pyatac.fragments.makeFragmentMat.

    Uses pysam to read tabix-indexed BED files. For more information, see:
    https://pysam.readthedocs.io/en/latest/usage.html#working-with-tabix-indexed-files

    Args:
    - fragments: path to the compressed, tabix-indexed BED file containing fragments. Tabix index
    expected at `fragments + '.tbi'`.
    - chrom: chromosome corresponding to chunk (genomic region) to process
    - start: start position for chunk to process
    - end: end position for chunk to process
    - lower: lower bound of fragment sizes to consider
    - upper: upper bound of fragment sizes to consider
    """
    
    # init matrix
    nrow = upper - lower
    ncol = end - start
    mat = np.zeros((nrow, ncol), dtype=np.float)

    # open tabix-indexed file
    tbx = pysam.TabixFile(fragments)

    # fetch to get an iterator over fragments in the specified chunk/region
    for row in tbx.fetch(chrom, start, end, parser=pysam.asBed()):

        # extract fragment info 
        fragment_start = int(row.start)
        fragment_end = int(row.end)

        # calculate the fragment size
        ilen = fragment_end - fragment_start

        # ensure the fragment size is within the specified bounds
        if lower <= ilen < upper:
            
            # calculate row and column idx in the matrix
            row = ilen - lower
            col = (ilen - 1) // 2 + fragment_start - start

            # check if the calculated column and row are within the matrix bounds
            if 0 <= col < ncol and 0 <= row < nrow:
                mat[row, col] += 1

    # close the tabix file
    tbx.close()

    return mat




def getInsertionsFromFragments(fragments, chrom, start, end, lower = 0, upper = 2000):
    """
    Calculate insertions from a tabix-indexed fragments file in BED format.
    This function generates an insertion profile (1D array) for a specified
    genomic region. Adapted from pyatac.fragments.getInsertions.

    Args:
    - fragments: path to the compressed, tabix-indexed BED file containing fragments.
    - chrom: chromosome corresponding to the region to process.
    - start: start position for the region to process.
    - end: end position for the region to process.
    - lower: lower bound for fragment sizes to consider.
    - upper: upper bound for fragment sizes to consider.

    Returns:
    - mat: NumPy array containing the count of insertions at each position within the specified region.
    """

    # initialize an array to hold the insertion counts within the specified region
    npos = end - start
    mat = np.zeros(npos, dtype=np.float)

    # open the tabix-indexed fragments file using pysam
    tbx = pysam.TabixFile(fragments)

    # fetch fragments in the specified region
    for row in tbx.fetch(chrom, start, end, parser=pysam.asBed()):
        
        # get coordinates of the fragment
        fragment_start = int(row.start)
        fragment_end = int(row.end)

        # calculate fragment size
        ilen = fragment_end - fragment_start

        # ensure the fragment size is within the specified bounds
        if lower <= ilen < upper:
            # calculate left (l_pos) and right (r_pos) insertion positions
            l_pos = fragment_start
            r_pos = fragment_end - 1

            # increment the count at the insertion positions if within bounds
            if start <= l_pos < end:
                mat[l_pos - start] += 1
            if start <= r_pos < end:
                mat[r_pos - start] += 1

    # close the tabix file
    tbx.close()

    return mat