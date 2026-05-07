"""
Classes for working with fragment distribution

@author: Alicia Schep, Greenleaf Lab, Stanford University. Updated by Selin Jessa,
Greenleaf Lab, Stanford University.
"""

import numpy as np
import pyximport
pyximport.install(setup_args={"include_dirs":np.get_include()})
from pyatac.fragments import getAllFragmentSizes, getFragmentSizesFromChunkList
from nucleoatac.fragments_handling import getAllFragmentSizesFromFragsFile,getAllFragmentSizesFromFragsFileFromChunkList


class FragmentSizes:
    """Class for storing fragment size distribution"""

    def __init__(self, lower, upper, atac = True, vals = None):
        self.lower = lower
        self.upper = upper
        self.vals = vals
        self.atac = atac
        self.raw = None

    def calculateSizes(self, input_file, input_type = "bam", chunks = None):
        """Calculate fragment size distribution from a BAM file or fragments file"""
        if input_type == "bam":
            if chunks is None:
                sizes = getAllFragmentSizes(input_file, self.lower, self.upper, atac = self.atac)
            else:
                sizes = getFragmentSizesFromChunkList(chunks, input_file, self.lower, self.upper, atac = self.atac)
        elif input_type == "fragments":
            if chunks is None:
                sizes = getAllFragmentSizesFromFragsFile(input_file, self.lower, self.upper)
            else:
                sizes = getAllFragmentSizesFromFragsFileFromChunkList(chunks, input_file, self.lower, self.upper)

        self.raw = np.asarray(sizes).astype(np.int64)
        self.vals = sizes / (np.sum(sizes) + (np.sum(sizes)==0))

    def get(self, lower = None, upper = None, size = None):
        if size:
            try:
                return self.vals[size - self.lower]
            except:
                raise Exception("Looks like size doesn't match FragmentSizes")
        else:
            if lower is None:
                lower = self.lower
            if upper is None:
                upper = self.upper
            y1 = lower - self.lower
            y2 = upper - self.lower
            try:
                return self.vals[y1:y2]
            except:
                raise Exception("Looks like dimensions from get probaby don't match FragmentSizes")
            
    def save(self, filename):
        """Save Fragment Distribution information"""
        with open(filename,"w") as f:
            f.write("#lower\n")
            f.write(str(self.lower)+"\n")
            f.write("#upper\n")
            f.write(str(self.upper)+"\n")
            f.write("#sizes\n")
            f.write("\t".join(map(str,self.get()))+"\n")
            if self.raw is not None:
                f.write("#raw_counts\n")
                f.write("\t".join(map(str, self.raw.astype(np.int64).tolist()))+"\n")

    @staticmethod
    def open(filename):
        """Create FragmentDistribution object from text descriptor file"""
        state = ''
        raw_counts = None
        with open(filename,'r') as infile:
            for line in infile:
                if '#lower' in line:
                    state = 'lower'
                elif '#upper' in line:
                    state = 'upper'
                elif '#sizes' in line:
                    state = 'sizes'
                elif '#raw_counts' in line:
                    state = 'raw_counts'
                elif '#' in line:
                    state = 'other'
                elif state == 'lower':
                    lower = int(line.strip('\n'))
                elif state == 'upper':
                    upper = int(line.strip('\n'))
                elif state == 'sizes':
                    fragmentsizes = np.array(list(map(float,line.rstrip("\n").split("\t"))))
                elif state == 'raw_counts':
                    raw_counts = np.array(list(map(int, line.rstrip("\n").split("\t"))), dtype=np.int64)
        try:
            new = FragmentSizes(lower, upper, vals = fragmentsizes)
        except NameError:
            raise Exception("FragmentDistribution decriptor file appeas to be missing some\
needed components")
        new.raw = raw_counts
        return new



