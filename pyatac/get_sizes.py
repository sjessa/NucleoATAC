"""
Script to get fragment size distribution

@author: Alicia Schep. Updated by Selin Jessa.
"""

##### IMPORT MODULES #####
# import necessary for python
import os
#import matplotlib as mpl
#mpl.use('PS')
import matplotlib.pyplot as plt
from pyatac.fragmentsizes import FragmentSizes
from pyatac.chunk import ChunkList

def get_sizes(args):
    """function to get fragment sizes

    """
    if args.out is None:
        args.out = '.'.join(os.path.basename(args.bam).split('.')[0:-1])

    sizes = FragmentSizes(lower = args.lower, upper = args.upper, atac = args.atac)

    # depending on the input, perform the appropriate frag size distribution calculation
    if args.bam:

        if args.bed:
            chunks = ChunkList.read(args.bed)
            chunks.merge()
            sizes.calculateSizes(input_file = args.bam, input_type = "bam", chunks = chunks)
        else:
            sizes.calculateSizes(input_file = args.bam, input_type = "bam")

    elif args.fragments:

        sizes.calculateSizes(input_file = args.fragments, input_type = "fragments")

    sizes.save(args.out+'.fragmentsizes.txt')

    if not args.no_plot:
        #make figure
        fig = plt.figure()
        plt.plot(range(sizes.lower,sizes.upper),sizes.get(sizes.lower,sizes.upper),label = args.out)
        plt.xlabel("Fragment Size")
        plt.ylabel("Frequency")
        fig.savefig(args.out+'.fragmentsizes.pdf')
        plt.close(fig)

