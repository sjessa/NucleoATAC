"""
Script to make nucleosome occupancy track!

@author: Alicia Schep
"""

##### IMPORT MODULES #####
# import necessary python modules
#import matplotlib as mpl
#mpl.use('PS')
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import traceback
import itertools
import pysam
from pyatac.utils import shell_command,read_chrom_sizes_from_bam,read_chrom_sizes_from_fasta
from pyatac.chunk import ChunkList
from nucleoatac.Occupancy import FragmentMixDistribution, OccupancyParameters, OccChunk
from pyatac.fragmentsizes import FragmentSizes
from pyatac.bias import PWM

def _occHelper(arg):
    """function to get occupancy for a set of bed regions

    """
    (chunk, params) = arg
    try:
        occ = OccChunk(chunk)
        occ.process(params)
        out = (occ.getNucDist(),
                occ.occ, [occ.peaks[i] for i in sorted(occ.peaks.keys())])
        occ.removeData()
    except Exception as e:
        print('Caught exception when processing:\n'+  chunk.asBed()+"\n")
        traceback.print_exc()
        print()
        raise e
    return out

def _writeOcc(track_queue, out):
    out_handle1 = open(out + '.occ.bedgraph','a')
    out_handle2 = open(out + '.occ.lower_bound.bedgraph','a')
    out_handle3 = open(out + '.occ.upper_bound.bedgraph','a')
    try:
        for track in iter(track_queue.get, 'STOP'):
            track.write_track(out_handle1, vals = track.smoothed_vals)
            track.write_track(out_handle2, vals = track.smoothed_lower)
            track.write_track(out_handle3, vals = track.smoothed_upper)
            track_queue.task_done()
    except Exception, e:
        print('Caught exception when writing occupancy track\n')
        traceback.print_exc()
        print()
        raise e
    out_handle1.close()
    out_handle2.close()
    out_handle3.close()
    return True


def _writePeaks(pos_queue, out):
    out_handle = open(out + '.occpeaks.bed','a')
    try:
        for poslist in iter(pos_queue.get, 'STOP'):
            for pos in poslist:
                pos.write(out_handle)
            pos_queue.task_done()
    except Exception, e:
        print('Caught exception when writing occupancy track\n')
        traceback.print_exc()
        print()
        raise e
    out_handle.close()
    return True

def run_occ(args):
    """run occupancy calling

    """

    # modified to deliberately get the chromosome sizes from the fasta file, 
    # to not use the BAM file.
    chrs = read_chrom_sizes_from_fasta(args.fasta)

    pwm = PWM.open(args.pwm)

    # modified to optionally keep only certain chromosomes
    if args.chroms_keep is not None:
        # parse comma separated list of chromosomes
        chroms_keep = args.chroms_keep.split(',')
        print "@ NOTE: restricting analysis to chromosomes: " + ", ".join(chroms_keep)
    else:
        chroms_keep = None

    # process peaks into chunks
    chunks = ChunkList.read(args.bed, chromDict = chrs, min_offset = args.flank + args.upper/2 + max(pwm.up,pwm.down) + args.nuc_sep/2, chroms_keep = chroms_keep)
    chunks.slop(chrs, up = args.nuc_sep/2, down = args.nuc_sep/2)
    chunks.merge()
    maxQueueSize = args.cores*10
    fragment_dist = FragmentMixDistribution(0, upper = args.upper)

    print "@ calculating fragment sizes..."

    # if sizes are provided, use them
    if args.sizes is not None:
        tmp = FragmentSizes.open(args.sizes)
        fragment_dist.fragmentsizes = FragmentSizes(0, args.upper, vals = tmp.get(0,args.upper))
    # otherwise, calculate them either from the BAM or fragments file
    else:
        if args.bam is not None:
            fragment_dist.getFragmentSizes(input_file = args.bam, input_type = "bam", chunklist = chunks)
        elif args.fragments is not None:
            fragment_dist.getFragmentSizes(input_file = args.fragments, input_type = "fragments", chunklist = chunks)
        
    # plot the fragment size distribution, as in pyatac.get_sizes
    fig = plt.figure()
    plt.plot(range(fragment_dist.fragmentsizes.lower,fragment_dist.fragmentsizes.upper),
             fragment_dist.fragmentsizes.get(fragment_dist.fragmentsizes.lower,fragment_dist.fragmentsizes.upper),
             label = args.out)
    # add title
    plt.title("Fragment Size Distribution")
    plt.xlabel("Fragment Size")
    plt.ylabel("Frequency")
    fig.savefig(args.out+'.fragmentsizes.pdf')
    plt.close(fig)

    # plot a second with wider domain, max size as in ArchR
    # fig = plt.figure()
    # plt.plot(range(0,750),
    #          fragment_dist.fragmentsizes.get(0, 750),
    #          label = args.out)
    # # add title
    # plt.title("Fragment Size Distribution")
    # plt.xlabel("Fragment Size")
    # plt.ylabel("Frequency")
    # fig.savefig(args.out+'.fragmentsizes_widerange.pdf')
    # plt.close(fig)

    print "@ fitting fragment size distribution..."
    fragment_dist.modelNFR()
    fragment_dist.plotFits(args.out + '.occ_fit.pdf')
    fragment_dist.fragmentsizes.save(args.out + '.fragmentsizes.txt')

    # set input_file/input_type based on whichever file is provided
    if args.bam is not None:
        input_file = args.bam
        input_type = "bam"
    elif args.fragments is not None:
        input_file = args.fragments
        input_type = "fragments"

    params = OccupancyParameters(fragment_dist, args.upper, args.fasta, args.pwm,
                                 sep = args.nuc_sep, min_occ = args.min_occ,
                                 flank = args.flank,
                                 # refactored to take in the input file and
                                 # corresponding file type (BAM or fragments)
                                 input_file = input_file,
                                 input_type = input_type,
                                 ci = args.confidence_interval, step = args.step)
    
    params.print_parameters()
    
    print "@ calculating occupancy..."
    sets = chunks.split(items = args.cores * 5)
    pool1 = mp.Pool(processes = max(1,args.cores-1))
    out_handle1 = open(args.out + '.occ.bedgraph','w')
    out_handle1.close()
    out_handle2 = open(args.out + '.occ.lower_bound.bedgraph','w')
    out_handle2.close()
    out_handle3 = open(args.out + '.occ.upper_bound.bedgraph','w')
    out_handle3.close()
    write_queue = mp.JoinableQueue(maxsize = maxQueueSize)
    write_process = mp.Process(target = _writeOcc, args=(write_queue, args.out))
    write_process.start()
    peaks_handle = open(args.out + '.occpeaks.bed','w')
    peaks_handle.close()
    peaks_queue = mp.JoinableQueue()
    peaks_process = mp.Process(target = _writePeaks, args=(peaks_queue, args.out))
    peaks_process.start()
    nuc_dist = np.zeros(args.upper)

    # NOTE: _occHlper takes a tuple as input, containing the chunk and the parameters.
    # This is the part that actually executes the occupancy calling. It doesn't need
    # to be modified directly, but calls occ.process in nucleoatac.Occupancy.OccChunk.process which we have modified slightly.
    for j in sets:
        tmp = pool1.map(_occHelper, zip(j,itertools.repeat(params)))
        for result in tmp:
            nuc_dist += result[0]
            write_queue.put(result[1])
            peaks_queue.put(result[2])
    pool1.close()
    pool1.join()
    write_queue.put('STOP')
    peaks_queue.put('STOP')
    write_process.join()
    peaks_process.join()

    print "@ compressing and indexing output files..."
    pysam.tabix_compress(args.out + '.occpeaks.bed', args.out + '.occpeaks.bed.gz',force = True)
    shell_command('rm ' + args.out + '.occpeaks.bed')
    pysam.tabix_index(args.out + '.occpeaks.bed.gz', preset = "bed", force = True)
    for i in ('occ','occ.lower_bound','occ.upper_bound'):
        pysam.tabix_compress(args.out + '.' + i + '.bedgraph', args.out + '.'+i+'.bedgraph.gz',force = True)
        shell_command('rm ' + args.out + '.' + i + '.bedgraph')
        pysam.tabix_index(args.out + '.' + i + '.bedgraph.gz', preset = "bed", force = True)

    dist_out = FragmentSizes(0, args.upper, vals = nuc_dist)
    dist_out.save(args.out + '.nuc_dist.txt')

    print "@ making figure"
    #make figure
    fig = plt.figure()
    plt.plot(range(0,args.upper),dist_out.get(0,args.upper),label = "Nucleosome Distribution")
    # add title
    plt.title("Nucleosome Distribution")
    plt.xlabel("Fragment Size")
    plt.ylabel("Frequency")
    fig.savefig(args.out+'.nuc_dist.pdf')
    plt.close(fig)

    print "@ done."











