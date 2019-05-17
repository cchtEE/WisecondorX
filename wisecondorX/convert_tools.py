# WisecondorX

import logging

import numpy as np
import pysam

from wisecondorX.props import ChromosomeMap

'''
Converts bam file to numpy array by transforming
individual reads to counts per bin.
'''

chromosome_map = ChromosomeMap()


def convert_bam(args):
    bins_per_chr = dict()
    for chr in range(1, chromosome_map.get_max() + 1):
        bins_per_chr[str(chr)] = None

    logging.info('Importing data ...')

    bam_file = pysam.AlignmentFile(args.infile, 'rb')

    reads_seen = 0
    reads_kept = 0
    reads_mapq = 0
    reads_rmdup = 0
    reads_pairf = 0
    larp = -1
    larp2 = -1

    logging.info('Converting bam ... This might take a while ...')

    for index, chr in enumerate(bam_file.references):

        chr_name = chr
        if chr_name[:3].lower() == 'chr':
            chr_name = chr_name[3:]

        x = chromosome_map.get_as_chr_name_str(chromosome_map.get_x_index())
        y = chromosome_map.get_as_chr_name_str(chromosome_map.get_y_index())

        if chr_name not in bins_per_chr and chr_name != x and chr_name != y:
            continue

        logging.info('Working at {}; processing {} bins'
                     .format(chr, int(bam_file.lengths[index] / float(args.binsize) + 1)))
        counts = np.zeros(int(bam_file.lengths[index] / float(args.binsize) + 1), dtype=np.int32)
        bam_chr = bam_file.fetch(chr)

        if chr_name == x:
            chr_name = str(chromosome_map.get_x_index())
        if chr_name == y:
            chr_name = str(chromosome_map.get_y_index())

        if args.paired:

            for read in bam_chr:
                if not (read.is_proper_pair and read.is_read1):
                    reads_pairf += 1
                    continue
                if larp == read.pos and larp2 == read.next_reference_start:
                    reads_rmdup += 1
                else:
                    if read.mapping_quality >= 1:
                        location = read.pos / args.binsize
                        counts[int(location)] += 1
                    else:
                        reads_mapq += 1

                larp2 = read.next_reference_start
                reads_seen += 1
                larp = read.pos
        else:

            for read in bam_chr:
                if larp == read.pos:
                    reads_rmdup += 1
                else:
                    if read.mapping_quality >= 1:
                        location = read.pos / args.binsize
                        counts[int(location)] += 1
                    else:
                        reads_mapq += 1

                reads_seen += 1
                larp = read.pos

        bins_per_chr[chr_name] = counts
        reads_kept += sum(counts)

    qual_info = {'mapped': bam_file.mapped,
                 'unmapped': bam_file.unmapped,
                 'no_coordinate': bam_file.nocoordinate,
                 'filter_rmdup': reads_rmdup,
                 'filter_mapq': reads_mapq,
                 'pre_retro': reads_seen,
                 'post_retro': reads_kept,
                 'pair_fail': reads_pairf}
    return bins_per_chr, qual_info
