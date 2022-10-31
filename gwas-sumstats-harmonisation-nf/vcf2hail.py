#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
VCF2Hail: a tool for converting a VCF file to Hail MatrixTable
"""

import hail as hl
import argparse

SCRIPT_DESCRIPTION = """
#==============================================================================
# VCF2Hail: a tool for converting a VCF file to Hail MatrixTable
#
# Author: David Piñeyro <davidp@lifebit.ai>
# Date: 2022-03-08
# Version: 0.0.1
#==============================================================================
"""

# Mapping for GRCh38 contigs
CONTIG_RECODING = {
    '1': 'chr1',
    '2': 'chr2',
    '3': 'chr3',
    '4': 'chr4',
    '5': 'chr5',
    '6': 'chr6',
    '7': 'chr7',
    '8': 'chr8',
    '9': 'chr9',
    '10': 'chr10',
    '11': 'chr11',
    '12': 'chr12',
    '13': 'chr13',
    '14': 'chr14',
    '15': 'chr15',
    '16': 'chr16',
    '17': 'chr17',
    '18': 'chr18',
    '19': 'chr19',
    '20': 'chr20',
    '21': 'chr21',
    '22': 'chr22'
}


def arguments():
    """This function uses argparse functionality to collect arguments."""
    parser = argparse.ArgumentParser(
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description=SCRIPT_DESCRIPTION)
    parser.add_argument('-v', '--vcf',
                        metavar='<str: input.vcf>',
                        type=str,
                        required=True,
                        help="""Path to the input VCF file.""")
    parser.add_argument('-g', '--genome_build',
                        metavar='<str: GRCh37 or GRCh38>',
                        type=str,
                        required=True,
                        help="""Either 'GRCh37' or 'GRCh38' indicating the
                                reference build of the input VCF.""")
    parser.add_argument('-o', '--out',
                        metavar='<str: output.mt',
                        type=str,
                        required=True,
                        help="""Path to the output Hail MatrixTable. It is
                                expected to be a folder.""")
    args = parser.parse_args()
    return args


###############################################################################
# MAIN
def main_program():
    """Main program."""
    args = arguments()

    if args.genome_build == 'GRCh38':
        mt = hl.import_vcf(args.vcf,
                           reference_genome=args.genome_build,
                           contig_recoding=CONTIG_RECODING)
    else:
        mt = hl.import_vcf(args.vcf,
                           reference_genome=args.genome_build)

    mt.write(args.out)


###############################################################################
# Conditional to run the script
if __name__ == '__main__':
    main_program()
