#! /usr/bin/env python3

# This script collapses rows of a sorted VCF with multialleic sites represented as multiple rows
# into a a single row for each multiallelic site. e.g.:
# This:
# 1	1097335	rs9442385	T	A	.	PASS	AF=0.013379	ES:SE:LP:AF:SS	0.003:0.015:0.075:0.013:1097
# 1	1097335	rs9442385	T	G	.	PASS	AF=0.834665	ES:SE:LP:AF:SS	0.002:0.004:0.187:0.930:1097
#
# Into:
# 1	1097335	rs9442385	T	A,G	.	PASS	AF=0.013379,0.834665	ES:SE:LP:AF:SS	0.003,0.002:0.015,0.004:0.075,0.187:0.013,0.930:1097,1097

import sys
import argparse

def merge_format_col(format_keys, r1_format_keys, r2_format_keys, r1_sample_col, r2_sample_col):
    '''merge the format values across two alleles for a particular sample column'''
    r1_format = { i[0]:i[1] for i in zip(r1_format_keys, r1_sample_col.split(':')) }
    r2_format = { i[0]:i[1] for i in zip(r2_format_keys, r2_sample_col.split(':')) }
    rnew_format = { k:f"{r1_format.get(k,'.')},{r2_format.get(k,'.')}" for k in format_keys }
    rnew_sample_col = ':'.join([ rnew_format[k] for k in format_keys ])
    return rnew_sample_col


def merge_allele_records(r1, r2):
    '''merge the two vcf lines into a single line (assuming different alleles of the same SNP)'''
    rnew = r1
    rnew[4] = f'{r1[4]},{r2[4]}'

    # combine info
    if r1[7] != '.':
        r1_info = { i.split('=')[0]:i.split('=')[1] for i in r1[7].split(';') }
    else:
        r1_info = dict()
    
    if r2[7] != '.':
        r2_info = { i.split('=')[0]:i.split('=')[1] for i in r2[7].split(';') }
    else:
        r2_info = dict()
    
    if len(r1_info) + len(r2_info) > 0:
        info_keys = list(set(r1_info.keys() + r2_info.keys()))
        rnew_info = { k:f"{r1_info.get(k,'.')},{r2_info.get(k,'.')}" for k in info_keys }
        rnew[7] = ';'.join([ f"{k}={v}" for k,v in rnew_info.items() ])
    else:
        rnew[7] = '.'


    # combine format keys
    r1_format_keys = r1[8].split(':')
    r2_format_keys = r2[8].split(':')
    # using set(r1_format_keys + r2_format_keys) would be more efficient
    # but we want to keep key order as stable as possible
    format_keys = r1_format_keys + [ k for k in r2_format_keys if k not in r1_format_keys]
    rnew[8] = ':'.join(format_keys)

    # combine format values in multiple sample columns
    for s_i in range(9, len(rnew)):
        rnew[s_i] = merge_format_col(format_keys, r1_format_keys, r2_format_keys, r1[s_i], r2[s_i])

    return rnew


def main(args):
    prev_id = ('.','.','.')
    prev = [""]

    with args.infile as f:
        for line in f:

            if line[0] == '#':
                if line[:2] == '##':
                    # print header as is
                    print(line.strip(), file=args.outfile)
                    continue
                else:
                    # this is the col names line so put put col names into prev
                    # to be printed in the next iteration
                    prev = line.strip().split('\t')
                    continue
            
            this = line.strip().split('\t')
            this_id = (this[0], this[1], this[3])

            if prev_id != this_id:
                # new line is not same SNP as previous so we can now print
                # previous record as we know there is nothing more to add to it
                print(*prev, sep='\t', file=args.outfile)

            else:
                # new line is same SNP as previous so we must combine 'this' with 'prev'
                this = merge_allele_records(prev, this)

            prev = this
            prev_id = this_id

        # print the final line that is now stored in prev
        print(*prev, sep='\t', file=args.outfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?',
        type=argparse.FileType('r'), default=sys.stdin,
        help="""Input VCF file.""")
    parser.add_argument('-o', '--outfile',
        type=argparse.FileType('w'), default=sys.stdout,
        help="""Output filename.""")
    args = parser.parse_args()

    main(args)
