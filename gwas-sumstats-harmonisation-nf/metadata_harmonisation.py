#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Metadata harmonisation
"""

import argparse
import pandas as pd

SCRIPT_DESCRIPTION = """
#==============================================================================
# metadata_harmonisation: a tool for harmonising metadata from EBI and IEU
#
# Author: David Piñeyro <davidp@lifebit.ai>
# Date: 2022-02-22
# Version: 0.0.1
#==============================================================================
"""

HARMONISED_FIELDS = [
    'pubmedId',
    'title',
    'author',
    'publication',
    'publicationDate',
    'associationCount',
    'note',
    'mr',
    'year',
    'group_name',
    'consortium',
    'sex',
    'priority',
    'population',
    'unit',
    'nsnp',
    'sample_size',
    'build',
    'id',
    'subcategory',
    'category',
    'ontology',
    'trait',
    'mappedLabel',
    'mappedUri'
]


def arguments():
    """This function uses argparse functionality to collect arguments."""
    parser = argparse.ArgumentParser(
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description=SCRIPT_DESCRIPTION)
    parser.add_argument('-i', '--input',
                        metavar='<str: input.tsv>',
                        type=str,
                        required=True,
                        help="""Path to the input TSV file with the original
                                metadata.""")
    parser.add_argument('-o', '--output',
                        metavar='<str: output.txt',
                        type=str,
                        required=True,
                        help="""Path to the output file with the metadata
                                comments to be added to the final VCF.""")
    args = parser.parse_args()
    return args


def harmonise(df, fields):
    """Harmonisation of GWAS metadata.

    Parameters
    ----------
    df : pandas.DataFrame
        The original metadata dataframe.
    fields : list
        A list of strings with the target harmonised fileds.
    Return
    ------
    out_l : list
        A list of metadata fields to be inserted in a VCF.
    """
    # Harmonise names
    df.rename(columns={'author_s': 'author',
                       'accessionId': 'id',
                       'traitName_s': 'trait'},
              inplace=True)
    out_l = []
    for c in fields:
        new_meta = '##meta=<' + c + '='
        if c in df:
            if pd.isna(df[c].values[0]):
                new_meta += 'NA>'
            else:
                new_meta += str(df[c].values[0]) + '>'
        else:
            new_meta += 'NA>'
        out_l.append(new_meta)
    return out_l


###############################################################################
# MAIN
def main_program():
    """Main program."""
    args = arguments()

    # Load input data
    input_df = pd.read_csv(args.input, sep='\t')

    # Harmonisation
    output_lines = harmonise(input_df, HARMONISED_FIELDS)

    # Write the output file
    with open(args.output, 'w') as out_f:
        for meta in output_lines:
            out_f.write(meta + '\n')

    print(f'Output harmonised metadata file written to: {args.output}')
    print('Script completed successfully!')


###############################################################################
# Conditional to run the script
if __name__ == '__main__':
    main_program()
