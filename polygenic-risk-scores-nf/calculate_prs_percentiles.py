#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Script for calculating PRS percentiles
"""

import argparse
import pandas as pd
from scipy import stats 

SCRIPT_DESCRIPTION = """
#==============================================================================
# script for calculating PRS percentiles for individual PRS scores
#
# Author: Eva Gradovich <eva@lifebit.ai>
# Date: 2022-03-02
# Version: 0.0.1
#==============================================================================
"""


def arguments():

    parser = argparse.ArgumentParser(
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description=SCRIPT_DESCRIPTION)
    parser.add_argument('-prs', '--input_prs_scores',
                        metavar='<str: prs_scores.tsv>',
                        type=str,
                        required=True,
                        help="""Path to the input file containing individual PRS scores.""")
    parser.add_argument('-out', '--output_prs_df',
                        metavar='<str: output_prs_df.tsv',
                        type=str,
                        required=True,
                        help="""Path to the output file containing PRS percentile info.""")
    parser.add_argument('-lower_pp', '--lower_prs_percentile',
                        metavar='<int: 0',
                        type=int,
                        required=True,
                        help="""Lower prs percentile for filtering.""")
    parser.add_argument('-upper_pp', '--upper_prs_percentile',
                        metavar='<int: 0',
                        type=int,
                        required=True,
                        help="""Upper prs percentile for filtering.""")
    parser.add_argument('-prs_col', '--prs_column_name',
                        metavar='<str: PRS',
                        type=str,
                        required=True,
                        help="""Name of column containing PRS scores.""")
    args = parser.parse_args()
    return args


def calculate_percentiles(df, prs_column_name):
    """Calculating PRS percentiles 
    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe containing individual PRS scores.
    prs_column_name : string
        String value corresponding to column name containing PRS.
    Return
    ------
    df : pandas.DataFrame
        Updated dataframe containing PRS percentiles as well as individual PRS scores.
    """

    prs_scores = df[prs_column_name].tolist()
    df['prs_percentile'] = df.apply(lambda row : stats.percentileofscore(prs_scores,row[prs_column_name]), axis = 1)
    return df

def filter_samples_by_percentiles(df, lower_percentile, upper_percentile):
    """Filtering cohort samples by PRS percentiles 

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe containing individual PRS scores.
    Return
    ------
    df : pandas.DataFrame
        Updated dataframe containing PRS percentiles as well as individual PRS scores.
    """
    filtered_df = df[df['prs_percentile'].between(lower_percentile, upper_percentile, inclusive=True)]
    filtered_df['FID'] = filtered_df['IID']
    filtered_samples_keep = filtered_df[['FID','IID']]
    return filtered_samples_keep

def main_program():

    args = arguments()

    # Load individual PRS scores. Using sep parameter to autodetect separator
    prs_score_df = pd.read_csv(args.input_prs_scores, sep=None)

    # Get PRS percentiles
    prs_full_df = calculate_percentiles(prs_score_df, args.prs_column_name)

    filtered_samples_keep_file = filter_samples_by_percentiles(prs_full_df, args.lower_prs_percentile, args.upper_prs_percentile)

    prs_full_df.to_csv(args.output_prs_df,index=False)
    filtered_samples_keep_file.to_csv("samples_filtered_percentile.tsv", header=None,index=False,sep='\t')

    print("PRS scores and percentiles written to: {}".format(args.output_prs_df))


if __name__ == '__main__':
    main_program()
