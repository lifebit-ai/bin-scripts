#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
convert_coeff: a tool for converting GWAS coefficients
"""

import argparse
import math
import numpy as np
import statistics

SCRIPT_DESCRIPTION = """
#==============================================================================
# convert_coeff: a tool for converting GWAS coefficients
#
# Author: David Piñeyro <davidp@lifebit.ai>
# Date: 2022-02-09
# Version: 0.0.1
#==============================================================================
"""


def arguments():
    """This function uses argparse functionality to collect arguments."""
    parser = argparse.ArgumentParser(
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description=SCRIPT_DESCRIPTION)
    parser.add_argument('-v', '--vcf',
                        metavar='<str: input.vcf>',
                        type=str,
                        required=True,
                        help="""Path to the input VCF file where the traits
                                comments will be added as part of
                                its header.""")
    parser.add_argument('-o', '--out_vcf',
                        metavar='<str: output.vcf',
                        type=str,
                        required=True,
                        help="""Path to the output VCF file with the traits
                                added as part of its header.""")
    parser.add_argument('--standardise',
                        action='store_true',
                        help="""Whether to standardise BETA and SE columns.""")
    parser.add_argument('--beta2or',
                        action='store_true',
                        help="""Whether to convert BETA to Odds Ratio.""")
    parser.add_argument('--or2beta',
                        action='store_true',
                        help="""Whether to convert Odds Ratio to BETA.""")
    parser.add_argument('--beta-col',
                        metavar='<str: COLNAME>',
                        type=str,
                        default='ES',
                        help="""BETA colname.""")
    parser.add_argument('--beta-se-col',
                        metavar='<str: COLNAME>',
                        type=str,
                        default='SE',
                        help="""BETA Standard Error colname.""")
    parser.add_argument('--or-col',
                        metavar='<str: COLNAME>',
                        type=str,
                        default='OR',
                        help="""Odds Ratio colname.""")
    parser.add_argument('--freq-col',
                        metavar='<str: COLNAME>',
                        type=str,
                        default='AF',
                        help="""Allele frequency colname. Required for BETA
                                standardisation.""")
    args = parser.parse_args()
    print('\nThe following operations will be performed in:', args.vcf)
    print('\tBETA standardisation ---->', args.standardise)
    print('\tBETA SE standardisation ->', args.standardise)
    print('\tBETA to OR conversion --->', args.beta2or)
    print('\tOR to BETA conversion --->', args.or2beta)
    print('Output VCF to:', args.out_vcf, '\n')
    return args


def modify_header(header, stnd, beta2or, or2beta):
    """Performs header additions when required.

    Parameters
    ----------
    header : list
        A list of strings from the original VCF. They still conserve
        all new line and space characters from the original.
    stnd : bool
        Whether to apply standardisation of BETAS and SEs.
    beta2or : bool
        Whether to convert BETA to Odds Ratio.
    or2beta : bool
        Whether to convert Odds Ratio to Beta.
    Return
    ------
    new_header : list
        The strings to be printed to the VCF header.
    """
    format_found = False
    format_written = False
    new_header = []
    for e in header:
        if e[:8] == '##FORMAT':
            format_found = True
        elif format_found and not format_written:
            if or2beta:
                new_e = '##FORMAT=<ID=ES,Number=A,Type=Float,Description="Effect size estimate relative to the alternative allele">\n'
                e = new_e + e
                format_written = True
            if beta2or:
                new_e = '##FORMAT=<ID=OR,Number=A,Type=Float,Description="Odds ratio of effect">\n'
                e = new_e + e
                format_written = True
            if stnd:
                new_e = '##FORMAT=<ID=SES,Number=A,Type=Float,Description="Standardised effect size">\n'
                new_e2 = '##FORMAT=<ID=SSE,Number=A,Type=Float,Description="Standardised standard error of effect size">\n'
                e = new_e + new_e2 + e
                format_written = True
        new_header.append(e)
    return new_header


def standardise_beta(beta, se, median_phe_se):
    """BETA and SE standardisation.

    Parameters
    ----------
    beta : str or float
        BETA (effect size).
    median_phe_se : float or None
        Median phenotypic standard error calculated from the input VCF.
    Return
    ------
    stnd_beta : float
        Standardised BETA.
    stnd_se : float
        Standardised BETA SE.
    """
    stnd_beta = float(beta) / median_phe_se
    stnd_se = float(se) / median_phe_se
    return round(stnd_beta, 6), round(stnd_se, 6)


def convert_beta2or(beta):
    """Converts BETA to OR.

    Based on documentation here:
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3154648/
    https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS9.html

    Parameters
    ----------
    beta : str or float
        BETA (effect size).
    Return
    ------
    or_ : float
        Odds ratio.
    """
    or_ = np.exp(float(beta))
    return round(or_, 6)


def convert_or2beta(oddsr):
    """Converts OR to BETA.

    Based on documentation here:
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3154648/
    https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS9.html

    Parameters
    ----------
    oddsr : str or float
        Odds ratio.
    Return
    ------
    beta_ : float
        BETA computed from OR.
    """
    beta_ = math.log(float(oddsr))
    return round(beta_, 6)


def modify_record(line, stnd, beta2or, or2beta, median_phe_se, beta_col,
                  beta_se_col, or_col):
    """Performs additions to the FORMAT and subsequent columns.

    Parameters
    ----------
    line : string
        A line record from the original VCF. It still conserves all
        new line and space characterds from the original.
    stnd : bool
        Whether to apply standardisation of BETAS and SEs.
    beta2or : bool
        Whether to convert BETA to Odds Ratio.
    or2beta : bool
        Whether to convert Odds Ratio to Beta.
    median_phe_se : float or None
        Median phenotypic standard error calculated from the input VCF.
    beta_col : str
        Column name of the BETA field.
    beta_se_col : str
        Column name of the BETA SE field.
    or_col : str
        Column name of the Odds Ration field.
    Return
    ------
    new_line : string
        A new line containing the required new information.
    """
    line_list = line.rstrip().split('\t')
    format_col = line_list[8]
    new_line = line_list
    if stnd:
        # Looking for BETA column
        f_col_list = format_col.split(':')
        try:
            beta_ix = f_col_list.index(beta_col)
        except ValueError as err:
            print(f'[ERROR] The input VCF does not contain {beta_col} in FORMAT column.')
            print('[ERROR] BETA standardisation not possible.')
            raise err
        try:
            se_ix = f_col_list.index(beta_se_col)
        except ValueError as err:
            print(f'[ERROR] The input VCF does not contain {beta_se_col} in FORMAT column.')
            print('[ERROR] BETA standardisation not possible.')
            raise err
        # Add extra values to FORMAT
        new_line[8] += ':SES:SSE'
        for i in range(9, len(line_list)):
            line_sample_list = line_list[i].split(':')
            sample_beta = line_sample_list[beta_ix]
            sample_se = line_sample_list[se_ix]
            [stnd_beta, stnd_se] = standardise_beta(sample_beta, sample_se, median_phe_se)
            new_line[i] += f':{stnd_beta}:{stnd_se}'
    if beta2or:
        # Looking for BETA column
        f_col_list = format_col.split(':')
        try:
            beta_ix = f_col_list.index(beta_col)
        except ValueError as err:
            print(f'[ERROR] The input VCF does not contain {beta_col} in FORMAT column.')
            print('[ERROR] BETA to OR conversion not possible.')
            raise err
        # Add extra values to FORMAT
        new_line[8] += ':OR'
        for i in range(9, len(line_list)):
            line_sample_list = line_list[i].split(':')
            sample_beta = line_sample_list[beta_ix]
            or_ = convert_beta2or(sample_beta)
            new_line[i] += f':{or_}'
    if or2beta:
        # Looking for OR column
        f_col_list = format_col.split(':')
        try:
            or_ix = f_col_list.index(or_col)
        except ValueError as err:
            print(f'[ERROR] The input VCF does not contain {or_col} in FORMAT column.')
            print('[ERROR] OR to BETA conversion not possible.')
            raise err
        # Add extra values to FORMAT
        new_line[8] += ':ES'
        for i in range(9, len(line_list)):
            line_sample_list = line_list[i].split(':')
            sample_or = line_sample_list[or_ix]
            beta_ = convert_or2beta(sample_or)
            new_line[i] += f':{beta_}'
    new_line[-1] += '\n'
    return ('\t').join(new_line)


def process_vcf(input_vcf, out_vcf, stnd, beta2or, or2beta, median_phe_se, beta_col, beta_se_col,
                or_col):
    """Function to read a VCF file, do the requested conversions and
    write an output VCF.

    Parameters
    ----------
    input_vcf : str
        Path to the input VCF.
    out_vcf : str
        Path to write the output VCF file.
    stnd : bool
        Whether to apply standardisation of BETAS and SEs.
    beta2or : bool
        Whether to convert BETA to Odds Ratio.
    or2beta : bool
        Whether to convert Odds Ratio to Beta.
    median_phe_se : float or None
        Median phenotypic standard error calculated from the input VCF.
    beta_col : str
        Column name of the BETA field.
    beta_se_col : str
        Column name of the BETA SE field.
    or_col : str
        Column name of the Odds Ration field.
    Return
    ------
    out_vcf : str
        Path to the output VCF file.
    """
    header = []
    final_vcf = open(out_vcf, 'w')
    header_written = False
    with open(input_vcf, 'r') as vcf:
        for line in vcf:
            # Capture all header.
            if line[0] == '#':
                header.append(line)
            else:
                if not header_written:
                    # Modify header with the new required lines and write it.
                    new_header = modify_header(header, stnd, beta2or, or2beta)
                    for r in new_header:
                        final_vcf.write(r)
                    header_written = True
                new_line = modify_record(line, stnd, beta2or, or2beta,
                                         median_phe_se, beta_col, beta_se_col,
                                         or_col)
                final_vcf.write(new_line)
    final_vcf.close()
    return out_vcf


def get_median_phe_se(input_vcf, beta_se_col, freq_col):
    """Computes the median phenotypic standard error from the input VCF.

    Parameters
    ----------
    input_vcf : str
        Path to the input VCF.
    beta_se_col : str
        Column name of the BETA SE field.
    freq_col : str
        Colname of the allele frequency column.
    Return
    ------
    median_phe_se : float
        Median phenotypic standard error.
    """
    phe_se = []
    # Calculate n_snps
    n_snps = 0
    with open(input_vcf, 'r') as vcf:
        for line in vcf:
            if line[0] != '#':
                n_snps += 1
    with open(input_vcf, 'r') as vcf:
        for line in vcf:
            if line[0] != '#':
                line_l = line.rstrip().split('\t')
                format_col = line_l[8].split(':')
                try:
                    se_ix = format_col.index(beta_se_col)
                except ValueError as err:
                    print('[ERROR] The input VCF does not contain SE in FORMAT column.')
                    print('[ERROR] BETA standardisation not possible.')
                    raise err
                try:
                    freq_ix = format_col.index(freq_col)
                except ValueError as err:
                    print(f'[ERROR] The input VCF does not contain {freq_col} in FORMAT column.')
                    print('[ERROR] BETA standardisation not possible.')
                    raise err
                gwas_col = line_l[9].split(':')
                sample_se = gwas_col[se_ix]
                sample_freq = gwas_col[freq_ix]
                sample_phe_se = float(sample_se) * math.sqrt(2 * n_snps * float(sample_freq)
                                                             * (1 - float(sample_freq)))
                phe_se.append(sample_phe_se)
    median_phe_se = statistics.median(phe_se)
    return median_phe_se


###############################################################################
# MAIN
def main_program():
    """Main program."""
    args = arguments()

    # First step calculating median phenotypic standard error
    if args.standardise:
        median_phe_se = get_median_phe_se(args.vcf, args.beta_se_col,
                                          args.freq_col)
    else:
        median_phe_se = None

    out = process_vcf(args.vcf, args.out_vcf, args.standardise, args.beta2or,
                      args.or2beta, median_phe_se, args.beta_col,
                      args.beta_se_col, args.or_col)

    print(f'\nOutput VCF file written to: {out}')
    print('Script completed successfully!')


###############################################################################
# Conditional to run the script
if __name__ == '__main__':
    main_program()
