#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script for obtaining score files from www.pgscatalog.org
"""

import argparse
import json
import urllib.request
import shutil
from contextlib import closing

SCRIPT_DESCRIPTION = """
#==============================================================================
# script for retrieving PGS score file from json retrieved from 
# www.pgscatalog.org/rest/score/search GET endpoint
#
# Author: Inês Mendes <ines@lifebit.ai>
# Date: 2022-03-15
# Version: 0.0.1
#==============================================================================
"""

def arguments():

    parser = argparse.ArgumentParser(
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description=SCRIPT_DESCRIPTION)
    parser.add_argument('-j', '--json',
                        type=str,
                        required=True,
                        help="""Path to the json input file.""")
    args = parser.parse_args()
    return args

def download_score_file(file_url):

    with closing(urllib.request.urlopen(file_url)) as r:
        with open(file_url.split('/')[-1], 'wb') as f:
            shutil.copyfileobj(r, f)

def main():

    args = arguments()

    # Opening JSON file
    f = open(args.json)

    # Load PGS json catalog
    pgs_catalog = json.load(f)

    score_file = ""
    sample_number = 0
    metadata = ""

    for result in pgs_catalog['results']:
        # Check whether there is sample number information
        n_present = len(result['samples_variants']) > 0
        # And if it's bigger than the previous one, update values
        if n_present and result['samples_variants'][0]['sample_number'] > sample_number:
            sample_number = result['samples_variants'][0]['sample_number']
            score_file = result['ftp_scoring_file']
            metadata = result['samples_variants'][0]
    # In case no sample number for any of the studies, pick the first one
    if score_file == "":
        score_file = pgs_catalog['results'][0]['ftp_scoring_file']
    # download score file
    download_score_file(score_file)

    # save metadata file
    metadata_file = score_file.split('/')[-1].split('.')[0] + '.json'
    with open(metadata_file, 'w') as outfile:
        json.dump(metadata, outfile)

if __name__ == '__main__':
    main()
