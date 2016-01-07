#!/usr/bin/env python3

# Program: detect_problem_scafs.py
# Purpose:
# Authors: Daniel Mapleson.


import sys
import os
import argparse
from Bio import SeqIO
from Bio import SeqRecord
from Bio import Seq
import re

class Problem:
    scaf = ""
    scaf_len = 0
    start = True
    prob_len = 0

parser=argparse.ArgumentParser("Detects scaffolds that contain annotations that are difficult to fix by looking at GFF entries on the modified scaffolds that overlap the modified regions")
parser.add_argument("bad_scafs", help="The bad scaffolds output from \"fasta_scaffolds2contigs\"")
parser.add_argument("gff", help="The GFF annotation")
args=parser.parse_args()

# Load problems into dict
problems = {}
with open(args.bad_scafs, "r") as bad_scafs:
    for line in bad_scafs:
        parts = line.split(sep="\t")
        p = Problem()
        p.scaf = parts[0]
        p.scaf_len = int(parts[1])
        p.start = parts[2].startswith("start")
        p.prob_len = int(parts[3])
        problems[parts[0]] = p

with open(args.gff, "r") as gff:
    for line in gff:
        parts = line.split(sep="\t")
        if parts[0] in problems:
            # This is a problem scaffold (now see if there is an annotation in the problem region
            p = problems[parts[0]]
            if p.start:
                start = 1
                end = p.prob_len
            else:
                start = p.scaf_len - p.prob_len
                end = p.scaf_len

            gff_start = int(parts[3])
            gff_end = int(parts[4])

            if gff_start >= start and gff_start <= end or gff_end >= start and gff_end <= end:
                print("Problem in scaffold: " + p.scaf)
