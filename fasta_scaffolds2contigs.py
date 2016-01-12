#!/usr/bin/env python3

# Program: fasta_scaffolds2contigs.py
# Purpose:
# Authors: Daniel Mapleson.  Derived from work by Richard Leggett, Mario Caccamo, which was originally based on abyss-fatoagp by Shaun Jackman


import sys
import os
import argparse
from Bio import SeqIO
from Bio import SeqRecord
from Bio import Seq
import re

# Minimum contig length in ENA submissions
ENA_MIN_LEN = 20

class AGP2:

    # Properties common to all component types
    object = ""
    object_beg = 0
    object_end = 0
    part_number = 0
    comp_type = ""

    # These properties assume comp_type is NOT set to 'N' or 'U'
    comp_id = ""
    comp_beg = 0
    comp_end = 0
    orientation = "?"

    # These properties assume comp_type IS set to 'N' or 'U'
    gap_length = 0
    gap_type = ""
    linkage = ""
    linkage_evidence = ""

    # Extra properties to enable sorting
    object_id = 0


    def __init__(self, *args):

        if len(args) == 0:
            self.object = ""
            self.object_beg = 0
            self.object_end = 0
            self.part_number = 0
            self.comp_type = ""

            self.gap_length = 0
            self.gap_type = ""
            self.linkage = ""
            self.linkage_evidence = ""

            self.comp_id = ""
            self.comp_beg = 0
            self.comp_end = 0
            self.orientation = ""

            self.object_id = 0
            return

        else:
            line = args[0]

        self.data = []
        parts = line.split("\t")

        self.object = parts[0]
        self.object_beg = int(parts[1])
        self.object_end = int(parts[2])
        self.part_number = int(parts[3])
        self.comp_type = parts[4]

        if (parts[4] == 'N' or parts[4] == 'U'):
            self.gap_length = int(parts[5])
            self.gap_type = parts[6]
            self.linkage = parts[7]
            self.linkage_evidence = parts[8].strip()    # Part 8 Will include newline
        else:
            self.comp_id = parts[5]
            self.comp_beg = int(parts[6])
            self.comp_end = int(parts[7])
            self.orientation = parts[8]

        self.object_id = 0

    def __str__(self):
        s = self.object + "\t" + str(self.object_beg) + "\t" + str(self.object_end) + "\t" + str(self.part_number) + "\t" + self.comp_type + "\t"
        if self.isGap():
            s += str(self.gap_length) + "\t" + self.gap_type + "\t" + self.linkage + "\t" + self.linkage_evidence
        else:
            s += self.comp_id + "\t" + str(self.comp_beg) + "\t" + str(self.comp_end) + "\t" + self.orientation

        return s

    def isGap(self):
        return self.comp_type == "N" or self.comp_type == "U"

    def __lt__(self, other):

        if self.object_id < other.object_id:
            return True
        elif self.object_id > other.object_id:
            return False
        else:
            if self.object_beg < other.object_beg:
                return True
            else:
                return False    # Think we can get away with this, we should have two AGP objects with same name with the same start coord.


def checkIfGap(lst):
    return lst[1:] == lst[:-1] and lst.startswith("N")

def createGap(length):
    newgap = ""
    for i in range(0,length):
        newgap += "N"

    return newgap

def trimNsFromStart(seq):
    newseq = ""
    start = True
    count = 0
    for b in seq:
        if start == True and b == "N":
            # skip
            c=1
            count += 1
        else:
            start = False
            newseq += b

    return (newseq, count)

def trimNsFromEnd(seq):
    newnewseq = ""
    end = True
    count = 0
    for b in reversed(seq):
        if end == True and b == "N":
            # skip
            c=1
            count += 1
        else:
            end = False
            newnewseq += b

    return (newnewseq[::-1], count)

# This might need modifying for each job
def breakScaffoldHeader(id):
    id_parts = str(id).split("_")
    stem = "_".join(id_parts[0:len(id_parts)-2])
    scaf_index = int(id_parts[4])
    return (stem, scaf_index)


parser=argparse.ArgumentParser("Create broken contigs and AGP file from a FASTA file of scaffolds. Assumes input scaffold file contains headers ending in \"_scaffold_<index>\"  Prefer this script to the original fasta2agp script (or RAMPART) when you must not break the defined naming convention in the input scaffold file.")
parser.add_argument("input", help="The input scaffolds FastA file")
parser.add_argument("-o", "--output", required=True, help="Output prefix for output files")
parser.add_argument("-n", "--min_n", type=int, default=10, help="minimum number of Ns to break sequence (default 10)")
parser.add_argument("-s", "--emit_scaffolds", default=False, help="Whether or not to emit scaffolds derived from contigs and AGP (should be the same as the original input)")
args=parser.parse_args()

contig_file = args.output + "_contigs.fa"
scaffold_file = args.output + "_scaffolds.fa"
agp_file = args.output + ".agp"
bad_scaf_file = args.output + ".bad_scafs"

contig_index = 0
actual_contig_index = 1
agp_list = list()
contig_list = list()

bad_scafs_out = open(bad_scaf_file, "w")

print ("Loading existing scaffold fasta...", end="")
scaffolds = ()
scaffold_index_map = {}
scaffold_indicies = list()
scafs_in = open(args.input, "rU")
assembly_size = 0
test1_file = args.output + "_test1.tsv"
with open(test1_file, "w") as test1_out:
    for scaf in SeqIO.parse(scafs_in, "fasta") :

        stem, scaf_index = breakScaffoldHeader(scaf.id)
        scaffold_index_map[scaf_index] = scaf
        scaffold_indicies.append(scaf_index)
        scaf_len = len(scaf.seq)
        assembly_size += scaf_len
        print(scaf.id + "\t" + str(scaf_len), file=test1_out)

    scafs_in.close()
print(" done.")
print("Loaded " + str(len(scaffold_indicies)) + " scaffolds")
print("Assembly size: " + str(assembly_size))
print("Sorting scaffolds by index...", end="")
scaffold_indicies.sort(key=int)
print(" done.")

contig_index = 1
other_index = 1

print ("Processing scaffolds...", end="")
for scaf_index in scaffold_indicies:

    scaf = scaffold_index_map[scaf_index]
    stem, scaf_index = breakScaffoldHeader(scaf.id)
    scaf_len = len(scaf.seq)

    if scaf_len < 20:
        print ("ERROR: Not outputting scaffold \"" + scaf.id + "\" as scaffold has length < 20, which is the minimum required by ENA")
        print (scaf.id + "\t" + str(scaf_len) + "\tScaffold too short", file=bad_scafs_out)
        exit(1)

    # Split scaffold by N's
    contig_seqs = re.split("([N|n]{" + str(args.min_n) + ",})", str(scaf.seq))

    part_number = 1
    position = 1

    ###############
    # Remove empty entries at beginning and ends to make sure that there are no issues there.  The split might result
    # in the first or last contig being empty if the last but one contig was made up of Ns
    if len(contig_seqs[0]) == 0:
        contig_seqs.pop(0)

    if len(contig_seqs[-1]) == 0:
        contig_seqs.pop()

    # Make sure this is done after the previous trimming
    num_contigs = len(contig_seqs)

    # Iterate over the cleaned contig_seqs (although we might still have short contigs in the middle
    i = 1
    gap_linkage = True
    for contig in contig_seqs:

        contig_len = len(contig)

        isgap = checkIfGap(contig) # Check if this contig is a gap

        # Initialise AGP object for this instance (always the same start for this scaffold)
        agp = AGP2()
        agp.object = scaf.id
        agp.object_beg = position
        agp.object_end = position + contig_len -1
        agp.part_number = part_number
        agp.object_id = scaf_index

        # Short contig
        if len(contig) < 20 and not isgap:
            if (num_contigs > 2):
                contig_id = stem + "_other_" + str(other_index)
                agp.comp_type = "O"
                agp.comp_id = contig_id
                agp.comp_beg = 1            # Can't be broken down further so always start at 1
                agp.comp_end = contig_len
                agp.orientation = "+"
                other_index += 1

                if i == 1:
                    gap_linkage = False
                    print ("WARNING: " + scaf.id + " starts with a short contig.  Using \"O\" component type for these contigs.  Length: " + str(contig_len))
                    print (scaf.id + "\t" + str(scaf_len) + "\tstarts with a short contig\t" + str(contig_len), file=bad_scafs_out)
                elif i == num_contigs:
                    # Need to modify previous gap
                    agp_list[-1].gap_type = "scaffold"
                    agp_list[-1].linkage = "no"
                    agp_list[-1].linkage_evidence = "na"
                    print ("WARNING: " + scaf.id + " ends with a short contig.  Using \"O\" component type for these contigs.  Length: " + str(contig_len))
                    print (scaf.id + "\t" + str(scaf_len) + "\tends with a short contig\t" + str(contig_len), file=bad_scafs_out)
                else:
                    print("ERROR: short contig in the middle of a scaffold.  Not sure what to do with this yet!")
                    exit(1)

            else:
                print ("ERROR: Scaffold " + scaf.id + " is too short!")
                exit(1)

        # Starts or ends with Ns
        if not isgap:

            if contig[0] == "N":

                if not i == 1:
                    print("ERROR: found contig starting with an N in the middle of scaffold: " + scaf.id + ".  Length: " + str(scaf_len) + ". Position: " + str(position) + ".  Not sure what to do with this yet!")
                    exit(1)

                mod_contig, gap_len = trimNsFromStart(contig)
                mod_contig_len = len(mod_contig)

                if contig_len < 20:
                    print("Modified contig length is shorter than that required by ENA.  Not sure what to do with this yet!")
                    exit(1)

                # Create a gap before actual contig
                agp.object_end = gap_len
                agp.comp_type = "N"
                agp.gap_length = gap_len
                agp.gap_type = "scaffold"
                agp.linkage = "no"
                agp.linkage_evidence = "na"
                agp_list.append(agp)

                part_number += 1
                position += gap_len

                agp = AGP2()
                agp.object = scaf.id
                agp.object_beg = position
                agp.object_end = position + mod_contig_len -1
                agp.part_number = part_number
                agp.object_id = scaf_index
                # Can carry on as normal now (contig has been modified)

            elif contig[-1] == "N":

                if not i == num_contigs:
                    print("ERROR: found contig ending with an N in the middle of scaffold: " + scaf.id + ".  Length: " + str(scaf_len) + ". Position: " + str(position) + ".  Not sure what to do with this yet!")
                    exit(1)

                mod_contig, gap_len = trimNsFromEnd(contig)
                mod_contig_len = len(mod_contig)
                contig_id = stem + "_contig_" + str(contig_index)

                if contig_len < 20:
                    print("Modified contig length is shorter than that required by ENA.  Not sure what to do with this yet!")
                    exit(1)

                agp.object_end = position + mod_contig_len -1
                agp.comp_type = "W"
                agp.comp_id = contig_id
                agp.comp_beg = 1            # Can't be broken down further so always start at 1
                agp.comp_end = mod_contig_len
                agp.orientation = "+"   # Always on the positive strand
                agp_list.append(agp)
                record = SeqRecord.SeqRecord(Seq.Seq(mod_contig), id=contig_id, description="")
                contig_list.append(record)
                contig_index += 1
                part_number += 1
                position += mod_contig_len

                # Now output the gap
                agp = AGP2()
                agp.object = scaf.id
                agp.object_beg = position
                agp.object_end = position + gap_len -1
                agp.part_number = part_number
                agp.object_id = scaf_index
                agp.comp_type = "N"
                agp.gap_length = gap_len
                agp.gap_type = "scaffold"
                agp.linkage = "no"
                agp.linkage_evidence = "na"
                agp_list.append(agp)

                part_number += 1
                position += gap_len
                break
            else:
                mod_contig = contig
                mod_contig_len = contig_len

            # Normal contig
            contig_id = stem + "_contig_" + str(contig_index)
            agp.comp_type = "W"
            agp.comp_id = contig_id
            agp.comp_beg = 1            # Can't be broken down further so always start at 1
            agp.comp_end = mod_contig_len
            agp.orientation = "+"   # Always on the positive strand
            agp_list.append(agp)
            gap_linkage = True

            record = SeqRecord.SeqRecord(Seq.Seq(mod_contig), id=contig_id, description="")
            contig_list.append(record)
            contig_index += 1
            part_number += 1
            position += mod_contig_len


        # Is gap
        else:

            # If we had a short contig previously, or this is a gap at the start or end of a scaffold, then assume this gap is probably a repeat
            if not gap_linkage or i == 1 or i == num_contigs:
                agp.comp_type = "N"
                agp.gap_length = contig_len
                agp.gap_type = "scaffold"
                agp.linkage = "no"
                agp.linkage_evidence = "na"
                gap_linkage = True
                print ("WARNING: " + scaf.id + " contains an unexpected gap.  Length: " + str(contig_len))
                print (scaf.id + "\t" + str(scaf_len) + "\tcontains repeat gap\t" + str(contig_len), file=bad_scafs_out)
            # Otherwise this is a normal gap caused by scaffolding
            else:
                agp.comp_type = "N"
                agp.gap_length = contig_len
                agp.gap_type = "scaffold"
                agp.linkage = "yes"
                agp.linkage_evidence = "paired-ends"
                gap_linkage = True

            agp_list.append(agp)
            part_number += 1
            position += contig_len

        i += 1


print ("Saving AGP...")
with open(agp_file, "w") as agp_out:
    for agp in agp_list:
        print(agp, file=agp_out)

print ("Saving contigs...")
with open(contig_file, "w") as contigs_out:
    for contig in contig_list:
        SeqIO.write(contig, contigs_out, "fasta") # Urgh!


# Mainly just for validation purposes
if (args.emit_scaffolds):
    # This is horrendous but will do for now (I have the memory!)
    # Make map of id to record first
    print ("Indexing contigs...", end="")
    fasta_map = {}
    id_list = list()
    for contig in contig_list:
        fasta_map[str(contig.id)] = contig
    print (" done.")

    print ("Saving scaffolds...")
    num_scafs = 0
    asm_size = 0
    test2_file = args.output + "_test2.tsv"
    with open(test2_file, "w") as test2_out:
        with open(scaffold_file, "w") as scaffolds_out:
            previous_id = agp_list[0].object
            current_id = agp_list[0].object
            scaffold = ""
            for agp in agp_list:

                if agp.object == current_id:
                    if (agp.isGap()):
                        scaffold += createGap(agp.gap_length)
                    else:
                        scaffold += str(fasta_map[agp.comp_id].seq)
                else:
                    record = SeqRecord.SeqRecord(Seq.Seq(scaffold), id=previous_id, description="")
                    SeqIO.write(record, scaffolds_out, "fasta") # Urgh!
                    previous_id = current_id
                    current_id = agp.object
                    num_scafs += 1
                    asm_size += len(scaffold)
                    print("" + previous_id + "\t" + str(len(scaffold)), file=test2_out)
                    scaffold = ""
                    if (agp.isGap()):
                        scaffold += createGap(agp.gap_length)
                    else:
                        scaffold += str(fasta_map[agp.comp_id].seq)


            # Output last scaffold
            record = SeqRecord.SeqRecord(Seq.Seq(scaffold), id=previous_id, description="")
            SeqIO.write(record, scaffolds_out, "fasta") # Urgh!
            num_scafs += 1
            asm_size += len(scaffold)
            print("" + previous_id + "\t" + str(len(scaffold)), file=test2_out)
            scaffold = ""

    print (" done.")

    print ("Wrote out " + str(num_scafs) + " scaffolds.  Assembly size: " + str(asm_size))

    print("Delta scafs: " + str(num_scafs - len(scaffold_indicies)) + ". Delta asm_size: " + str(asm_size - assembly_size))
