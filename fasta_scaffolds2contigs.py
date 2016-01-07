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


parser=argparse.ArgumentParser("Create broken contigs and AGP file from a FASTA file of scaffolds. Assumes input scaffold file contains headers ending in \"_scaffold_<index>\"  Prefer this script to the original fasta2agp script (or RAMPART) when you must not break the defined naming convention in the input scaffold file.")
parser.add_argument("input", help="The input scaffolds FastA file")
parser.add_argument("-o", "--output", required=True, help="Output prefix for output files")
parser.add_argument("-n", "--min_n", type=int, default=10, help="minimum number of Ns to break sequence (default 10)")
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

print ("Loading data...")
scafs_in = open(args.input, "rU")
for scaf in SeqIO.parse(scafs_in, "fasta") :

    scaf_len = len(scaf.seq)

    if scaf_len < 20:
        print ("Warning: Not outputting scaffold \"" + scaf.id + "\" as scaffold has length < 20, which is the minimum required by ENA")
        print (scaf.id + "\t" + str(scaf_len) + "\tScaffold too short", file=bad_scafs_out)
    else:
        # Get stem of scaffold id
        id_parts = str(scaf.id).split("_")
        stem = "_".join(id_parts[0:len(id_parts)-2])
        scaf_index = int(id_parts[4])

        # Split scaffolds by N's
        contig_seqs = re.split('([N|n]{10,})', str(scaf.seq))

        part_number = 1
        position = 1

        contig_index = 1

        ###############
        ##  Fix beginning and ends to make sure that there are no issues there.

        # The split might result in the first or last contig being empty if the last but one contig was made up of Ns
        if len(contig_seqs[0]) == 0:
            contig_seqs.pop(0)

        if len(contig_seqs[-1]) == 0:
            contig_seqs.pop()

        if len(contig_seqs[0]) < 20 and not checkIfGap(contig_seqs[0]):
            if len(contig_seqs) > 2:
                # remove first short contig and next gap
                gap_len = len(contig_seqs[0]) + len(contig_seqs[1])
                contig_seqs.pop(0)
                contig_seqs.pop(0)
                print ("WARNING: " + scaf.id + " starts with a short contig.  Removing first contig and gap.  Length: " + str(gap_len))
                print (scaf.id + "\t" + str(scaf_len) + "\tstarts with a short contig\t" + str(gap_len), file=bad_scafs_out)
            else:
                print ("ERROR: Found a scaffold that's too small for ENA!")
                exit(1)

        if len(contig_seqs[-1]) < 20 and not checkIfGap(contig_seqs[-1]):
            if len(contig_seqs) > 2:
                gap_len = len(contig_seqs[-1]) + len(contig_seqs[-2])
                contig_seqs.pop()
                contig_seqs.pop()
                print ("WARNING: " + scaf.id + " ends with a short contig.  Removing last contig and adjacent gap.  Length: " + str(gap_len))
                print (scaf.id + "\t" + str(scaf_len) + "\tends with a short contig\t" + str(gap_len), file=bad_scafs_out)
            else:
                print ("ERROR: Found a scaffold that's too small for ENA!")
                exit(1)

        if checkIfGap(contig_seqs[0]):
            if len(contig_seqs) > 1:
                gap_len = len(contig_seqs[0])
                contig_seqs.pop(0)
                print ("WARNING: " + scaf.id + " starts with a gap.  Removing gap. Length: " + str(gap_len))
                print (scaf.id + "\t" + str(scaf_len) + "\tstarts with a gap\t" + str(gap_len), file=bad_scafs_out)
            else:
                print ("ERROR: Found a scaffold that just contains Ns!")
                exit(1)

        if checkIfGap(contig_seqs[-1]):
            if len(contig_seqs) > 1:
                gap_len = len(contig_seqs[-1])
                contig_seqs.pop()
                print ("WARNING: " + scaf.id + " end with a gap.  Removing gap. Length: " + str(gap_len))
                print (scaf.id + "\t" + str(scaf_len) + "\tends with a short contig\t" + str(gap_len), file=bad_scafs_out)
            else:
                print ("ERROR: Found a scaffold that just contains Ns!")
                exit(1)


        #Trim first and last contigs if they starts with an N
        if contig_seqs[0][0] == "N":
            contig_seqs[0], gap_len = trimNsFromStart(contig_seqs[0])
            print ("WARNING: " + scaf.id + " starts with a short run of Ns.  Trimming. Length: " + str(gap_len))
            print (scaf.id + "\t" + str(scaf_len) + "\tstarts with a short run of Ns\t" + str(gap_len), file=bad_scafs_out)

        if contig_seqs[-1][-1] == "N":
            contig_seqs[-1], gap_len = trimNsFromEnd(contig_seqs[-1])
            print ("WARNING: " + scaf.id + " ends with a short run of Ns.  Trimming. Length: " + str(gap_len))
            print (scaf.id + "\t" + str(scaf_len) + "\tends with a short run of Ns\t" + str(gap_len), file=bad_scafs_out)



        # Iterate over the cleaned contig_seqs (although we might still have short contigs in the middle
        i = 1
        for contig in contig_seqs:

            contig_len = len(contig)

            isgap = checkIfGap(contig) # Check if this contig is a gap

            isGood = False

            # Initialise AGP object for this instance
            agp = AGP2()
            agp.object = scaf.id
            agp.object_beg = position
            agp.object_end = position + contig_len -1
            agp.part_number = part_number
            agp.object_id = scaf_index

            if isgap:
                # component_type gap_length gap_type linkage
                agp.comp_type = "N"
                agp.gap_length = contig_len
                agp.gap_type = "scaffold"
                agp.linkage = "yes"
                agp.linkage_evidence = "paired-ends"
            else:

                contig_id = stem + "_contig_" + str(actual_contig_index)

                agp.comp_type = "W"
                agp.comp_id = contig_id
                agp.comp_beg = 1            # Can't be broken down further so always start at 1
                agp.comp_end = contig_len
                agp.orientation = "+"   # Always on the positive strand

            #  Check to make sure there are no short contigs in the middle of a scaffold
            if (isgap and part_number > 1 and agp_list[-1].isGap()) or (not isgap and contig_len < 20):

                # This case will end up converting this to a gap and merging this with the previous and later gaps
                # into one
                agp_list[-1].object_end = agp.object_end    # Adjust the end of the current gap to include this
                agp_list[-1].gap_length += contig_len
                print ("Warning: found scaffold with short contig inside... merging region")
                print (scaf.id + "\t" + str(scaf_len) + "\tScaffold contains short contig", file=bad_scafs_out)
            else:
                # Assume this is a good contig or gap
                isgood = True

            # Only output AGP if good
            if isgood:
                agp_list.append(agp)
                # Only output contig if this isn't a gap
                if not isgap:
                    record = SeqRecord.SeqRecord(Seq.Seq(contig), id=contig_id, description="")
                    contig_list.append(record)
                    contig_index += 1
                    actual_contig_index += 1
                part_number += 1
                position += contig_len

            i += 1


scafs_in.close()


print ("Sorting AGP...")
agp_list.sort()

print ("Sorting Fasta...")
# This is horrendous but will do for now (I have the memory!)
# Make map of id to record first
fasta_map = {}
contig_id_map = {}
id_list = list()
for contig in contig_list:
    #print (str(contig.id))
    fasta_map[str(contig.id)] = contig
    parts = str(contig.id).split("_")   # This is custom to my specific assembly
    contig_id_map[parts[4]] = contig.id
    id_list.append(int(parts[4]))

id_list.sort(key=int)

print ("Saving AGP...")
with open(agp_file, "w") as agp_out:
    for agp in agp_list:
        print(agp, file=agp_out)

print ("Saving contigs...")
with open(contig_file, "w") as contigs_out:
    for index in id_list:
        SeqIO.write(fasta_map[contig_id_map[str(index)]], contigs_out, "fasta") # Urgh!

print ("Saving scaffolds...")
with open(scaffold_file, "w") as scaffolds_out:
    current_id = agp_list[0].object
    scaffold = ""
    for agp in agp_list:
        if current_id == agp.object:
            if (agp.isGap()):
                scaffold += createGap(agp.gap_length)
            else:
                scaffold += str(fasta_map[agp.comp_id].seq)
        else:
            record = SeqRecord.SeqRecord(Seq.Seq(scaffold), id=current_id, description="")
            SeqIO.write(record, scaffolds_out, "fasta") # Urgh!
            # Reset scaffold sequence and id
            scaffold = str(fasta_map[agp.comp_id].seq)  # Must be a contig, not a gap
            current_id = agp.object

    # Output last scaffold
    record = SeqRecord.SeqRecord(Seq.Seq(scaffold), id=current_id, description="")
    SeqIO.write(record, scaffolds_out, "fasta") # Urgh!
