#!/usr/bin/python3
import os
import sys

Usage = """
formats the GFF3 to final version - similar to NAM gff3 file
including the canonical transcript assignment tag
The file is acceptable for MaizeGDB submission

Usage:

  maizegdb_gff3_formatter.py <renamed_agat_formatted.gff3> <canonical_transcript_ids.txt> 

Input:

  renamed_agat_formatted.gff3 : The gff3 file sanitized using the agat_sp_manage_IDs.pl script
                                typically, you should first obtain gene id prefix from maizeGDB and run this as follows
                                
                                agat_sp_manage_IDs.pl --gff input.gff3 --prefix Ab00001aa --tair  --output prefinal.gff3


  canonical_transcript_ids.txt : list of transcript ids that are considered as primary transcript. You can run the TRaCE program to determine
                                the canonical transcript and create a list of mRNA ids (one per line). The number should be equal to the gene count in GFF3


eg command:

  maizegdb_gff3_formatter.py renamed_agat_formatted.gff3 canonical_transcript_ids.txt > maizeGDB_specifications.gff3 


Arun Seetharam
arnstrm@iastate.edu
09/15/2022
"""
if len(sys.argv) < 2:
    print(Usage)
    exit()
else:
   cmdargs = str(sys.argv)
   file1 = sys.argv[1]
   file2 = sys.argv[2]


def remap_gene(oldid, oldparent):
    spPrefix = oldid[0:9]
    geneNum = oldid[9:]
    newid = spPrefix + geneNum.zfill(6)
    newparent = spPrefix + geneNum.zfill(6)
    return(newid, newparent)


def remap_mrna(oldid, oldparent):
    spPrefix = oldid[0:9]
    tid = oldid.split(".")[1]
    geneNum = oldid.split(".")[0][9:]
    newid = spPrefix + geneNum.zfill(6) + "_T" + tid.zfill(3)
    newparent = spPrefix + geneNum.zfill(6)
    return(newid, newparent)


def remap_exon(oldid, oldparent):
    spPrefix = oldid[0:9]
    exonNum = oldid.split("-")[1].strip("exon")
    tid = oldid.replace("-", ".").split(".")[1]
    geneNum = oldid.split(".")[0][9:]
    newid = spPrefix + geneNum.zfill(6) + \
        "_T" + tid.zfill(3) + ".exon." + exonNum
    newparent = spPrefix + geneNum.zfill(6) + "_T" + tid.zfill(3)
    return(newid, newparent)


def remap_cds(oldid, oldparent):
    spPrefix = oldid[0:9]
    cdsNum = oldid.split("-")[1].strip("cds")
    tid = oldid.replace("-", ".").split(".")[1]
    geneNum = oldid.split(".")[0][9:]
    newid = spPrefix + geneNum.zfill(6) + \
        "_P" + tid.zfill(3)
    newparent = spPrefix + geneNum.zfill(6) + "_T" + tid.zfill(3)
    return(newid, newparent)


def remap_5utr(oldid, oldparent):
    spPrefix = oldid[0:9]
    utr = oldid.split("-")[1].strip("five_prime_utr")
    tid = oldid.replace("-", ".").split(".")[1]
    geneNum = oldid.split(".")[0][9:]
    newid = spPrefix + geneNum.zfill(6) + \
        "_T" + tid.zfill(3) + ".five_prime_utr." + utr
    newparent = spPrefix + geneNum.zfill(6) + "_T" + tid.zfill(3)
    return(newid, newparent)


def remap_3utr(oldid, oldparent):
    spPrefix = oldid[0:9]
    utr = oldid.split("-")[1].strip("three_prime_utr")
    tid = oldid.replace("-", ".").split(".")[1]
    geneNum = oldid.split(".")[0][9:]
    newid = spPrefix + geneNum.zfill(6) + \
        "_T" + tid.zfill(3) + ".three_prime_utr." + utr
    newparent = spPrefix + geneNum.zfill(6) + "_T" + tid.zfill(3)
    return(newid, newparent)


# read GFF file
my_gff = []
my_comments = []
with open(file1, 'r') as gff_file:
    for line in gff_file:
        if not line.startswith("#"):
            features = line.rstrip("\n").split("\t")
            my_gff.append(features)
with open(file1, 'r') as gff_file:
    for line in gff_file:
        if line.startswith("##"):
            comments = line.rstrip("\n")
            my_comments.append(comments)
# read primary ids
txt_file = open(file2, "r")
file_content = txt_file.read()
my_primary_ids = file_content.split("\n")
# process gff3 lines
for comments in my_comments:
    print(comments)
for features in my_gff:
    seqid, source, type, start, end, score, strand, phase, attributes = features
    d = dict(x.split("=") for x in attributes.split(";"))
    
    if type == "gene":
        d["ID"], d["Name"] = remap_gene(d["ID"], d["Name"])
        ID = 'ID='+d["ID"]
        req_attr = [ID, 'biotype=protein_coding', 'logic_name=BIND_gene']
    if type == "mRNA":
        d["ID"], d["Parent"] = remap_mrna(d["ID"], d["Parent"])
        ID = 'ID='+d["ID"]
        Parent = 'Parent='+d["Parent"]
        biotype = 'biotype=protein_coding'
        transcript_id = 'transcript_id='+d["ID"]
        if d["ID"] in my_primary_ids:
            primary = 'canonical_transcript=1'
            req_attr = [ID, Parent, biotype, transcript_id, primary]
        else:
            req_attr = [ID, Parent, biotype, transcript_id]
    if type == "exon":
        d["ID"], d["Parent"] = remap_exon(d["ID"], d["Parent"])
        ID = 'ID='+d["ID"]
        Parent = 'Parent='+d["Parent"]
        Name = 'ID='+d["ID"]
        exonID = 'ID='+d["ID"]
        strsplit = exonID.split('.')
        rank = 'rank='+strsplit[-1]
        req_attr = [Parent, Name, exonID, rank]
    if type == "CDS":
        d["ID"], d["Parent"] = remap_cds(d["ID"], d["Parent"])
        ID = 'ID='+d["ID"]
        Parent = 'Parent='+d["Parent"]
        protein_id = 'protein_id='+d["ID"]
        req_attr = [ID, Parent, protein_id]
    if type == "five_prime_UTR":
        d["ID"], d["Parent"] = remap_5utr(d["ID"], d["Parent"])
        ID = 'ID='+d["ID"]
        Parent = 'Parent='+d["Parent"]
        req_attr = [ID, Parent]
    if type == "three_prime_UTR":
        d["ID"], d["Parent"] = remap_3utr(d["ID"], d["Parent"])
        ID = 'ID='+d["ID"]
        Parent = 'Parent='+d["Parent"]
        req_attr = [ID, Parent]
    attributes = ";".join(req_attr)
    source = "PanAnd"
    print(seqid, source, type, start, end, score,
          strand, phase, attributes, sep='\t')


