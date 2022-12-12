#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 13:07:14 2021

@author: Avery Davis Bell
Generates an EMASE-format genes-to-transcripts mapping file from a GTF (first column = gene ID; subsequent columns = transcript IDs; tab-delimited)
"""
import argparse
#%%
def readtranscripts(strGTF):
    """
    Reads all transcript lines of GTF, sorts them into dictionary keyed by gene id

    Parameters
    ----------
    strGTF : str
        Path to GTF to process.

    Returns
    -------
    hGns : dictionary
        Keys are gene_ids, values = list of transcript_ids for that gene

    """
    hGns = {}
    istm = open(strGTF, "r")
    
    for strLine in istm:
        if strLine[0] == "#":
            continue
        astrLine = strLine.strip().split("\t")
        if astrLine[2] == "transcript":
            astrInfo = astrLine[8].split("; ")
            hInfo = {}
            for strInfo in astrInfo:
                astrOne = strInfo.split(" ")
                hInfo[astrOne[0]] = astrOne[1]
            strGeneID = hInfo["gene_id"].strip('"')
            strTrID = hInfo["transcript_id"].strip('"')
            ############################print([strGeneID, strTrID])
            if strGeneID in hGns.keys():
                hGns[strGeneID].append(strTrID)
            else:
                hGns[strGeneID] = [strTrID]
    
    return(hGns)
#%%

def writeemasemapping(hGns, strOut):
    """
    Writes genes-to-transcripts mapping in EMASE format: first column is gene ID, subsequent = transcript ID(s); all tab-delimited
    Sorted by gene ID!

    Parameters
    ----------
    hGns : dictionary
        Keys are gene_ids, values = list of transcript_ids for that gene
    strOut : string
        output filepath to write to

    Returns
    -------
    None.

    """
    ostm = open(strOut, "w")
    for strGene in sorted(hGns.keys()):
        astrOut = [strGene] + sorted(hGns[strGene])
        ostm.write("\t".join(astrOut) + "\n")
    ostm.close()
    
#%%
def main():
    # Parse arguments
    argp = argparse.ArgumentParser(prog = "gtftoemasegenemappings.py",
                                   description = "Generates an EMASE-format genes-to-transcripts mapping file from a GTF (first column = gene ID; subsequent columns = transcript IDs; tab-delimited)")
    argp.add_argument("-gtf", dest = "strGTF", required = True, type = str,
                      metavar = "GTF file",
                      help = "Path to (uncompressed) GTF file")
    argp.add_argument("-o", dest = "strOut", required = True, type = str,
                      metavar = "Output file",
                      help = "Output filepath for writing EMASE-format genes-to-transcripts mapping file (tsv)")

    args = argp.parse_args()
    
    # Run
    hGns = readtranscripts(args.strGTF)
    writeemasemapping(hGns, args.strOut)


if __name__=="__main__":
	main()
     
