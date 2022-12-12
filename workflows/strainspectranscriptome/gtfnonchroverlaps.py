#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 09:37:18 2021

@author: Avery Davis Bell
Splits a GTF file into records whose ends are before chromosome ends and those that aren't (for intermediate step in creating diploid transcriptome)
"""
import argparse

#%%
def readchrlns(strChrlns):
    """
    Reads a tab-delimited chromosome length file into a dict

    Parameters
    ----------
    strChrlns : string
        Path to tab-delimited no-header file with first column having chromosome ID, second column having length of chromosome

    Returns
    -------
    hChrs: dict
        Keys are chromosome IDs, values are chromosome lengths (integer)

    """
    istm = open(strChrlns, "r")
    hChrs = {}
    
    for strLine in istm:
        strID, strLn = strLine.strip().split("\t")
        hChrs[strID] = int(strLn)
    istm.close()
    
    return(hChrs)
#%%

def splitgtf(hChrs, strGTF, strOutInclude, strOutExclude):
    """
    Splits a GTF file into records whose ends are before chromosome ends and those that aren't (for intermediate step in creating diploid transcriptome)

    Parameters
    ----------
    hChrs : dict
        Keys are chromosome IDs, values are chromosome lengths (integer)
    strGTF : string
        Path to full GTF to read. unzipped.
    strOutInclude : str
        Output filepath for writing GTF records that fit within chromosomes (majority!!)
    strOutExclude : str
        Output filepath for writing GTF records that DON'T fit within chromosomes
        (ends after chromosome ends). Could be none.

    Returns
    -------
    None.

    """
    # Open input, output files
    istm = open(strGTF, "r")
    ostmIn = open(strOutInclude, "w")
    ostmEx = open(strOutExclude, "w")
    
    # Line by line process GTF
    for strLine in istm:
        if strLine[0] == "#":
            continue # some GTFs start with some comments
        astrLine = strLine.strip().split("\t")
        strChr = astrLine[0]
        iEnd = int(astrLine[4])
        if hChrs[strChr] < iEnd:
            ostmEx.write("\t".join(astrLine) + "\n") # Write to exclude

        else:
            ostmIn.write("\t".join(astrLine) + "\n") # Write to include
       
    # Close all filestreams
    istm.close()
    ostmIn.close()
    ostmEx.close()
    
#%%
def main():
    # Parse arguments
    argp = argparse.ArgumentParser(prog = "gtfnonchroverlaps.py",
                                   description = "Splits a GTF file into records whose ends are before chromosome ends and those that aren't (for intermediate step in creating diploid transcriptome)")
    argp.add_argument("-gtf", dest = "strGTF", required = True, type = str,
                      metavar = "GTF file",
                      help = "Path to (uncompressed) GTF file")
    argp.add_argument("-lens", dest = "strChrlns", required = True, type = str, 
                      metavar = "input chromosome lengths file",
                      help = "Path to e Path to tab-delimited no-header file with first column having chromosome ID, second column having length of chromosome")
    argp.add_argument("-okeep", dest = "strOutInclude", required = True, type = str,
                      metavar = "Output path for included/kept records",
                      help = "Output filepath for writing GTF records that fit within chromosomes (majority!!)")
    argp.add_argument("-oexcl", dest = "strOutExclude", required = True, type = str,
                      metavar = "Output path for excluded records",
                      help = "Output filepath for writing GTF records that DON'T fit within chromosomes (could be none)")
    args = argp.parse_args()
    
    # Run
    hChrs = readchrlns(args.strChrlns)
    splitgtf(hChrs, args.strGTF, args.strOutInclude, args.strOutExclude)

if __name__=="__main__":
	main()
