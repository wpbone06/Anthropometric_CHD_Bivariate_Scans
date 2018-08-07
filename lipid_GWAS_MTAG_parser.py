#!/usr/bin/python2.7

import sys
import getopt
import pdb
import os
import re
import pandas as pd


def main():
    argsList = sys.argv[1:] ## grab arguments from command line
    helpstatement = "Usage -i lipid_GWAS file -o output (MTAG ready) file. Generates a MTAG compatible file given a lipid GWAS file"

    if len(argsList) != 4:
        print "\nYou have supplied too many or too few arguments!"
        print helpstatement
        exit(1)

    # Parsing arguments from command line
    else:
        try:
            opts, args = getopt.getopt(argsList,"i:o:")

        except:
            print("Usage: -i lipid_GWAS file  -o output file" % sys.argv[0])
            print helpstatement
            exit(1)

        for (opt, arg) in opts:

            if opt == '-i':
                lipidGWAS_str = arg
                print "\nlipid_GWAS file found %s \n" % lipidGWAS_str


            elif opt == '-o':
                output_str = arg
                #output_file = open(output_str, "w")
                print "The output file name will be %s \n" % output_str

            else:
                print "\n%s is not an accepted option for this script" % opt
                print helpstatement
                exit(1)

    lipid_DF = pd.read_table(lipidGWAS_str ,sep="\t", low_memory=False)

    #split the "SNP_hg19 column on colons
    lipid_DF["chr"], lipid_DF["bpos"] = zip(*lipid_DF["SNP_hg19"].map(lambda x: x.split(":")))

    lipid_DF = lipid_DF.drop(columns=["SNP_hg18","SNP_hg19"])

    # change the required MTAG column headers to the MTAG standard column name
    lipid_DF = lipid_DF.rename(columns={"rsid":"snpid","A1":"a1","A2":"a2","N":"n","P-value":"pval","Freq.A1.1000G.EUR":"freq"})

    #remove the NaNs from the freq column
    lipid_DF = lipid_DF[lipid_DF.freq.notnull()]

    # remove the rows with a freq of 0 or 1
    badMAFsList = [0.0,1.0]
    lipid_DF = lipid_DF[~lipid_DF.freq.isin(badMAFsList)]

    #randomly select half of the  entries to flip the betas and alleles
    lipid_DF_NegativeBsDF = lipid_DF.sample(frac=0.5, random_state=10)

    # multiply the betas by -1
    lipid_DF_NegativeBsDF.beta *=-1

    #rename subsets allele columns switching the column names a1 now a2 and vice versa
    lipid_DF_NegativeBsDF = lipid_DF_NegativeBsDF.rename(columns={"a1":"a2","a2":"a1"})


    #flip the frequencies to match the alleles that were flipped
    lipid_DF_NegativeBsDF["flipFreq"] = 1 - lipid_DF_NegativeBsDF["freq"]

    lipid_DF_NegativeBsDF = lipid_DF_NegativeBsDF.drop(columns=["freq"])

    lipid_DF_NegativeBsDF = lipid_DF_NegativeBsDF.rename(columns={"flipFreq":"freq"})

    #make a list of the rownames/ index numbers that were flipped to negative betas
    NegativeRowNames = list(lipid_DF_NegativeBsDF.index)

    # make a data frame of all the rows that did not have their betas negatived
    lipid_DF_PositiveBsDF = lipid_DF[~lipid_DF.index.isin(NegativeRowNames)]

    lipid_DF_PosNegBs = lipid_DF_PositiveBsDF.append(lipid_DF_NegativeBsDF)

    lipid_DF_PosNegBs = lipid_DF_PosNegBs.sort_index(axis=0,ascending=True)

    lipid_DF_PosNegBs["z"] = lipid_DF_PosNegBs["beta"] / lipid_DF_PosNegBs["se"]


    #remove the beta and se columns since not used by MTAG
    lipid_DF_PosNegBs = lipid_DF_PosNegBs.drop(columns=["beta","se"])

    # sort the columns to be the same as those in MTAG Tutorial
    lipid_DF_PosNegBs =  lipid_DF_PosNegBs[["snpid", "chr", "bpos", "a1", "a2", "freq", "z", "pval", "n"]]

    lipid_DF_PosNegBs.to_csv(output_str,sep="\t",index=False,header=True)

if __name__ == "__main__":
    main()
