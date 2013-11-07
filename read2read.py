#!/usr/bin/env python

#Run the pipeline for creating graph of NGS data.

import os, glob, commands
import string
import sys, optparse

if __name__ == "__main__":

        scriptDir = os.path.dirname(os.path.realpath(__file__))
        print "Script directory is " + scriptDir

        parser = optparse.OptionParser("usage: %prog [options] fastaFile outputDir")
        parser.add_option("-P", "--protein", dest="protein",default="F", help="specify formatdb input, protein or nucleotide (F for nucleotide). Default=F")
        parser.add_option("-p", "--percentage", dest="similarity", default=85, type="int", help="specify sequence similarity percentage. Default=85.")
        parser.add_option("-W", "--wordsize", dest="wordsize", default=18, type="int", help="specify magablast word size. Default=18.")
        parser.add_option("-a", "--parallel", dest="numthreads", default=1, type="int", help="specify magablast number of threads. Default=1.")
        parser.add_option("-l", "--coverage", dest="coverage", default=55, type="int", help="specify percentage of length coverage for creating the graph. Default=55.")
        parser.add_option("-c", "--contigassembly", action="store_true", dest="contig", default=False, help="create contings using CAP3 and create BioLayout class files. Default=False.")
        parser.add_option("-q", "--verbose", action="store_false", dest="verbose", default=True, help="quit verbose. Default=True.")
        (options, args) = parser.parse_args()

        if (len(args) != 2):
                parser.error("incorrect number of arguments")

        fastaFile = args[0]
        outputDir = args[1]
        protein = options.protein
        percentage = options.similarity
        wordsize = options.wordsize
        numthreads = options.numthreads
        coverage = options.coverage
        contig = options.contig
        verbose=options.verbose

        if not os.path.isdir(outputDir):
                os.makedirs(outputDir)
        else:
                print "Output directory already exists!!"
                sys.exit(1)

        os.chdir(outputDir)
        nameDB = os.path.basename(string.split(fastaFile,".")[0])
        cmd = scriptDir + "/formatdb -t " + nameDB + " -i " + fastaFile + " -p " +  protein + " -n " + nameDB
        if verbose:
                print "Creating database..."
                print cmd
        execCmd = os.popen(cmd, "r")
        outCmd = execCmd.read()
        execCmd.close()

        cmd = scriptDir + "/megablast -d " + nameDB + " -i " + fastaFile + " -p " + str(percentage) + " -W " + str(wordsize) + " -UT -X40 -JF -F mD -v90000000 -b90000000 -D3 -a " +  str(numthreads) + " -o " + outputDir + "/" + nameDB + "_megablast.txt"
        if verbose:
                print "Making alignment..."
                print cmd
        execCmd = os.popen(cmd, "r")
        outCmd = execCmd.read()
        execCmd.close()

        cmd = "python " + scriptDir + "/megablast2ncol.py " + fastaFile + " " + outputDir + "/" + nameDB + "_megablast.txt" + " " + outputDir + "/" + nameDB +"_pairwise.txt" + " " + str(percentage/float(100)) + " " + str(coverage/float(100))
        if verbose:
                print "Creating graph..."
                print cmd

        execCmd = os.popen(cmd, "r")
        outCmd = execCmd.read()
        execCmd.close()

