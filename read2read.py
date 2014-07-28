#!/usr/bin/env python

#Run the pipeline for creating graph of NGS data.

import os, glob, commands
import string
import sys, optparse

if __name__ == "__main__":

        scriptDir = os.path.dirname(os.path.realpath(__file__))
        print "Script directory is " + scriptDir

        parser = optparse.OptionParser("usage: %prog [options] fastaFile outputDir")
        parser.add_option("-p", "--percentage", dest="similarity", default=85, type="int", help="specify sequence similarity percentage. Default=85.")
        parser.add_option("-W", "--wordsize", dest="wordsize", default=18, type="int", help="specify megablast word size. Default=18.")
        parser.add_option("-a", "--parallel", dest="numthreads", default=1, type="int", help="specify megablast number of threads. Default=1.")
        parser.add_option("-l", "--coverage", dest="coverage", default=55, type="int", help="specify percentage of length coverage for creating the graph. Default=55.")
        parser.add_option("-c", "--contigassembly", action="store_true", dest="contig", default=False, help="create contings using CAP3 and create BioLayout class files. Default=False.")
        parser.add_option("-k", "--keepoutput", action="store_true", dest="keepoutput", default=False, help="keep intermediate data files. Default=False.")
        parser.add_option("-q", "--verbose", action="store_false", dest="verbose", default=True, help="quit verbose. Default=True.")
        (options, args) = parser.parse_args()

        if (len(args) != 2):
                parser.error("incorrect number of arguments")

        fastaFile = args[0]
        outputDir = args[1]
        percentage = options.similarity
        wordsize = options.wordsize
        numthreads = options.numthreads
        coverage = options.coverage
        contig = options.contig
        keepoutput = options.keepoutput
        verbose=options.verbose

        if not os.path.isdir(outputDir):
                os.makedirs(outputDir)
        else:
                print "Output directory already exists!"
                sys.exit(1)

        os.chdir(outputDir)
        nameDB = os.path.basename(string.split(fastaFile,".")[0])
        cmd = scriptDir + "/makeblastdb -dbtype nucl -out " + nameDB + " -input_type fasta -in " + fastaFile + " -title " + nameDB + " -max_file_sz 2GB"
        if verbose:
                print "Creating database..."
                print cmd
        execCmd = os.popen(cmd, "r")
        outCmd = execCmd.read()
        exitCode = execCmd.close()

        if exitCode != None:
          print outCmd
          print "...failed"
          exit(exitCode)

        megablast_txtFile = outputDir + "/" + nameDB + "_megablast.txt"
        pairwise_txtFile = outputDir + "/" + nameDB + "_pairwise.txt"

        cmd = scriptDir + "/blastn -db " + nameDB + " -query " + fastaFile + " " \
            "-perc_identity " + str(percentage) + " -word_size " + str(wordsize) + " " \
            "-xdrop_gap 40 -num_alignments 90000000 -dust yes -soft_masking true -lcase_masking " + \
            " -outfmt \"7 qacc sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\" " + \
            "-num_threads " + str(numthreads) + " -out " + megablast_txtFile
        if verbose:
                print "Making alignment..."
                print cmd
        execCmd = os.popen(cmd, "r")
        outCmd = execCmd.read()
        exitCode = execCmd.close()

        if exitCode != None:
          print outCmd
          print "...failed"
          exit(exitCode)

        cmd = "python " + scriptDir + "/megablast2ncol.py " + fastaFile + " " + megablast_txtFile + " " + \
            pairwise_txtFile + " " + str(percentage/float(100)) + " " + str(coverage/float(100))
        if verbose:
                print "Creating graph..."
                print cmd

        execCmd = os.popen(cmd, "r")
        outCmd = execCmd.read()
        execCmd.close()

        if not keepoutput:
                print "Removing " + megablast_txtFile + "..."
                os.remove(megablast_txtFile)
