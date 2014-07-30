#!/usr/bin/env python

#Run the pipeline for creating graph of NGS data.

import os, glob, commands
import string
import sys, optparse
import math
import shlex
import subprocess

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

  fastaFile = os.path.abspath(args[0])
  outputDir = os.path.abspath(args[1])
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
  nameDB = os.path.splitext(os.path.basename(fastaFile))[0]
  cmd = scriptDir + "/makeblastdb -dbtype nucl -out " + nameDB + " -input_type fasta -in " + fastaFile + " -title " + nameDB + " -max_file_sz 2GB"
  if verbose:
    print cmd
  args = shlex.split(cmd)
  p = subprocess.Popen(args, stdout=subprocess.PIPE)

  output = p.communicate()[0]
  if p.returncode != 0:
    print output
    print "...failed"
    exit(p.returncode)

  numLines = sum(1 for line in open(fastaFile))
  linesPerThread = -(-numLines // numthreads)

  # We need an even number of lines because fasta format is pairs of lines
  if linesPerThread % 2 != 0:
    linesPerThread += 1

  remainingLines = 0
  thread = 0
  fastaThreadFile = None
  for line in open(fastaFile):
    fastaThreadFilename = outputDir + "/" + nameDB + ".t" + `int(thread)` + ".fasta"
    if remainingLines == 0:
      if fastaThreadFile != None:
        fastaThreadFile.close()

      fastaThreadFile = open(fastaThreadFilename, 'a+')
      remainingLines = linesPerThread
      thread += 1

    fastaThreadFile.write(line)
    remainingLines -= 1

  if fastaThreadFile != None:
    fastaThreadFile.close()

  # May not be able to divide into requested number of threads
  numthreads = thread

  megablastProcesses = []
  for thread in range(0, numthreads):
    fastaThreadFilename = outputDir + "/" + nameDB + ".t" + `int(thread)` + ".fasta"
    megablastThreadFilename = outputDir + "/" + nameDB + ".t" + `int(thread)` + "_megablast.txt"

    cmd = scriptDir + "/blastn -db " + nameDB + " -query " + fastaThreadFilename + " " \
        "-perc_identity " + str(percentage) + " -word_size " + str(wordsize) + " " \
        "-xdrop_gap 40 -num_alignments 90000000 -dust yes -soft_masking true -lcase_masking " + \
        "-outfmt \"7 qacc sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\" " + \
        "-out " + megablastThreadFilename
    if verbose:
      print cmd

    args = shlex.split(cmd)
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    megablastProcesses.append({'process' : p, 'blastFilename' : megablastThreadFilename, 'fastaFilename' : fastaThreadFilename})

  megablastFilename = outputDir + "/" + nameDB + "_megablast.txt"
  with open(megablastFilename, 'w') as megablastFile:
    for megablastProcess in megablastProcesses:
      p = megablastProcess['process']
      megablastThreadFilename = megablastProcess['blastFilename']
      fastaThreadFilename = megablastProcess['fastaFilename']

      output = p.communicate()[0]
      if p.returncode != 0:
        print output
        print "...failed"
        exit(p.returncode)

      with open(megablastThreadFilename) as megablastThreadFile:
        for line in megablastThreadFile:
          megablastFile.write(line)
      if not keepoutput:
        os.remove(megablastThreadFilename)
        os.remove(fastaThreadFilename)

  pairwise_txtFile = outputDir + "/" + nameDB + "_pairwise.txt"
  cmd = "python " + scriptDir + "/megablast2ncol.py " + fastaFile + " " + megablastFilename + " " + \
      pairwise_txtFile + " " + str(percentage/float(100)) + " " + str(coverage/float(100))
  if verbose:
    print cmd

  execCmd = os.popen(cmd, "r")
  outCmd = execCmd.read()
  execCmd.close()

  if not keepoutput:
    os.remove(megablastFilename)
