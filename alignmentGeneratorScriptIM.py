#!/usr/bin/env python3
# Name: Eric Berg (edberg@ucsc.edu)
# Jurica Laboratory (mjurica@ucsc.edu)

import csv
import subprocess
import re, locale
import os

'''
Utilizes the Stretcher program from the EBMOSS suite to 
generate text files that contain global alignment 
information.
'''

class alignmentGenerator:

   def __init__(self, compFile, index):

      self.compFile = compFile
      self.index = [int(i) for i in index.split(',')]   #Badass one-liner.

   def inputParse(self):
      '''Gets human and yeast protein names (Uniprot).'''

      with open(self.compFile, newline = '') as f:
        reader = csv.reader(f)
        for row in reader:
          protNames = [row[x] for x in self.index]
          self.callStretcher(protNames)
   
   def callStretcher(self, protNames):
      '''Calls and manipulates Stretcher program via command line'''
      
      nameList = []
      subprocess.call(["stretcher", "-aformat3",
                       "markx3", "uniprot:" + protNames[0],
                       "uniprot:" + protNames[1],
                       "renameFile.txt"])
      try:                                  #checks if alignment was valid
         f = open('renameFile.txt', 'r')
         f.close()
      except FileNotFoundError:
         return

      with open('renameFile.txt', 'r') as rename:  #rename output
         for line in rename:
            if ">" in line:
               line.replace('\n','')
               lineClean = re.sub(r'[^a-zA-Z0-9\_]','',line)
               nameList.append(lineClean)
      nameString = ','.join(nameList) + '.txt'
      os.rename('renameFile.txt', nameString)
                       
      return      
      

compFile = input('Enter .csv file with protein names organized by column: ')
index = input('Enter column indices: ')
protCompare = alignmentGenerator(compFile, index)
protCompare.inputParse()
 

