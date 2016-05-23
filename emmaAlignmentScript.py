#!/usr/bin/env python3
# Name: Eric Berg (edberg@ucsc.edu)
# Jurica Laboratory (mjurica@ucsc.edu)

import csv
import subprocess
import re
import os

'''
Utilizes the Emma program and ClustalW from the EBMOSS suite to 
generate text files that contain global alignment 
information.
'''

class alignmentGenerator:

  def __init__(self, compFile, index):

      self.compFile = compFile
      self.index = [int(i) for i in index.split(',')]   #Badass one-liner.

  def inputParse(self):
      '''Gets human and yeast protein names (Uniprot).'''

      tempFile = 'tempFile.txt'   
      numberOfProteins = 0      #Keeps
      with open(self.compFile, newline = '') as f:
        reader = csv.reader(f)
        for row in reader:
          humanProt = [row[x] for x in self.index]
          print(humanProt)
          with open(tempFile, 'w') as emmaInput:
            for x in humanProt:
              if len(x) > 0:
                 emmaInput.write('uniprot:' + x + '\n')
                 numberOfProteins = numberOfProteins + 1
              
          if numberOfProteins >= 3:  #Only evaluates 3 or more proteins.
             self.callEmma(tempFile) 
             numberOfProteins = 0   #resets counter.
          else: 
             numberOfProteins = 0   #resets counter.

      return

  
  def callEmma(self, tempFile):
      '''Calls and manipulates Emma program via command line'''
      
      nameList = []
      subprocess.call(['emma', '@' + tempFile,
                      'renameFile.txt',
                      'removeWhenDone.txt'])
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

  def cleanDir(self):
      '''gets rid of temporary files.'''
      
      os.remove('removeWhenDone.txt')
      os.remove('tempFile.txt.')
       
      return
     
      

compFile = input('Enter .csv file with protein names organized by column: ')
index = input('Enter column indices (comma-separated): ')
protCompare = alignmentGenerator(compFile, index)
protCompare.inputParse()
protCompare.cleanDir()


