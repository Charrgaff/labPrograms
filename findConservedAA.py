#!/usr/bin/env python3
# Name: Eric Berg (edberg@ucsc.edu)
# Jurica Laboratory (mjurica@ucsc.edu)

'''Quick little program to find conserved animo acids given a pair of 
   sequences.
'''

import csv

class findConservedAA:

   def __init__(self, inputFile):
      self.inputFile = inputFile
  
   def parseFile(self):
      '''basic file parser using .csv module.'''

      with open(self.inputFile, 'r') as file:
         conservedList = []
         reader = csv.reader(file)
         for line in reader:
           if '-' not in line[1]:
             if '-' not in line[2]:
               try:
                 conservedList.append(str(line[2][1:]))
               except ValueError:
                 pass
      print(conservedList)
      conservedString = ','.join(conservedList)
      print(conservedString)
          
#############################################################
#Control flow.
#############################################################

inputFile = input("Enter name of alignment file (.csv): ")
AAFinder = findConservedAA(inputFile)
AAFinder.parseFile()
