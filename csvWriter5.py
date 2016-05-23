#!/usr/bin/env python3
# Name: Eric Berg (edberg@ucsc.edu)
# Jurica Laboratory (mjurica@ucsc.edu)

'''
Reformats protein alignment data into a spreadsheet-readable format.
Reads either a single text file or a directory of text files.
'''

#############################################################
#Pseudocode
# 1. Asks for either a directory or single text file (.txt).
# 2. If single text file is read, program directly parses it. If a
#    directory is read, the program lists the directory files and finds/
#    parses the specified text files individually.
# 3. ">" identifier is used to determine where sequences start. The 
#    protein names are then formatted and placed into arrays.
# 4. "#" identifier is used to determine when last sequence has been 
#    reached.
# 5. Protein name and sequence data is stored in multidimensional lists.
# 6. Numpy arrays are used to specify specific indices
#    when generating .csv files.
#
# Note: All data from single text files is placed in array before
#       generation of .csv file. This may cause a stack overflow
#       when attempting to reformat very large alignments.
# Acceptable formats:
#    1. markx3.txt (emboss/stretcher)
#    2. FASTA-type (ex. >Q15393)
############################################################ 

'''
Example usage:
   >>>python3 csvWriter3.py
   Enter directory or alignment file (.txt): textFile.txt

Example usage2:
   >>>python3 csvWriter3.py
   Enter directory or alignment file (.txt): /home/alignmentFiles
'''

###############################################################
#Beginning of program.
###############################################################
    

import csv
import numpy
import re, locale
import os

class textParser:
   '''Creates spreadsheet-readable .csv file from protein alignment data.
   
   Args:
      textfile(str): Text file name.
      separateResidues(bool): Determines whether user wants specific index.
      outPutLoc(str): Location of .csv output
   Attributes:
      protFileName(str): Text file name.
      protNames(list): List of all protein names for a single text file.
      initProtSeq(list): Temporary list that contains unformatted sequences.
      allSeqs(list): Final list of all formatted sequences from single file.
   '''  

   protNames = []
   protFileName = []
   initProtSeq = []
   allSeqs = []

   def __init__(self, textfile, outPutLoc, separateResidues):
      self.textfile = textfile
      self.separateResidues = separateResidues
      self.outPutLoc = outPutLoc
      
   def findProtein(self):
      '''Determines location of protein names and sequences.
    
      This method goes through all lines in a single text file and sends
      protein information to protCleaner and listCleaner for formatting.

      Note:
         ">" indicates start of new sequence.
         "#" indicates end of final sequence.        
      Args:
         None.
      Returns:
         None.
      '''
     
      with open(self.textfile, 'r') as alignmentOutput: 
         foundProt = 0       
         for line in alignmentOutput:  
            if ">" in line:          
               self.listCleaner()       #Formats sequence.
               self.protCleaner(line)   #Formats protein name.
               foundProt = 1
            elif foundProt == 1:        #Checks if at end of final sequence.
               if "#" in line:
                  self.listCleaner()    #Formats final seqence.
                  break
               else:
                  textParser.initProtSeq.extend(line) #Adds to initProtSeq.
     
      self.listCleaner()    #Formats final sequence.
             
      return

   def protCleaner(self, line):
      '''Formats protein name data from line and appends it to protName list.

      Once a line is formatted, it is appended to the protNames list.
      Control is then returned to the findProtein method.

      Args:
         line(str): A single line from a text file.
      Returns:
         None.
      '''
  
      line.replace('\n','')    
      cleanProtName = re.sub(r'[^a-zA-Z0-9\_]','',line)
      textParser.protFileName.append(cleanProtName)
      textParser.protNames.append(cleanProtName)
      if separateResidues == True:
         textParser.protNames.append("..........")
      return      
   
   def listCleaner(self):
      '''Formats sequence data from initProtSeq and appends it to allSeqs. 
       
      Once initProtSeq is formatted, it is appended to the allSeqs list.
      Control is then returned to the findProtein method.
      '''

      if self.initProtSeq == []:  #Checks if there is data to parse.
         return
      else:
         cleanSeq = []    #Temporary list.      
         for x in range(0, len(textParser.initProtSeq)):
            if (textParser.initProtSeq[x]) != '\n':
               cleanSeq.append(textParser.initProtSeq[x])
            else:
               continue
         ticker = 1       #Keeps track of protein index.
         for x in range(0, len(cleanSeq)):
            if cleanSeq[x] != "-":
               stringTicker = str(ticker)  #Converts to string.
               cleanSeq[x] = cleanSeq[x] + stringTicker
               ticker = ticker + 1
            else:
               continue
            
         textParser.allSeqs.append(cleanSeq)  #Appends formatted list.
         textParser.initProtSeq = []          #Clears initProtSeq.
         
      return
   
   def writeCSV(self):
      '''Uses protNames and allSeqs lists to generate .csv file.

      Numpy multidimensional arrays are used to keep track of residue
      type and location in individual proteins, as well as their relation
      to other proteins.
      '''

      outputFileName = self.outPutLoc+'/'+','.join(textParser.protFileName)+".csv"
      textParser.protNames = ["Index"] + textParser.protNames
      
      numpyAllSeqs = numpy.array(textParser.allSeqs) #Creates numpy array.
      alignmentLength = int(len(numpyAllSeqs[0]))  #Longest len b/n all read.
      
      with open(outputFileName, 'a', newline='') as f: #Generation of output.
         writer = csv.writer(f)
         writer.writerow(textParser.protNames)
         for x in range(0, alignmentLength):
           currentRow = numpyAllSeqs[:, x]
           indexedRow = currentRow.tolist()
           if separateResidues == True:
              outPutRow = []
              for x in range(0, len(indexedRow)):
                 outPutRow.append(indexedRow[x][0])
                 outPutRow.append(indexedRow[x][1:])
              outPutRow = [x + 1] + outPutRow
              writer.writerow(outPutRow)
           elif separateResidues == False:
              indexedRow = [x + 1] + indexedRow
              writer.writerow(indexedRow)
           else:
              print("Unable to set index for row.")
              return
              
      textParser.protNames = []    #Clears lists
      textParser.initProtSeq = []
      textParser.allSeqs = []
      textParser.protFileName = []
      return
      

def inputFileParse(textDataString, outPutLoc, separateResidues):
   '''Parses input as either a single file or a directory of list files.

   If a directory is given, this method reads and formats files one at a time

   Args:
      textDataString(str): User input (Either text file or directory).
   Returns:
      None.'''

   if ".txt" in textDataString:  #text file case.
      parseText = textParser(textDataString, outPutLoc, separateResidues) 
      parseText.findProtein()
      parseText.writeCSV()  
   else:                         #Directory case.
      for fileName in os.listdir(textDataString):
         direcFile = textDataString + '/' + fileName
         print(direcFile)
         parseText = textParser(direcFile, outPutLoc, separateResidues) 
         parseText.findProtein()
         parseText.writeCSV()  
   return 
       
####################################################################
#End of program.
####################################################################

separateResidues = False
textDataString = input("Enter directory path or alignment file (.txt): ")
outPutLoc = input("Enter output directory: ")
choice = input("Place residue type and number in separate columns (y/n)? ")
if choice == 'y':
   separateResidues = True
elif choice == 'n':
   separateResidues = False
else:
   print("Unable to read choice, resorting to default.")
inputFileParse(textDataString, outPutLoc, separateResidues)





