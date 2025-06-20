# Code from https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/napari-stitching/napari-cosmx-stitching.html

import numpy
import pandas as pd
import os
from os import listdir
from glob import iglob
import sys
import subprocess

# Arguments from the command line
# Directory containing multiple slides
BatchFolder = sys.argv[1]
# Output path
OutputFolder = sys.argv[2]

### Functions

def check_folder(path: str) -> bool:
  """ Checks if a selected folder is a valid slide.

      Description: 
          A folder is a valid slide if it has the following folders:
          CellStatsDir, AnalysisResults/*, and RunSummary
        
      Returns:
          bool: True if valid, False if not valid
  """
  isValid = True
  if not os.path.isdir(path + '/CellStatsDir'):
    print("No valid CellStatsDir")
    isValid = False
  if not os.path.isdir(path + '/RunSummary'):
    print("No valid RunSummary")
    isValid = False 
  if not os.path.isdir(path + '/AnalysisResults'):
    print("No valid AnalysisResults Parent folder")
    isValid = False
  else: 
      # check if /AnalysisResults/<random_subfolder_name> exists
      analysis_sub_dir = [i for i in listdir(path + '/AnalysisResults') if not i.startswith('.')]
      if(len(analysis_sub_dir)!=1):
        print("No valid AnalysisResults subfolder")
        isValid = False
  return isValid

### Processing

for dir in listdir(BatchFolder):
  ParentFolder = os.path.join(BatchFolder, dir)
  if os.path.isfile(ParentFolder):
    print("Skipping file " + dir)
    continue
  else: 
    print('Processing Parent Folder: ' + dir)
    # Process slides within the Parent folder
    for SlideDir in listdir(ParentFolder):
        if os.path.isfile(os.path.join(ParentFolder, SlideDir)):
          print("Skipping file " + SlideDir)
          continue
        else:
          # check that it is a valid slide
          if not check_folder(os.path.join(ParentFolder, SlideDir)):
            print("Skipping folder " + SlideDir)
          else:
            SlideOutputDir = os.path.join(OutputFolder, dir)
            CellStatsDir = os.path.join(ParentFolder, SlideDir, "CellStatsDir")
            RunSummaryDir = os.path.join(ParentFolder, SlideDir, "RunSummary")
            AnalysisDirParent = os.path.join(ParentFolder, SlideDir, 'AnalysisResults')
            AnalysisDirSubBasename = [i for i in listdir(AnalysisDirParent) if not i.startswith('.')]
            AnalysisDir = os.path.join(AnalysisDirParent, AnalysisDirSubBasename[0])
            if os.path.exists(SlideOutputDir):
              print("Skipping output " + SlideOutputDir + ". Folder exists already")
              continue
            else:
              os.makedirs(SlideOutputDir)
              cmd_stitch = 'stitch-images -i "' + CellStatsDir + '" -f "' + RunSummaryDir + '" -o ' + SlideOutputDir
              print(cmd_stitch)
              subprocess.run(cmd_stitch, shell=True, check=True)
              cmd_read_targets = 'read-targets "' + AnalysisDir + '"' + ' -o ' + SlideOutputDir
              print(cmd_read_targets)
              subprocess.run(cmd_read_targets, shell=True, check=True)
              print("\n")