# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************


"""
This protocol is used to perform a pocket search on a protein structure using the P2Rank software

"""
import os, gzip

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pyworkflow.object import String
from pwem.protocols import EMProtocol
import pwem.convert as emconv

from pwchem.objects import SetOfStructROIs, PredictStructROIsOutput, StructROI
from pwchem.utils import writePDBLine, splitPDBLine, runOpenBabel, cifFromASFile, writeCIFLine, gunzipFile, \
  getBaseName, performBatchThreading
from pwchem.constants import CIF_DEF_COLS, CIF_DEF_HEADER

from p2rank import Plugin


class P2RankFindPockets(EMProtocol):
    """
    Executes the p2rank software to look for protein pockets.
    """
    _label = 'Find pockets'
    _possibleOutputs = PredictStructROIsOutput
    stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputAtomStruct', params.PointerParam,
                       pointerClass='AtomStruct', allowsNull=False,
                       label="Input atom structure",
                       help='Select the atom structure to be fitted in the volume')
        form.addParallelSection(threads=4)

    def _getP2RankArgs(self):
      args = ['-f', os.path.abspath(self._getCifFile())]
      args += ['-o', os.path.abspath(self._getExtraPath())]
      args += ['-threads', self.getScipionThreads()]

      return args

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.P2RankStep)
        self._insertFunctionStep(self.createOutputStep)

    def convertInputStep(self):
        inpFile = self.inputAtomStruct.get().getFileName()
        cifFromASFile(inpFile, self._getCifFile())

    def P2RankStep(self):
        Plugin.runP2Rank(self, 'predict', args=self._getP2RankArgs(), cwd=self._getExtraPath())

    def createOutputStep(self):
        nt = self.numberOfThreads.get()
        inpStruct = self.inputAtomStruct.get()
        outASPath = os.path.relpath(self._getCifFile())
        pocketFiles = self._divideOutputPockets()

        outSet = SetOfStructROIs(filename=self._getExtraPath('StructROIs.sqlite'))
        outputPocks = performBatchThreading(self.performOutputCreation, pocketFiles, nt, cloneItem=False,
                                            inpStruct=inpStruct, propFile=self.getPropertiesFile(), asFile=outASPath)
        for i, pock in enumerate(outputPocks):
          outSet.append(pock)

        if len(outSet) > 0:
          outSet.buildPDBhetatmFile()
        self._defineOutputs(**{self._possibleOutputs.outputStructROIs.name: outSet})

    def performOutputCreation(self, pocketFiles, molLists, it, propFile, inpStruct, asFile):
      outPocks = []
      for pFile in pocketFiles:
        pock = StructROI(pFile, asFile, propFile, pClass='P2Rank')
        if len(pock.getPointsCoords()) > 2:  # minimum size for building pocket. cannot calculate volume otherwise
          pock.setVolume(pock.getPocketVolume())
          if str(type(inpStruct).__name__) == 'SchrodingerAtomStruct':
            pock._maeFile = String(inpStruct.getFileName())
          outPocks.append(pock)

      molLists[it] = outPocks

    # --------------------------- Utils functions --------------------
    def _getInputName(self):
        return os.path.splitext(os.path.basename(self.inputAtomStruct.get().getFileName()))[0]

    def _getCifFile(self):
        return os.path.abspath(self._getExtraPath(self._getInputName() + '.cif'))

    def getPdbInputStructName(self):
      return self._getCifFile().split('/')[-1]

    def getPropertiesFile(self):
        return self._getExtraPath(self.getPdbInputStructName()+'_predictions.csv')

    def _divideOutputPockets(self):
      '''Creates individiual pocket files'''
      gzfile = self._getExtraPath('visualizations/data/{}_points.pdb.gz'.format(
        self.getPdbInputStructName()))
      pocketDic = self.getPocketDic(gzfile)

      oDir = self._getExtraPath('pocketFiles')
      if not os.path.exists(oDir):
          os.mkdir(oDir)

      pFiles = []
      for pocketK in sorted(pocketDic):
          pFile = os.path.join(oDir, f'pocketFile_{pocketK}.cif')
          with open(pFile, 'w') as f:
              f.write(self.formatPocketStr(pocketDic[pocketK], pocketK))
          pFiles.append(pFile)
      return pFiles

    def formatPocketStr(self, pocketLines, pocketK):
      cifCols = '\n'.join(CIF_DEF_COLS)
      outStr = CIF_DEF_HEADER.format(cifCols)

      for i, pLine in enumerate(pocketLines):
          pLine = self.splitP2RankPDBLine(pLine)
          coords = [float(c) for c in pLine[6:9]]
          replacements = [str(i+1), f'C{i+1}', 'STP', 'C', 1, pocketK, *coords]
          cifLine = writeCIFLine(*replacements)
          outStr += cifLine
      return outStr

    def getPocketDic(self, pointsFile):
      dic={}
      with gzip.open(pointsFile) as f:
        for line in f:
          line = line.decode('utf-8')
          splittedLine = self.splitP2RankPDBLine(line)
          pocketId = int(splittedLine[5])
          if pocketId != 0:
            if pocketId in dic:
              dic[pocketId] += [line]
            else:
              dic[pocketId] = [line]
      return dic

    def _countNumberOfChains(self, inpFile):
      structureHandler = emconv.AtomicStructHandler()
      structureHandler.read(inpFile)
      structureHandler.getStructure()
      listOfChains, listOfResidues = structureHandler.getModelsChains()
      return len(listOfChains)

    def _countNumberOfAtoms(self, inpFile):
      with open(inpFile) as f:
        fileStr = f.read()
      return fileStr.count('ATOM')

    def splitP2RankPDBLine(self, line):
        '''Split lines taking into account the multiple exceptions found in P2Rank pdbs'''
        lenElem = len(line.split())
        if lenElem == 11:
            return line.split()
        else:
            lenLine = len(line.strip())
            #This happens when there are more than 9999 points (atom number collides with HETAM)
            if lenLine != 66:
                # This happens when the pocket number is higher than 99 (coordinates and later are displaced right)
                line = line[:28] + line[29:]
            return splitPDBLine(line)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

