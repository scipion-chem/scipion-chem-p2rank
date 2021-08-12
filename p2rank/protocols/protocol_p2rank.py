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
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pyworkflow.utils import Message
import pyworkflow.utils as pwutils
from pwem.objects.data import AtomStruct
import pwem.convert as emconv
from pwem.convert.atom_struct import toPdb

import os, gzip
from p2rank import Plugin
from pwchem.objects import SetOfPockets
from ..objects import P2RankPocket
from pwchem.constants import PML_STR
from pwchem.utils import writeRawPDB, writePDBLine, writeSurfPML

class P2RankFindPockets(EMProtocol):
    """
    Executes the p2rank software to look for protein pockets.
    """
    _label = 'Find pockets'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputAtomStruct', params.PointerParam,
                       pointerClass='AtomStruct', allowsNull=False,
                       label="Input atom structure",
                       help='Select the atom structure to be fitted in the volume')
        form.addParam('numThreads', params.IntParam, default=1,
                       label="Number of threads",
                       help='Select the number of threads to perform the random forest')

    def _getP2RankArgs(self):
      args = ['-f', os.path.abspath(self.pdbFile)]
      args += ['-o', os.path.abspath(self._getExtraPath())]
      args += ['-threads', self.numThreads.get()]

      return args

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('P2RankStep')
        self._insertFunctionStep('createOutputStep')

    def convertInputStep(self):
      self.pdbFile = self._getPdbInputStruct()

    def P2RankStep(self):
        Plugin.runP2Rank(self, 'predict', args=self._getP2RankArgs(), cwd=self._getExtraPath())

    def createOutputStep(self):
        outFile = os.path.abspath(self.getOutFileName())
        outStruct = AtomStruct(outFile)
        self._defineOutputs(outputAtomStruct=outStruct)

        self.createOutputProteinFile()
        self.createPML()

        pocketFiles = self._divideoutputPockets()
        outPockets = SetOfPockets(filename=self._getExtraPath('pockets.sqlite'))
        for pFile in pocketFiles:
            pock = P2RankPocket(os.path.abspath(pFile), outFile)
            outPockets.append(pock)

        pmlFileName = '{}/{}_surf.pml'.format(self._getExtraPath(), self.getPDBName())
        writeSurfPML(outPockets, pmlFileName)
        self._defineOutputs(outputPockets=outPockets)


    # --------------------------- Utils functions --------------------
    def _getPdbInputStruct(self):
      inpStruct = self.inputAtomStruct.get()
      name, ext = os.path.splitext(inpStruct.getFileName())
      if ext == '.cif':
          cifFile = inpStruct.getFileName()
          pdbFile = self._getExtraPath(pwutils.replaceBaseExt(cifFile, 'pdb'))
          toPdb(cifFile, pdbFile)

      else:
        pdbFile = inpStruct.getFileName()
      return os.path.abspath(pdbFile)

    def getPdbInputStructName(self):
      return self._getPdbInputStruct().split('/')[-1]

    def getPDBName(self):
      return self.getPdbInputStructName().split('.')[0]
    
    def getOutFileName(self):
      pdbName = self.getPDBName()
      return self._getExtraPath('{}_out.pdb'.format(pdbName))

    def _divideoutputPockets(self):
      '''Creates individiual pocket files'''
      pocketDir = self._getExtraPath('pocketFiles')
      os.mkdir(pocketDir)

      pFiles = []
      pocketsFile = self._getExtraPath(self.getPdbInputStructName()) + '_predictions.csv'
      with open(pocketsFile) as fAll:
          headerLine = fAll.readline()
          for i, line in enumerate(fAll):
              pFile = os.path.join(pocketDir, self.getPdbInputStructName() + '_pocket{}.csv'.format(i+1))
              with open(pFile, 'w') as f:
                  f.write(headerLine)
                  f.write(line)
              pFiles.append(pFile)
      return pFiles

    def createOutputProteinFile(self):
      pdbName = self.getPDBName()
      oripdbFile = self._getPdbInputStruct()
      pdbFile = self._getExtraPath('{}_raw.pdb'.format(pdbName))
      writeRawPDB(oripdbFile, pdbFile, ter=False)

      gzfile = self._getExtraPath('visualizations/data/{}_points.pdb.gz'.format(
        self.getPdbInputStructName()))
      outFile = self.getOutFileName()
      with open(pdbFile) as fpdb:
        outStr = '\n'.join(fpdb.read().split('\n')[:-1])

      # Creates a pdb(qt) with the HETATOM corresponding to pocket points
      pocketDic = self.getPocketDic(gzfile)
      for pocketK in sorted(pocketDic):
        outStr += self.formatPocketStr(pocketDic[pocketK], pocketK)
      with open(outFile, 'w') as f:
        f.write(outStr)

    def createPML(self):
      pmlFile = self._getExtraPath('{}.pml'.format(self.getPDBName()))
      # Creates the pml for pymol visualization
      with open(pmlFile, 'w') as f:
        f.write(PML_STR.format(self.getOutFileName().split('/')[-1]))

    def createSurfPML(self, pockets):
      pmlFile2 = self._getExtraPath('{}_surf.pml'.format(self.getPDBName()))
      with open(pmlFile2, 'w') as f:
        f.write(createSurfacePml(pockets))

    def formatPocketStr(self, pocketLines, pocketK):
      outStr=''
      for i, pLine in enumerate(pocketLines):
          pLine = pLine.split()
          replacements = ['HETATM', str(i+1), 'APOL', 'STP', 'C', str(pocketK), *pLine[6:], '', 'Ve']
          pdbLine = writePDBLine(replacements)
          outStr += pdbLine
      return outStr

    def getPocketDic(self, pointsFile):
      dic={}
      with gzip.open(pointsFile) as f:
        for line in f:
          pocketId = int(line.split()[5].decode('utf-8'))
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

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def validate(self):
        """ Try to find errors on define params. """
        errors = []
        inpFile = os.path.abspath(self.inputAtomStruct.get().getFileName())
        if not 'pdb' in inpFile:
            nChains, nAtoms = self._countNumberOfChains(inpFile), self._countNumberOfAtoms(inpFile)
            if nChains > 62:
              errors.append('The atom structure file {} is too big for converting to pdb, '
                            'which is needed for running imodfit. Number of chains ({}) > 62'
                            .format(inpFile.split('/')[-1], nChains))
            elif nAtoms > 99999:
              errors.append('The atom structure file {} is too big for converting to pdb, '
                            'which is needed for running imodfit. Number of atoms ({}) > 99999'.
                            format(inpFile.split('/')[-1], nAtoms))

        return errors
