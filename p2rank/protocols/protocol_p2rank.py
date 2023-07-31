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
import os, gzip, shutil

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pyworkflow.object import String
from pwem.protocols import EMProtocol
import pwem.convert as emconv
from pwem.convert.atom_struct import toPdb

from pwchem.objects import SetOfStructROIs, PredictStructROIsOutput, StructROI
from pwchem.utils import writePDBLine, splitPDBLine, runOpenBabel

from p2rank import Plugin


class P2RankFindPockets(EMProtocol):
    """
    Executes the p2rank software to look for protein pockets.
    """
    _label = 'Find pockets'
    _possibleOutputs = PredictStructROIsOutput

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputAtomStruct', params.PointerParam,
                       pointerClass='AtomStruct', allowsNull=False,
                       label="Input atom structure",
                       help='Select the atom structure to be fitted in the volume')
        form.addParallelSection(threads=4, mpi=1)

    def _getP2RankArgs(self):
      args = ['-f', os.path.abspath(self.pdbFile)]
      args += ['-o', os.path.abspath(self._getExtraPath())]
      args += ['-threads', self.numberOfThreads.get()]

      return args

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('P2RankStep')
        self._insertFunctionStep('createOutputStep')

    def convertInputStep(self):
      self.pdbFile = self._convertInputPDB()

    def P2RankStep(self):
        Plugin.runP2Rank(self, 'predict', args=self._getP2RankArgs(), cwd=self._getExtraPath())

    def createOutputStep(self):
        inpStruct = self.inputAtomStruct.get()
        outASPath = os.path.relpath(self._getPDBFile())
        pocketFiles = self._divideOutputPockets()

        outPockets = SetOfStructROIs(filename=self._getExtraPath('StructROIs.sqlite'))
        for pFile in pocketFiles:
            pock = StructROI(pFile, outASPath, self.getPropertiesFile(), pClass='P2Rank')
            if len(pock.getPointsCoords()) > 2: #minimum size for building pocket. cannot calculate volume otherwise
                pock.setVolume(pock.getPocketVolume())
                if str(type(inpStruct).__name__) == 'SchrodingerAtomStruct':
                    pock._maeFile = String(inpStruct.getFileName())
                outPockets.append(pock)

        outPockets.buildPDBhetatmFile()
        self._defineOutputs(**{self._possibleOutputs.outputStructROIs.name: outPockets})


    # --------------------------- Utils functions --------------------
    def _getInputName(self):
        return os.path.splitext(os.path.basename(self.inputAtomStruct.get().getFileName()))[0]

    def _getPDBFile(self):
        return os.path.abspath(self._getExtraPath(self._getInputName() + '.pdb'))

    def _convertInputPDB(self):
      inpStruct = self.inputAtomStruct.get()
      name, ext = os.path.splitext(inpStruct.getFileName())
      if ext == '.cif':
          cifFile = inpStruct.getFileName()
          toPdb(cifFile, self._getPDBFile())

      elif str(type(inpStruct).__name__) == 'SchrodingerAtomStruct':
          inpStruct.convert2PDB(outPDB=self._getPDBFile())

      elif ext == '.pdbqt':
          pdbFile = os.path.abspath(self._getPDBFile())
          args = ' -ipdbqt {} -opdb -O {}'.format(os.path.abspath(inpStruct.getFileName()), pdbFile)
          runOpenBabel(protocol=self, args=args, cwd=self._getExtraPath())

      else:
          shutil.copy(inpStruct.getFileName(), self._getPDBFile())
      return self._getPDBFile()

    def getPdbInputStructName(self):
      return self._getPDBFile().split('/')[-1]

    def getPropertiesFile(self):
        return self._getExtraPath(self.getPdbInputStructName()+'_predictions.csv')

    def _divideOutputPockets(self):
      '''Creates individiual pocket files'''
      gzfile = self._getExtraPath('visualizations/data/{}_points.pdb.gz'.format(
        self.getPdbInputStructName()))
      pocketDic = self.getPocketDic(gzfile)

      os.mkdir(self._getExtraPath('pocketFiles'))
      pFiles = []
      for pocketK in sorted(pocketDic):
          pFile = self._getExtraPath('pocketFiles/pocketFile_{}.pdb'.format(pocketK))
          with open(pFile, 'w') as f:
              f.write(self.formatPocketStr(pocketDic[pocketK], pocketK))
          pFiles.append(pFile)
      return pFiles

    def formatPocketStr(self, pocketLines, pocketK):
      outStr=''
      for i, pLine in enumerate(pocketLines):
          pLine = self.splitP2RankPDBLine(pLine)
          replacements = ['HETATM', str(i+1), 'APOL', 'STP', 'C', str(pocketK), *pLine[6:], '', 'Ve']
          pdbLine = writePDBLine(replacements)
          outStr += pdbLine
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

    def validate(self):
        """ Try to find errors on define params. """
        errors = []
        inpStruct = self.inputAtomStruct.get()
        inpFile = os.path.abspath(inpStruct.getFileName())

        if str(type(inpStruct).__name__) == 'SchrodingerAtomStruct':
            inpFile = inpStruct.convert2PDB()

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
