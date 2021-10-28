# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import os

from ..protocols import P2RankFindPockets
import pyworkflow.viewer as pwviewer
from pwchem.viewers import PyMolViewer, PocketPointsViewer, ContactSurfaceViewer
import pyworkflow.protocol.params as params

P2RANK_PYMOL, VOLUME_PYMOL, VOLUME_PYMOL_SURF = 0, 1, 2

class viewerP2Rank(pwviewer.ProtocolViewer):
  _label = 'Viewer pockets'
  _targets = [P2RankFindPockets]

  def __init__(self, **kwargs):
    pwviewer.ProtocolViewer.__init__(self, **kwargs)

  def _defineParams(self, form):
    form.addSection(label='Pymol visualization of predicted pockets')
    form.addParam('displayAtomStruct', params.EnumParam,
                  choices=['P2Rank Viewer', 'PyMol (Pocket Points)', 'PyMol (Contact Surface)'],
                  default=VOLUME_PYMOL,
                  display=params.EnumParam.DISPLAY_HLIST,
                  label='Display output AtomStruct with',
                  help='*Default Viewer*: display general pocket visualization in pymol\n'
                       '*P2Rank Viwer*: Display pocket with own P2Rank visualization in pymol'
                  )
    form.addParam('displayBBoxes', params.BooleanParam,
                  default=False, label='Display pocket bounding boxes',
                  help='Display the bounding boxes in pymol to check the size for the localized docking')
    form.addParam('pocketRadiusN', params.FloatParam, label='Grid radius vs pocket radius: ',
                  default=1.1, condition='displayBBoxes',
                  help='The radius * n of each pocket will be used as grid radius')

  def _getVisualizeDict(self):
    return {
      'displayAtomStruct': self._showAtomStruct,
    }

  def _showAtomStruct(self, paramName=None):
    if self.displayAtomStruct == VOLUME_PYMOL:
      return self._showAtomStructPyMolPoints()

    elif self.displayAtomStruct == VOLUME_PYMOL_SURF:
      return self._showAtomStructPyMolSurf()

    elif self.displayAtomStruct == P2RANK_PYMOL:
      pdbFileName = self.protocol.getPdbInputStructName()
      outDir = os.path.abspath(self.protocol._getExtraPath('visualizations'))
      pmlFile = outDir + '/' + pdbFileName + '.pml'
      return self._showAtomStructPyMol(pmlFile, outDir)

  def _showAtomStructPyMolPoints(self):
    bBox = self.displayBBoxes.get()
    if bBox:
      bBox = self.pocketRadiusN.get()

    pymolV = PocketPointsViewer(project=self.getProject())
    pymolV._visualize(self.protocol.outputPockets, bBox=bBox)

  def _showAtomStructPyMolSurf(self):
    bBox = self.displayBBoxes.get()
    if bBox:
      bBox = self.pocketRadiusN.get()

    pymolV = ContactSurfaceViewer(project=self.getProject())
    pymolV._visualize(self.protocol.outputPockets, bBox=bBox)

  def _showAtomStructPyMol(self, pmlFile, outDir):
    pymolV = PyMolViewer(project=self.getProject())
    pymolV.visualize(pmlFile, cwd=outDir)


class P2RankViewer(pwviewer.Viewer):
  _label = 'Viewer pockets'
  _environments = [pwviewer.DESKTOP_TKINTER]
  #_targets = [P2RankFindPockets]

  def _validate(self):
    return []

  # =========================================================================
  # ShowAtomStructs
  # =========================================================================

  def getInputAtomStructFile(self):
    return os.path.abspath(self.protocol.inputAtomStruct.get().getFileName())

  def _visualize(self, obj, **kwargs):
    pdbFileName = self.getInputAtomStructFile().split('/')[-1]
    outDir = os.path.abspath(self.protocol._getExtraPath('visualizations'))
    pymolFile = outDir + '/' + pdbFileName + '.pml'

    pymolV = PyMolViewer(project=self.getProject())
    pymolV.visualize(pymolFile, cwd=outDir)
