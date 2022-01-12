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
from pwchem.viewers import ViewerGeneralPockets
import pyworkflow.protocol.params as params

class viewerP2Rank(ViewerGeneralPockets):
  _label = 'Viewer P2Rank pockets'
  _targets = [P2RankFindPockets]

  def __init__(self, **kwargs):
    super().__init__(**kwargs)

  def _defineParams(self, form):
    super()._defineParams(form)
    form.addSection(label='P2Rank Pymol visualization')
    form.addParam('displayP2Rank', params.LabelParam,
                  label='Display with P2Rank viewer',
                  help='Display pocket with own P2Rank visualization in pymol'
                  )

  def _getVisualizeDict(self):
    dispDic = super()._getVisualizeDict()
    dispDic.update({'displayP2Rank': self._showP2Rank})
    return dispDic

  def _showP2Rank(self, paramName=None):
      pdbFileName = self.protocol.getPdbInputStructName()
      outDir = os.path.abspath(self.protocol._getExtraPath('visualizations'))
      pmlFile = outDir + '/' + pdbFileName + '.pml'
      return self._showAtomStructPyMol(pmlFile, outDir)
