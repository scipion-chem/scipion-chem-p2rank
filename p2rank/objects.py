# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Alberto Manuel Parra Pérez (amparraperez@gmail.com)
# *
# * Biocomputing Unit, CNB-CSIC
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************


from pwchem.objects import ProteinPocket
from .constants import ATTRIBUTES_MAPPING as AM

class P2RankPocket(ProteinPocket):
  """ Represent a pocket file from p2rank"""
  def __init__(self, filename=None, proteinFile=None, **kwargs):
    self.properties, self.pocketId = self.parseFile(filename)
    kwargs.update(self.getKwargs(self.properties, AM))
    super().__init__(filename, proteinFile, **kwargs)
    self.setObjId(self.pocketId)

  def __str__(self):
    s = 'P2Rank pocket {}\nFile: {}'.format(self.pocketId, self.getFileName())
    return s

  def getVolume(self):
    return self.properties['volume']

  def getPocketScore(self):
    return self.properties['score']

  def getNumberOfPoints(self):
    return self.properties['sas_points']

  def getCenterOfMass(self):
    return self.properties['center_x'], self.properties['center_y'], self.properties['center_z']

  def parseFile(self, filename):
    props = {}
    with open(filename) as f:
      keys = f.readline().split(',')
      values = f.readline().split(',')
    for i, k in enumerate(keys):
      props[k.strip()] = values[i]

    props['residue_ids'] = props['residue_ids'].strip().replace(' ', '-')
    props['surf_atom_ids'] = props['surf_atom_ids'].strip().replace(' ', '-')
    props['class'] = 'P2Rank'

    return props, int(values[1])

