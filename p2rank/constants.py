# **************************************************************************
# *
# * Authors: Daniel Del Hoyo Gomez
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************


P2RANK_HOME = 'P2RANK_HOME'

# Supported Versions
V2_3 = '2.3'
P2RANK_DEFAULT_VERSION = V2_3

P2RANK = 'p2rank'

PYMOL_HOME = 'PYMOL_HOME'

PML_STR = '''from pymol import cmd,stored
load {}
#select pockets, resn STP
stored.list=[]
cmd.iterate("(resn STP)","stored.list.append(resi)")	#read info about residues STP
#print stored.list
lastSTP=stored.list[-1]	#get the index of the last residu
hide lines, resn STP

#show spheres, resn STP
for my_index in range(1,int(lastSTP)+1): cmd.select("pocket"+str(my_index), "resn STP and resi "+str(my_index))
for my_index in range(2,int(lastSTP)+2): cmd.color(my_index,"pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.show("spheres","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_scale","0.3","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_transparency","0.1","pocket"+str(my_index))'''

ATTRIBUTES_MAPPING = {'score': 'score', 'sas_points': 'nPoints', 'class': 'class',
                      'surf_atom_ids': 'contactAtoms', 'residue_ids': 'contactResidues'}
