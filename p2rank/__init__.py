# **************************************************************************
# *
# * Authors:  Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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

from os.path import join

import pwem
from scipion.install.funcs import InstallHelper

from pwchem import Plugin as pwchemPlugin

from .constants import *

_version_ = '0.1'
_logo = "p2rank_logo.png"
_references = ['']

P2RANK_DIC = {'name': 'p2rank', 'version': '2.5.1', 'home': 'P2RANK_HOME', 'java': "17.0"}


class Plugin(pwchemPlugin):
    _homeVar = P2RANK_DIC['home']
    _pathVars = [P2RANK_DIC['home']]
    _supportedVersions = [P2RANK_DIC['version']]

    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file.
        """
        cls._defineEmVar(P2RANK_DIC['home'], P2RANK_DIC['name'] + '-' + P2RANK_DIC['version'])

    @classmethod
    def defineBinaries(cls, env):
        installationCmd = 'wget %s -O %s && ' % (cls._getP2RankDownloadUrl(), cls._getP2RankTar())
        installationCmd += 'tar -xf %s --strip-components 1 && ' % cls._getP2RankTar()
        installationCmd += 'rm %s ' % cls._getP2RankTar()

        # Instantiating the install helper
        installer = InstallHelper(P2RANK_DIC['name'], packageHome=cls.getVar(P2RANK_DIC['home']),
                                  packageVersion=P2RANK_DIC['version'])

        # Generating AutoSite installation commands
        installer.addCommand(installationCmd, 'P2RANK_DOWNLOADED') \
            .getCondaEnvCommand(P2RANK_DIC['name'], binaryVersion=P2RANK_DIC['version'], pythonVersion='3.10') \
            .addCondaPackages([f'openjdk={P2RANK_DIC["java"]}'], channel='conda-forge', targetName='JAVA_CONDA')\
            .addPackage(env, ['conda'])

    @classmethod
    def runP2Rank(cls, protocol, program, args, cwd=None):
        """ Run P2Rank command from a given protocol. """
        actEnv = f'{cls.getEnvActivationCommand(P2RANK_DIC)} && '
        p2RankCommand = actEnv + join(cls.getVar(P2RANK_DIC['home']), f'prank {program}')
        protocol.runJob(p2RankCommand, args, cwd=cwd)

    # ---------------------------------- Utils functions  -----------------------
    @classmethod
    def _getP2RankDownloadUrl(cls):
        return "\'https://github.com/rdk/p2rank/releases/download/{}/p2rank_{}.tar.gz\'".\
            format(P2RANK_DIC['version'], P2RANK_DIC['version'])

    @classmethod
    def _getP2RankTar(cls):
        pluginHome = join(pwem.Config.EM_ROOT, P2RANK_DIC['name'] + '-' + P2RANK_DIC['version'])
        return pluginHome + '/' + P2RANK_DIC['name'] + '-' + P2RANK_DIC['version'] + '.tar.gz'

