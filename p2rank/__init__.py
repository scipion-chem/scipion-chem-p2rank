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

import pwem
from os.path import join, exists
from .constants import *

_version_ = '0.1'
_logo = "p2rank_logo.png"
_references = ['']


class Plugin(pwem.Plugin):
    _homeVar = P2RANK_HOME
    _pathVars = [P2RANK_HOME]
    _supportedVersions = [V2_3]
    _pluginHome = join(pwem.Config.EM_ROOT, P2RANK + '-' + P2RANK_DEFAULT_VERSION)

    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file.
        """
        cls._defineEmVar(P2RANK_HOME, P2RANK + '-' + P2RANK_DEFAULT_VERSION)

    @classmethod
    def defineBinaries(cls, env):
        installationCmd = 'wget %s -O %s && ' % (cls._getP2RankDownloadUrl(), cls._getP2RankTar())
        installationCmd += 'tar -xf %s --strip-components 1 && ' % cls._getP2RankTar()
        installationCmd += 'rm %s && ' % cls._getP2RankTar()

        # Creating validation file
        P2RANK_INSTALLED = '%s_installed' % P2RANK
        installationCmd += 'touch %s' % P2RANK_INSTALLED  # Flag installation finished

        env.addPackage(P2RANK,
                       version=P2RANK_DEFAULT_VERSION,
                       tar='void.tgz',
                       commands=[(installationCmd, P2RANK_INSTALLED)],
                       neededProgs=["conda"],
                       default=True)

    @classmethod
    def runP2Rank(cls, protocol, program, args, cwd=None):
        """ Run P2Rank command from a given protocol. """
        print(protocol)
        protocol.runJob(join(cls._pluginHome, 'prank {}'.format(program)), args, cwd=cwd)

    @classmethod  #  Test that
    def getEnviron(cls):
        pass

    # ---------------------------------- Utils functions  -----------------------
    @classmethod
    def _getP2RankDownloadUrl(cls):
        return "\'https://github.com/rdk/p2rank/releases/download/2.3/p2rank_2.3.tar.gz\'"

    @classmethod
    def _getP2RankTar(cls):
        pluginHome = join(pwem.Config.EM_ROOT, P2RANK + '-' + P2RANK_DEFAULT_VERSION)
        return pluginHome + '/' + P2RANK + '-' + P2RANK_DEFAULT_VERSION + '.tar.gz'

