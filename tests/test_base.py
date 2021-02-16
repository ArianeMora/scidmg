###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import os
import shutil
import tempfile
import unittest

from scidmg import SciDMG


class TestClass(unittest.TestCase):

    @classmethod
    def setup_class(self):
        local = True
        # Create a base object since it will be the same for all the tests
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))

        self.data_dir = os.path.join(THIS_DIR, 'data/')
        if local:
            self.tmp_dir = os.path.join(THIS_DIR, 'data/tmp/')
            if os.path.exists(self.tmp_dir):
                shutil.rmtree(self.tmp_dir)
            os.mkdir(self.tmp_dir)
        else:
            self.tmp_dir = tempfile.mkdtemp(prefix='scircm_tmp_')
        # Setup the default data for each of the tests
        self.dmc_file = os.path.join(self.data_dir, 'methylKit_DMC.csv')
        self.percent_meth = os.path.join(self.data_dir, 'methylKit_percentMeth.csv')
        self.dmr_file = os.path.join(self.data_dir, 'methylSig_prom.csv')

    @classmethod
    def teardown_class(self):
        shutil.rmtree(self.tmp_dir)


class TestSciRCM(TestClass):

    def test_run_dmg(self):
        dmg = SciDMG(self.dmr_file, self.percent_meth, self.dmc_file, 'meth.diff', 'meth_diff', 'qvalue', 'fdr', 'external_gene_name',
                     dmc_uid='uid', top_cpg_method='meth_diff', dmr_uid='uid', min_perc_agreement=0.5, min_meth_diff=0.01
                     )
        dmg.run()
        dmg.print_stats()
        stats = dmg.get_stats()
        assert stats['Length grouped by DMRs'] == 6
        assert stats['Number of CpGs to keep from grouped DMRs'] == 4  # Since we have two DMRs with multiple CpGs
        assert stats['Length of merged DMR and DMC'] == 16
        assert stats['Length grouped by Genes'] == 3
        assert stats['Number of Genes with DMRs that disagree in direction'] == 1
