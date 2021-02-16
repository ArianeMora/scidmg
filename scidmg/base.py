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

from sciviso import Volcanoplot
import pandas as pd
import matplotlib.pyplot as plt
import os
from sciutil import SciUtil, SciException


class SciDMGException(SciException):
    def __init__(self, message=''):
        Exception.__init__(self, message)


class SciDMG:

    def __init__(self, dmr_file: str, perc_meth_file: str, dmc_file: str, dmc_methdiff='meth.diff', dmr_methdiff='meth_diff',
                 dmc_padj='qvalue', dmr_padj='fdr', gene_name='external_gene_name', dmc_uid='uid', top_cpg_method='meth_diff',
                 dmr_uid='uid', min_perc_agreement=0.6, min_meth_diff=0.01, sciutil=None, plot=False, output_dir='.',
                 dmg_filename='DMG_output.csv', dmc_padj_cutoff=0.1, dmr_padj_cutoff=0.1):
        self.u = SciUtil() if sciutil is None else sciutil
        self.dmr_df = pd.read_csv(dmr_file)
        self.dmc_df = pd.read_csv(dmc_file)
        self.dmg_df = None
        if perc_meth_file is not None:
            self.perc_meth_df = pd.read_csv(perc_meth_file)
            self.merge_perc_meth = True
        else:
            self.merge_perc_meth = False
        self.dmc_uid = dmc_uid
        self.dmc_methdiff = dmc_methdiff
        self.dmr_methdiff = dmr_methdiff
        self.dmc_padj = dmc_padj
        self.dmr_padj = dmr_padj
        self.gene_name = gene_name
        self.dmr_uid = dmr_uid
        self.uid = 'gene_dmc-uid'
        self.plot = plot
        self.min_meth_diff = min_meth_diff
        self.min_perc_agreement = min_perc_agreement
        self.top_cpg_method = top_cpg_method if top_cpg_method in ['meth_diff', 'pvalue'] else None
        self.output_dir = output_dir
        self.dmg_filename = dmg_filename
        self.dmc_padj_cutoff = dmc_padj_cutoff
        self.dmr_padj_cutoff = dmr_padj_cutoff
        self.stats = {}
        if not self.top_cpg_method:
            msg = f'ERR: you passed a non-accepted value for the top cpg method, passed {top_cpg_method},' \
                  f' expected: meth_diff or pvalue'
            self.u.err_p([msg])
            raise SciDMGException(msg)

    def run(self):
        # Generate unique IDs
        self.__gen_uid()
        # Plot Volcanos
        if self.plot:
            self.plot_dmc_volcano()
            self.plot_dmr_volcano()
        # Merge the percentage methylation results with the DMC analysis
        if self.merge_perc_meth:
            self.merge_dmc_perc_meth()
        # Filter significant genes from both DMC and DMRs
        self.filter_sig()
        # Merge the DMR with the DMC
        self.merge_dmc_dmr()
        # Group by DMRs and only keep the DMRs that have at least 60% of CpGs that agree with the direction
        self.__grp_by_dmr(self.min_perc_agreement, self.min_meth_diff)
        self.__grp_by_gene(self.min_meth_diff)
        # Keep only the top CpG
        self.__filter_top_cpg()
        # Lastly save to DF
        self.dmg_df.to_csv(os.path.join(self.output_dir, self.dmg_filename))

    def get_stats(self):
        return self.stats

    def print_stats(self):
        self.u.dp(["Printing stats: "])
        for s in self.stats:
            print(s, self.stats[s])

    def __grp_by_gene(self, cutoff=0.01):
        """
        Group by Gene.
        Once we have filtered on the DMRs (previous step) we want to ensure that we only get one CpG per gene,
        which we denote the "driver CpG".
        """
        meth_all_df_filtered_grouped = self.dmg_df.groupby(self.gene_name)

        self.stats['Length grouped by Genes'] = len(meth_all_df_filtered_grouped)  # Keep track of some of the stats

        self.u.dp(["Length of filtered methylation dataframe: ", len(self.dmg_df),
                   "\nNumber of genes with Methylation:", len(meth_all_df_filtered_grouped)])

        cpgs_to_keep = {}
        meths_to_drop = []
        for m in meth_all_df_filtered_grouped:
            region_id = m[0]
            dmc_meth_diffs = m[1][self.dmc_methdiff].values  # Want the largest value from this that agrees with the
            dmr_meth_diffs = m[1][self.dmc_methdiff].values  # DMR methdiff --> note if we have conflicting DMR meth diffs we drop it

            meth_ids = m[1][self.dmc_uid].values
            qvals = m[1][self.dmc_padj].values
            count_pos = 0
            count_neg = 0
            count_pos_dmr, count_neg_dmr = 0, 0
            cutoff_pos, cutoff_neg = cutoff, -1 * cutoff  # i.e. has non-zero data
            max_cpg_idx, max_cpg = None, 0
            for i, meth_diff in enumerate(dmc_meth_diffs):
                if qvals[i] <= self.dmc_padj_cutoff:  # Ensure we're only looking at the DMCs
                    if meth_diff < cutoff_neg:
                        count_neg += 1
                    elif meth_diff > cutoff_pos:
                        count_pos += 1
                    if abs(meth_diff) > abs(max_cpg):
                        max_cpg = meth_diff
                        max_cpg_idx = i
                    # Also keep track of teh DMR direction
                    if dmr_meth_diffs[i] > 0:
                        count_pos_dmr += 1
                    elif dmr_meth_diffs[i] < 0:
                        count_neg_dmr += 1
            # if we have any disagreement between the regions we drop it
            if count_pos_dmr > 0 and count_neg_dmr > 0:
                meths_to_drop.append(region_id)
            elif count_pos > 0:  # Let's just keep the max cpg idx
                cpgs_to_keep[meth_ids[max_cpg_idx]] = count_pos
            elif count_neg > 0:
                cpgs_to_keep[meth_ids[max_cpg_idx]] = count_neg
        # Create a new column with the meth keep idxs
        meths_to_keep = []
        for meth_id in self.dmg_df[self.dmc_uid].values:
            if cpgs_to_keep.get(meth_id):
                meths_to_keep.append(1)
            else:
                meths_to_keep.append(0)

        self.dmg_df['CpGsToKeep'] = meths_to_keep
        self.dmg_df = self.dmg_df[self.dmg_df['CpGsToKeep'] == 1]
        self.u.dp(["Dropping any genes with disagreeing DMRs: ", len(meths_to_drop)])

        self.stats['Number of Genes with DMRs that disagree in direction'] = len(meths_to_drop)
        self.dmg_df.to_csv('DMG_after_grouped_by_gene.csv', index=False)

    def __filter_top_cpg(self):
        # Keep only the first CpG i.e. the one with the largest change
        if self.top_cpg_method == 'meth_diff':
            self.dmg_df['abs_logfc'] = abs(self.dmg_df[self.dmc_methdiff].values)
            self.dmg_df = self.dmg_df.sort_values([self.gene_name, 'abs_logfc'], ascending=False).drop_duplicates(subset=[self.gene_name], keep="first")
        elif self.top_cpg_method == 'pvalue':
            # Select the smallest p-value from the DMC
            self.dmg_df = self.dmg_df.sort_values([self.gene_name, self.dmc_padj], ascending=True).drop_duplicates(
                subset=[self.gene_name], keep="first")
        self.u.dp(["Length of dataframe filtered to only keep top DMC mapped to genes:", len(self.dmg_df)])
        self.dmg_df.to_csv('DMG_after_selecting_top_cpg.csv', index=False)

    def __grp_by_dmr(self, min_agreement=0.6, cutoff=0.01):
        """
        Group by DMR ID.
        Here we count whether there are at least X percent of CpGs that agree with our call of the differentially
        methylated region.
        """
        meth_grouped = self.dmg_df.groupby(self.dmr_uid)
        self.stats['Length grouped by DMRs'] = len(meth_grouped)  # Keep track of some of the stats

        cpg_to_keep = {}
        dmrs_to_keep = {}
        non_sigs = 0
        for m in meth_grouped:
            region_qvalue = m[1][self.dmr_padj].values[0]
            dmr_id = m[1][self.dmr_uid].values[0]
            if region_qvalue <= self.dmr_padj_cutoff:  # Only keep regions that meet the qvalue threshold
                meth_diffs = m[1][self.dmc_methdiff].values
                meth_ids = m[1][self.dmc_uid].values
                qvals = m[1][self.dmc_padj].values

                count_pos = 0
                count_neg = 0
                cutoff_pos, cutoff_neg = cutoff, -1 * cutoff
                max_cpg_idx, max_cpg = None, 0
                # Get the direction of the DMRseq region
                region_direction = m[1][self.dmr_methdiff].values[0]
                for i, meth_diff in enumerate(meth_diffs):
                    if qvals[i] <= self.dmc_padj_cutoff:  # Only keep CpGs that meet the qvalue threshold
                        if meth_diff < cutoff_neg:
                            count_neg += 1
                        elif meth_diff > cutoff_pos:
                            count_pos += 1
                        # Keep max CpG in that direction
                        if region_direction > 0 and meth_diff > max_cpg:
                            max_cpg = meth_diff
                            max_cpg_idx = i
                        elif region_direction < 0 and meth_diff < max_cpg:
                            max_cpg = meth_diff
                            max_cpg_idx = i

                if count_pos == 0 and count_neg == 0:
                    non_sigs += 1
                elif region_direction > 0:
                    if count_pos > 0 and count_neg == 0:
                        cpg_to_keep[meth_ids[max_cpg_idx]] = count_pos
                        dmrs_to_keep[dmr_id] = True
                    elif count_pos / (count_neg * 1.0) >= min_agreement:
                        # Check the max cpg is in the same direction as the dmrseq region
                        if max_cpg > 0:
                            cpg_to_keep[meth_ids[max_cpg_idx]] = True
                            dmrs_to_keep[dmr_id] = True
                        else:  # We need to find the max that agrees with the direction > 0
                            max_cpg_idx, max_cpg = None, 0
                            for i, meth_diff in enumerate(meth_diffs):
                                if meth_diff > max_cpg:
                                    max_cpg = meth_diff
                                    max_cpg_idx = i
                            cpg_to_keep[meth_ids[max_cpg_idx]] = True
                            dmrs_to_keep[dmr_id] = True
                else:
                    if count_pos == 0 and count_neg > 0:
                        cpg_to_keep[meth_ids[max_cpg_idx]] = count_neg
                        dmrs_to_keep[dmr_id] = True
                    elif count_neg / (count_pos * 1.0) >= min_agreement:
                        # Check the max cpg is in the same direction as the dmrseq region
                        if max_cpg < 0:
                            cpg_to_keep[meth_ids[max_cpg_idx]] = True
                            dmrs_to_keep[dmr_id] = True
                        else:  # We need to find the max that agrees with the direction > 0
                            max_cpg_idx, max_cpg = None, 0
                            for i, meth_diff in enumerate(meth_diffs):
                                if meth_diff < max_cpg:
                                    max_cpg = meth_diff
                                    max_cpg_idx = i
                            cpg_to_keep[meth_ids[max_cpg_idx]] = True
                            dmrs_to_keep[dmr_id] = True
        # Print out the results
        self.u.dp(["Length of all merged methylation data: ", len(self.dmg_df)])
        self.u.dp(["Length of merged methylation data grouped by region: ", len(meth_grouped)])
        self.u.dp(["Number of CpGs to keep based on the regions: ", len(cpg_to_keep)])

        self.stats['Number of CpGs to keep from grouped DMRs'] = len(cpg_to_keep)  # Keep track of some of the stats
        self.stats['Length of merged DMR and DMC'] = len(self.dmg_df)

        # Create a new column with the meth keep idxs
        meths_to_keep = []
        dmr_ids = self.dmg_df[self.dmr_uid].values
        for i, meth_id in enumerate(self.dmg_df[self.dmc_uid].values):
            if cpg_to_keep.get(meth_id) and dmrs_to_keep.get(dmr_ids[i]):
                meths_to_keep.append(1)
            else:
                meths_to_keep.append(0)

        self.dmg_df['CpGsToKeep'] = meths_to_keep
        self.dmg_df = self.dmg_df[self.dmg_df['CpGsToKeep'] == 1]
        self.dmg_df.to_csv('DMG_after_grouped_by_DMR.csv', index=False)

    def plot_dmc_volcano(self, p=0.05, methdiff=25, output_dir='.'):
        # Volcano plot of the fold change nad padjusted
        if self.dmg_df is None:
            self.u.warn_p(["Please run dmg.run() before plotting the volcanos."])
            return
        volcanoplot = Volcanoplot(self.dmg_df, self.dmc_methdiff, self.dmc_padj, self.dmc_uid,
                                  'DNA methylation change and significance (DMC)',
                                  'Methylation difference', '-log10(p adj)',
                                  p_val_cutoff=p,
                                  label_big_sig=True, log_fc_cuttoff=methdiff)
        volcanoplot.plot()
        plt.savefig(os.path.join(self.output_dir, 'volcano_meth_DMC-unfiltered.pdf'))
        plt.show()

    def plot_dmr_volcano(self, p=0.05, methdiff=25):
        # Volcano plot of the fold change nad padjusted
        if self.dmg_df is None:
            self.u.warn_p(["Please run dmg.run() before plotting the volcanos."])
            return
        volcanoplot = Volcanoplot(self.dmg_df, self.dmr_methdiff, self.dmr_padj, self.dmr_uid,
                                  'DNA methylation change and significance (DMR)',
                                  'Methylation difference', '-log10(p adj)',
                                  p_val_cutoff=p,
                                  label_big_sig=True, log_fc_cuttoff=methdiff)
        volcanoplot.plot()
        plt.savefig(os.path.join(self.output_dir, 'volcano_meth_DMR-unfiltered.pdf'))
        plt.show()

    def merge_dmc_perc_meth(self):
        """ Assumes both are the unaltered output of MethylSig. """
        self.u.warn_p(['WARNING: Running merge_dmc_perc_meth which assumes that you have run the DMC '
                       'analysis using MethylKit.'])
        genes = self.perc_meth_df[self.gene_name].values
        genes_cpgs = self.dmc_df[self.gene_name].values

        if genes.all() != genes_cpgs.all():
            print('UNABLE TO JOIN')

        # Since they match perfectly, we can just do a normal join
        for c in self.perc_meth_df.columns:
            self.dmc_df[c] = self.perc_meth_df[c].values

    def merge_dmc_dmr(self):
        """ merges the two datasets """

        self.dmg_df = self.dmc_df.merge(self.dmr_df, on=self.gene_name,
                                         how='inner', suffixes=['_dmc', '_dmr']).drop_duplicates()
        # Update the iDs if they were the same
        if self.dmc_uid == self.dmr_uid:
            self.dmc_uid = self.dmc_uid + '_dmc'
            self.dmr_uid = self.dmr_uid + '_dmr'

    def filter_sig(self, dmc_padj_cutoff=0.1, dmr_padj_cutoff=0.1):
        self.dmc_df = self.dmc_df[self.dmc_df[self.dmc_padj] <= dmc_padj_cutoff]
        self.dmr_df = self.dmr_df[self.dmr_df[self.dmr_padj] <= dmr_padj_cutoff]

    def __gen_uid(self):
        gene_uid = []
        genes = self.dmc_df[self.gene_name].values
        for i, u_id in enumerate(self.dmc_df[self.dmc_uid].values):
            gene_uid.append(f'{genes[i]} {u_id}')
        self.dmc_df[self.uid] = gene_uid

        gene_uid = []
        genes = self.dmr_df[self.gene_name].values
        for i, u_id in enumerate(self.dmr_df[self.dmr_uid].values):
            gene_uid.append(f'{genes[i]} {u_id}')
        self.dmr_df[self.uid] = gene_uid




