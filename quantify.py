"""Read data from prot.xml file and extract
quantification information"""

from pyteomics import protxml
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
from functools import lru_cache


class Quantify:
    """class for quantification of proteomics data using
    SIn or NSAF"""

    def __init__(self):
        pass


    def generate_master_dict(self, lists, impute_method='average'):
        """Generate master dictionary"""
        master_dict = {}
        rep_no = len(lists)

        for idx, file_name in enumerate(lists):
            cur_list = self.get_protein_list(file_name)
            for prot, cur_sin, cur_nsaf in cur_list:
                if prot not in master_dict:
                    master_dict[prot] = {'sin': rep_no*[None,], 'nsaf': rep_no*[None,]}
                master_dict[prot]['sin'][idx] = float(cur_sin)
                master_dict[prot]['nsaf'][idx]  = float(cur_nsaf)

        master_dict = self.impute_missing_values(master_dict, rep_no, impute_method)

        return master_dict


    def get_protein_list(self, file_name):
        """Get the protein list"""

        data = protxml.read(file_name)
        protein_info = [self.get_quant_data(x) for x in data]
        protein_info = [x for x in protein_info if x != ()]

        return protein_info


    @staticmethod
    def get_quant_data(item):
        """extract protein data from an entry in protXML file"""
        protein = item['protein'][0]
        if 'analysis_result' in protein.keys():
            protein_name = protein['protein_name'].split('|')[-1]
            quant_data = protein['analysis_result'][0]['StPeterQuant']
            si_n = np.float32(quant_data['SIn'])
            nsaf = np.float32(quant_data['NSAF'])
            return (protein_name, si_n, nsaf)
        return ()


    @staticmethod
    def impute_missing_values(all_dicts, no_rep, method='average'):
        """Impute the missing values"""
        def impute_ave(scores):
            sin = scores['sin']
            ave_sin = np.mean([x for x in sin if x is not None])
            nsaf = scores['nsaf']
            ave_nsaf = np.mean([x for x in nsaf if x is not None])

            new_sin = [x if isinstance(x, float) else ave_sin + np.random.rand()*0.001 for x in sin ]
            new_nsaf = [x if isinstance(x, float) else ave_nsaf + np.random.rand()*0.001 for x in nsaf]

            return {'sin': new_sin, 'nsaf': new_nsaf}

        def get_minima(all_dicts, no_rep):
            min_sins = []
            min_nsafs = []
            for idx in range(no_rep):
                min_sin = min((all_dicts[x]['sin'][idx] for x in all_dicts.keys() if all_dicts[x]['sin'][idx] is not None))
                min_nsaf = min((all_dicts[x]['nsaf'][idx] for x in all_dicts.keys() if all_dicts[x]['sin'][idx] is not None))
                min_sins.append(min_sin)
                min_nsafs.append(min_nsaf)
            return min_sins, min_nsafs
        
        def impute_min(scores, min_sins, min_nsafs):

            new_sin = [x if isinstance(x, float) else min_sins[idx] + np.random.rand()*0.001 for idx, x in enumerate(scores['sin']) ]
            new_nsaf = [x if isinstance(x, float) else min_nsafs[idx] + np.random.rand()*0.001 for idx, x in enumerate(scores['nsaf'])]

            return {'sin': new_sin, 'nsaf': new_nsaf}


        if method == 'average':
            new_dict = dict((k, impute_ave(v)) if ((None in v['sin']) | (None in v['nsaf'])) else (k, v) for k, v in all_dicts.items())

        if method == 'minimum':
            min_sins, min_nsafs = get_minima(all_dicts, no_rep)
            new_dict = dict((k, impute_min(v, min_sins, min_nsafs)) if ((None in v['sin']) | (None in v['nsaf'])) else (k, v) for k, v in all_dicts.items())

        return new_dict


    @staticmethod
    def welch_test(conts, treats, mode='sin'):
        """get p-values for Welch t-test between control and treatment"""
        # extract list of all proteins shared by the control and treatments
        all_proteins = set.union(set(conts.keys()), set(treats.keys()))

        output = {} # protein: (p_value, log2(fold change)
        for prot in all_proteins:
            if prot not in conts:
                output[prot] = (1e-16, 5) # protein definitely upregulated
            elif prot not in treats:
                output[prot] = (1e-16, -5) # protein definitely downregulated
            else:
                cont_vals = conts[prot][mode]
                treat_vals = treats[prot][mode]
                p_val = st.ttest_ind(cont_vals, treat_vals, equal_var=False).pvalue
                output[prot] = (p_val, np.mean(treat_vals)-np.mean(cont_vals))

        return output


    @staticmethod
    def volcano_plot(data, fdr_lim, fold_lim):
        """generate volcano plot for the quantification data with
        return information about down- and up-regulated proteins"""
        #pval_cutoff = -np.log(self.get_fdr_threshold(data, fdr_lim))

        fold_vals = [x[1] for x in data.values()]
        log_pvals = -np.log10([x[0] for x in data.values()])

        c_accept = 'green'
        c_reject = (0.90196078, 0.81568627, 0.85882353, 0.8)
        cutoff_style = 'dotted'
        cutoff_color = 'grey'

        x_lims = (min(fold_vals)-1, max(fold_vals)+1)
        y_lims = (0, max(log_pvals) + 1)


        colors = [c_accept if abs(x) > fold_lim and abs(y) > fdr_lim else c_reject for (x, y) in zip(fold_vals, log_pvals)]
        
        fig, axes = plt.subplots(figsize=(5,5))
        axes.scatter(fold_vals, log_pvals, s=5, color=colors)
        axes.hlines(y=fdr_lim, xmin=x_lims[0], xmax=x_lims[1], linestyles=cutoff_style, color=cutoff_color)
        axes.vlines(x=[fold_lim, -fold_lim], ymin=y_lims[0], ymax=y_lims[1], linestyles=cutoff_style, color=cutoff_color)
        axes.set_xlabel("log2(fold change)")
        axes.set_ylabel("-log(p-value)")
        axes.set_ylim(y_lims[0], y_lims[1])
        axes.set_xlim(x_lims[0], x_lims[1])
        fig.tight_layout()



    def get_stats(self, data, fdr_lim, fold_lim):
        """get the p-value cutoff for the given FDR threshold
        and output the list of up and downregulated proteins"""
        pval_th = self.get_fdr_threshold(data, fdr_lim)
        sel_data = [x for x in data.items() if (x[1][0] <= pval_th and  abs(x[1][1]) >= fold_lim)]
        #up_reg = [x for x in data.items() if (x[1][0] <=pval_cutoff and  x[1][1] >= fold_lim)]
        #down_reg = [x for x in data.items() if (x[1][0] <=pval_cutoff and  -x[1][1] >= fold_lim)]
        return sel_data


    @staticmethod
    def get_fdr_threshold(data, fdr_lim):
        """find the adjusted p-value cutoff for the specified
        FDR threshold"""
        q_idx = 1 # index of quantification data in the tuple
        p_idx = 0 # index of p-value in the quantification data tuple
        sorted_d = sorted(data.items(), key=lambda x: x[q_idx][p_idx]) # sort by p-values
        adj_ps = [(idx+1)*val[1][0]/len(sorted_d) for idx, val in enumerate(sorted_d)]
        cutoff = sorted_d[len([x for x in adj_ps if x <= fdr_lim])-1][q_idx][p_idx]
        return cutoff
