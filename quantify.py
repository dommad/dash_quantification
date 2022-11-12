"""Read data from prot.xml file and extract
quantification information"""

import random
from pyteomics import protxml
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt


class Quantify:
    """class for quantification of proteomics data using
    SIn or NSAF"""

    def __init__(self):
        self._all_proteins = {}

    @staticmethod
    def combine_replicates(files):
        """combine dictionaries from all files into
        one master dictionary"""
        master_list = [None]*len(files)
        for idx, file in enumerate(files):
            master_list[idx] = ExtractData(file).get_stpeter_data()
        return master_list

    def imput_missing(self, master_list, option=1):
        """replace zeros with values based on
        the imputation method selected, 0 for average,
        1 for minimum value"""
        self._all_proteins = self.get_all_peptides(master_list)[0]
        if option:
            self.imput_minimum(master_list)
        else:
            self.imput_average(master_list)

        return master_list

    def imput_average(self, master_list):
        """replace zeros with average for the selected replicate"""
        for item in master_list:
            missing_proteins = self._all_proteins.difference(set(item.keys()))
            for protein in missing_proteins:
                item[protein] = self.get_average(master_list, protein)

    @staticmethod
    def get_average(master_list, key):
        """get average values of SIn and NSAF for selected protein (key) """
        sin = []
        nsaf = []
        for data in master_list:
            value = data.get(key)
            if value is not None:
                sin.append(2**value[0])
                nsaf.append(2**value[1])
        ave_sin = np.mean(sin)
        ave_nsaf = np.mean(nsaf)
        return (np.log2(ave_sin), np.log2(ave_nsaf))


    def imput_minimum(self, master_list):
        """replace zeros with minnimum for the selected replicate"""
        for item in master_list:
            min_si = min(x[0] for x in item)
            min_nsaf = min(x[1] for x in item)
            missing_proteins = self._all_proteins.difference(set(item.keys()))
            for protein in missing_proteins:
                item[protein] = (min_si, min_nsaf)


    def welch_test(self, conts, treats, mode='SIn'):
        """get p-values for Welch t-test between control and treatment"""
        # extract list of all proteins shared by the control and treatments
        cont_p, treat_p = self.get_all_peptides(conts, treats)
        all_proteins = set.union(cont_p, treat_p)
        mode_idx = 0
        if mode == 'NSAF':
            mode_idx = 1

        output = {} # protein: (p_value, log2(fold change)
        for prot in all_proteins:
            if prot not in cont_p:
                output[prot] = (1e-16, 5) # protein definitely upregulated
            elif prot not in treat_p:
                output[prot] = (1e-16, -5) # protein definitely downregulated
            else:
                cont_vals = self.fetch_quant_data(conts, prot, mode_idx)
                treat_vals = self.fetch_quant_data(treats, prot, mode_idx)
                p_val = st.ttest_ind(cont_vals, treat_vals, equal_var=False).pvalue
                output[prot] = (p_val, np.log2(np.mean(treat_vals)/np.mean(cont_vals)))
        return output

    @staticmethod
    def fetch_quant_data(master_list, prot, mode_idx):
        """retrieve data for given protein;
        to prevent failure of Welch's t test in case of nearly identical data,
        we may need to add some noise to the data"""
        data = [2**float(x[prot][mode_idx]) for x in master_list]
        if len(set(data)) == 1:
            return [x + random.random()*min(data)*0.0001 for x in data]
        return data

    @staticmethod
    def get_all_peptides(*args):
        """find set of all peptides in the replicates"""
        return list(map(lambda item: set.union(*[set(x.keys()) for x in item]), args))

    # plotting
    def volcano_plot(self, data, fdr_lim, fold_lim):
        """generate volcano plot for the quantification data with
        return information about down- and up-regulated proteins"""
        pval_cutoff = -np.log(self.get_fdr_threshold(data, fdr_lim))
        fold_idx = 1
        pval_idx = 0
        x_s = [x[fold_idx] for x in data.values()]
        y_s = [-np.log(x[pval_idx]) for x in data.values()]
        fig, axes = plt.subplots(figsize=(5,5))
        axes.scatter(x_s, y_s, s=1)
        axes.hlines(y=pval_cutoff, xmin=min(x_s), xmax=max(x_s), linestyles='dashed')
        axes.vlines(x=[fold_lim, -fold_lim], ymin=min(y_s), ymax=max(y_s), linestyles='dashed')
        axes.set_xlabel("log2(fold change)")
        axes.set_ylabel("-log(p-value)")
        fig.tight_layout()
        fig.savefig('volcano_plot.png', dpi=500)

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

class ExtractData:
    """extract quantification data from StPeter's prot.xml file"""

    def __init__(self, file) -> None:
        self.file = file

    def get_stpeter_data(self):
        """generate dictionary of data from StPeter's protXML output"""
        items = protxml.read(self.file)
        data = filter(lambda p: p != (), map(self.get_quant_data, items))
        dicts = {x[0]: (x[1], x[2]) for x in data}
        return dicts

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
