from pyteomics import protxml
import numpy as np

class Quantify:

    def __init__(self) -> None:
        self._all_peptides = {}
    
    @staticmethod
    def combine_replicates(files):
        """combine dictionaries from all files into 
        one master dictionary"""
        master_list = [None]*len(files)
        for idx, file in enumerate(files):
            master_list[idx] = ExtractData(file).get_stpeter_data()
        return master_list

    def imput_missing(self, master_list, option=1):
        """replace zeros with something based
        on the imputation method selected, 0 for average,
        1 for minimum value"""
        self.all_peptides = self.get_all_peptides(master_list)
        if option == 0:
            self.imput_average(master_list)
        else:
            self.imput_minimum(master_list)

        return master_list

    def imput_average(self, master_list):
        for item in master_list:
            ave_si = sum(map(lambda x: item[x][0], item))/len(item)
            ave_nsaf = sum(map(lambda x: item[x][1], item))/len(item)
            missing_peptides = self.all_peptides.difference(set(item.keys()))
            for peptide in missing_peptides:
                item[peptide] = (ave_si, ave_nsaf)

    def imput_minimum(self, master_list):
        for item in master_list:
            ave_si = min(map(lambda x: item[x][0], item))
            ave_nsaf = min(map(lambda x: item[x][1], item))
            missing_peptides = self.all_peptides.difference(set(item.keys()))
            for peptide in missing_peptides:
                item[peptide] = (ave_si, ave_nsaf)
                
    @staticmethod
    def get_all_peptides(master_list):
        all_peptides = set.union(*[set(x.keys()) for x in master_list])
        return all_peptides

class ExtractData:
    """quantification of proteomics data"""

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

    
if __name__ == "__main__":
    FILES = ("interact-f10.ipro.prot.xml", "interact-f07.ipro.prot.xml")
    results = Quantify().combine_replicates(FILES)
    new = Quantify().imput_missing(results, option=0)

