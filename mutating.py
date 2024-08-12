## This code addes mutation found by CoronApp in epitopes. 


import pandas as pd
import numpy as np
import seaborn as sns

# Создается словарь вида:
#
# белок : 
#
#      { [i] : 
#      
#            {начало:, конец:, 
#            эпитоп до:[], эпитоп текущий:[] (обновлять после каждого пациента до исх.), эпитопы после:set(),
#            число мутаций за каждый из проходов:[], мутаций за всё время:[], hlas:[]}, ...}
#            
# Такой словарь надо инициировать заново для каждого пациента, в конце сохранять. 


# словарь для перевода белков из таблички с эпитопами в такие, как у coronapp
prot_translation = {'leader':'NSP1', 'nsp2':'NSP2', 'nsp3':'NSP3', 'nsp4':'NSP4', '3C':'NSP5',
                    'nsp6':'NSP6', 'nsp7':'NSP7', 'nsp8':'NSP8', 'nsp9':'NSP9', 'nsp10':'NSP10',
                    'RdRp':'NSP12', 'helicase':'NSP13', 'exonuclease':'NSP14', 'endornase':'NSP15',
                    'methyltransferase':'NSP16', 'Spike':'S', 'ORF3a':'ORF3a', 'E':'E', 'M':'M',
                    'ORF6':'ORF6', 'ORF7a':'ORF7a', 'ORF7b':'ORF7b', 'ORF8':'ORF8', 'N':'N', 'ORF10':'ORF10'}


def create_fresh_epitops_dict():
    epits = pd.read_csv('immun_epit_all2_HLA_Wuhan2.csv')
    epitops = {}
    for i, row in epits.iterrows():
        pept_name = prot_translation[row[0]]
        if pept_name not in epitops:
            epitops[pept_name] = {}
        start, end = row[3], row[4]
        epitop, hla = row[1], row[5]
        epitops[pept_name][i] = {}
        epitops[pept_name][i]['start'] = start
        epitops[pept_name][i]['end'] = end
        epitops[pept_name][i]['epi_before'] = epitop
        epitops[pept_name][i]['epi_curr'] = epitop
        epitops[pept_name][i]['epi_after'] = set()
        epitops[pept_name][i]['muts_per_patient'] = 0
        epitops[pept_name][i]['muts_all_time'] = 0
        epitops[pept_name][i]['hla'] = hla
        
    return epitops

# Функция ниже мутирует эпитопы отдельно для каждой временной точки
# словарь epitops передается от временной точки к временной точке (но он обнуляется для новых пациентов)
# если мутация была в день 0, а в день n её нет, то считается как будто она всегда есть. 

# Краткий тутор как этим пользоваться:
# создание объекта: obj = Mutations_to_epitops(file, epitops, status) status = [ancestal, in patient
# to "mutate" ancestal obj: obj.epitops_mutations_ancestral()
# to "mutate" in_patient obj: obj.epitops_mutations_in_patient
# to get датафрейм с мутированными эпитопами: obj.return_result_dataframe()
# to get словарь с мутированными эпитопами: obj.return_epitops_dict()
# to get число мутаций (любых): obj.return_number_of_mutations()
# to get число мутаций, попавших в иммуногенные эпитопы: obj.return_number_of_immunogenic_mutations()


class Mutations_to_epitops:
    
    def __init__(self, filename, epitops, status):
        self.filename = filename
        self.epitops = epitops
        self.status = status
        self.cnt_mutations = 0
        self.cnt_immunogenic_mutations = 0
        self.result_df = None 
        
    # геттер для датафрейма с мутированными эпитопами
    def return_result_dataframe(self):
        return self.result_df
    
    # геттер для словаря с мутированными эпитопами
    def return_epitops_dict(self):
        return self.epitops
    
    # геттер для числа мутаций
    def return_number_of_mutations(self):
        return self.cnt_mutations
    
    # геттер для числа мутаций, попавших в иммуногенные эпитопы
    def return_number_of_immunogenic_mutations(self):
        return self.cnt_immunogenic_mutations
    
    # сеттер для мутирования эпитопов ancestral
    def epitops_mutations_ancestral(self):
        data = self.filename
        epitops = self.epitops
        df = pd.DataFrame({'sample':[], 'status':[],
                             'protein':[], 'mutations':[], 'ref_peptide':[], 'mut_peptide':[], 'Allele Name':[]})

        cnt_mutations = 0
        cnt_immunogenic_mutations = 0

        for i, row in data.iterrows():
            if row['refAA'] != row['qAA'] and row['qAA'] != '.':
                cnt_mutations += 1
                refpos = int(''.join(filter(str.isdigit, row['variant'])))
                #refpos = int(row['variant'][1:-1])
                pept_name = row['protein']
                if pept_name == 'NSP12a' or pept_name == 'NSP12b':
                    pept_name = 'NSP12'

                for epi in epitops[pept_name]:
                    if epitops[pept_name][epi]['start'] <= refpos <= epitops[pept_name][epi]['end']:
                        cnt_immunogenic_mutations += 1
                        change = refpos - epitops[pept_name][epi]['start']
                        mutated_epi = epitops[pept_name][epi]['epi_curr'][:change]+row['qAA']+epitops[pept_name][epi]['epi_curr'][change+1:]
                        epitops[pept_name][epi]['epi_curr'] = mutated_epi
    #                     if epitops[pept_name][i]['muts_per_patient'] == []:
    #                         epitops[pept_name][i]['muts_per_patient'].append(1)
                        epitops[pept_name][epi]['muts_per_patient'] = epitops[pept_name][epi]['muts_per_patient']+1

        for pept_name in epitops:
            for i in epitops[pept_name]:
                if epitops[pept_name][i]['epi_curr'] != epitops[pept_name][i]['epi_before']:

    # если мутация была в день 0, а в день n её нет, то считается как будто она всегда есть

                    # epitops[pept_name][i]['epi_after'].add(epitops[pept_name][i]['epi_curr'])
                    # epitops[pept_name][i]['muts_all_time'] += 1

                    new_line = {'sample':row['sample'], 'status':row['status'], 'protein':pept_name,
                            'mutations':[],
                            'ref_peptide':epitops[pept_name][i]['epi_before'],
                            'mut_peptide':epitops[pept_name][i]['epi_curr'], 
                            'Allele Name':epitops[pept_name][i]['hla']}

                    epitops[pept_name][i]['epi_before'] = epitops[pept_name][i]['epi_curr']
                    epitops[pept_name][i]['muts_per_patient'] = 0

                    df = pd.concat([df, pd.DataFrame.from_records([new_line])], ignore_index=True)
        
        self.result_df = df
        self.epitops = epitops
        self.cnt_mutations = cnt_mutations
        self.cnt_immunogenic_mutations = cnt_immunogenic_mutations

        
    # сеттер для мутирования эпитопов in_patient
    def epitops_mutations_in_patient(self):
        data = self.filename
        epitops = self.epitops
        df = pd.DataFrame({'sample':[], 'status':[],
                             'protein':[], 'mutations':[], 'ref_peptide':[], 'mut_peptide':[], 'Allele Name':[]})

        cnt_mutations = 0
        cnt_immunogenic_mutations = 0

        for i, row in data.iterrows():
            if row['refAA'] != row['qAA'] and row['qAA'] != '.':
                cnt_mutations += 1
                refpos = int(''.join(filter(str.isdigit, row['variant'])))
                pept_name = row['protein']
                if pept_name == 'NSP12a' or pept_name == 'NSP12b':
                    pept_name = 'NSP12'

                for epi in epitops[pept_name]:
                    if epitops[pept_name][epi]['start'] <= refpos <= epitops[pept_name][epi]['end']:
                        cnt_immunogenic_mutations += 1
                        change = refpos - epitops[pept_name][epi]['start']
                        mutated_epi = epitops[pept_name][epi]['epi_curr'][:change]+row['qAA']+epitops[pept_name][epi]['epi_curr'][change+1:]
                        epitops[pept_name][epi]['epi_curr'] = mutated_epi
    #                     if epitops[pept_name][i]['muts_per_patient'] == []:
    #                         epitops[pept_name][i]['muts_per_patient'].append(1)
                        epitops[pept_name][epi]['muts_per_patient'] += 1

        for pept_name in epitops:
            for i in epitops[pept_name]:
                if epitops[pept_name][i]['epi_curr'] != epitops[pept_name][i]['epi_before']:

    # если мутация была в день 0, а в день n её нет, то считается как будто она всегда есть

                    epitops[pept_name][i]['epi_after'].add(epitops[pept_name][i]['epi_curr'])
                    epitops[pept_name][i]['muts_all_time'] += 1

                    new_line = {'sample':row['sample'], 'status':row['status'], 'protein':pept_name,
                            'mutations':[],
                            'ref_peptide':epitops[pept_name][i]['epi_before'],
                            'mut_peptide':epitops[pept_name][i]['epi_curr'], 
                            'Allele Name':epitops[pept_name][i]['hla']}

                    epitops[pept_name][i]['epi_curr'] = epitops[pept_name][i]['epi_before']
                    epitops[pept_name][i]['muts_per_patient'] = 0

                    df = pd.concat([df, pd.DataFrame.from_records([new_line])], ignore_index=True)
                
        self.result_df = df
        self.epitops = epitops
        self.cnt_mutations = cnt_mutations
        self.cnt_immunogenic_mutations = cnt_immunogenic_mutations        





