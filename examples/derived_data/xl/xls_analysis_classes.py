import pandas as pd
import numpy as np
import math
import itertools
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)
import math
import copy
import re
from collections import OrderedDict 
import Bio
from Bio import SeqIO

def get_unique_links(DB):
    unique_links = set(DB['Linkage Info'])
    return len(unique_links)    

def get_files_names(DB):
    files = [f for f in set((list(DB['Search Filename']))) if not (isinstance(f, float) and math.isnan(float(f)))]
    return files
    
def get_protein_names(DB):
    '''
    '''
    prots_conv = {}
    SA = DB[['Common/Gene Name A','Subunit Name A']]
    for index, row in SA.iterrows():
        name = row['Common/Gene Name A']
        t = row['Subunit Name A']
    
        if name not in prots_conv.keys():
            if isinstance(t, float) and math.isnan(float(t)):
                continue
            else:
                prots_conv[name] = t
    
    SB = DB[['Common/Gene Name B','Subunit Name B']]
    for index, row in SB.iterrows():
        name = row['Common/Gene Name B']
        t = row['Subunit Name B']
    
        if name not in prots_conv.keys():
            if isinstance(t, float) and math.isnan(float(t)):
                continue
            else:
                prots_conv[name] = t
    '''
    elif 'Common/Gene Name' in DB.columns.values:
        SA = DB[['Common/Gene Name','Subunit Name']]
        for index, row in SA.iterrows():
            name = row['Common/Gene Name B']
            t = row['Subunit Name B']
    
            if name not in prots_conv.keys():
                if isinstance(t, float) and math.isnan(float(t)):
                    continue
                else:
                    prots_conv[name] = t
    '''
                
    return prots_conv

def process_interlinks_DB(DB, files, prots_conv, out_name):
    
    names = ['id','Subunit1','Protein_name1','Subunit2','Protein_name2','Residue1','Residue2', \
             'N_count','Redundant_count','Unique_count','Best_Av_score','Best_Min_score','full_name1', \
             'full_name2','id_all'] + files
    files_0 = list(np.zeros(len(files)))
    
    unique_links = set(DB['Linkage Info'])
    print('Number of unique links:', len(unique_links))
    
    DE_unique = pd.DataFrame(columns = names, dtype=object)
    
    id_num = 1
    for u in unique_links:
        dimer = 0
        if isinstance(u, float) and math.isnan(float(u)):
            continue
        else:
            T = DB[DB['Linkage Info'] == u]
            # Protein 1
            prot1_full = T['Common/Gene Name A'].iloc[0]
            prot1_code = T['Subunit A'].iloc[0]
            pp1 = prots_conv[prot1_full]
        
            # Protein 2
            prot2_full = T['Common/Gene Name B'].iloc[0]
            prot2_code = T['Subunit B'].iloc[0]
            pp2 = prots_conv[prot2_full]
       
            # Check that there is a unique XL'ed residue
            XL1 = T['XL A']
            if len(pd.unique(XL1)) != 1:
                print('>WARNING: Non-unique XL for linkage '+ u)
            XL2 = T['XL B']
            if len(pd.unique(XL2)) != 1:
                print('>WARNING: Non-unique XL for linkage '+ u)
        
            #  Check if dimer
            if (pd.unique(XL1) == pd.unique(XL2)) and  (prot1_full ==  prot2_full):
                dimer = 1
        
            # Check score
            if 'Ave Score' in T.columns.values:
                av_score = np.array(T['Ave Score'].values).astype('float')
            elif 'Average Score' in T.columns.values:
                av_score = np.array(T['Average Score'].values).astype('float')
            else:
                v1 = np.array(T['Score A'].values).astype('float')
                v2 = np.array(T['Score B'].values).astype('float')
                av_score = (v1+v2)/2.
                print('No Average Score column')
            
            if 'Min Score' in T.columns.values:
                min_score = np.array(T['Min Score'].values).astype('float')
            elif 'Minimum Score' in T.columns.values:
                min_score = np.array(T['Minimum Score'].values).astype('float')
            else:
                v1 = np.array(T['Score A'].values).astype('float')
                v2 = np.array(T['Score B'].values).astype('float')
                min_score = np.min([v1,v2])
                print('No Minimum Score column')

            
            if av_score.max()>=1.0 and min_score.max()>=1:
        
                # Gather data for new dataframe
                XL1 = pd.unique(XL1)[0] 
                XL2 = pd.unique(XL2)[0] 
                
                temp = T[['Sequence A','Sequence B']].drop_duplicates()
                redundant_count = len(T)
                unique_count = len(temp)
                
                N = pd.DataFrame([[id_num, prot1_code, pp1, prot2_code, pp2, XL1, XL2, len(T), redundant_count, unique_count,av_score.max(),min_score.max(),prot1_full,prot2_full,id_num]+files_0],columns=names, dtype=object)
                DE_unique = DE_unique.append(N, ignore_index=True)
                # Add Id to original DB
                for k in T.index:
                    file_id = DB.loc[k,'Search Filename']
                    DE_unique.loc[(DE_unique['id_all']==id_num),file_id] = DE_unique.loc[(DE_unique['id_all']==id_num),file_id]+1
                    DB.loc[k,'UniqueID'] = id_num
                id_num += 1
            else:
                for k in T.index:
                    DB.loc[k,'UniqueID'] = np.nan
    print('Number of unique links, number in DE_unique ', len(unique_links), len(DE_unique))
    DE_unique.to_csv(out_name+'.csv',index=False) 
    DB.to_csv(out_name+'_with_ID.csv',index=False)
    return DE_unique

# Intralinks
def process_intralinks_DB(IN, out):
    unique_links = pd.unique(IN['Linkage Info'])

    names = ['id','Subunit Name', 'Protein Name', 'Protein Code', 'XL', 'N_count','Redundant_count','Unique_count','Score']
    DB_intralinks_unique = pd.DataFrame(columns = names, dtype=object)

    id_num = 1
    for u in unique_links:
        if isinstance(u, float) and math.isnan(float(u)):
            continue
        else:
            T = IN[IN['Linkage Info'] == u]
            # Protein
            protein_full = T['Common/Gene Name'].iloc[0]
            subunit_name = T['Subunit Name'].iloc[0]
            if subunit_name in prot_names.keys():
                protein_name = prot_names[subunit_name]
            else:
                protein_name = subunit_name
        
            # Check that there is a unique XL'ed residue
            XL = T['XL']
            if len(pd.unique(XL)) != 1:
                print('>WARNING: Non-unique XL for linkage '+ u)
            XL = pd.unique(XL)[0] 
        
            # Check score
            score = T['Score']
            if score.max()>=1.0:    
                temp = T['Sequence'].drop_duplicates()
                redundant_count = len(T)
                unique_count = len(temp) 
                N = pd.DataFrame([[id_num, subunit_name, protein_name, protein_full, XL, len(T), redundant_count,unique_count,score.max()]],columns=names)
                DB_intralinks_unique = DB_intralinks_unique.append(N, ignore_index=True)
                id_num += 1
                for k in T.index:
                    IN.loc[k,'UniqueID'] = id_num
            else:
                for k in T.index:
                    IN.loc[k,'UniqueID'] = 'NaN'
    
    DB_intralinks_unique.to_csv(out+'.csv', index=False)
    return DB_intralinks_unique 

def process_deadends_DB(DE, prot_names, out):
    unique_links = pd.unique(DE['Linkage Info'])

    names = ['id','Subunit Name', 'Protein Name', 'Protein Code', 'XL', 'N_count','Redundant_count','Unique_count','Score']
    DE_deadends_unique = pd.DataFrame(columns = names, dtype=object)

    id_num = 1
    for u in unique_links:
        if isinstance(u, float) and math.isnan(float(u)):
            continue
        else:
            T = DE[DE['Linkage Info'] == u]
            # Protein
            protein_full = T['Common/Gene Name'].iloc[0]
            subunit_name = T['Subunit Name'].iloc[0]
            if subunit_name in prot_names.keys():
                protein_name = prot_names[subunit_name]
            else:
                protein_name = subunit_name
        
            # Check that there is a unique XL'ed residue
            XL = T['XL']
            if len(pd.unique(XL)) != 1:
                print('>WARNING: Non-unique XL for linkage '+ u)
            XL = pd.unique(XL)[0] 
        
            # Check score
            score = T['Score']
            if score.max()>=1.0:    
                temp = T['Sequence'].drop_duplicates()
                redundant_count = len(T)
                unique_count = len(temp) 
                N = pd.DataFrame([[id_num, subunit_name, protein_name, protein_full, XL, len(T), redundant_count,unique_count,score.max()]],columns=names)
                DE_deadends_unique = DE_deadends_unique.append(N, ignore_index=True)
                id_num += 1
                for k in T.index:
                    DE.loc[k,'UniqueID'] = id_num
            else:
                for k in T.index:
                    DE.loc[k,'UniqueID'] = 'NaN'
    
    DE_deadends_unique.to_csv(out+'.csv', index=False)
    return DE_deadends_unique    

def remove_ambiguity(DE_unique, out_name):
    '''
    Check ambiguos XLs (;&| separator)
    '''
    
    names2 = ['id','Subunit1','Protein_name1','Subunit2','Protein_name2','Residue1','Residue2','N_count','Unique_count','Best_Av_score','Best_Min_score','full_name1','full_name2','id_all']

    DE_unique_nonambiguos = pd.DataFrame(columns = names2, dtype=object)

    added_rows = []
    amb = 0
    not_amb = 0
    for index, row in DE_unique.iterrows():
        p1 = row['Protein_name1']
        p2 = row['Protein_name2']
        r1 = row['Residue1']
        r2 = row['Residue2']
        if ';' in r1 and '|' not in r1:
            r1s = r1.split(';')
        elif '|' in r1 and ';' not in r1:
            r1s = r1.split('|')
        elif '|' in r1 and ';' in r1:
            r1 = r1.replace('|',';')
            r1s = r1.split(';')
        else:
            r1s = [r1]
            
        if ';' in r2 and '|' not in r2:
            r2s = r2.split(';')
        elif '|' in r2 and ';' not in r2:
            r2s = r2.split('|')
        elif '|' in r2 and ';' in r2:
            r2 = r2.replace('|',';')
            r2s = r2.split(';')
        else:
            r2s = [r2]
        
        if len(r1s)>1 or len(r2s)>1:
            amb += 1
            T = DE_unique[(DE_unique['Protein_name1']==p1) & 
                              (DE_unique['Residue1'].isin(r1s)) &
                              (DE_unique['Protein_name2']==p2) & 
                              (DE_unique['Residue2'].isin(r2s))]
            sum_counts = T['N_count'].sum()
            
            if len(T)> 0:
                for temp_index, temp_row in T.iterrows(): 
                    frac = temp_row['N_count']/float(sum_counts)
                    to_add = round(row.N_count*frac)
                    if to_add > 0:
                        id_sel = temp_row['id']
                        if id_sel in added_rows:
                            continue
                        else:
                            added_rows.append(id_sel)
                            sel_DE = DE_unique_nonambiguos[DE_unique_nonambiguos['id_all'].isin([id_sel])]
                            if len(sel_DE)==0:
                                row_sel = DE_unique[DE_unique['id']==id_sel]
                                DE_unique_nonambiguos = DE_unique_nonambiguos.append(row_sel, ignore_index=True)
                            ind_o = DE_unique_nonambiguos[DE_unique_nonambiguos['id']==id_sel].index
                            av_score = np.max([DE_unique_nonambiguos['Best_Av_score'].loc[ind_o].values[0],row['Best_Av_score']])
                            min_score = np.max([DE_unique_nonambiguos['Best_Min_score'].loc[ind_o].values[0],row['Best_Min_score']])
                    
                            DE_unique_nonambiguos.loc[DE_unique_nonambiguos['id']==id_sel, 'Best_Av_score'] = np.max([av_score, row['Best_Av_score']])
                            DE_unique_nonambiguos.loc[DE_unique_nonambiguos['id']==id_sel, 'Best_Min_score'] = np.max([min_score, row['Best_Min_score']])
                    
                            DE_unique_nonambiguos.loc[DE_unique_nonambiguos['id']==id_sel,'N_count'] += round(row.N_count*frac)
                            DE_unique_nonambiguos.loc[DE_unique_nonambiguos['id']==id_sel,'Redundant_count'] += round(row.Redundant_count*frac) 
                            #DE_unique_nonambiguos.loc[DE_unique_nonambiguos['id']==id_sel,'Updated N_count'] += round(row.N_count*frac)
                            DE_unique_nonambiguos.loc[DE_unique_nonambiguos['id']==id_sel,'id_all'] = str(id_sel)+'|' + str(row['id'])
                    
            else:
                DE_unique_nonambiguos = DE_unique_nonambiguos.append(row, ignore_index=True)
            
        else:
            not_amb += 1
            id_sel = row['id']
            
            if id_sel not in DE_unique_nonambiguos['id'].values:
                # Check for flipped XL
                V = DE_unique_nonambiguos[(DE_unique_nonambiguos['Protein_name1']==p2) & 
                                         (DE_unique_nonambiguos['Residue1'].isin(r2s)) &
                                         (DE_unique_nonambiguos['Protein_name2']==p1) & 
                                         (DE_unique_nonambiguos['Residue2'].isin(r1s))]
                
                
                if len(V)> 0:

                    id_sel = V['id'].values[0]
                    ind_o = DE_unique_nonambiguos[DE_unique_nonambiguos['id']==id_sel].index
                    av_score = np.max(V['Best_Av_score'].values)
                    min_score = np.max(V['Best_Min_score'].values)
  
    
                    DE_unique_nonambiguos.loc[DE_unique_nonambiguos['id']==id_sel, 'Best_Av_score'] = np.max([av_score, row['Best_Av_score']])
                    DE_unique_nonambiguos.loc[DE_unique_nonambiguos['id']==id_sel, 'Best_Min_score'] = np.max([min_score, row['Best_Min_score']])
                    
                    DE_unique_nonambiguos.loc[DE_unique_nonambiguos['id']==id_sel,'N_count'] += round(row.N_count)
                    DE_unique_nonambiguos.loc[DE_unique_nonambiguos['id']==id_sel,'Redundant_count'] += round(row.Redundant_count) 
                    DE_unique_nonambiguos.loc[DE_unique_nonambiguos['id']==id_sel,'Unique_count'] += round(row.Unique_count)
                    DE_unique_nonambiguos.loc[DE_unique_nonambiguos['id']==id_sel,'id_all'] = str(id_sel)+'|' + str(row['id'])
                    
                    
                else:
                    added_rows.append(id_sel)
                    # Add to DE_unique_nonambiguos
                    DE_unique_nonambiguos = DE_unique_nonambiguos.append(row, ignore_index=True)
            
    print('Ambiguity', amb, not_amb)
    DE_unique_nonambiguos.to_csv(out_name+'_nonambiguos.csv',index=False)   
    return DE_unique_nonambiguos

    
def XLs_for_modeling(DE_unique_nonambiguos, out_name):
    # File for Integrative modeling
    names3 = ['Protein1','Protein2','AbsPos1','AbsPos2','Id','Num','Score']
    names4=['Protein1','AbsPos1','Protein2','AbsPos2']
    
    DE_unique_mod = pd.DataFrame(columns = names3, dtype=object)
    DE_unique_circos = pd.DataFrame(columns = ['Protein1','AbsPos1','Protein2','AbsPos2'],dtype=object)
    for index, row in DE_unique_nonambiguos.iterrows():
        p1 = row['Protein_name1']
        p2 = row['Protein_name2']
        r1 = row['Residue1']
        r2 = row['Residue2']
        
        idn = row['id']
        num = row['N_count']
        sc = row['Best_Min_score']
        if sc <15:
            if num <15:
                score = 0.5
            else:
                score = 0.1
        if sc >= 15 and sc<25:
            if num <15:
                score = 0.1
            else:
                score = 0.01
        if sc >=25:
            score = 0.01
        
        # All proteins
        p1_name = p1
        p2_name = p2
        # Check for Nterm linkage
        if 'N-Term' in str(r1):
            r1 = str(1)
        if 'N-Term' in str(r2):
            r2 = str(1)
            
        if ('|' not in r1) and (';' not in r1):
            r1 = r1.lstrip('K')
        elif (';' in r1):
            rs1 = r1.split(';')
            if len(rs1) >= 3:
                r1 = rs1[1].lstrip('K')
            else:
                r1 = rs1[0].lstrip('K')
        elif ('|' in r1):
            rs1 = r1.split('|')
            if len(rs1) == 3:
                r1 = rs1[1].lstrip('K')
            else:
                r1 = rs1[0].lstrip('K')
            
        if ('|' not in r2) and (';' not in r2):
            r2 = r2.lstrip('K')
        elif (';' in r2):
            rs2 = r2.split(';')
            if len(rs2) >= 3:
                r2 = rs2[1].lstrip('K')
            else:
                r2 = rs2[0].lstrip('K')
        elif ('|' in r2):
            rs2 = r2.split('|')
            if len(rs2) >= 3:
                r2 = rs2[1].lstrip('K')
            else:
                r2 = rs2[0].lstrip('K')  
        N = pd.DataFrame([[p1_name,p2_name,r1, r2, idn,num, score]],columns=names3)
        DE_unique_mod = DE_unique_mod.append(N, ignore_index=True)
        
        N_circos = pd.DataFrame([[p1_name,r1,p2_name,r2]],columns=names4)
        DE_unique_circos = DE_unique_circos.append(N_circos, ignore_index=True)
    DE_unique_mod.to_csv(out_name+'_unique_modeling.csv',index=False)    
    DE_unique_circos.to_csv(out_name+'_unique_circos.csv',index=False) 
    
    return DE_unique_mod

class get_sequences_coverage(object):
    def __init__(self,
                 dataset,
                 sequences,
                 pmids,
                 prots = ['Subunit A','Subunit B'],
                 frags = ['Sequence A','Sequence B']):

        self.dataset = dataset
        self.sequences = sequences
        self.pmids = pmids
        self.prot_fields = prots
        self.frag_fields = frags
        
        # Run
        self.build_coverage_dict()
        self.get_all_matches()
        
    def get_all_matches(self):
        for index,row in self.dataset.iterrows():
            prots = row[self.prot_fields].values
            frags = row[self.frag_fields].values
            # Match sequences 
            if all(elem in self.pmids.keys() for elem in prots):
                sequences = [self.sequences[self.pmids[p]] for p in prots]
                for k, p in enumerate(prots):
                    for i, s in enumerate(sequences):
                        m = self.match_seq(frags[i], s)
                        if len(m)>0:
                            self.coverage_sequences[self.pmids[p]][m[0]:(m[0]+len(frags[i]))] += 1
                            continue
            else:
                print('Missing proteins:', prots)
        
    def build_coverage_dict(self):
        self.coverage_sequences = {}
        for k, v in self.sequences.items():
            self.coverage_sequences[k] = np.zeros((len(v),))

    def match_seq(self, fragment, sequence):
        match = [m.start() for m in re.finditer(fragment, str(sequence))]
        return match
    
    def coverage_output(self):
        return self.coverage_sequences

'''
class get_occupancy_XLs(object):
    def __init__(self,
                 df,
                 sequences):
        self.sequences = sequences

        self.build_occupancy_dict()
        self.occupancy_residues()

    def occupacy_residues(self):
        for i,row in D1_unique.iterrows():
            protein1 = row['Protein_name1']
            protein2 = row['Protein_name2']
            XL1 = row['Residue1']
            XL2 = row['Residue2']
            N_count = row['N_count']
        
        

    def build_occupancy_dict(self):
        self.occupancy_XLs = {}
        for k, v in self.sequences.items():
            self.occupancy_XLs[k] = np.zeros((len(v),))
'''

def process_XLinkX(DB):
    DB_grouped = DB.groupby(['Protein Accession A',
                             'Protein Accession B',
                             'Leading Protein Position A',
                             'Leading Protein Position B']).agg({'XlinkX Score':['mean', 'min', 'max', 'size']})
    DB_grouped.columns = ['XlinkX_Score_mean', 'XlinkX_Score_min', 'XlinkX_Score_max','size']
    DB_grouped = DB_grouped.reset_index()

    for index, row in DB_grouped.iterrows():
        temp = DB_grouped[(DB_grouped['Protein Accession A']==row['Protein Accession B']) &
                          (DB_grouped['Protein Accession B']==row['Protein Accession A']) &
                          (DB_grouped['Leading Protein Position A']==row['Leading Protein Position B']) &
                          (DB_grouped['Leading Protein Position B']==row['Leading Protein Position A'])]
        if len(temp)>0:
            DB_grouped.loc[index,'size'] = row['size'] + temp['size'].values 
            t_mean = (row['size']*row['XlinkX_Score_mean'] + temp['size']*temp['XlinkX_Score_mean'])/(row['size']+temp['size'])
            DB_grouped.loc[index,'XlinkX_Score_mean'] = t_mean.values
            DB_grouped.loc[index,'XlinkX_Score_min'] = np.min([row['XlinkX_Score_min'],temp['XlinkX_Score_min'].values])
            DB_grouped.loc[index,'XlinkX_Score_max'] = np.min([row['XlinkX_Score_max'],temp['XlinkX_Score_max'].values])
            DB_grouped = DB_grouped.drop(temp.index)

    return DB_grouped


def order_DF(DB):
    for index, row in DB.iterrows():
        if row['Protein Accession A']==row['Protein Accession B']:
            oo = np.argsort(row[['Leading Protein Position A','Leading Protein Position B']]) 
            vals_r = list(row[['Leading Protein Position A','Leading Protein Position B']][oo])
            DB.loc[index,['Leading Protein Position A','Leading Protein Position B']] = vals_r
        else:
            op = np.argsort(row[['Protein Accession A','Protein Accession B']])
            vals_p = list(row[['Protein Accession A','Protein Accession B']][op])
            vals_r = list(row[['Leading Protein Position A','Leading Protein Position B']][op])
            DB.loc[index,['Protein Accession A','Protein Accession B']] = vals_p
            DB.loc[index,['Leading Protein Position A','Leading Protein Position B']] = vals_r

    return DB

def get_occupancy_XLs(DB, seqs, pmids, sel_field=['XlinkX, BS3'], count_field='size_BS3'):
    occupancy_XLs = {}
    for k, v in seqs.items():
        occupancy_XLs[k] = np.zeros((len(v),))

    df_XLs = DB[(DB[sel_field] == 1).all(1)]
    
    for prot, v in seqs.items():
        for field in ['A','B']:
            prot_pmid = [k for k,v in pmids.items() if v==prot][0]
            sel = df_XLs[df_XLs[f'Protein Accession {field}']==prot_pmid]
            for i, row in sel.iterrows():
                resi = int(row[f'Leading Protein Position {field}'])-1
                count = row[count_field]
                if count>1:
                    #print('---',count, row)
                    count = math.log2(count)
                    #print(count)
                occupancy_XLs[prot][resi] += count
    return occupancy_XLs

def get_occupancy_DEs(DB, seqs, pmids, sel_field='XlinkX, BS3', count_field='size_BS3'):
    occupancy_DEs = {}
    for k, v in seqs.items():
        occupancy_DEs[k] = np.zeros((len(v),))
    
    df_XLs = DB[DB[sel_field]==1]
    
    for prot, v in seqs.items():
        prot_pmid = [k for k,v in pmids.items() if v==prot][0]
        sel = df_XLs[df_XLs[f'Protein Accession']==prot_pmid]
        for i, row in sel.iterrows():
            resi = int(row[f'Leading Protein Position'])-1
            count = row[count_field]
            if count>1:
                count = math.log2(count)
            occupancy_DEs[prot][resi] += count
    return occupancy_DEs

def fragment_occupancy_pdb(prot, pdb_in, pdb_out, occupancies, chains):
    out = open(pdb_out,'w')
    
    for line in open(pdb_in):
        if line[0:4]=='ATOM' and line[21] in chains.keys():
            vals = line.split()
            prot = chains[vals[4]]
            resi = int(vals[5])
            o = occupancies[prot][resi-1]
            out.write(f'{line[0:56]} round(o,2) {line[61:-1]}')
    out.close()
            
    



    
'''
if isinstance(row['Subunit Name A'], str):
        try: 
            prot1 = prot_names[row['Subunit Name A']]
            prot2 = prot_names[row['Subunit Name B']]
        except:
            continue
        seq1 = row['Sequence A']
        seq2 = row['Sequence B']
        XL1 = row['XL A']
        XL2 = row['XL B']
        
        matches1 = [m.start() for m in re.finditer(seq1, str(seqs[prot1]))] 
        matches2 = [m.start() for m in re.finditer(seq2, str(seqs[prot2]))] 
'''
