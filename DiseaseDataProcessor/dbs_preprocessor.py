import time

import pandas as pd
import requests as req
import json
import regex as re
import obonet
import os
from threading import Thread

data_folder='Data/'
def preprocess_BioMuta():
    print('preprocessing BioMuta')
    folder=f'{data_folder}/BioMuta'
    pd.read_csv(f'{folder}/BioMuta_diseases.tsv',sep='\t').drop_duplicates(subset=['DiseaseID']).dropna(subset=['DiseaseID'])\
        .set_index('DiseaseID').sort_index().to_csv(f'{folder}/processed_BioMuta.tsv',sep='\t')
def preprocess_ClinGen():
    print('preprocessing ClinGen')
    folder=f'{data_folder}/ClinGen/'
    pd.read_csv(f'{folder}/Clingen-Curation-Activity-Summary-Report-2022-09-01.csv',sep=',')[['mondo_id','disease_label']]\
        .drop_duplicates(subset=['mondo_id']).dropna(subset=['mondo_id']).rename(columns={'mondo_id':'DiseaseID','disease_label':'DiseaseName'})\
        .set_index('DiseaseID').sort_index().to_csv(f'{folder}/processed_ClinGen.tsv',sep='\t')
def preprocess_ClinVar():
    print('preprocessing ClinVar')
    folder=f'{data_folder}/ClinVar'
    df1=pd.read_csv(f'{folder}/disease_names.txt',sep='\t')[['ConceptID','DiseaseName','SourceName','SourceID']].\
        drop_duplicates(subset=['ConceptID']).dropna(subset=['ConceptID']).rename(columns={'ConceptID':'DiseaseID'})
    df1['DiseaseID']='UMLS:'+df1['DiseaseID']
    df1.set_index('DiseaseID').sort_index().to_csv(f'{folder}/processed_ClinVar.tsv',sep='\t')

    # df2 = pd.read_csv(f'{folder}/gene_condition_source_id.txt', sep='\t')[
    #     ['ConceptID', 'DiseaseName', 'SourceName','SourceID']].drop_duplicates(subset=['ConceptID']) \
    #     .dropna(subset=['ConceptID']).rename(columns={'ConceptID': 'DiseaseID'})
    # df2['DiseaseID'] = 'UMLS:' + df2['DiseaseID']
    # unified_df=pd.concat([df1,df2],axis=0).drop_duplicates(subset='DiseaseID').set_index('DiseaseID').sort_index()
    # df2.set_index('DiseaseID').sort_index().to_csv(f'{folder}/processed_ClinVar2.tsv', sep='\t')
    # unified_df.to_csv(f'{folder}/processed_ClinVar.tsv',sep='\t')

def preprocess_COSMIC():
    print('preprocessing COSMIC')
    folder = f'{data_folder}/COSMIC'
    df = pd.read_csv(f'{folder}/classification.txt', sep=',')[
        ['NCI', 'HISTOLOGY','HISTOLOGY_COSMIC']].drop_duplicates(subset=['NCI']) \
        .dropna(subset=['NCI']).rename(columns={'NCI': 'DiseaseID','HISTOLOGY':'DiseaseName'})
    df['DiseaseID'] = 'NCIT:' + df['DiseaseID']
    df['DiseaseName']=df['DiseaseName']+'/'+df['HISTOLOGY_COSMIC']
    df=df[['DiseaseID','DiseaseName']]
    df.set_index('DiseaseID').sort_index().to_csv(f'{folder}/processed_COSMIC.tsv', sep='\t')
def preprocess_CTDbase():
    print('preprocessing CTDbase')
    folder = f'{data_folder}/CTDbase'
    pd.read_csv(f'{folder}/CTD_diseases.tsv', sep='\t')[
        ['DiseaseID','DiseaseName']].drop_duplicates(subset=['DiseaseID']) \
        .dropna(subset=['DiseaseID']).set_index('DiseaseID').sort_index().to_csv(f'{folder}/processed_CTDbase.tsv',sep='\t')

def preprocess_DISEASES():
    print('preprocessing DISEASES')
    folder=f'{data_folder}/DISEASES'
    # this database is too large to read with pandas
    df=pd.DataFrame({'DiseaseID':[],'DiseaseName':[]})
    disease_ids=[]
    disease_names=[]
    #Columns: 'GeneID','GeneName','DiseaseID','DiseaseName','Z-Score','ConfidenceScore','URL'
    with open (f'{folder}/human_disease_textmining_full.tsv','r') as f:
        for line in f:
            disease_ids.append(line.split('\t')[2])
            disease_names.append(line.split('\t')[3])
    df['DiseaseID']=disease_ids
    df['DiseaseName']=disease_names
    df.drop_duplicates(subset=['DiseaseID']) \
        .dropna(subset=['DiseaseID']).set_index('DiseaseID').sort_index().to_csv(f'{folder}/processed_DISEASES.tsv',
                                                                                 sep='\t')
def preprocess_DisGeNET():
    print('preprocessing DisGeNET')
    folder = f'{data_folder}/DisGeNET/'
    pd.read_csv(f'{folder}/disease_mappings.tsv', sep='\t')[
        ['diseaseId','name','vocabulary','code']].drop_duplicates(subset=['diseaseId']).dropna(subset=['diseaseId'])\
        .rename(columns={'diseaseId':'DiseaseID','name':'DiseaseName'}).set_index('DiseaseID').sort_index()\
        .to_csv(f'{folder}/processed_DisGeNET.tsv',sep='\t')
def preprocess_DOID():
    print('preprocessing DOID')
    folder=f'{data_folder}/DOID'
    graph = obonet.read_obo(f"{folder}/doid.obo.txt")
    disease_ids=[]
    disease_names=[]
    for id_,data in graph.nodes(data=True):
        disease_names.append(data["name"])
        disease_ids.append(id_)
    disease_ids=pd.Series(disease_ids)
    disease_names=pd.Series(disease_names)
    df=pd.DataFrame({"DiseaseID":disease_ids,"DiseaseName":disease_names})
    df.drop_duplicates(subset=['DiseaseID']).dropna(subset=['DiseaseID']).set_index('DiseaseID').sort_index()\
        .to_csv(f'{folder}/processed_DOID.tsv',sep='\t')
def preprocess_EFO():
    print('preprocessing EFO')
    folder = f'{data_folder}/EFO'
    with open(f"{folder}/efo.obo.txt", 'r', encoding='utf-8') as f:
        graph=obonet.read_obo(f)
        disease_ids=[]
        disease_names=[]
        for node in graph.nodes(data=True):
            try:
                name=node[1].get("name")
                id_=node[0]
                id_type=id_[:id_.index(":")]
                if name is not None and id_type in ["EFO","DOID","HP","MONDO","GO","NCIT","Orphanet","UMLS"]:
                    if id_type=="EFO":
                        is_disease = False
                        for property_value in node[1].get("property_value"):
                            if property_value.startswith("exactMatch") or property_value.startswith("closeMatch"):
                                is_disease=True
                                break
                        if is_disease:
                            disease_names.append(name)
                            disease_ids.append(id_)
                    else:
                        disease_names.append(name)
                        disease_ids.append(id_)
            except:
                continue
        df=pd.DataFrame({"DiseaseID":disease_ids, "DiseaseName":disease_names})
        df.drop_duplicates(subset=['DiseaseID']).dropna(subset=['DiseaseID']).set_index('DiseaseID').sort_index()\
            .to_csv(f'{folder}/processed_EFO.tsv',sep='\t')
def preprocess_G2P():
    print('preprocessing G2P')
    folder = f'{data_folder}/G2P'
    df_final=pd.DataFrame()
    for file in os.listdir(folder):
        if file.endswith('.csv'):
            df=pd.read_csv(f"{folder}/{file}")[['disease mim','disease name']]\
                .rename(columns={'disease mim':'DiseaseID','disease name':'DiseaseName'})
            df["DiseaseID"]="OMIM:"+df["DiseaseID"]
            df=df.set_index('DiseaseID')
            df=df.filter(regex="OMIM:[0-9]+",axis=0)
            df=df.reset_index()
            df_final=df_final.append(df)
    df_final.drop_duplicates(subset=['DiseaseID']).dropna(subset=['DiseaseID']).set_index('DiseaseID').sort_index()\
        .to_csv(f'{folder}/processed_G2P.tsv', sep='\t')
def preprocess_GEO():
    print('preprocessing GEO')
    folder = f'{data_folder}/GEO'
    disease_ids = []
    disease_names=[]
    df=pd.DataFrame({"DiseaseID":[],"DiseaseName":[]})
    with open(f'{folder}/Disease_Perturbations_from_GEO_down.txt', 'r') as f:
        for disease in f:
            disease=disease.split('\t')[0]
            # Only considering human data
            if 'human' in disease:
                x=re.search(r"\b[DOCUI-]+[0-9]+",disease)
                mapping=x.group()
                if mapping[0]=="C":
                    if mapping.startswith("CUI"):
                        mapping="UMLS:"+mapping[4:]
                    else:
                        mapping="UMLS:"+mapping
                elif mapping[0]=="D":
                    mapping=mapping.replace("-",":")
                disease_ids.append(mapping)
                disease_names.append(disease[:int(x.regs[0][0])])
    df["DiseaseID"]=pd.Series(disease_ids)
    df["DiseaseName"]=pd.Series(disease_names)
    df.drop_duplicates(subset=['DiseaseID']).dropna(subset=['DiseaseID']).set_index('DiseaseID').sort_index()\
        .to_csv(f'{folder}/processed_GEO.tsv',sep='\t')

def preprocess_GWAS():
    print('preprocessing GWAS')
    folder = f'{data_folder}/GWAS'
    df=pd.read_csv(f'{folder}/gwas-catalog-associations_ontology-annotated.tsv',sep='\t',dtype=str)[['MAPPED_TRAIT_URI','MAPPED_TRAIT']]\
        .rename(columns={'MAPPED_TRAIT':'DiseaseName','MAPPED_TRAIT_URI':'DiseaseID'})
    df['DiseaseID']=df.DiseaseID.str.extract(r'\b/(\w+)_', expand=True)+':'+df.DiseaseID.str.extract(r'_(\d+)$',expand=True)
    df.drop_duplicates(subset=['DiseaseID']).dropna(subset=['DiseaseID']).set_index('DiseaseID').sort_index()\
        .to_csv(f'{folder}/processed_GWAS.tsv',sep='\t')
def preprocess_HPO():
    print('preprocessing HPO')
    folder = f'{data_folder}/HPO'
    pd.read_csv(f'{folder}/hp_kgx_tsv_nodes.tsv.txt',sep='\t')[['DiseaseID','DiseaseName','has_db_xref']]\
        .drop_duplicates(subset=['DiseaseID']).dropna(subset=['DiseaseID']).set_index('DiseaseID')\
        .filter(regex='HP:.*|GO:.*',axis=0).sort_index().to_csv(f'{folder}/processed_HPO.tsv',sep='\t')
def preprocess_KEGG():
    def my_thread(entries:list,disease_pairs:list):
        for entry in entries:
            disease_id=''
            disease_name=''
            print(f"Checking {entry}")
            try:
                request=req.get(f'{api_address}/{entry}')
                # if not request.ok:
                #     print('ERROR')
                #     print(request)
                #     exit()
            except:
                print(f"request failed for {entry}")
                try:
                    request = req.get(f'{api_address}/{entry}')
                except:
                    print(f'request failed for {entry} again!')
                    continue
                # break
            try:
                disease_id='MESH:'+re.findall("[C-D][0-9]+",request.text[request.text.find("MeSH"):])[0]
            except:
                try:
                    disease_id='OMIM:'+re.findall("entry/([0-9]+)",request.text[request.text.find("OMIM:"):])[0]
                except:
                    try:
                        disease_id="ICD10:"+re.findall(">([0-9.]*[A-Z][0-9.]+[A-Z]*)",request.text[request.text.find("ICD-10:"):])[0]
                    except:
                        try:
                            disease_id="ICD11:" + re.findall(">([0-9.]*[A-Z][0-9.]+[A-Z]*)", request.text[request.text.find("ICD-11:"):])[0]
                        except:
                            error_entries.append(entry)
                            continue
            try:
                disease_name=re.findall("DISEASE: ([A-Za-z0-9,\- ]+)",request.text)[0]
                disease_pairs.append((disease_id,disease_name))
            except:
                print("error in name")
                print(entry)

    print('preprocessing KEGG')
    folder = f'{data_folder}/KEGG'
    api_address = 'https://www.genome.jp/entry/'
    disease_pairs=[]
    error_entries=[]
    all_threads = list()
    with open(f'{folder}/entries.txt','r') as f:
        entries=f.read().split(',')
        query_size=300
        for i in range(len(entries)//query_size+1):
            thread=Thread(target=my_thread,args=(entries[i*query_size:(i+1)*query_size],disease_pairs))
            thread.start()
            all_threads.append(thread)
        while len(list(filter(lambda x: x.is_alive(), all_threads))) > 0:
            time.sleep(10)
            print("*****",len(disease_pairs)," entries checked *****")
            continue
        disease_ids,disease_names=list(zip(*disease_pairs))
        pd.DataFrame({'DiseaseID': disease_ids, 'DiseaseName': disease_names}).drop_duplicates(subset='DiseaseID') \
            .dropna(subset=['DiseaseID']).set_index('DiseaseID').sort_index().to_csv(f'{folder}/processed_KEGG.tsv',sep='\t')

def preprocess_MedGen():
    print('preprocessing MedGen')
    folder = f'{data_folder}/MedGen'
    df=pd.DataFrame({'DiseaseID':[],'DiseaseName':[]})
    disease_names=[]
    disease_ids=[]
    with open (f'{folder}/medgen_result.txt','r',encoding='utf-8') as f:
        for line in f:
            if line.startswith('Title'):
                disease_names.append(line[line.find(':')+2:-1])
            elif line.startswith('Concept ID'):
                disease_ids.append('UMLS:'+line[line.find(':')+2:-1])
    df['DiseaseID']=disease_ids
    df['DiseaseName']=disease_names
    df.drop_duplicates(subset=['DiseaseID']).dropna(subset=['DiseaseID']).set_index('DiseaseID').sort_index()\
        .to_csv(f'{folder}/processed_MedGen.tsv',sep='\t')
def preprocess_OMIM():
    print('preprocessing OMIM')
    folder = f'{data_folder}/OMIM'
def preprocess_OncoMX():
    print('preprocessing OncoMX')
    folder = f'{data_folder}/OncoMX'
    pd.read_csv(f'{folder}/OncoMX_diseases.txt',sep='\t').drop_duplicates(subset=['DiseaseID'])\
        .dropna(subset=['DiseaseID']).set_index('DiseaseID').sort_index().to_csv(f'{folder}/processed_OncoMX.tsv',sep='\t')
def preprocess_OpenTargetsPlatform():
    def find_last_slash(string):
        try:
            return string.rindex('/')
        except:
            return -1
    print('preprocessing OpenTargetsPlatform')
    folder = f'{data_folder}/OpenTargetsPlatform'
    diseaseNames=[]
    diseaseIDs=[]
    for file_name in os.listdir(folder):
        if file_name.endswith(".json"):
            print(f'Getting IDS from {file_name}')
            with open(f'{folder}/{file_name}','r',encoding='utf-8') as f:
                for json_obj in f:
                    disease=json.loads(json_obj,)
                    id_=disease['id'].replace('_',':')[find_last_slash(disease['id'])+1:].upper()
                    diseaseIDs.append(id_)
                    diseaseNames.append(disease['name'])

    df=pd.DataFrame({"DiseaseName":diseaseNames,"DiseaseID":diseaseIDs})
    df.drop_duplicates(subset=['DiseaseID']).dropna(subset=['DiseaseID']).set_index('DiseaseID').sort_index()\
        .to_csv(f'{folder}/processed_OpenTargetsPlatform.tsv',sep='\t')
def preprocess_Pharos():
    print('preprocessing Pharos')
    folder = f'{data_folder}/Pharos'
    pd.read_csv(f'{folder}/query results.csv',sep=',')[['Mondo ID','Associated Disease']]\
        .rename(columns={'Mondo ID':'DiseaseID','Associated Disease':'DiseaseName'})\
        .drop_duplicates(subset=['DiseaseID']).dropna(subset=['DiseaseID']).set_index('DiseaseID').sort_index()\
        .to_csv(f'{folder}/processed_Pharos.tsv',sep='\t')
def preprocess_RepoDB():
    print('preprocessing RepoDB')
    folder = f'{data_folder}/RepoDB'
    df=pd.read_csv(f'{folder}/full.csv',sep=',')[['ind_id','ind_name']]\
        .rename(columns={'ind_name':'DiseaseName','ind_id':'DiseaseID'})
    df['DiseaseID']='UMLS:'+df['DiseaseID']
    df.drop_duplicates(subset=['DiseaseID']).dropna(subset=['DiseaseID']).set_index('DiseaseID').sort_index()\
        .to_csv(f'{folder}/processed_RepoDB.tsv',sep='\t')

def preprocess_SIDER():
    print('preprocessing SIDER')
    folder = f'{data_folder}/SIDER'
    df_final=pd.DataFrame({'DiseaseID':[],'DiseaseName':[]})
    df1 = pd.read_csv(f'{folder}/meddra_all_indications.tsv', sep='\t',header=None)
    df2 = pd.read_csv(f'{folder}/meddra_all_se.tsv', sep='\t', header=None)
    df3 = pd.read_csv(f'{folder}/meddra_freq.tsv', sep='\t', header=None)
    df_final=df_final.append(df1[[1,3]].rename(columns={1:'DiseaseID',3:'DiseaseName'}))\
        .append(df1[[5,6]].rename(columns={5:'DiseaseID',6:'DiseaseName'}))\
        .append(df2[[4,5]].rename(columns={4:'DiseaseID',5:'DiseaseName'}))\
        .append(df3[[8,9]].rename(columns={8:'DiseaseID',9:'DiseaseName'}))
    df_final['DiseaseID']='UMLS:'+df_final['DiseaseID']
    df_final.drop_duplicates(subset=['DiseaseID']).dropna(subset=['DiseaseID']).set_index('DiseaseID').sort_index()\
        .to_csv(f'{folder}/processed_SIDER.tsv',sep='\t')
def preprocess_TTD():
    print('preprocessing TTD')
    folder = f'{data_folder}/TTD'
    disease_ids=[]
    disease_names=[]
    with open(f'{folder}/P1-05-Drug_disease.txt',encoding="utf-8") as file:
        for line in file:
            if line.startswith("INDICATI"):
                disease_names.append(line[9:line.find("[")-1])
                disease_ids.append(line[line.find("[")+1:line.find("]")].replace(" ",''))

    with open(f'{folder}/P1-06-Target_disease.txt',encoding="utf-8") as file:
        for line in file:
            if line.__contains__("INDICATI"):
                disease_names.append(line[line.rfind('\t')+1:line.find("[")-1])
                disease_ids.append(line[line.find("[")+1:line.find("]")].replace(" ",''))
    with open(f'{folder}/P1-08-Biomarker_disease.txt',encoding="utf-8") as file:
        for line in file:
            if line.__contains__("INDICATI"):
                disease_names.append(line[line.rfind('\t')+1:line.find("[")-1])
                disease_ids.append(line[line.find("[")+1:line.find("]")].replace(" ",''))
    df=pd.DataFrame({'DiseaseID':disease_ids,'DiseaseName':disease_names})
    df.drop_duplicates(subset=['DiseaseID']).dropna(subset=['DiseaseID']).set_index('DiseaseID').sort_index()\
        .to_csv(f'{folder}/processed_TTD.tsv',sep='\t')



if __name__ == '__main__':
    databases=['BioMuta','ClinGen','ClinVar','COSMIC','CTDbase','DISEASES','DisGeNET','DOID','EFO','G2P',
               'GEO','GWAS','HPO','KEGG','MedGen','OMIM','OncoMX','OpenTargetsPlatform','Pharos','RepoDB','SIDER','TTD']
    databases_preprocess_functions=[preprocess_BioMuta,preprocess_ClinGen,preprocess_ClinVar,preprocess_COSMIC,
                                    preprocess_CTDbase,preprocess_DISEASES,preprocess_DisGeNET,preprocess_DOID,
                                    preprocess_EFO,preprocess_G2P,preprocess_GEO,preprocess_GWAS,preprocess_HPO,
                                    preprocess_KEGG,preprocess_MedGen,preprocess_OMIM,preprocess_OncoMX,preprocess_OpenTargetsPlatform,
                                    preprocess_Pharos,preprocess_RepoDB,preprocess_SIDER,preprocess_TTD]

    while True:
        print('\n'.join([f'{i}.{databases[i-1]}' for i in range(1,len(databases)+1)]))
        num=input('Enter a database number to preprocess:')
        if num.isnumeric() and 0<int(num)< len(databases)+1:
            databases_preprocess_functions[int(num)-1]()
        elif num=='Exit':
            exit()
        else:
            print('Try again!')
        print()


