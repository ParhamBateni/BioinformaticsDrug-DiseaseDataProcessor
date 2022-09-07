import pandas as pd
import numpy as np
import json
import os
import requests
import time
import winsound
from threading import Thread
import math

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


# Connect VPN first!
data_folder_name='Data'
id_converter_file= 'UMLS_ID_converter.json'
id_history_file_address=f'{data_folder_name}/ClinVar/ConceptID_history.txt'

result_file='disease_results.tsv'

global oxo_api_address, headers
oxo_api_address = "https://www.ebi.ac.uk/spot/oxo/api/search"
headers = {
    'Content-Type': 'application/json',
    'Accept': 'application/json'
}

global unmapped_ids_df,ids_history
class File:
    def __init__(self,folder,name):
        self.folder=folder
        self.name=name
        self.address=data_folder_name+'/'+folder+'/'+name

def clean_id(id: str):
    id_number = id[id.find(":") + 1:]
    id_type = id[:id.find(":")].upper()
    if id_type in ["MSH", "MESH"]:
        return "MESH:" + id_number
    elif id_type in ["NCIT", "NCI"]:
        return "NCIT:" + id_number
    elif id_type in ["ORPHANET", "ORDO"] or id_number.upper().startswith('ORPHANET_'):
        if id_number.upper().startswith('ORPHANET_'):
            return 'ORPHANET:' + id_number[9:]
        else:
            return "ORPHANET:" + id_number
    elif id_type in ["UMLS", "UMLS_CUI"]:
        return "UMLS:" + id_number
    elif id_type in ["SNOMEDCT_US", "SNOMEDCT", "SCTID"] or id_type.startswith("SNOMEDCT") or id_type.startswith(
            "SCTID"):
        return "SCTID:" + id_number
    elif id_type in ["HPO", "HP"]:
        return "HP:" + id_number
    elif id_type in ["DO", "DOID"]:
        return "DOID:" + id_number
    elif id_type in ["UMLS_US", "UMLS"]:
        return "UMLS:" + id_number
    elif id_type in ["OMIMPS", "OMIM"] or (id_number.startswith('PS') and id_number[2:].isnumeric()):
        if id_number.startswith('PS'):
            return 'OMIM:' + id_number[2:]
        else:
            return "OMIM:" + id_number
    elif id_type in ["ICD9", "ICD-9", "ICD9CM"]:
        return "ICD9:" + id_number
    elif id_type in ["ICD10", "ICD-10", "ICD10CM"]:
        return "ICD10:" + id_number
    elif id_type in ["ICD11", "ICD-11", "ICD11CM"]:
        return "ICD11:" + id_number
    elif id_type in ["ICDO", "ICD-O"]:
        return "ICDO:" + id_number
    else:
        return id_type + ':' + id_number


def get_data_files():
    files=[]
    folder_names=os.listdir(data_folder_name)
    for folder in folder_names:
        for file in os.listdir(data_folder_name+'/'+folder):
            files.append(File(folder,file))
    return files

def update_df(new_df,old_df):
    #Removing duplicate indexes in the new dataframe
    is_duplicate = new_df.index.duplicated(keep="first")
    not_duplicate = ~is_duplicate
    new_df = new_df[not_duplicate]


    updated_df=pd.concat([old_df,new_df],axis=0)
    updated_df=updated_df.groupby(updated_df.index,sort=False)[updated_df.columns.tolist()].first()

    return updated_df
def get_old_df():
    global unmapped_ids_df
    print("Loading the previous stored database")
    try:
        df= pd.read_csv(result_file, '\t').drop(columns=['Sum','Combinations'],errors='ignore').set_index('DiseaseID')
        unmapped_ids_df=pd.read_csv('unmapped_ids.tsv',sep='\t').set_index('DiseaseID')
        return df
    except:
        return pd.DataFrame()
def get_id_converter():
    with open(id_converter_file,'r') as f:
        data=f.read()
    converter=json.loads(data)
    return converter

def load_ids_history_data():
    global ids_history
    ids_history=dict()
    df = pd.read_csv(id_history_file_address, sep='\t')[['Previous ConceptID', 'Current ConceptID']]\
        .rename(columns={'Previous ConceptID':'PreviousUMLS','Current ConceptID':'CurrentUMLS'})
    df['PreviousUMLS'],df['CurrentUMLS']='UMLS:'+df['PreviousUMLS'],'UMLS:'+df['CurrentUMLS']
    df=df.drop_duplicates(subset=['PreviousUMLS','CurrentUMLS'])
    for index,row in df.iterrows():
        if ids_history.get(row['PreviousUMLS']) is not None:
            if row['CurrentUMLS']=='UMLS:No longer reported':
                ids_history[str(row['PreviousUMLS']).strip()]=[str(row['PreviousUMLS']).strip()+'(No longer reported)']
            else:
                ids_history[str(row['PreviousUMLS']).strip()].append(row['CurrentUMLS'])
        else:
            if row['CurrentUMLS']=='UMLS:No longer reported':
                ids_history[str(row['PreviousUMLS']).strip()]=[str(row['PreviousUMLS']).strip()+'(No longer reported)']
            else:
                ids_history[str(row['PreviousUMLS']).strip()] = [row['CurrentUMLS']]

def my_thread(ids,mappings,start_i):
    global oxo_api_address, headers
    query_limit = 20
    time_out=120
    t0 = time.time()
    for i in range(start_i,min(len(ids),start_i+query_limit*20),query_limit):
        data = '{"ids" : [ "%s" ],"inputSource" : null, "mappingTarget" : ["UMLS"],"mappingSource" : [],"distance" : 1}' % "\" , \"".join(
            ids[i:min(i + query_limit, len(ids))])
        response_state=0
        while response_state!=200:
            response = requests.post(oxo_api_address, headers=headers, data=data, timeout=time_out)
            response_content = response.content.strip()
            response_state = response.status_code
        search_result_json = json.loads(response_content)
        for query in search_result_json['_embedded']['searchResults']:
            mapped_ids = list()
            for mapped_id in query['mappingResponseList']:
                mapped_ids.append(clean_id(mapped_id['curie']))
            mappings[clean_id(query['queryId'])] = mapped_ids
        del (search_result_json)
    print(f"thread number ({start_i//(query_limit*20)}/{len(ids)//(query_limit*20)}) took {time.time() - t0:0.3} seconds to send and get request")
    return mappings
def check_oxo(ids:list):
    global oxo_api_address,headers
    query_limit=20
    mappings=dict()
    all_threads=list()
    for i in range(0,len(ids),query_limit*20):
        thread=Thread(target=my_thread,args=(ids,mappings,i,))
        thread.start()
        all_threads.append(thread)

    while len(list(filter(lambda x:x.is_alive(), all_threads)))>0:
        time.sleep(1)
        continue
    del(all_threads)
    return mappings

def process(file,id_converter,seperator):
    global unmapped_ids_df,ids_history
    print(f"Reading {file.name} from {file.folder} database")
    df=pd.read_csv(file.address,sep=seperator)
    is_duplicate = df.index.duplicated(keep="first")
    not_duplicate = ~is_duplicate
    df = df[not_duplicate]
    converted_df=pd.DataFrame({'DiseaseID':[],'DiseaseName':[]},dtype=str)
    converted_ids=[]
    converted_names=[]
    #Converting IDs
    print('Converting IDs started')
    count_mapped=0
    i=-1
    unmapped_disease_ids=[]
    unmapped_disease_names=[]
    oxo_mappings=check_oxo(list(filter(lambda x:not x.startswith("UMLS"),df['DiseaseID'].tolist())))
    for id_ in df['DiseaseID']:
        i+=1
        if str(id_).startswith("UMLS"):
            count_mapped+=1
            if ids_history.get(id_) is not None:
                for current_id in ids_history[id_]:
                    converted_ids.append(current_id)
                    converted_names.append(df['DiseaseName'].iloc[i])
            else:
                converted_ids.append(str(id_))
                converted_names.append(df['DiseaseName'].iloc[i])
        else:
            id_=clean_id(id_)
            try:
                mapped_ids=list(set(oxo_mappings[id_]+id_converter[id_]))
            except:
                mapped_ids=list(set(oxo_mappings[id_]))
            if mapped_ids !=[]:
                count_mapped+=1
                for xref in mapped_ids:
                    if ids_history.get(xref) is not None:
                        for current_id in ids_history[xref]:
                            converted_ids.append(current_id)
                            converted_names.append(df['DiseaseName'].iloc[i])
                    else:
                        converted_ids.append(xref)
                        converted_names.append(df['DiseaseName'].iloc[i])
            else:
                unmapped_disease_names.append(df['DiseaseName'].iloc[i])
                unmapped_disease_ids.append(id_)
    unmapped_df=pd.DataFrame({"DiseaseID":unmapped_disease_ids,"DiseaseName":unmapped_disease_names})
    unmapped_df=unmapped_df.set_index("DiseaseID")
    unmapped_ids_df=pd.concat([unmapped_ids_df,unmapped_df],axis=0)
    unmapped_ids_df = unmapped_ids_df.groupby(unmapped_ids_df.index, sort=False)[unmapped_ids_df.columns.tolist()].first()
    print("Converting IDs finished")


    converted_df['DiseaseID']=pd.Series(converted_ids)
    converted_df['DiseaseName']=pd.Series(converted_names)

    converted_df=converted_df.set_index('DiseaseID')
    converted_df[file.folder] = 1
    converted_df.loc['CountOfMappedIDs', file.folder] = count_mapped
    converted_df.loc['CountOfTotalIDs', file.folder] = int(len(df.index))
    return converted_df

if __name__ == '__main__':
    global unmapped_ids_df
    unmapped_ids_df = pd.DataFrame({"DiseaseID": [], "DiseaseName": []})
    unmapped_ids_df = unmapped_ids_df.set_index("DiseaseID")
    old_df=get_old_df()
    all_data_files=get_data_files()
    id_converter=get_id_converter()
    load_ids_history_data()
    old_databases=[]
    for column in old_df.columns:
        old_databases.append(str(column))
    data_base_modified=False
    counter=0
    for file in all_data_files:
        if file.folder not in old_databases:
            if file.name.startswith('processed'):
                data_base_modified=True
                if counter==15:
                    break
                if file.name.endswith('.tsv'):
                    counter += 1
                    old_df=update_df(process(file,id_converter,seperator='\t'),old_df)
                elif file.name.endswith('.csv'):
                    counter += 1
                    old_df = update_df(process(file, id_converter, seperator=','), old_df)
    print(f"Writing the result into {result_file}")
    if data_base_modified:
        old_df = old_df.fillna(0)
        old_df['Sum'] = old_df.loc[:, old_df.columns != 'DiseaseName'].sum(axis=1)
        old_df=old_df.sort_index(axis=0)
        i=0
        combinations=[]
        for ind_,row in old_df.iterrows():
            combination=""
            j=0
            for column in old_df.columns:
                if column!="DiseaseName" and column!='Sum':
                    if str(old_df[column].iloc[i])=='1.0':
                        combination+=str(j)+","
                j+=1
            i+=1
            combinations.append(combination[:-1])
        old_df['Combinations'] = combinations
        old_df.to_csv(result_file,sep='\t',index=True)
        unmapped_ids_df=unmapped_ids_df.sort_index()
        unmapped_ids_df.to_csv("unmapped_ids.tsv",sep="\t")
    print(f"The result is written in {result_file}")