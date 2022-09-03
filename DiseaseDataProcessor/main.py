import pandas as pd
import numpy as np
import json
import os
import requests
import time
import winsound
from threading import Thread
import math


# Connect VPN first!
data_folder_name='Data'
id_converter_file= 'UMLS_ID_converter.json'

global oxo_api_address, headers
oxo_api_address = "https://www.ebi.ac.uk/spot/oxo/api/search"
headers = {
    'Content-Type': 'application/json',
    'Accept': 'application/json'
}

global unmapped_ids_df
class File:
    def __init__(self,folder,name):
        self.folder=folder
        self.name=name
        self.address=data_folder_name+'/'+folder+'/'+name

def clean_id(id: str):
    id_number = id[id.find(":"):]
    id_type = id[:id.find(":")].upper()
    if id_type in ["MSH", "MESH"]:
        return "MESH" + id_number
    elif id_type in ["NCIT","NCI"]:
        return "NCIT" + id_number
    elif id_type in ["ORPHANET","ORDO"]:
        return "ORPHANET" + id_number
    elif id_type in ["UMLS", "UMLS_CUI"]:
        return "UMLS" + id_number
    elif id_type in ["SNOMEDCT_US", "SNOMEDCT", "SCTID"] or id_type.startswith("SNOMEDCT") or id_type.startswith("SCTID"):
        return "SCTID" + id_number
    elif id_type in ["HPO","HP"]:
        return "HP"+id_number
    elif id_type in ["DO","DOID"]:
        return "DOID"+id_number
    elif id_type in ["UMLS_US", "UMLS"]:
        return "UMLS" + id_number
    elif id_type in ["OMIMPS","OMIM"]:
        return "OMIM"+id_number
    elif id_type in ["ICD9","ICD-9","ICD9CM"]:
        return "ICD9"+id_number
    elif id_type in ["ICD10","ICD-10","ICD10CM"]:
        return "ICD10"+id_number
    elif id_type in ["ICD11","ICD-11","ICD11CM"]:
        return "ICD11"+id_number
    elif id_type in ["ICDO","ICD-O"]:
        return "ICDO"+id_number
    else:
        return id_type + id_number

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
    print("Loading the previous stored database")
    try:
        df= pd.read_csv('disease_results.tsv', '\t')
        df=df.drop('Sum',axis=1)
        df=df.drop('Combinations',axis=1)
        df=df.set_index('DiseaseID')
        return df
    except:
        return pd.DataFrame()
def get_id_converter():
    with open(id_converter_file,'r') as f:
        data=f.read()
    converter=json.loads(data)
    return converter

def my_thread(ids,mappings,start_i):
    global oxo_api_address, headers
    query_limit = 20
    time_out=120
    t0 = time.time()
    for i in range(start_i,min(len(ids),start_i+query_limit*20),query_limit):
        data = '{"ids" : [ "%s" ],"inputSource" : null, "mappingTarget" : ["UMLS"],"mappingSource" : [],"distance" : 1}' % "\" , \"".join(
            ids[i:min(i + query_limit, len(ids))])
        search_state=0
        while search_state!=200:
            search_result = requests.post(oxo_api_address, headers=headers, data=data,timeout=time_out)
            search_result_content = search_result.content.strip()
            search_state = search_result.status_code
        search_result_json = json.loads(search_result_content)
        #check for distance 2 and 3
        # print(search_result_json)
        # exit(0)
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
        # print(list(filter(lambda x:x.is_alive(), all_threads)))
        time.sleep(1)
        continue
    del(all_threads)
    return mappings

def process(file,id_converter,seperator):
    global unmapped_ids_df
    print(f"Reading {file.name} from {file.folder} database")
    df=pd.read_csv(file.address,sep=seperator)
    is_duplicate = df.index.duplicated(keep="first")
    not_duplicate = ~is_duplicate
    df = df[not_duplicate]
    converted_df=pd.DataFrame({'DiseaseID':[],'DiseaseName':[]})
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
    old_databases=[]
    for column in old_df.columns:
        old_databases.append(str(column))
    data_base_modified=False
    counter=0
    for file in all_data_files:
        if file.folder not in old_databases:
            data_base_modified=True
            if counter==3:
                break
            if file.name.endswith('.tsv'):
                counter += 1
                old_df=update_df(process(file,id_converter,seperator='\t'),old_df)
            elif file.name.endswith('.csv'):
                counter += 1
                old_df = update_df(process(file, id_converter, seperator=','), old_df)
    print("Writing the result into result.tsv")
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
        # indexes_except_sum_and_count = [i for i in old_df.index if i != 'Sum' and i != 'DiseaseCounts']
        # old_df = old_df.loc[['CountOfMappedIDs'] + ['CountOfTotalIDs']+indexes_except_sum_and_count ]
        old_df.to_csv("disease_results.tsv",sep='\t',index=True)
        unmapped_ids_df=unmapped_ids_df.sort_index()
        unmapped_ids_df.to_csv("unmapped_ids.tsv",sep="\t")
    print("The result is written in result.tsv")