import os
import pandas as pd
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
    elif id_type in ["SNOMEDCT_US", "SNOMEDCT", "SCTID"] or id_type.startswith("SNOMEDCT"):
        return "SCTID" + id_number
    elif id_type in ["HPO","HP"]:
        return "HP"+id_number
    elif id_type in ["DO","DOID"]:
        return "DOID"+id_number
    elif id_type in ["UMLS_US", "UMLS"]:
        return "UMLS" + id_number
    elif id_type in ["OMIMPS","OMIM"]:
        return "OMIM"+id_number
    else:
        return id_type + id_number
if __name__ == '__main__':
    ID_stat_df=pd.DataFrame({'database':[], 'IDs':[]})
    for folder in os.listdir('../Data'):
        for file in os.listdir(f'../Data/{folder}'):
            if file.endswith('.tsv'):
                ID_types=dict()
                df=pd.read_csv(f'../Data/{folder}/{file}',sep='\t')
                for index,row in df.iterrows():
                    try:
                        cleaned_ID=clean_id(row['DiseaseID'])
                        ID_type=cleaned_ID[:cleaned_ID.index(":")]
                        if ID_types.get(ID_type) is None:
                            ID_types[ID_type]=1
                        else: ID_types[ID_type]+=1
                    except:
                        continue
                print(f"{folder}:")
                ID_types_listed=list(ID_types)
                ID_types_listed.sort(key=lambda x:int(ID_types[x]),reverse=True)

                print({id_:ID_types[id_] for id_ in ID_types_listed})
                print()


