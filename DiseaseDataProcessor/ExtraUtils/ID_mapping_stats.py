import pandas as pd
import obonet
import os
import json

mapping_stats = dict()


def clean_id(id: str):
    id_number = id[id.find(":")+1:]
    id_type = id[:id.find(":")].upper()
    if id_type in ["MSH", "MESH"]:
        return "MESH:" + id_number
    elif id_type in ["NCIT", "NCI"]:
        return "NCIT:" + id_number
    elif id_type in ["ORPHANET", "ORDO"] or id_number.upper().startswith('ORPHANET_'):
        if id_number.upper().startswith('ORPHANET_'):return 'ORPHANET:'+id_number[9:]
        else:return "ORPHANET:" + id_number
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
            return 'OMIM:'+id_number[2:]
        else: return "OMIM:" + id_number
    elif id_type in ["ICD9", "ICD-9", "ICD9CM"]:
        return "ICD9:" + id_number
    elif id_type in ["ICD10", "ICD-10", "ICD10CM"]:
        return "ICD10:" + id_number
    elif id_type in ["ICD11", "ICD-11", "ICD11CM"]:
        return "ICD11:" + id_number
    elif id_type in ["ICDO", "ICD-O"]:
        return "ICDO:" + id_number
    else:
        return id_type + ':'+id_number

def check_DisGeNET():
    print('Checking DisGeNET')
    df = pd.read_csv("../Data/DisGeNET/disease_mappings.tsv", sep='\t')
    mapping_stat = dict()
    for ind, row in df.iterrows():
        cleaned_id = clean_id(row['vocabulary'] + ":" + row['code'])
        try:
            cleaned_id_type = cleaned_id[:cleaned_id.index(":")]
        except:
            continue
        if mapping_stat.get(cleaned_id_type) is not None:
            mapping_stat[cleaned_id_type] += 1
        else:
            mapping_stat[cleaned_id_type] = 1
    mapping_stats["DisGeNET"] = mapping_stat


def check_CLINVAR():
    print('Checking CLINVAR')
    mapping_stat = dict()
    df = pd.read_csv("../Data/ClinVar/processed_ClinVar.tsv", sep='\t', dtype=str).dropna(subset=['SourceID'])
    for ind, row in df.iterrows():
        if str(row['SourceID']).startswith('ORPHA') and row['SourceName'] == 'Orphanet':
            disease_id = 'ORPHANET:' + disease_id[5:]
        elif str(row['SourceID']).isnumeric() and row['SourceName'] == 'OMIM':
            disease_id = 'OMIM:' + row['SourceID']
        elif str(row['SourceName']) == 'OMIM phenotypic series':
            disease_id = 'OMIM:' + row['SourceID'][2:]
        elif str(row['SourceName'])=='NCBI curation':disease_id='NCBI:'+row['SourceID']
        elif str(row['SourceName']) == 'PharmGKB':
            disease_id = 'PharmGKB:' + str(row['SourceID'])
        elif str(row['SourceName'])=='GeneReviews':disease_id='GeneReviews:'+row['SourceID']
        else:
            disease_id = str(row['SourceID'])
        cleaned_id = clean_id(disease_id)
        try:
            cleaned_id_type = cleaned_id[:cleaned_id.index(":")]
        except:
            continue
        if mapping_stat.get(cleaned_id_type) is not None:
            mapping_stat[cleaned_id_type] += 1
        else:
            mapping_stat[cleaned_id_type] = 1
    mapping_stats["CLINVAR"] = mapping_stat


def check_DOID():
    print('Checking DOID')
    graph = obonet.read_obo("../Data/DOID/doid.obo.txt")
    mapping_stat = dict()
    for id_, data in graph.nodes(data=True):
        try:
            reference_id_xrefs = list(filter(lambda x: x.startswith("UMLS"), data["xref"]))
        except:continue
        other_ids = list(filter(lambda x: not x.startswith("UMLS"), data["xref"]))
        other_ids.append(id_)
        for other_id in other_ids:
            for reference_id in reference_id_xrefs:
                other_id = clean_id(other_id)
                try:
                    other_id_type = other_id[:other_id.index(":")]
                except:
                    continue
                if mapping_stat.get(other_id_type) is None:
                    mapping_stat[other_id_type] = 1
                else:
                    mapping_stat[other_id_type] += 1
    mapping_stats["DOID"] = mapping_stat


def check_EFO():
    print('Checking EFO')
    with open("../Data/EFO/efo.obo.txt", 'r', encoding='utf-8') as f:
        graph = obonet.read_obo(f)
        mapping_stat = dict()
        for id_, data in graph.nodes(data=True):
            try:
                reference_id_xrefs = list(filter(lambda x: x.startswith("UMLS"), data["xref"]))
            except:continue
            other_ids = list(filter(lambda x: not x.startswith("UMLS"), data["xref"]))
            other_ids.append(id_)
            for other_id in other_ids:
                for reference_id in reference_id_xrefs:
                    other_id = clean_id(other_id)
                    try:
                        other_id_type = other_id[:other_id.index(":")]
                    except:
                        continue
                    if mapping_stat.get(other_id_type) is None:
                        mapping_stat[other_id_type] = 1
                    else:
                        mapping_stat[other_id_type] += 1
    mapping_stats["EFO"] = mapping_stat
    pass


def check_HPO():
    df = pd.read_csv("../Data/HPO/hp_diseases_results.tsv", sep='\t').dropna(subset=['has_db_xref'])
    mapping_stat = dict()
    print('Checking HPO')
    for ind, row in df.iterrows():
        reference_id_xrefs = list(filter(lambda x: x.startswith("UMLS"), row["has_db_xref"].split("|")))
        other_ids = list(filter(lambda x: not x.startswith("UMLS"), row["has_db_xref"].split("|")))
        other_ids.append(row["DiseaseID"])
        for other_id in other_ids:
            for reference_id in reference_id_xrefs:
                other_id = clean_id(other_id)
                try:
                    other_id_type = other_id[:other_id.index(":")]
                except:
                    continue
                if mapping_stat.get(other_id_type) is None:
                    mapping_stat[other_id_type] = 1
                else:
                    mapping_stat[other_id_type] += 1

    mapping_stats["HPO"] = mapping_stat
    pass


def check_OTP():
    print('Checking OTP')
    disease_files_name = os.listdir("../Data/OpenTargetsPlatform")
    mapping_stat = dict()
    for file_name in disease_files_name:
        if file_name.endswith("json"):
            print(f'Checking {file_name}')
            with open(f"../Data/OpenTargetsPlatform/{file_name}", 'r', encoding='utf-8') as f:
                for json_obj in f:
                    disease = json.loads(json_obj, )
                    ids = [disease['id'].replace('_', ':')] + disease['dbXRefs']
                    disease_ref_ids = list(filter(lambda x: x.startswith("UMLS"), ids))
                    other_ids = list(filter(lambda x: not x.startswith("UMLS"), ids))
                    for other_id in other_ids:
                        for ref_id in disease_ref_ids:
                            other_id = clean_id(other_id)
                            try:
                                other_id_type = other_id[:other_id.index(":")]
                            except:
                                continue
                            if mapping_stat.get(other_id_type) is None:
                                mapping_stat[other_id_type] = 1
                            else:
                                mapping_stat[other_id_type] += 1
    mapping_stats["OTP"] = mapping_stat


def check_MedGen():
    print('Checking MedGen')
    mapping_stat=dict()
    df = pd.read_csv('../Data/MedGen/MedGenIDMappings.txt', sep='|')
    for index, row in df.iterrows():
        if row['source'] == 'MONDO':
            disease_id = row['source_id']
        elif row['source'] == 'HPO':
            disease_id = row['source_id']
        elif row['source'] == 'OMIM Phenotypic Series':
            disease_id = 'OMIM:' + row['source_id'][2:]
        elif row['source'] == 'OMIM included' or row['source'] == 'OMIM allelic variant':
            disease_id = 'OMIM:' + row['source_id']
        else:
            disease_id = row['source'] + ':' + row['source_id']
        try:
            disease_id_type = disease_id[:disease_id.find(':')]
        except:
            continue
        if mapping_stat.get(disease_id_type) is None:
            mapping_stat[disease_id_type] = 1
        else:
            mapping_stat[disease_id_type] += 1



def command1():
    print(1)
    check_DisGeNET()
    print(2)
    check_OTP()
    print(3)
    check_CLINVAR()
    print(4)
    check_HPO()
    print(5)
    check_EFO()
    print(6)
    check_DOID()
    print(7)
    check_MedGen()

    for key in mapping_stats:
        mapping_stats_listed = list(mapping_stats[key])
        mapping_stats_listed.sort(key=lambda x: int(mapping_stats[key][x]), reverse=True)
        print(f"{key}:{dict({id_: mapping_stats[key][id_] for id_ in mapping_stats_listed})}")
def command2():
    mapping_stats = dict()
    with open("../UMLS_ID_converter_new.json", "r") as f:
        data = f.read()
        UMLS_id_converter = json.loads(data)
    for id_ in UMLS_id_converter:
        try:
            id_ = clean_id(id_)
            id_type=id_[:id_.find(":")]
            if mapping_stats.get(id_type) is None:
                mapping_stats[id_type] = 1
            else:
                mapping_stats[id_type] += 1
        except:
            continue
    print("Count of each ID type mapped in UMLS ID converter:")
    mapping_stats_listed = list(mapping_stats)
    mapping_stats_listed.sort(key=lambda x: int(mapping_stats[x]), reverse=True)
    print(dict({id_: mapping_stats[id_] for id_ in mapping_stats_listed}))
def command3():
    df=pd.read_csv('../unmapped_ids.tsv',sep='\t')
    print(df['DiseaseID'].str.extract('(\w+):.+').value_counts())
if __name__ == '__main__':
    while True:
        print("Enter your command number:\n1- Get each database UMLS mapping stats\n2- Get UMLS ID converter stats\n3- Get unmapped IDs stats\n4- Exit")
        command=input()
        if command=="1":
            command1()
        elif command=="2":
            command2()
        elif command=="3":
            command3()
        elif command=="4":
            exit()
        else:
            print('Try again!')


