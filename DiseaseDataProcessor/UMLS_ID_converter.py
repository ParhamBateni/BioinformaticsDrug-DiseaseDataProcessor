import os,json,obonet,pandas as pd

global mapping
mapping=dict()

result_file_name = 'UMLS_ID_converter.json'
def update_mapping(possible_xref,reference_id,ignore_redundant_id=True):
    if reference_id=="UMLS:nan" or 'nan' in possible_xref:return
    clean_possible_xref = clean_id(possible_xref)
    try:
        if reference_id not in mapping[clean_possible_xref]:
            if not ignore_redundant_id:
                mapping[clean_possible_xref].add(clean_id(reference_id))
    except:
        mapping[clean_possible_xref] = set()
        mapping[clean_possible_xref].add(clean_id(reference_id))
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


def find_colon(string):
    try:
        return string.index(':')
    except:
        return 0
def convert_DisGeNET():
    df=pd.read_csv("Data/DisGeNET/processed_DisGeNET.tsv",sep='\t').dropna(subset=['code'])
    for ind,row in df.iterrows():
        update_mapping(row["vocabulary"] + ":" + row["code"], row["DiseaseID"])
def convert_CLINVAR():
    df=pd.read_csv("Data/ClinVar/processed_ClinVar.tsv",sep='\t',dtype=str).dropna(subset=['SourceID'])
    for ind, row in df.iterrows():
        if str(row['SourceID']).startswith('ORPHA') and row['SourceName'] == 'Orphanet':
            disease_id = 'ORPHANET:' + disease_id[5:]
        elif str(row['SourceID']).isnumeric() and row['SourceName'] == 'OMIM':
            disease_id = 'OMIM:' + row['SourceID']
        elif str(row['SourceName']) == 'OMIM phenotypic series':
            disease_id = 'OMIM:' + row['SourceID'][2:]
        elif str(row['SourceName']) == 'NCBI curation':
            disease_id = 'NCBI:' + str(row['SourceID'])
        elif str(row['SourceName']) == 'PharmGKB':
            disease_id = 'PharmGKB:' + str(row['SourceID'])
        elif str(row['SourceName']) == 'GeneReviews':
            disease_id = 'GeneReviews:' + row['SourceID']
        else:
            disease_id = str(row['SourceID'])
        update_mapping(disease_id, row["DiseaseID"])
def is_exact_match(id_,property_values):
    id_=clean_id(id_)
    id_type=id_[:id_.index(":")]
    id_number=id_[id_.index(":")+1:]
    if id_type in ["ICD9","ICD10","ICD11","ICDO"]:
        return True
    elif id_type=="SCTID":
        for property_value in property_values:
            if not property_value.startswith("exactMatch") and not property_value.startswith("closeMatch"):
                continue
            property_value=str(property_value).upper()
            if id_number in property_value and "SNOMEDCT" in property_value:
                return True
        return False
    else:
        for property_value in property_values:
            if not property_value.startswith("exactMatch") and not property_value.startswith("closeMatch"):
                continue
            property_value=str(property_value).upper()
            if id_number in property_value and id_type in property_value:
                return True
        return False

def convert_DOID():
    graph=obonet.read_obo("Data/DOID/doid.obo.txt")
    for id_,data in graph.nodes(data=True):
        try:
            reference_id_xrefs=list(filter(lambda x:x.startswith("UMLS") ,data["xref"]))
            other_ids=list(filter(lambda x:not x.startswith("UMLS"),data["xref"]))
        except:
            continue
        other_ids.append(id_)
        for other_id in other_ids:
            for reference_id in reference_id_xrefs:
                update_mapping(other_id,reference_id)


def convert_EFO():
    with open("Data/EFO/efo.obo.txt", 'r', encoding='utf-8') as f:
        graph = obonet.read_obo(f)
        for id_, data in graph.nodes(data=True):
            try:
                reference_id_xrefs = list(filter(lambda x: x.startswith("UMLS")and is_exact_match(x,data.get("property_value")), data["xref"]))
                other_ids = list(filter(lambda x: not x.startswith("UMLS")and is_exact_match(x,data.get("property_value")), data["xref"]))
            except:
                continue
            other_ids.append(id_)
            for other_id in other_ids:
                for reference_id in reference_id_xrefs:
                    update_mapping(other_id, reference_id,ignore_redundant_id=False)
def convert_HPO():
    df = pd.read_csv("Data/HPO/processed_HPO.tsv", sep='\t').dropna(subset=['has_db_xref'])
    for ind, row in df.iterrows():
        try:
            reference_id_xrefs = list(filter(lambda x: x.startswith("UMLS"), row["has_db_xref"].split("|")))
        except:
            continue
        other_ids = list(filter(lambda x: not x.startswith("UMLS"), row["has_db_xref"].split("|")))
        other_ids.append(row["DiseaseID"])
        for other_id in other_ids:
            for reference_id in reference_id_xrefs:
                update_mapping(other_id, reference_id)
    pass
def convert_OTP():
    global mapping
    disease_files_name = os.listdir("Data/OpenTargetsPlatform")
    for file_name in disease_files_name:
        if file_name.endswith("json"):
            print(f'Getting IDS from {file_name}')
            with open(f"Data/OpenTargetsPlatform/{file_name}", 'r', encoding='utf-8') as f:
                for json_obj in f:
                    disease = json.loads(json_obj, )
                    ids = [disease['id'].replace('_', ':')] + disease['dbXRefs']
                    disease_ref_ids = list(filter(lambda x: x.startswith("UMLS"), ids))
                    other_ids=list(filter(lambda x:not x.startswith("UMLS"),ids))
                    for possible_xref in other_ids:
                        for ref_id in disease_ref_ids:
                            update_mapping(possible_xref, ref_id)


def convert_MedGen():
    global mapping
    df = pd.read_csv('Data/MedGen/MedGenIDMappings.txt', sep='|').dropna(subset=['source_id'])
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
        update_mapping(disease_id,'UMLS:'+row['CUI'])

def serialize_sets(obj):
    if isinstance(obj, set):
        return list(obj)
    raise TypeError ("Type %s is not serializable" % type(obj))
if __name__ == '__main__':
    # convert_MedGen()
    # print(mapping)
    # exit()
    candidate_databases=["EFO","ClinVar","DisGeNET","HPO","DOID","OpenTargetsPlatform","MedGen"]
    for candidate in candidate_databases:
        print(f"Converting {candidate}")
        if candidate=="DisGeNET":
            convert_DisGeNET()
        elif candidate=="ClinVar":
            convert_CLINVAR()
        elif candidate=="DOID":
            convert_DOID()
        elif candidate=="EFO":
            convert_EFO()
        elif candidate=="HPO":
            convert_HPO()
        elif candidate=="OpenTargetsPlatform":
            convert_OTP()
        elif candidate=="MedGen":
            convert_MedGen()
    print(f"Writing the results in {result_file_name}")
    with open(result_file_name, 'w') as f:
        json.dump(mapping, f,default=serialize_sets,sort_keys=True)
