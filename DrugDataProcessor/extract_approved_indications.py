import json

import pandas as pd
from tqdm import tqdm

from DiseaseDataProcessor.UMLS_ID_converter import clean_id
from DiseaseDataProcessor.main import load_ids_history_data


def extract_indications_from_RepoDB(disease_IDs: list, drug_CIDs: list):
    raw_df = pd.read_csv('Data/RepoDB/repodb.txt', sep=',', dtype=(str, str))
    parsed_df = pd.read_csv('Data/RepoDB/parsed_RepoDB.tsv', sep='\t', dtype=str)
    indication_pairs = dict()
    ids_history = load_ids_history_data()
    for ind, row in raw_df.iterrows():
        if row['status'] == 'Approved':
            try:
                indication_pairs[row['drug_name']].append(row['ind_id'])
            except:
                indication_pairs[row['drug_name']] = [row['ind_id']]
    for drug_name in indication_pairs:
        query = parsed_df.query(f"DrugName=='{drug_name}'")
        for CID in query['CID']:
            for disease_ID in indication_pairs[drug_name]:
                if ids_history.get("UMLS:" + disease_ID) is not None:
                    for id_history in ids_history.get('UMLS:' + disease_ID):
                        disease_IDs.append(id_history)
                        drug_CIDs.append(CID)
                else:
                    disease_IDs.append("UMLS:" + disease_ID)
                    drug_CIDs.append(CID)


def extract_indications_from_DrugRepurposingHub(disease_IDs: list, drug_CIDs: list):
    # Extracting indications from DrugRepurposingHub

    # response = req.request(
    #     url='http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.06/output/etl/json/indication/', method='get')
    #
    # files = re.findall('href="(part[A-Za-z0-9-]*.json)', response.text)
    # for i in tqdm(range(len(files)),'Extracting approved indications from DrugRepurposingHub'):
    #     response_i = req.request(
    #         url=f'http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.06/output/etl/json/indication/{files[i]}',
    #         method='get')
    #     text_i = '[' + response_i.text + ']'
    #     text_i = text_i.replace('\n', ',')
    #     text_i = text_i[:-2] + ']'
    #     drugs_json=json.loads(text_i)
    #     for drug in drugs_json:
    #         for approved_indication in drug['approvedIndications']:
    #             diseases.append(approved_indication)
    #             drugs.append(drug['id'])
    # pd.DataFrame({'DiseaseID':disease_IDs,'DrugID':drugs_IDs}).sort_values(by='DiseaseID').to_csv('DrugRepurposingHub_indications.tsv',sep='\t')
    df = pd.read_csv('Approved Indication data/DrugRepurposingHub_indications.tsv', sep='\t', dtype=str)
    UMLS_ID_converter = json.loads(open('../DiseaseDataProcessor/UMLS_ID_converter.json').read())
    parsed_Chembl_to_CID = pd.read_csv('Approved Indication data/parsed_DrugRepurposingHub_Chembl_to_CID.txt', sep='\t',
                                       index_col=None, dtype=str)
    ids_history = load_ids_history_data()
    for i in tqdm(range(len(df))):
        row = df.iloc[i]
        id_ = row['DiseaseID'].replace('_', ':')
        mapped_ids = UMLS_ID_converter.get(clean_id(id_))
        if mapped_ids is None:
            # print(clean_id(id_))
            continue
        for mapped_id in mapped_ids:
            if ids_history.get(mapped_id) is not None:
                for id_history in ids_history.get(mapped_id):
                    disease_IDs.append(id_history)
                    drug_CIDs.append(parsed_Chembl_to_CID.query(f'DrugID=="{row["DrugID"]}"')['CID'].values[0])
            else:
                disease_IDs.append(mapped_id)
                drug_CIDs.append(parsed_Chembl_to_CID.query(f'DrugID=="{row["DrugID"]}"')['CID'].values[0])


def find_indication_pairs():
    disease_IDs = []
    drug_CIDs = []

    extract_indications_from_RepoDB(disease_IDs, drug_CIDs)
    extract_indications_from_DrugRepurposingHub(disease_IDs, drug_CIDs)

    indication_df = pd.DataFrame({'DiseaseID': disease_IDs, 'DrugCID': drug_CIDs})
    # print(len(set(indication_df['DiseaseID'])))
    # print(len(set(indication_df['DrugCID'])))

    indication_df = indication_df.drop_duplicates().set_index('DiseaseID').sort_index().dropna()
    indication_df.to_csv("Approved Indication data/indication_pairs.tsv", sep='\t')


def save_indication_diseases():
    disease_df = pd.read_csv('../DiseaseDataProcessor/disease_results.tsv', sep='\t', dtype=str).set_index('DiseaseID')
    disease_df.loc[
    pd.read_csv('../DrugDataProcessor/Approved Indication data/indication_pairs.tsv', sep='\t', dtype=str)[
        'DiseaseID'].tolist(), :].to_csv('../IndicationAnalyser/Indications_diseases.tsv', sep='\t')


def save_indication_drugs():
    drug_df = pd.read_csv('../DrugDataProcessor/drug_results.tsv', sep='\t', dtype=str).set_index('CID')
    drug_df.loc[pd.read_csv('../DrugDataProcessor/Approved Indication data/indication_pairs.tsv', sep='\t', dtype=str)[
                    'DrugCID'].tolist(),
    :].to_csv('../IndicationAnalyser/Indications_drugs.tsv', sep='\t')


if __name__ == '__main__':
    find_indication_pairs()
    save_indication_diseases()
    save_indication_drugs()
