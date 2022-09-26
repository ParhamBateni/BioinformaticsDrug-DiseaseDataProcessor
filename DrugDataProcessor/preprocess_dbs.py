import json

import pandas as pd

data_folder = 'Data/'


def process_BindingDB():
    df = pd.read_csv(data_folder + 'BindingDB/BindingDB_CID.txt', sep='\t', dtype='str', index_col=0, header=None)[1]
    CIDs = df.drop_duplicates()
    processed_df = pd.DataFrame({'CID': CIDs}).set_index('CID').sort_index()

    # pubchem id exchange limit for count of ids is 500000
    processed_df[:500000].to_csv(data_folder + 'BindingDB/processed_BindingDB1.txt', header=None, sep='\t')
    processed_df[500000:].to_csv(data_folder + 'BindingDB/processed_BindingDB2.txt', header=None, sep='\t')


def process_BioGRID():
    df = pd.read_csv(data_folder + "BIOGRID/BIOGRID-CHEMICALS-4.4.212.chemtab.txt", sep='\t')["Chemical Name"]
    names = df.drop_duplicates().tolist()
    processed_df = pd.DataFrame({"DrugName": []})
    processed_df['DrugName'] = names
    processed_df = processed_df.set_index("DrugName")
    processed_df = processed_df.sort_index()
    processed_df.to_csv(data_folder + "BioGRID/processed_BioGRID.txt", sep='\t')


def process_Cancerrxgene():
    df = pd.read_csv(data_folder + "Cancerrxgene/cancerrxgene.txt", sep='\t', header=0)
    names = df.drop_duplicates()
    processed_df = pd.DataFrame({"DrugName": []})
    processed_df['DrugName'] = names
    processed_df = processed_df.set_index("DrugName")
    processed_df = processed_df.sort_index()
    processed_df.to_csv(data_folder + "Cancerrxgene/processed_Cancerrxgene.txt", sep='\t')


def process_CLUE():
    df = pd.read_csv(data_folder + 'CLUE/compoundinfo_beta.txt', sep='\t')
    df['canonical_smiles'].dropna().drop_duplicates().to_csv(data_folder + 'CLUE/processed_CLUE_smiles.txt', sep='\t',
                                                             index=0, header=None)
    df[df['canonical_smiles'].isna()]['cmap_name'].dropna().drop_duplicates() \
        .to_csv(data_folder + 'CLUE/processed_CLUE_nan_smiles_names.txt', sep='\t', index=0, header=None)


def process_CTDbase():
    df = pd.read_csv(data_folder + "CTDbase/CTD_chemicals.txt", sep='\t')["ChemicalName"]
    names = df.drop_duplicates()
    processed_df = pd.DataFrame({"DrugName": []})
    processed_df['DrugName'] = names
    processed_df = processed_df.set_index("DrugName")
    processed_df = processed_df.sort_index()
    processed_df.to_csv(data_folder + "CTDbase/processed_CTDbase.txt", sep='\t')


def process_DGIdb():
    df = pd.read_csv(data_folder + "DGIdb/dgidb.txt", sep='\t')["drug_name"]
    names = df.drop_duplicates()
    processed_df = pd.DataFrame({"DrugName": []})
    processed_df['DrugName'] = names
    processed_df = processed_df.set_index("DrugName")
    processed_df = processed_df.sort_index()
    processed_df.to_csv(data_folder + "DGIdb/processed_DGIdb.txt", sep='\t')


def process_DrugBank():
    drug_names = []
    with open(data_folder + 'DrugBank/full database.xml', 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('  <name>'):
                drug_names.append(line[8:line.find('</name>')])
    pd.DataFrame({'DrugName': drug_names}).drop_duplicates().set_index('DrugName').sort_index() \
        .to_csv(data_folder + 'DrugBank/processed_DrugBank2.txt', sep='\t')


def process_DrugRepurposingHub():
    df = pd.read_csv(data_folder + "DrugRepurposingHub/drug_repurposing-hub.txt", sep='\t')["pert_iname"]
    names = df.drop_duplicates()
    processed_df = pd.DataFrame({"DrugName": []})
    processed_df['DrugName'] = names
    processed_df = processed_df.set_index("DrugName")
    processed_df = processed_df.sort_index()
    processed_df.to_csv(data_folder + "DrugRepurposingHub/processed_DrugRepurposingHub.txt", sep='\t')


def process_DrugCentral():
    df = pd.read_csv(data_folder + "DrugCentral/drugcentral.txt", sep='\t')["DRUG_NAME"]
    names = df.drop_duplicates()
    processed_df = pd.DataFrame({"DrugName": []})
    processed_df['DrugName'] = names
    processed_df = processed_df.set_index("DrugName")
    processed_df = processed_df.sort_index()
    processed_df.to_csv(data_folder + "DrugCentral/processed_DrugCentral.txt", sep='\t')


def process_DrugComb():
    drug_names = []
    cids = []
    with open(data_folder + 'DrugComb/drugs.json', 'r', encoding='utf-8') as f:
        entities = [e for e in f.read()[1:-1].split('},{')]
        for entity in entities:
            if entity[0] != '{': entity = '{' + entity
            if entity[-1] != '}': entity = entity + '}'
            drug = json.loads(entity)
            drug_names.append(drug['dname'])
            cids.append(str(drug['cid']))
    df = pd.DataFrame({'DrugName': drug_names, 'CID': cids})
    df.replace(['0', 'None'], '').set_index('DrugName').sort_index().to_csv(
        data_folder + 'DrugComb/parsed_DrugComb.tsv', sep='\t')


def process_DrugCombDB():
    # with open(data_folder + "DrugCombDB/drug_chemical_info.txt",'r',encoding='utf-8',errors='ignore') as f:
    #     f.read()
    df = pd.read_csv(data_folder + "DrugCombDB/drug_chemical_info.txt", sep=',', encoding_errors='ignore')["drugName"]
    names = df.drop_duplicates()
    processed_df = pd.DataFrame({"DrugName": []})
    processed_df['DrugName'] = names
    processed_df = processed_df.set_index("DrugName")
    processed_df = processed_df.sort_index()
    processed_df.to_csv(data_folder + "DrugCombDB/processed_DrugCombDB.txt", sep='\t')


def process_DrugPath():
    df = pd.read_csv(data_folder + "DrugPath/drugpath.txt", sep='\t', header=0)
    df = df[df.columns[1]]
    names = df.drop_duplicates()
    processed_df = pd.DataFrame({"DrugName": []})
    processed_df['DrugName'] = names
    processed_df = processed_df.set_index("DrugName")
    processed_df = processed_df.sort_index()
    processed_df.to_csv(data_folder + "DrugPath/processed_DrugPath.txt", sep='\t')


def process_DsigDB():
    drugs = set()
    with open(f"{data_folder}/DsigDB/dsigdb.txt", 'r')as f:
        for line in f:
            if line.startswith("Compound"):
                drugs.add(line[line.index(":") + 2:-1])
    processed_df = pd.DataFrame({"DrugName": list(drugs)})
    processed_df = processed_df.set_index("DrugName")
    processed_df = processed_df.sort_index()
    processed_df.to_csv(data_folder + "DsigDB/processed_DsigDB.txt", sep='\t')


def process_GtopDB():
    df = pd.read_csv(data_folder + "GtopDB/gtopdb.txt", sep=',')["Name"]
    names = df.drop_duplicates()
    processed_df = pd.DataFrame({"DrugName": []})
    processed_df['DrugName'] = names
    processed_df = processed_df.set_index("DrugName")
    processed_df = processed_df.sort_index()
    processed_df.to_csv(data_folder + "GtopDB/processed_GtopDB.txt", sep='\t')


def process_Hetionet():
    df = pd.read_csv(data_folder + "Hetionet/hetionet-v1.0-nodes.txt", sep='\t').query("kind=='Compound'")['name']
    names = df.drop_duplicates()
    processed_df = pd.DataFrame({"DrugName": []})
    processed_df['DrugName'] = names
    processed_df = processed_df.set_index("DrugName")
    processed_df = processed_df.sort_index()
    processed_df.to_csv(data_folder + "Hetionet/processed_Hetionet.txt", sep='\t')


def process_PDSPKIDB():
    df = pd.read_csv(data_folder + "PDSPKiDB/kidb.txt", sep=',', encoding_errors='ignore', header=0)
    names = df.drop_duplicates()
    processed_df = pd.DataFrame({"DrugName": []})
    processed_df['DrugName'] = names
    processed_df = processed_df.set_index("DrugName")
    processed_df = processed_df.sort_index()
    processed_df.to_csv(data_folder + "PDSPKiDB/processed_PDSPKiDB.txt", sep='\t')


def process_pharmGKB():
    df = pd.read_csv(data_folder + "pharmGKB/pharmgkb.txt", sep='\t')["Chemicals"]
    names = df.drop_duplicates()
    processed_df = pd.DataFrame({"DrugName": []})
    processed_df['DrugName'] = names
    processed_df = processed_df.set_index("DrugName")
    processed_df = processed_df.sort_index()
    processed_df.to_csv(data_folder + "pharmGKB/processed_pharmGKB.txt", sep='\t')


def process_RepoDB():
    df = pd.read_csv(data_folder + "RepoDB/repodb.txt", sep=',')["drug_name"]
    names = df.drop_duplicates()
    processed_df = pd.DataFrame({"DrugName": []})
    processed_df['DrugName'] = names
    processed_df = processed_df.set_index("DrugName")
    processed_df = processed_df.sort_index()
    processed_df.to_csv(data_folder + "RepoDB/processed_RepoDB.txt", sep='\t')


def process_SIDER():
    # this database has already CID info in it but the CIDs are different from PubChem!!
    # todo: check SIDER database
    df = pd.read_csv(data_folder + "SIDER/SIDER_drug_names.txt", sep='\t')["DrugName"]
    names = df.drop_duplicates()
    processed_df = pd.DataFrame({"DrugName": []})
    processed_df['DrugName'] = names
    processed_df = processed_df.set_index("DrugName")
    processed_df = processed_df.sort_index()
    processed_df.to_csv(data_folder + "SIDER/processed_SIDER.txt", sep='\t')


def process_SMPDB():
    df = pd.read_csv(data_folder + "SMPDB/smpdb_inchi_keys.txt", sep='\t', header=0)
    keys = df.drop_duplicates()
    processed_df = pd.DataFrame({"InchiKey": []})
    processed_df['InchiKey'] = keys
    processed_df = processed_df.set_index("InchiKey")
    processed_df = processed_df.sort_index()
    processed_df.to_csv(data_folder + "SMPDB/processed_SMPDB.txt", sep='\t')


def process_SWEETLEAD():
    df = pd.read_csv(data_folder + 'SWEETLEAD/SWEETLEAD.txt', sep='\t', dtype=str)
    CIDs = []
    ChEBIs = []
    for index, row in df.iterrows():
        cid = str(row['PubChem IDs'])
        if cid == 'nan':
            chebi = str(row['ChEBI IDs'])
            if chebi == 'nan':
                continue
            ChEBIs.extend(str(row['ChEBI IDs']).split(','))
        else:
            CIDs.extend(cid.split(','))
    pd.DataFrame({'CID': CIDs}).dropna().drop_duplicates().set_index('CID').sort_index() \
        .to_csv(data_folder + 'SWEETLEAD/processed_SWEETLEAD_CIDs.txt', sep='\t', header=None)
    pd.DataFrame({'ChEBI': ChEBIs}).dropna().drop_duplicates().set_index('ChEBI').sort_index() \
        .to_csv(data_folder + 'SWEETLEAD/processed_SWEETLEAD_ChEBIs.txt', sep='\t', header=None)


def process_TTD():
    drug_ids = []
    with open(data_folder + 'TTD/P1-02-TTD_drug_download.txt', 'r', encoding='utf-8') as f:
        start = 0
        for line in f:
            if start == 0:
                if line.strip() == '________________________________________________________________________':
                    start = 1
            else:
                if 'DRUG__ID' in line:
                    drug_ids.append(line.split()[-1])
    pd.DataFrame({'TTD_ID': drug_ids}).drop_duplicates().set_index('TTD_ID').sort_index().to_csv(
        data_folder + 'TTD/processed_TTD.txt', sep='\t')


if __name__ == '__main__':
    databases_processors = {'BindingDB': process_BindingDB, 'BioGRID': process_BioGRID,
                            'Cancerrxgene': process_Cancerrxgene,
                            'CLUE': process_CLUE, 'CTDbase': process_CTDbase, 'DGIdb': process_DGIdb,
                            'DrugBank': process_DrugBank,
                            'DrugCentral': process_DrugCentral, 'DrugComb': process_DrugComb,
                            'DrugCombDB': process_DrugCombDB,
                            'DrugPath': process_DrugPath, 'DrugRepurposingHub': process_DrugRepurposingHub,
                            'DsigDB': process_DsigDB,
                            'GtopDB': process_GtopDB, 'Hetionet': process_Hetionet, 'PDSPKiDB': process_PDSPKIDB,
                            'pharmGKB': process_pharmGKB, 'RepoDB': process_RepoDB, 'SIDER': process_SIDER,
                            'SMPDB': process_SMPDB,
                            'SWEETLEAD': process_SWEETLEAD, 'TTD': process_TTD}
    while True:
        print('\n'.join([f'{i + 1}.{database_name}' for i, database_name in enumerate(databases_processors.keys())]))
        print("Select a database to preprocess:", end='')
        command = input()
        if command.isdigit() and int(command) - 1 < len(databases_processors):
            databases_processors[list(databases_processors.keys())[int(command) - 1]]()
        elif command == 'exit':
            exit()
        else:
            print("Try again! (enter a number)")
        print()
