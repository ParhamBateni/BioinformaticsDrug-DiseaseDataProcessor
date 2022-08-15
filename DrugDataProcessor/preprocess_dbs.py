import pandas as pd
import requests
import json
data_folder='Data/'
def process_BioGRID():
    df=pd.read_csv(data_folder+"BIOGRID/BIOGRID-CHEMICALS-4.4.212.chemtab.txt",sep='\t')["Chemical Name"]
    names=df.drop_duplicates().tolist()
    processed_df=pd.DataFrame({"DrugName":[]})
    processed_df['DrugName']=names
    processed_df=processed_df.set_index("DrugName")
    processed_df = processed_df.sort_index()
    processed_df.to_csv(data_folder + "BioGRID/processed_BioGRID.txt", sep='\t')
def process_Cancerrxgene():
    df = pd.read_csv(data_folder + "Cancerrxgene/cancerrxgene.txt", sep='\t',header=0)
    names = df.drop_duplicates()
    processed_df=pd.DataFrame({"DrugName":[]})
    processed_df['DrugName'] = names
    processed_df = processed_df.set_index("DrugName")
    processed_df = processed_df.sort_index()
    processed_df.to_csv(data_folder + "Cancerrxgene/processed_Cancerrxgene.txt", sep='\t')
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
def process_DrugCombDB():
    # with open(data_folder + "DrugCombDB/drug_chemical_info.txt",'r',encoding='utf-8',errors='ignore') as f:
    #     f.read()
    df = pd.read_csv(data_folder + "DrugCombDB/drug_chemical_info.txt", sep=',',encoding_errors='ignore')["drugName"]
    names = df.drop_duplicates()
    processed_df = pd.DataFrame({"DrugName": []})
    processed_df['DrugName'] = names
    processed_df = processed_df.set_index("DrugName")
    processed_df = processed_df.sort_index()
    processed_df.to_csv(data_folder + "DrugCombDB/processed_DrugCombDB.txt", sep='\t')
def process_DrugPath():
    df = pd.read_csv(data_folder + "DrugPath/drugpath.txt", sep='\t',header=0)
    df=df[df.columns[1]]
    names = df.drop_duplicates()
    processed_df = pd.DataFrame({"DrugName": []})
    processed_df['DrugName'] = names
    processed_df = processed_df.set_index("DrugName")
    processed_df = processed_df.sort_index()
    processed_df.to_csv(data_folder + "DrugPath/processed_DrugPath.txt", sep='\t')
def process_DsigDB():
    drugs=set()
    with open(f"{data_folder}/DsigDB/dsigdb.txt",'r')as f:
        for line in f:
            if line.startswith("Compound"):
                drugs.add(line[line.index(":")+2:-1])
    processed_df=pd.DataFrame({"DrugName":list(drugs)})
    processed_df=processed_df.set_index("DrugName")
    processed_df=processed_df.sort_index()
    processed_df.to_csv(data_folder+"DsigDB/processed_DsigDB.txt",sep='\t')
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
    df = pd.read_csv(data_folder + "PDSPKiDB/kidb.txt", sep=',',encoding_errors='ignore',header=0)
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
    #todo: check SIDER database
    df = pd.read_csv(data_folder + "SIDER/SIDER_drug_names.txt", sep='\t')["DrugName"]
    names = df.drop_duplicates()
    processed_df = pd.DataFrame({"DrugName": []})
    processed_df['DrugName'] = names
    processed_df = processed_df.set_index("DrugName")
    processed_df = processed_df.sort_index()
    processed_df.to_csv(data_folder + "SIDER/processed_SIDER.txt", sep='\t')
def process_SMPDB():
    df = pd.read_csv(data_folder + "SMPDB/smpdb_inchi_keys.txt", sep='\t',header=0)
    names = df.drop_duplicates()
    processed_df = pd.DataFrame({"DrugName": []})
    processed_df['DrugName'] = names
    processed_df = processed_df.set_index("DrugName")
    processed_df = processed_df.sort_index()
    processed_df.to_csv(data_folder + "SMPDB/processed_SMPDB.txt", sep='\t')

if __name__ == '__main__':
    command=""
    while True:
        print("1.BioGRID\n2.CTDbase\n3.DGIdb\n4.DrugRepurposingHub\n5.DrugCentral\n6.DrugCombDB\n7.DrugPath\n"
              "8.DsigDB\n9.GtopDB\n10.Hetionet\n11.PDSPKIDB\n12.pharmGKB\n13.RepoDB\n14.SIDER\n15.SMPDB\n16.Cancerrxgene\n17.exit")
        print("Select a database to preprocess:",end='')
        command=input()
        if command=='1':
            process_BioGRID()
        elif command=='2':
            process_CTDbase()
        elif command=='3':
            process_DGIdb()
        elif command=='4':
            process_DrugRepurposingHub()
        elif command=='5':
            process_DrugCentral()
        elif command=='6':
            process_DrugCombDB()
        elif command=='7':
            process_DrugPath()
        elif command=='8':
            process_DsigDB()
        elif command=='9':
            process_GtopDB()
        elif command=='10':
            process_Hetionet()
        elif command=='11':
            process_PDSPKIDB()
        elif command=='12':
            process_pharmGKB()
        elif command=='13':
            process_RepoDB()
        elif command=='14':
            process_SIDER()
        elif command=='15':
            process_SMPDB()
        elif command=='16':
            process_Cancerrxgene()
        elif command=='17':
            exit()
        else:
            print("Try again! (enter a number)")
        print()