#todo: remove this file :)


import pandas as pd

def load_dataset(dataset_address):
    df=pd.read_csv(dataset_address,sep='\t')
    return df
def get_diseases_common(combination,df):
    disease_names=[]
    disease_ids=[]
    for ind,row in df.iterrows():
        try:
            if combination.issubset(set([int(num) for num in str(row["Combinations"]).split(",")])):
                disease_names.append(row['DiseaseName'])
                disease_ids.append(row['DiseaseID'])
        except:
            continue
    return pd.DataFrame({"DiseaseID":disease_ids,"DiseaseName":disease_names})
def get_drugs_common(combination,df):
    drug_names = []
    drug_ids = []
    for ind, row in df.iterrows():
        try:
            if combination.issubset(set([int(num) for num in str(row["Combinations"])[:-1].split(",")])):
                drug_names.append(row["Name"])
                drug_ids.append(row["CID"])
        except:
            continue
    return pd.DataFrame({"CID": drug_ids, "Name": drug_names})

def process_drugs():
    dataset_address = "DrugRepurposing - DrugDBsResults.tsv"
    df = load_dataset(dataset_address)
    best_combinations = [(1, 3, 4, 12), (1, 3, 4, 10, 12), (1, 3, 4, 5, 10, 12), (1, 3, 4, 5, 6, 10, 12),
                          (1, 3, 4, 5, 6, 9, 10, 12), (1, 3, 4, 5, 6, 9, 10, 11, 12),
                          (1, 3, 4, 5, 6, 8, 9, 10, 11, 12)]
    result_file_name = "Common_drugs/Best_combinations_drugs_"
    for k in range(4, 11):
        drugs = get_drugs_common(set(best_combinations[k - 4]), df)
        drugs = drugs.set_index("CID")
        drugs.to_csv(result_file_name + str(k) + ".tsv", sep='\t')
def process_diseases():
    dataset_address = "DrugRepurposing - DiseaseDBsResults.tsv"
    df = load_dataset(dataset_address)
    best_combinations = [(2, 3, 7, 20), (2, 7, 9, 18, 20), (2, 3, 7, 9, 18, 20), (2, 3, 5, 7, 9, 18, 20),
                             (2, 3, 5, 7, 8, 9, 18, 20), (2, 3, 5, 7, 8, 9, 15, 18, 20),
                             (2, 3, 5, 7, 8, 9, 15, 18, 19, 20)]
    result_file_name = "Common_diseases/Best_combinations_diseases_"
    for k in range(4, 11):
        diseases = get_diseases_common(set(best_combinations[k - 4]), df)
        diseases = diseases.set_index("DiseaseID")
        diseases.to_csv(result_file_name + str(k) + ".tsv", sep='\t')
if __name__ == '__main__':
    while True:
        print("Drug or Disease?")
        inp=input()
        if inp=="Drug":
           process_drugs()
           break
        elif inp=="Disease":
            process_diseases()
            break
        else:
            print("Try again!")

