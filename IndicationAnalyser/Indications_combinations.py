import pandas as pd

global disease_best_combinations, drug_best_combinations
global indication_diseases_df, indication_drugs_df
disease_best_combinations = [(3, 7, 20, 21), (7, 9, 18, 20, 21), (3, 7, 9, 18, 20, 21), (3, 5, 7, 9, 18, 20, 21),
                             (3, 5, 7, 8, 9, 18, 20, 21), (3, 5, 7, 8, 9, 15, 18, 20, 21),
                             (3, 5, 7, 8, 9, 15, 18, 19, 20, 21)]
drug_best_combinations = [(3, 4, 12, 15), (3, 4, 10, 12, 15), (3, 4, 5, 10, 12, 15), (3, 4, 5, 6, 10, 12, 15),
                          (3, 4, 5, 6, 9, 10, 12, 15), (1, 3, 4, 5, 6, 9, 10, 12, 15),
                          (1, 3, 4, 5, 6, 9, 10, 11, 12, 15)]


def get_count_indications_in_best_combinations(disease_combination: set, drug_combination: set,i,j):
    global disease_best_combinations, drug_best_combinations
    global indication_diseases_df, indication_drugs_df
    count = 0
    best_indexes = []
    best_indications=pd.DataFrame({"DiseaseID":[],"DrugCID":[]})
    disease_IDs=[]
    drug_CIDs=[]
    for index in range(len(indication_drugs_df)):
        if disease_combination.issubset(
                set(list(map(int, indication_diseases_df.loc[index, "Combinations"].split(","))))):
            if drug_combination.issubset(
                    set(list(map(int, indication_drugs_df.loc[index, "Combinations"].split(","))))):
                count += 1
                best_indexes.append(index)
                disease_IDs.append(indication_diseases_df.loc[index,"DiseaseID"])
                drug_CIDs.append(indication_drugs_df.loc[index,"CID"])
    # print(best_indexes[:10])
    best_indications['DiseaseID']=disease_IDs
    best_indications['DrugCID']=drug_CIDs
    best_indications.to_csv(f"Best Drug-Disease Indications/Disease{j}-Drug{i} indications.tsv",sep='\t',index=0)
    return count


if __name__ == '__main__':
    # rows are for diseases and columns are for drugs
    df = pd.DataFrame({"K=4": [], "K=5": [], "K=6": [], "K=7": [], "K=8": [], "K=9": [], "K=10": []})

    indication_diseases_df = pd.read_csv("Indications_diseases.tsv", sep="\t")
    indication_drugs_df = pd.read_csv("Indications_drugs.tsv", sep="\t")

    for i in range(4, 11):
        col_counts = []
        for j in range(4, 11):
            col_counts.append(get_count_indications_in_best_combinations(set(disease_best_combinations[j - 4]),
                                                                         set(drug_best_combinations[i - 4]),i,j))
        df[f"K={i}"] = col_counts
    df.index = ["K=4", "K=5", "K=6", "K=7", "K=8", "K=9", "K=10"]
    df.to_csv("Indications_combinations_results.tsv", sep="\t")
    print(df)
