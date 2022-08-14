import pandas as pd
from itertools import combinations

global disease_best_combinations, drug_best_combinations
global indication_diseases_df, indication_drugs_df
disease_best_combinations=[None]*7
drug_best_combinations=[None]*7
# disease_best_combinations = [(7, 9, 18, 21), (3, 7, 9, 18, 21), (3, 5, 7, 9, 18, 21), (3, 5, 7, 8, 9, 18, 21),
#                              (3, 5, 7, 8, 9, 15, 18, 21),(3, 5, 7, 8, 9, 15, 18, 19, 21),
#                              (3, 5, 6, 7, 8, 9, 15, 18, 19, 21)]
# drug_best_combinations = [(3, 4, 10, 12), (3, 4, 5, 10, 12), (3, 4, 5, 6, 10, 12), (3, 4, 5, 6, 9, 10, 12),
#                           (1, 3, 4, 5, 6, 9, 10, 12), (1, 3, 4, 5, 6, 9, 10, 11, 12),
#                           (1, 3, 4, 5, 6, 8, 9, 10, 11, 12)]


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
    best_indications['DiseaseID']=disease_IDs
    best_indications['DrugCID']=drug_CIDs
    best_indications.to_csv(f"Best Drug-Disease Indications/Disease{j}-Drug{i} indications.tsv",sep='\t',index=0)
    return len(set(drug_CIDs)),len(set(disease_IDs)),count

global analysis_result,analysis_output_file_address
analysis_result=""

def get_combinations(file_address):
    # First exclude all rows with nan value for Combinations column(like Sum row) and then call this function!
    df=pd.read_csv(file_address, sep='\t')
    return [set(list(map(int,str(comb).split(",")))) for comb in list(df["Combinations"])]

def get_count_common_in_databases(comb,database_combinations):
    count=sum([1 for d in database_combinations if comb.issubset(d)])
    return count


def analyse_drugs():
    global analysis_result,analysis_output_file_address
    file_address = "Indications_drugs.tsv"
    analysis_output_file_address = "Indications drugs combination analysis.txt"
    database_combinations = get_combinations(file_address)
    #Databases with small count of drugs mapped(2,14) and RepoDB(15) are ignored
    databases_to_inspect=[1,3,4,5,6,7,8,9,10,11,12,13,16,17]
    #
    analysis_result=""
    for k in range(4, 11):
        combs_stats = dict()
        combs = combinations(databases_to_inspect, k)
        for comb in combs:
            combs_stats[comb] = get_count_common_in_databases(set(comb), database_combinations)
        sorted_combs_stats = {k: v for k, v in sorted(combs_stats.items(), key=lambda item: item[1], reverse=True)}

        analysis_result += f"k={k}: " + str(list(sorted_combs_stats.items())[:20]) + "\n"
        drug_best_combinations[k-4]=list(sorted_combs_stats.items())[0][0]
        print(f"k={k}: ", end="")
        print(list(sorted_combs_stats.items())[:10])

    with open(analysis_output_file_address, 'w') as f:
        f.write(analysis_result)

def analyse_diseases():
    global analysis_result,analysis_output_file_address
    file_address = "Indications_diseases.tsv"
    analysis_output_file_address = "Indications diseases combination analysis.txt"
    database_combinations = get_combinations(file_address)
    #Databases with small count of diseases mapped(1,11,17,22) and RepoDB(20) are ignored
    databases_to_inspect=[2,3,4,5,6,7,8,9,10,12,13,14,15,16,18,19,21]

    analysis_result=""
    for k in range(4, 11):
        combs_stats = dict()
        combs = combinations(databases_to_inspect, k)
        for comb in combs:
            combs_stats[comb] = get_count_common_in_databases(set(comb), database_combinations)
        sorted_combs_stats = {k: v for k, v in sorted(combs_stats.items(), key=lambda item: item[1], reverse=True)}
        analysis_result += f"k={k}: " + str(list(sorted_combs_stats.items())[:20]) + "\n"
        disease_best_combinations[k-4]=list(sorted_combs_stats.items())[0][0]
        print(f"k={k}: ", end="")
        print(list(sorted_combs_stats.items())[:10])
    with open(analysis_output_file_address, 'w') as f:
        f.write(analysis_result)
if __name__ == '__main__':
    print("Disease:")
    analyse_diseases()
    print("Drug:")
    analyse_drugs()

    # rows are for diseases and columns are for drugs
    df = pd.DataFrame({"K=4": [], "K=5": [], "K=6": [], "K=7": [], "K=8": [], "K=9": [], "K=10": []})

    indication_diseases_df = pd.read_csv("Indications_diseases.tsv", sep="\t")
    indication_drugs_df = pd.read_csv("Indications_drugs.tsv", sep="\t")

    for i in range(4, 11):
        col_counts = []
        for j in range(4, 11):
            count_unique_drugs,count_unique_diseases,count=get_count_indications_in_best_combinations(set(disease_best_combinations[j - 4]),
                                                                         set(drug_best_combinations[i - 4]),i,j)
            col_counts.append(f"{count} ({count_unique_diseases},{count_unique_drugs})")
        df[f"K={i}"] = col_counts
    df.index = ["K=4", "K=5", "K=6", "K=7", "K=8", "K=9", "K=10"]
    df.to_csv("Indications_combinations_results.tsv", sep="\t")
    print(df)
