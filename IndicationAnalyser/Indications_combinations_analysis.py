from itertools import combinations

import pandas as pd

global disease_best_combinations, drug_best_combinations
global indication_diseases_df, indication_drugs_df
disease_best_combinations = [None] * 7
drug_best_combinations = [None] * 7


# disease_best_combinations = [(7, 9, 18, 21), (3, 7, 9, 18, 21), (3, 5, 7, 9, 18, 21), (3, 5, 7, 8, 9, 18, 21),
#                              (3, 5, 7, 8, 9, 15, 18, 21),(3, 5, 7, 8, 9, 15, 18, 19, 21),
#                              (3, 5, 6, 7, 8, 9, 15, 18, 19, 21)]
# drug_best_combinations = [(3, 4, 10, 12), (3, 4, 5, 10, 12), (3, 4, 5, 6, 10, 12), (3, 4, 5, 6, 9, 10, 12),
#                           (1, 3, 4, 5, 6, 9, 10, 12), (1, 3, 4, 5, 6, 9, 10, 11, 12),
#                           (1, 3, 4, 5, 6, 8, 9, 10, 11, 12)]


def get_count_indications_in_best_combinations(disease_combination: set, drug_combination: set, i, j):
    global disease_best_combinations, drug_best_combinations
    global indication_diseases_df, indication_drugs_df
    count = 0
    best_indexes = []
    best_indications = pd.DataFrame({"DiseaseID": [], "DrugCID": []})
    disease_IDs = []
    drug_CIDs = []
    for index in range(len(indication_drugs_df)):
        if disease_combination.issubset(
                set(list(map(int, indication_diseases_df.loc[index, "Combinations"].split(","))))):
            if drug_combination.issubset(
                    set(list(map(int, indication_drugs_df.loc[index, "Combinations"].split(","))))):
                count += 1
                best_indexes.append(index)
                disease_IDs.append(indication_diseases_df.loc[index, "DiseaseID"])
                drug_CIDs.append(indication_drugs_df.loc[index, "CID"])
    best_indications['DiseaseID'] = disease_IDs
    best_indications['DrugCID'] = drug_CIDs
    best_indications.to_csv(f"Best Drug-Disease Indications/Disease{j}-Drug{i} indications.tsv", sep='\t', index=0)
    return len(set(drug_CIDs)), len(set(disease_IDs)), count


global analysis_result, analysis_output_file_address
analysis_result = ""


def get_combinations(file_address):
    # First exclude all rows with nan value for Combinations column(like Sum row) and then call this function!
    df = pd.read_csv(file_address, sep='\t')
    return [set(list(map(int, str(comb).split(",")))) for comb in
            list(df[df['Sum'].astype(float) >= 4]["Combinations"])]


def get_count_common_in_databases(comb, database_combinations):
    count = sum([1 for d in database_combinations if comb.issubset(d)])
    return count


def analyse_drugs():
    global analysis_result, analysis_output_file_address
    file_address = "Indications_drugs.tsv"
    analysis_output_file_address = "Indications drugs combination analysis.txt"
    database_combinations = get_combinations(file_address)
    # Databases with small count of drugs mapped(3,11,14,17,19) and RepoDB(18) are ignored
    databases_to_inspect = [1, 2, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 16, 20, 21, 22]

    analysis_result = ""
    analysis_table = pd.DataFrame({'K': [], 'Best Combination': [], 'Count of Indications Drugs': []})
    for k in range(4, 11):
        combs_stats = dict()
        combs = combinations(databases_to_inspect, k)
        for comb in combs:
            combs_stats[comb] = get_count_common_in_databases(set(comb), database_combinations)
        sorted_combs_stats = {k: v for k, v in sorted(combs_stats.items(), key=lambda item: item[1], reverse=True)}

        analysis_result += f"k={k}: " + str(list(sorted_combs_stats.items())[:20]) + "\n"
        drug_best_combinations[k - 4] = list(sorted_combs_stats.items())[0][0]
        print(f"k={k}: ", end="")
        print(list(sorted_combs_stats.items())[:10])
        analysis_table = analysis_table.append(
            {'K': str(k), 'Best Combination': str(list(sorted_combs_stats.items())[0][0]),
             'Count of Indications Drugs': str(list(sorted_combs_stats.items())[0][1])},
            ignore_index=True)

    with open(analysis_output_file_address, 'w') as f:
        f.write(analysis_result)
    analysis_table.to_csv('Indications drugs combination analysis table.tsv', sep='\t', index=None)


def analyse_diseases():
    global analysis_result, analysis_output_file_address
    file_address = "Indications_diseases.tsv"
    analysis_output_file_address = "Indications diseases combination analysis.txt"
    database_combinations = get_combinations(file_address)
    # Databases with small count of diseases mapped(2,13,19,24) and DOID(9), EFO(10), HP(15), RepoDB(22) are ignored
    databases_to_inspect = [1, 3, 4, 5, 6, 7, 8, 11, 12, 14, 16, 17, 18, 20, 21, 23]

    analysis_result = ""
    analysis_table = pd.DataFrame({'K': [], 'Best Combination': [], 'Count of Indications Diseases': []})
    for k in range(4, 11):
        combs_stats = dict()
        combs = combinations(databases_to_inspect, k)
        for comb in combs:
            combs_stats[comb] = get_count_common_in_databases(set(comb), database_combinations)
        sorted_combs_stats = {k: v for k, v in sorted(combs_stats.items(), key=lambda item: item[1], reverse=True)}
        analysis_result += f"k={k}: " + str(list(sorted_combs_stats.items())[:20]) + "\n"
        disease_best_combinations[k - 4] = list(sorted_combs_stats.items())[0][0]
        print(f"k={k}: ", end="")
        print(list(sorted_combs_stats.items())[:10])
        analysis_table = analysis_table.append(
            {'K': str(k), 'Best Combination': str(list(sorted_combs_stats.items())[0][0]),
             'Count of Indications Diseases': str(list(sorted_combs_stats.items())[0][1])},
            ignore_index=True)
    with open(analysis_output_file_address, 'w') as f:
        f.write(analysis_result)
    analysis_table.to_csv('Indications diseases combination analysis table.tsv', sep='\t', index=None)


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
            count_unique_drugs, count_unique_diseases, count = get_count_indications_in_best_combinations(
                set(disease_best_combinations[j - 4]),
                set(drug_best_combinations[i - 4]), i, j)
            col_counts.append(f"{count} ({count_unique_diseases},{count_unique_drugs})")
        df[f"K={i}"] = col_counts
    df.index = ["K=4", "K=5", "K=6", "K=7", "K=8", "K=9", "K=10"]
    df.to_csv("Indications_combinations_results.tsv", sep="\t")
    print(df)
