from itertools import combinations

import pandas as pd

global analysis_result, analysis_output_file_address
analysis_result = ""


def get_combinations(file_address):
    # First exclude all rows with nan value for Combinations column(like Sum row) and then call this function!
    df = pd.read_csv(file_address, sep='\t', dtype=str).fillna(0)
    return [set(list(map(int, str(comb).split(",")))) for comb in
            list(df[df['Sum'].astype(float) >= 4]["Combinations"])]


def get_count_common_in_databases(comb, database_combinations):
    count = sum([1 for d in database_combinations if comb.issubset(d)])
    return count


def analyse_drugs():
    global analysis_result, analysis_output_file_address
    file_address = "../DrugDataProcessor/drug_results.tsv"
    analysis_output_file_address = "Drugs combination analysis.txt"
    database_combinations = get_combinations(file_address)
    # Databases with small count of drugs mapped(3,11,14,17,19) and RepoDB(18) are ignored
    databases_to_inspect = [1, 2, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 16, 20, 21, 22]

    analysis_table = pd.DataFrame({'K': [], 'Best Combination': [], 'Count of Drugs': []})
    for k in range(4, 11):
        combs_stats = dict()
        combs = combinations(databases_to_inspect, k)
        for comb in combs:
            combs_stats[comb] = get_count_common_in_databases(set(comb), database_combinations)
        sorted_combs_stats = {k: v for k, v in sorted(combs_stats.items(), key=lambda item: item[1], reverse=True)}
        analysis_result += f"k={k}: " + str(list(sorted_combs_stats.items())[:20]) + "\n"
        print(f"k={k}: ", end="")
        print(list(sorted_combs_stats.items())[:10])
        analysis_table = analysis_table.append(
            {'K': str(k), 'Best Combination': str(list(sorted_combs_stats.items())[0][0]),
             'Count of Drugs': str(list(sorted_combs_stats.items())[0][1])},
            ignore_index=True)
    analysis_table.to_csv('Drugs combination analysis table.tsv', sep='\t', index=None)


def analyse_diseases():
    global analysis_result, analysis_output_file_address
    file_address = "../DiseaseDataProcessor/disease_results.tsv"
    analysis_output_file_address = "Diseases combination analysis.txt"
    database_combinations = get_combinations(file_address)
    # Databases with small count of diseases mapped(2,13,19,24) and DOID(9), EFO(10), HP(15), RepoDB(22) are ignored
    databases_to_inspect = [1, 3, 4, 5, 6, 7, 8, 11, 12, 14, 16, 17, 18, 20, 21, 23]

    analysis_table = pd.DataFrame({'K': [], 'Best Combination': [], 'Count of Diseases': []})
    for k in range(4, 11):
        combs_stats = dict()
        combs = combinations(databases_to_inspect, k)
        for comb in combs:
            combs_stats[comb] = get_count_common_in_databases(set(comb), database_combinations)
        sorted_combs_stats = {k: v for k, v in sorted(combs_stats.items(), key=lambda item: item[1], reverse=True)}
        analysis_result += f"k={k}: " + str(list(sorted_combs_stats.items())[:20]) + "\n"
        print(f"k={k}: ", end="")
        print(list(sorted_combs_stats.items())[:10])
        analysis_table = analysis_table.append(
            {'K': str(k), 'Best Combination': str(list(sorted_combs_stats.items())[0][0]),
             'Count of Diseases': str(list(sorted_combs_stats.items())[0][1])},
            ignore_index=True)
    analysis_table.to_csv('Diseases combination analysis table.tsv', sep='\t', index=None)


if __name__ == '__main__':
    while True:
        print("Drug or Disease?")
        inp = input()
        if inp == "Drug":
            analyse_drugs()
            break
        elif inp == "Disease":
            analyse_diseases()
            break
        else:
            print("Try again!")
    with open(analysis_output_file_address, 'w') as f:
        f.write(analysis_result)
