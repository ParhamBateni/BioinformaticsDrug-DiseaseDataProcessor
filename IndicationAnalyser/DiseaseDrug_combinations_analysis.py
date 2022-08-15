import pandas as pd
from itertools import combinations

global analysis_result,analysis_output_file_address
analysis_result=""

def get_combinations(file_address):
    # First exclude all rows with nan value for Combinations column(like Sum row) and then call this function!
    df=pd.read_csv(file_address, sep='\t',dtype=str).fillna(0)
    return [set(list(map(int,str(comb).split(",")))) for comb in list(df["Combinations"])]

def get_count_common_in_databases(comb,database_combinations):
    count=sum([1 for d in database_combinations if comb.issubset(d)])
    return count


def analyse_drugs():
    global analysis_result,analysis_output_file_address
    file_address = "../DrugDataProcessor/drug_results.tsv"
    analysis_output_file_address = "Drugs combination analysis.txt"
    database_combinations = get_combinations(file_address)
    #Databases with small count of drugs mapped(2,14) and RepoDB(15) are ignored
    databases_to_inspect=[1,3,4,5,6,7,8,9,10,11,12,13,16,17]

    for k in range(4, 11):
        combs_stats = dict()
        combs = combinations(databases_to_inspect, k)
        for comb in combs:
            combs_stats[comb] = get_count_common_in_databases(set(comb), database_combinations)
        sorted_combs_stats = {k: v for k, v in sorted(combs_stats.items(), key=lambda item: item[1], reverse=True)}
        analysis_result += f"k={k}: " + str(list(sorted_combs_stats.items())[:20]) + "\n"
        print(f"k={k}: ", end="")
        print(list(sorted_combs_stats.items())[:10])

def analyse_diseases():
    global analysis_result,analysis_output_file_address
    file_address = "../DiseaseDataProcessor/disease_results.tsv"
    analysis_output_file_address = "Diseases combination analysis.txt"
    database_combinations = get_combinations(file_address)
    #Databases with small count of diseases mapped(1,11,17,22) and RepoDB(20) are ignored
    databases_to_inspect=[2,3,4,5,6,7,8,9,10,12,13,14,15,16,18,19,21]


    for k in range(4, 11):
        combs_stats = dict()
        combs = combinations(databases_to_inspect, k)
        for comb in combs:
            combs_stats[comb] = get_count_common_in_databases(set(comb), database_combinations)
        sorted_combs_stats = {k: v for k, v in sorted(combs_stats.items(), key=lambda item: item[1], reverse=True)}
        analysis_result += f"k={k}: " + str(list(sorted_combs_stats.items())[:20]) + "\n"
        print(f"k={k}: ", end="")
        print(list(sorted_combs_stats.items())[:10])



if __name__ == '__main__':
    while True:
        print("Drug or Disease?")
        inp=input()
        if inp=="Drug":
            analyse_drugs()
            break
        elif inp=="Disease":
            analyse_diseases()
            break
        else:
            print("Try again!")
    with open(analysis_output_file_address,'w') as f:
        f.write(analysis_result)