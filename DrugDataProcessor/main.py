import os

import pandas as pd

data_folder_name = "Data"
global result_df
result_df = pd.DataFrame({})
unmapped_drugs = set()


def update_df(new_df, result_df):
    # Removing duplicate indexes in the new dataframe
    is_duplicate = new_df.index.duplicated(keep="first")
    not_duplicate = ~is_duplicate
    new_df = new_df[not_duplicate]

    updated_df = pd.concat([result_df, new_df], axis=0)
    updated_df = updated_df.groupby(updated_df.index, sort=False)[updated_df.columns.tolist()].first()

    return updated_df


def process(file_name, folder_name, file_address):
    global unmapped_drugs
    print(f"Reading {file_name} from {folder_name} database")
    df = pd.read_csv(file_address, sep='\t', dtype=('str', 'str'))[["CID", "DrugName"]]
    count_total_names = len(df.index)

    unmapped_names = df[df["CID"].isna()]["DrugName"].tolist()
    unmapped_drugs = unmapped_drugs.union(set(unmapped_names))

    df.loc[df['DrugName'].isna(), 'DrugName'] = 'Not available'
    df.dropna()
    df = df.set_index('CID')
    df[folder_name] = 1
    df.loc['CountOfMappedNames', folder_name] = count_total_names - len(unmapped_names)
    df.loc['CountOfTotalNames', folder_name] = count_total_names
    return df


if __name__ == '__main__':
    for folder in os.listdir(data_folder_name):
        for file in os.listdir(f"{data_folder_name}/{folder}"):
            if file.startswith("parsed") and file.endswith(".tsv"):
                file_address = f"{data_folder_name}/{folder}/{file}"
                result_df = update_df(process(file, folder, file_address), result_df)
    result_df = result_df.fillna(0)
    result_df['Sum'] = result_df.loc[:, result_df.columns != 'DrugName'].sum(axis=1)
    result_df = result_df.sort_index()
    comb_df = result_df.reset_index()[:-2]
    for i in range(len(comb_df.columns) - 3):
        if i == 0:
            comb_df['Combinations'] = comb_df[comb_df.columns[2]].astype(int) * pd.Series([',1'] * len(comb_df.index))
        else:
            comb_df['Combinations'] = comb_df['Combinations'] + comb_df[comb_df.columns[2 + i]].astype(int) * pd.Series(
                [f',{i + 1}'] * len(comb_df.index))
    comb_df['Combinations'] = comb_df['Combinations'].str.extract(',(.*)')
    result_df = comb_df.append(result_df.iloc[-2:].reset_index()).set_index('CID')
    result_df = result_df.iloc[-2:].append(result_df.iloc[:-2])
    print('Writing the result into drug_results.tsv')
    result_df.to_csv("drug_results.tsv", sep='\t')
    count_parts = 5
    part_size = len(result_df.index) // count_parts
    for i in range(count_parts):
        if i != count_parts - 1:
            result_df.iloc[i * part_size:(i + 1) * part_size].to_csv(f'drug_results_part{i + 1}.tsv', sep='\t')
        else:
            result_df.iloc[i * part_size:].to_csv(f'drug_results_part{i + 1}.tsv', sep='\t')
    pd.DataFrame({"DrugName": list(unmapped_drugs)}).set_index("DrugName").to_csv("unmapped_drugs.tsv", sep='\t')
