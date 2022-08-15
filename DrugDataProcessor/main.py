import os
import pandas as pd
data_folder_name="Data"
global result_df
result_df=pd.DataFrame({})
unmapped_drugs=set()
def update_df(new_df,result_df):
    #Removing duplicate indexes in the new dataframe
    is_duplicate = new_df.index.duplicated(keep="first")
    not_duplicate = ~is_duplicate
    new_df = new_df[not_duplicate]


    updated_df=pd.concat([result_df,new_df],axis=0)
    updated_df=updated_df.groupby(updated_df.index,sort=False)[updated_df.columns.tolist()].first()

    return updated_df

def process(file_name,folder_name,file_address):
    global unmapped_drugs
    print(f"Reading {file_name} from {folder_name} database")
    df=pd.read_csv(file_address,sep='\t',dtype=('str','str'))[["CID","DrugName"]]
    count_total_names=len(df.index)

    unmapped_names = df[df["CID"].isna()]["DrugName"].tolist()
    unmapped_drugs=unmapped_drugs.union(set(unmapped_names))


    df=df.dropna()
    df=df.set_index('CID')
    df[folder_name] = 1
    df.loc['CountOfMappedNames', folder_name] = count_total_names-len(unmapped_names)
    df.loc['CountOfTotalNames', folder_name] = count_total_names
    return df
if __name__ == '__main__':
    counter=0
    for folder in os.listdir(data_folder_name):
        for file in os.listdir(f"{data_folder_name}/{folder}"):
            # if counter==6:continue
            if file.startswith("parsed") and file.endswith(".tsv"):
                counter+=1
                file_address=f"{data_folder_name}/{folder}/{file}"
                result_df = update_df(process(file,folder,file_address), result_df)
    result_df = result_df.fillna(0)
    result_df['Sum'] = result_df.loc[:, result_df.columns != 'DrugName'].sum(axis=1)
    result_df = result_df.sort_index()
    combinations = []
    i=0
    for ind_, row in result_df.iterrows():
        combination = ""
        j = 0
        for column in result_df.columns:
            if column != "DrugName" and column != 'Sum' and column!='CID':
                if str(result_df[column].iloc[i]) == '1.0':
                    combination += str(j) + ","
            j += 1
        i += 1
        combinations.append(combination[:-1])
    result_df['Combinations'] = combinations
    indexes_except_sum_and_count = [i for i in result_df.index if  i != 'CountOfMappedNames' and i !="CountOfTotalNames"]
    result_df = result_df.loc[['CountOfMappedNames'] + ['CountOfTotalNames']+indexes_except_sum_and_count ]
    result_df.to_csv("drug_results.tsv", sep='\t')
    pd.DataFrame({"DrugName": list(unmapped_drugs)}).set_index("DrugName").to_csv("unmapped_drugs.tsv", sep='\t')