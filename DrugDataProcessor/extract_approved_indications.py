import pandas as pd


def find_indication_pairs():
    raw_df=pd.read_csv('Data/RepoDB/repodb.txt',sep=',',dtype=(str,str))
    parsed_df=pd.read_csv('Data/RepoDB/parsed_RepoDB.tsv',sep='\t',dtype=str)
    indication_df=pd.DataFrame({"DiseaseID":[],"DrugCID":[]})
    disease_IDs=[]
    drug_CIDs=[]

    indication_pairs=dict()
    for ind,row in raw_df.iterrows():
        if row['status']=='Approved':
            try:
                indication_pairs[row['drug_name']].append(row['ind_id'])
            except:
                indication_pairs[row['drug_name']]=[row['ind_id']]
    for drug_name in indication_pairs:
        # print(drug_name)
        query=parsed_df.query(f"DrugName=='{drug_name}'")
        for CID in query['CID']:
            for disease_ID in indication_pairs[drug_name]:
                disease_IDs.append(disease_ID)
                drug_CIDs.append(CID)
    indication_df['DiseaseID']=disease_IDs
    indication_df['DrugCID']=drug_CIDs
    # print(len(set(indication_df['DiseaseID'])))
    # print(len(set(indication_df['DrugCID'])))
    indication_df=indication_df.set_index('DiseaseID').sort_index().dropna()

    indication_df.to_csv("indication_pairs.tsv",sep='\t')
    # print(indication_pairs)
    # print(parsed_df)



if __name__ == '__main__':
    find_indication_pairs()