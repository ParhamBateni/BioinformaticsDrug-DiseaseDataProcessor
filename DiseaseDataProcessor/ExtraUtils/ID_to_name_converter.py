import pandas as pd
import json
df=pd.read_csv("diseases.csv", index_col=[0])
unmapped_ids_df=pd.read_csv("../unmapped_ids.tsv", sep='\t')
unmapped_ids_df=unmapped_ids_df.set_index("DiseaseID")
with open("../UMLS_ID_converter.json", 'r') as f:
    data=f.read()
UMLS_id_converter=json.loads(data)
results_df=pd.read_csv("../results_final.tsv", sep="\t")
results_df=results_df.set_index("DiseaseID")
disease_names=[]
indexes_not_supported=[]
for index,row in df.iterrows():
    disease_id=row["disease"][9:]
    disease_name=''
    try:
        names_hash=dict()
        names=[]
        for ref_id in UMLS_id_converter[disease_id]:
            name=results_df.loc[ref_id,"DiseaseName"]
            if names_hash.get(results_df.loc[ref_id,"DiseaseName"].upper()) is None:
                names_hash[name.upper()]=1
                names.append(name)
        disease_name="/ ".join(names)
    except Exception as e:
        try:
            disease_name=unmapped_ids_df.loc[disease_id]["DiseaseName"]
        except:
            indexes_not_supported.append(index)
    disease_names.append(disease_name.strip())
df["disease_name"]=pd.Series(disease_names)
df.to_csv("updated_diseases.csv")
print(indexes_not_supported)
print(len(indexes_not_supported))
