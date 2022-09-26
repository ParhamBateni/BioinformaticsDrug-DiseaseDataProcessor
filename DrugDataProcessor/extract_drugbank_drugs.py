import pandas as pd

drug_names = []
df = pd.DataFrame({"DrugName": []})
with open("Data/DrugBank/full database.xml", 'r', encoding='utf-8') as f:
    for line in f:
        if line.startswith("  <name>"):
            drug_names.append(line[line.index("  <name>") + 8:line.rindex("</name>")])
df['DrugName'] = drug_names
df.set_index('DrugName').sort_index().to_csv('Data/DrugBank/processed_DrugBank.txt', sep='\t')
