#<codecell1>
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import re
import csv
import numpy as np

df= pd.read_csv('/Users/romina/Documentos/UNIVERSIDAD/tesis/data_codes/gene_description.csv', header=0, names=['genes', 'description'])
genes_list= dict(zip(list(df.genes), list(df.description)))
df.drop('description', axis='columns', inplace= True)
print(df)

def data_separator(data_string):
    data_list=[]
    res_1=data_string.split(' ')
    data_list.append(res_1[0])
    res_2=re.findall(r'\[.*?\]', data_string)
    for value in res_2:
        data_list.append(value)
    return (data_list)

for key, value in genes_list.items():
    new_value=data_separator(value)
    genes_list[key]=new_value

#<codecell2>
def value_to_sub(df, dict, column_tag, row_tag):
    row_index=list(dict.keys())
    key_name=row_index[row_tag]
    if column_tag=='gene_id':
        substring= 'lcl'
    else:
        substring = column_tag
    for value in dict[key_name]:
        if substring in value:
            if substring == 'protein':
                if substring in value and 'protein_id' not in value:
                    value_to_sub = value
                else:
                    pass
            else:
                value_to_sub=value
        else:
            value_to_sub= 'None'
        if value_to_sub != 'None':
            break
    return value_to_sub


#df['gene_id']= np.NaN

#for i, row in df.iterrows():
    #sub_with= value_to_sub(df, genes_list, 'gene_id', i)
    #df.at[i, 'gene_id']= str(sub_with)


#<codecell3>

column_names=['gene_id', 'gene', 'locus_tag', 'db_xref','protein', 'protein_id', 'location', 'gbkey']

for name in column_names:
    for i, row in df.iterrows():
        sub_with= value_to_sub(df, genes_list,name, i)
        df.at[i, name]= str(sub_with)

print(df)

df.to_csv('/Users/romina/Documentos/UNIVERSIDAD/tesis/data_codes/gene_table_clean_2.csv')
