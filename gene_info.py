#<cell1>
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd

df=pd.read_csv('/Users/romina/Documentos/UNIVERSIDAD/tesis/data_codes/conditions_no_outliers/all_conditions_pivot_FDR_by3.csv')
#df.info()
gene_names=list(df['genes'])

gene_description_dict={}
with open('/Users/romina/Documentos/UNIVERSIDAD/tesis/data_codes/LT708304_1.fasta') as filename:
    for record in SeqIO.parse(filename,'fasta'):
        for gene in gene_names:
            if gene in record.description:
                gene_description_dict[gene]=record.description

#for key, value in gene_description_dict.items():
    #d= '] '
    #s=[e+d for e in value.split(d) if e]
    #gene_description_dict[key]=s

description_df=pd.DataFrame.from_dict(gene_description_dict, orient='index', columns=['description'])

#print(description_df)

#clean_description_df=description_df['description'].str.split(('['), expand = True)
#description_df[['id,','locus_tag', 'protein', 'protein_id','location','gb_key']] =description_df.description.str.split(n=1, expand=True)

description_df.to_csv('/Users/romina/Documentos/UNIVERSIDAD/tesis/data_codes/gene_description.csv')
