#<codecell1>
import pandas as pd
import numpy as np

class Experiment_condition:

    def __init__(self, path, condition, cut_point):
        self.path= path
        self.condition = condition
        self.cut_point= cut_point
        self.df =pd.read_csv(self.path, index_col='genes')
        self.df.drop(self.df.columns[self.df.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
        self.__significant_values = self.df[(self.df.FDR <0.05)]

    def significant_values(self):
        return self.__significant_values

    def cutoff_by_fold_change(self):
        by_fold_change=self.__significant_values[(self.__significant_values.logFC >= self.cut_point) | (self.__significant_values.logFC <= (- self.cut_point))].reset_index()
        by_fold_change['exp_condition']=self.condition
        return by_fold_change

    def x_regulated (self, ur_or_dr):
        if ur_or_dr == 'ur':
            x_regulated_genes= self.__significant_values[(self.__significant_values.logFC >0)].reset_index()
            by_cut_point= x_regulated_genes[(x_regulated_genes.logFC > self.cut_point)].reset_index()
        elif ur_or_dr == 'dr':
            x_regulated_genes = self.__significant_values[(self.__significant_values.logFC <0)].reset_index()
            by_cut_point= x_regulated_genes[(x_regulated_genes.logFC < -self.cut_point)].reset_index()

        return by_cut_point



#<condition_3>
df_3='/Users/romina/Documentos/UNIVERSIDAD/tesis/data_codes/conditions_no_outliers/3_wt-pza0_pnca-pza0.csv'

condition_3=Experiment_condition(df_3, '3_wt-pza0_pnca-pza0', 3)
condition_3_cutoff=condition_3.cutoff_by_fold_change()
print(condition_3_cutoff)
ur_genes_3=condition_3.x_regulated('ur')
#print(ur_genes_3)
dr_genes_3=condition_3.x_regulated('dr')
#print(dr_genes_3)

#<condition_4>
df_4='/Users/romina/Documentos/UNIVERSIDAD/tesis/data_codes/conditions_no_outliers/4_wt-pza0_pnca-pza50.csv'

condition_4=Experiment_condition(df_4, '4_wt-pza0_pnca-pza50', 3)
condition_4_cutoff=condition_4.cutoff_by_fold_change()
print(condition_4_cutoff)
ur_genes_4=condition_4.x_regulated('ur')
#print(ur_genes_4)
dr_genes_4=condition_4.x_regulated('dr')
#print(dr_genes_4)

#<condition_5>
df_5='/Users/romina/Documentos/UNIVERSIDAD/tesis/data_codes/conditions_no_outliers/5_pnca-pza0_wt-pza50.csv'

condition_5=Experiment_condition(df_5, '5_pnca-pza0_wt-pza50', 3)
condition_5_cutoff=condition_5.cutoff_by_fold_change()
print(condition_5_cutoff)
ur_genes_5=condition_5.x_regulated('ur')
#print(ur_genes_5)
dr_genes_5=condition_5.x_regulated('dr')
#print(dr_genes_5)

#<condition_6>
df_6='/Users/romina/Documentos/UNIVERSIDAD/tesis/data_codes/conditions_no_outliers/6_wt-pza50_pnca-pza50_no_outliers.csv'

condition_6=Experiment_condition(df_6, '6_wt-pza50_pnca-pza50', 3)
condition_6_cutoff=condition_6.cutoff_by_fold_change()
print(condition_6_cutoff)
ur_genes_6=condition_6.x_regulated('ur')
#print(ur_genes_6)
dr_genes_6=condition_6.x_regulated('dr')
#print(dr_genes_6)

#<all_conditions>
all_conditions= pd.concat([condition_3_cutoff,condition_4_cutoff,condition_5_cutoff,condition_6_cutoff])
all_conditions_clean= all_conditions.drop(columns=['logCPM', 'F', 'PValue'])
#print(all_conditions_clean)
all_conditions_pivot= all_conditions_clean.pivot(index='genes', columns='exp_condition', values='logFC')
print(all_conditions_pivot)

all_conditions_pivot.to_csv('/Users/romina/Documentos/UNIVERSIDAD/tesis/data_codes/conditions_no_outliers/all_conditions_pivot_FDR_by3.csv')
