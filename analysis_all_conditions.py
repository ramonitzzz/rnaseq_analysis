#<codecell1>
import pandas as pd
import numpy as np

class Experiment_condition:

    def __init__(self, path, condition, cut_point):
        self.path= path
        self.condition = condition
        self.cut_point= cut_point
        self.df =pd.read_csv(self.path)
        self.__significant_values = self.df[(self.df.FDR <0.05)]

    def significant_values(self):
        return self.__significant_values
        #df=pd.read_csv(self.path)
        #significant_df=df[(df.FDR <0.05)]
        #return significant_df

#class Cutoff(Experiment_condition):
    def cutoff_by_fold_change(self):

        by_fold_change=self.__significant_values[(self.__significant_values.logFC >= self.cut_point) | (self.__significant_values.logFC <= (- self.cut_point))].reset_index()
        return by_fold_change

    def x_regulated (self, ur_or_dr):
        if ur_or_dr == 'ur':
            x_regulated_genes= self.__significant_values[(self.__significant_values.logFC >0)].reset_index()
            by_cut_point= x_regulated_genes[(x_regulated_genes.logFC > self.cut_point)].reset_index()
        elif ur_or_dr == 'dr':
            x_regulated_genes = self.__significant_values[(self.__significant_values.logFC <0)].reset_index()
            by_cut_point= x_regulated_genes[(x_regulated_genes.logFC < -self.cut_point)].reset_index()

        return by_cut_point



#<condition_6>
df_6='/Users/romina/Documentos/UNIVERSIDAD/tesis/data_codes/conditions_no_outlier/6_wt-pza50_pnca-pza50_no_outliers.csv'

condition_6=Experiment_condition(df_6, '6_wt-pza50_pnca-pza50', 2.5)
condition_6_cutoff=condition_6.cutoff_by_fold_change()
print(condition_6_cutoff)
ur_genes_6=condition_6.x_regulated('ur')
print(ur_genes_6)
dr_genes_6=condition_6.x_regulated('dr')
print(dr_genes_6)
