# 这是一个示例 Python 脚本。

# 按 Shift+F10 执行或将其替换为您的代码。
# 按 双击 Shift 在所有地方搜索类、文件、工具窗口、操作和设置。




#help(TCRrep)
# from Levenshtein import distance as lev_distance
# from Bio.Align import PairwiseAligner
# from Bio.Align import substitution_matrices

# alpha_tcr = TCRrep(cell_df=alpha_sequences_HomoSapiens,
#                    organism='human',
#                    chains=['alpha'])
                  #  # 对于alpha链，使用cdr3_a_aa_col指定CDR3α序列列的列名
                  #v_a_gene_col='v.segm')  # 对于alpha链，使用v_a_gene_col指定V基因列的列名
# try:
#     with open(r'D:\pycharm\py_project\py_project1\mp_data_processing\combo_xcr_2024.tsv', 'r') as file:
#         print("File is accessible.")
# except Exception as e:
#     print("Error accessing the file:", e)
import pandas as pd
import numpy as np
# import pytest
# import tcrdist as td
# from tcrdist import mappers
from tcrdist.repertoire import TCRrep
import Levenshtein as lev
import seaborn as sns
import matplotlib.pyplot as plt
data = pd.read_excel('miniproject_002.xlsx')
#data.dropna(inplace=True)   #commented out, this would cause error

data_human = data[data['species'] == 'HomoSapiens']
data_mouse = data[data['species'] == 'MusMusculus']
data_rhesus = data[data['species'] == 'MacacaMulatta']
#print(data_rhesus)
alpha_data_human = data_human[data_human['gene'] == 'TRA']
beta_data_human = data_human[data_human['gene'] == 'TRB']

alpha_data_mouse = data_mouse[data_mouse['gene'] == 'TRA']
beta_data_mouse = data_mouse[data_mouse['gene'] == 'TRB']

alpha_data_rhesus = data_rhesus[data_rhesus['gene'] == 'TRA']  #focus on this line, it caused a problem, solved in this code
beta_data_rhesus = data_rhesus[data_rhesus['gene'] == 'TRB']
# print(alpha_data_rhesus)
# print(beta_data_rhesus)
new_cols_alpha = {'CDR3':'cdr3_a_aa', 'V': 'v_a_gene', 'J': 'j_a_gene'}
alpha_data_human = alpha_data_human.rename(columns = new_cols_alpha)
new_cols_beta = {'CDR3':'cdr3_b_aa', 'V': 'v_b_gene', 'J': 'j_b_gene'}
beta_data_human = beta_data_human.rename(columns = new_cols_beta)

new_cols_alpha = {'CDR3':'cdr3_a_aa', 'V': 'v_a_gene', 'J': 'j_a_gene'}
alpha_data_mouse = alpha_data_mouse.rename(columns = new_cols_alpha)
new_cols_beta = {'CDR3':'cdr3_b_aa', 'V': 'v_b_gene', 'J': 'j_b_gene'}
beta_data_mouse = beta_data_mouse.rename(columns = new_cols_beta)

new_cols_alpha = {'CDR3':'cdr3_a_aa', 'V': 'v_a_gene', 'J': 'j_a_gene'}
alpha_data_rhesus = alpha_data_rhesus.rename(columns = new_cols_alpha)
new_cols_beta = {'CDR3':'cdr3_b_aa', 'V': 'v_b_gene', 'J': 'j_b_gene'}
beta_data_rhesus = beta_data_rhesus.rename(columns = new_cols_beta)
print(alpha_data_rhesus) # for testing, no J and V information
datasets = {
    'human': {
        'alpha': alpha_data_human,
        'beta': beta_data_human,
    },
    'mouse': {
        'alpha': alpha_data_mouse,
        'beta': beta_data_mouse,
    },
    'rhesus': {
        #'alpha': alpha_data_rhesus, # this line must be commented out of the code, cuz in rhesus TRA data, V and J sequences are empty
        'beta': beta_data_rhesus,
    }
    # 其他物种和链类型
}

for species, chains_data in datasets.items():
    for chain_type, data in chains_data.items():
        # 根据链类型设置参数
        if chain_type == 'alpha':
            rep = TCRrep(cell_df=data,
                         organism=species,
                         chains=['alpha'],
                         db_file= r'D:\pycharm\py_project\py_project1\mp_data_processing\combo_xcr_2024.tsv')
        elif chain_type == 'beta':
            rep = TCRrep(cell_df=data,
                         organism=species,
                         chains=['beta'],
                         db_file= r'D:\pycharm\py_project\py_project1\mp_data_processing\combo_xcr_2024.tsv')

        # 在这里进行后续的处理，比如计算距离矩阵
        rep.compute_distances()

        # 导出距离矩阵或进行其他需要的分析
        distance_matrix = rep.pw_alpha if chain_type == 'alpha' else rep.pw_beta
        file_name = f"{species}_{chain_type}_distance_matrix.txt"
        np.savetxt(file_name, distance_matrix, fmt="%s")


        #pd.DataFrame(distance_matrix).to_csv(file_name)
#print(alpha_data_rhesus)
## np.savetxt('alpha_data_rhesus.txt', alpha_data_rhesus, fmt="%s") # for testing
new_cols_alpha = {'cdr3_a_aa':'cdr3_alpha_aa'}
alpha_data_human = alpha_data_human.rename(columns = new_cols_alpha)
new_cols_beta = {'cdr3_b_aa':'cdr3_beta_aa'}
beta_data_human = beta_data_human.rename(columns = new_cols_beta)
#print(sequences)
new_cols_alpha = {'cdr3_a_aa':'cdr3_alpha_aa'}
alpha_data_mouse = alpha_data_mouse.rename(columns = new_cols_alpha)
new_cols_beta = {'cdr3_b_aa':'cdr3_beta_aa'}
beta_data_mouse = beta_data_mouse.rename(columns = new_cols_beta)

new_cols_beta = {'cdr3_b_aa':'cdr3_beta_aa'}
beta_data_rhesus = beta_data_rhesus.rename(columns = new_cols_beta)

datasets_1 = {
    'human': {
        'alpha': alpha_data_human,
        'beta': beta_data_human,
    },
    'mouse': {
        'alpha': alpha_data_mouse,
        'beta': beta_data_mouse,
    },
    'rhesus': {
        #'alpha': alpha_data_rhesus, # this line must be commented out of the code, cuz in rhesus TRA data, V and J sequences are empty
        'beta': beta_data_rhesus,
    }
    # 其他物种和链类型
}
print(beta_data_human.cdr3_beta_aa)

for species, chains_data in datasets_1.items():
    for chain_type, data in chains_data.items():  # Use a different variable name here
        cdr3_column = f'cdr3_{chain_type.lower()}_aa'  # Dynamically construct column name
        if cdr3_column in data.columns:
            cdr3_lengths = data[cdr3_column].apply(len)

            plt.figure(figsize=(10, 6))
            sns.histplot(cdr3_lengths, kde=True)
            plt.title(f'Distribution of CDR3 Lengths - {species} {chain_type.capitalize()} Chains')
            plt.xlabel('CDR3 Length')
            plt.ylabel('Frequency')
            plt.show()
        else:
            print(f'Column {cdr3_column} not found in the dataset for {species} {chain_type} chains.')
# 计算距离矩阵

sequences = alpha_data_rhesus['cdr3_a_aa']
n = len(sequences)
distance_matrix_TRA_rhesus = np.zeros((n, n), dtype=int)
sequences = sequences.tolist()
for i in range(n):
    for j in range(n):
        distance_matrix_TRA_rhesus[i, j] = lev.distance(sequences[i], sequences[j])
print(distance_matrix_TRA_rhesus)
np.savetxt('distance_matrix_TRA_rhesus.txt', distance_matrix_TRA_rhesus, fmt="%s")
print(alpha_data_human)

# df_alpha_HomoSapiens.dropna(inplace=True)  # to remove null
# df_beta_HomoSapiens.dropna(inplace=True)  # to remove null

# alpha_rep = TCRrep(cell_df=alpha_data_human,
#                    organism='human',
#                    chains=['alpha'], # Use 'beta' here because tcrdist uses 'alpha' and 'beta' to distinguish chain types, not gene names
#                    db_file= r'D:\pycharm\py_project\py_project1\mp_data_processing\combo_xcr_2024.tsv')
#
# alpha_rep.compute_distances()
#
# beta_rep = TCRrep(cell_df=beta_data_human,
#                    organism='human',
#                    chains=['beta'], # Use 'beta' here because tcrdist uses 'alpha' and 'beta' to distinguish chain types, not gene names
#                    db_file=r'D:\pycharm\py_project\py_project1\mp_data_processing\combo_xcr_2024.tsv')
#
# beta_rep.compute_distances()
#
# distance_matrix_alpha_homosapiens = alpha_rep.pw_alpha
# distance_matrix_beta_homosapiens = beta_rep.pw_beta
# # ids = data['reference.id'].tolist()
# # distance_matrix_df_alpha_homosapiens = pd.DataFrame(distance_matrix_alpha_homosapiens, index=ids, columns=ids)
# np.savetxt("distance_matrix_alpha_homosapiens.txt", distance_matrix_alpha_homosapiens, fmt="%s")
# np.savetxt("distance_matrix_beta_homosapiens.txt", distance_matrix_beta_homosapiens, fmt="%s")
# alpha_sequences_HomoSapiens = pd.DataFrame(df_alpha_HomoSapiens[['cdr3_a_aa', 'v_a_gene']])
# #alpha_sequences = ["CASSLGTQYIY", "CASRAGSSYEQY", "CASSFGGGTDTQY"] # an example
# beta_sequences_HomoSapiens = pd.DataFrame(df_beta_HomoSapiens[['cdr3_a_aa', 'v_a_gene']])
# beta_sequences_HomoSapiens.rename(columns={'cdr3_a_aa': 'cdr3_b_aa', 'v_a_gene': 'v_b_gene'}, inplace=True)
# from tcrdist.vis_tools import bostock_cat_colors, cluster_viz
# def calculate_distance(pair):
#     seq1, seq2 = pair
#     return levenshtein_distance(seq1, seq2)

# Load the data to see the first few rows and understand its structure
#data = pd.read_excel('mini_project_data.xlsx')
# pd_df = pd.read_csv('vdjDB_PMID28636592.tsv', sep = "\t")
#
# t_df = td.mappers.vdjdb_to_tcrdist2(pd_df=pd_df)
# print(t_df)
# species_counts = t_df.organism.value_counts()
# print(species_counts)
# unique_species = t_df.organism.unique()
# print(unique_species)
# for species in unique_species:
#     # 过滤出当前物种的数据
#     index_species = t_df.organism == species
#     species_df = t_df.loc[index_species, :].copy()
#
#     # 将物种名称映射到'human'或'mouse'
#     if species == 'HomoSapiens':
#         organism_for_tcrrep = 'human'
#     elif species == 'MacacaMulatta':
#         organism_for_tcrrep = 'human'
#     elif species == 'MusMusculus':
#         organism_for_tcrrep = 'mouse'
#     else:
#         raise ValueError(f"Unknown species: {species}")
#
#     # 为当前物种创建TCRrep实例
#    # tr = TCRrep(cell_df=species_df, organism=organism_for_tcrrep)


# alpha_data = data[data['gene'] == 'TRA']
# beta_data = data[data['gene'] == 'TRB']
# seqs_alpha = alpha_data.tolist()
# print(seqs_alpha)
# print(seqs_alpha[1])
# print(seqs_alpha[0])
# seqs_beta = beta_data.tolist()
# seq_pairs_alpha = [(seqs_alpha[i], seqs_alpha[j]) for i in range(len(seqs_alpha)) for j in range(i+1, len(seqs_alpha))]
# seq_pairs_beta = [(seqs_beta[i], seqs_beta[j]) for i in range(len(seqs_beta)) for j in range(i+1, len(seqs_beta))]
#
# with Pool(processes=4) as pool:  # 根据机器性能调整进程数
#     distances_alpha = pool.map(calculate_distance, seq_pairs_alpha)
#
# with Pool(processes=4) as pool:  # 根据机器性能调整进程数
#     distances_beta = pool.map(calculate_distance, seq_pairs_beta)
# new_cols_alpha = {'cdr3':'cdr3_a_aa', 'v_gene': 'v_a_gene', 'j.segm': 'j_a_gene'}
# alpha_data = alpha_data.rename(columns = new_cols_alpha)
# # Setup for alpha chains
# alpha_rep = TCRrep(cell_df=alpha_data,
#                    organism='human',
#                    chains=['alpha'], # Use 'beta' here because tcrdist uses 'alpha' and 'beta' to distinguish chain types, not gene names
#                    db_file='combo_xcr_2024-03-05.tsv')
#
# # customize = lambda a : {new_cols_alpha.get(k, k):v for k,v in a.items()} #(5)
# # alpha_rep.metrics_a = customize(alpha_rep.metrics_a)
# # alpha_rep.weights_a = customize(alpha_rep.weights_a)
# # alpha_rep.kargs_a = customize(alpha_rep.kargs_a)
# # Compute the distance matrix for alpha chains
# alpha_rep.compute_distances()
#
# new_cols_beta = {'cdr3':'cdr3_b_aa', 'v_gene': 'v_b_gene', 'j.segm': 'j_b_gene'}
# # Setup for beta chains
# beta_data = beta_data.rename(columns = new_cols_beta)
# beta_rep = TCRrep(cell_df=beta_data,
#                   organism='human',
#                   chains=['beta'],
#                   # infer_all_genes=True,
#                   # infer_cdrs=False,  # (2)s
#                   # compute_distances=False,  # (3)
#                   # deduplicate=False,  # (4)
#                   db_file='alphabeta_gammadelta_db.tsv')
#
# # customize = lambda b : {new_cols_beta.get(k, k):v for k,v in b.items()} #(5)
# # beta_rep.metrics_b = customize(beta_rep.metrics_b)
# # beta_rep.weights_b = customize(beta_rep.weights_b)
# # beta_rep.kargs_b = customize(beta_rep.kargs_b)
# # Compute the distance matrix for beta chains
# beta_rep.compute_distances()
#
# # # 打印或保存距离矩阵
# # print("Alpha Chain Distance Matrix:")
# # print(alpha_distances)
# #
# # print("Beta Chain Distance Matrix:")
# # print(beta_distances)
#
# # 如果需要，可以将距离矩阵保存到文件
# # alpha_distances_HomoSapiens.to_csv('alpha_chain_distances_HomoSapiens.csv', index=False)
# # beta_distances_HomoSapiens.to_csv('beta_chain_distances_HomoSapiens.csv', index=False)
#
# # # 计算距离/相似性矩阵
# # alpha_matrix = create_distance_matrix(alpha_sequences)
# # beta_matrix = create_distance_matrix(beta_sequences)
# #
# # alpha_distance_matrix = calculate_distance_matrix(alpha_sequences)
# # x = alpha_matrix
# # np.savetxt(r'mat1.txt', x, fmt='%d', delimiter=',')
# # print(type(alpha_sequences))