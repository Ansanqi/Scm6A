import pandas as pd
import numpy as np
import pickle
import joblib
import os
from pathlib import Path
import shutil
from pandas.core.frame import DataFrame
import argparse
from sklearn.linear_model import LinearRegression
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.svm import SVR
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings("ignore", category=UserWarning)


#inputpath = '/Users/mac/Documents/An/Scm6A/GPB返修/3、与scm6A-seq比较/input.csv'
#outputpath = '/Users/mac/Documents/An/单细胞m6A数据库/Scm6A_database'

a = [i for i in range(2,43) if i not in [10, 24, 38]]
loc = os.path.dirname(os.path.abspath('Scm6A.py'))
#loc = '/Users/mac/Documents/An/单细胞m6A数据库/Scm6A'
filesavepath = loc + "/inputdata"

def read_file(inputpath):
    file_extension = inputpath.split(".")[-1].lower()  # 获取文件扩展名并转换为小写

    if file_extension == "csv":
        return pd.read_csv(inputpath, header=0, index_col=0)
    elif file_extension in ["xls", "xlsx"]:
        return pd.read_excel(inputpath, header=0, index_col=0)
    elif file_extension == "tsv":
        return pd.read_csv(inputpath, header=0, index_col=0, sep='\t')
    elif file_extension == "txt":
        return pd.read_csv(inputpath, header=0, index_col=0, sep='\t')
    else:
        print("不支持的文件格式，请提供csv、xls、xlsx、tsv或txt文件。")
        return None

def preprocess():
    RBPlist = pd.read_csv(loc +"/RBP.csv",header=0,index_col=0)
    #sample = pd.read_csv(inputpath , header = 0,index_col=0)
    sample = read_file(inputpath)

    filepath,fullflname = os.path.split(inputpath)
    fname,ext = os.path.splitext(fullflname)

    sample_RBP = pd.merge(RBPlist, sample, left_index=True, right_index=True,how='left')
    sort = sample_RBP.loc[RBPlist.index]
    sort.to_csv(loc + '/' + fname + '_RBP.csv', index = True, header = True)
    module = pd.read_csv(loc + "/RBP-module.csv")
    A = pd.merge(module, sort, on='RBP', how='left')

    if os.path.exists(filesavepath):
        shutil.rmtree(filesavepath)
        os.makedirs(filesavepath)
    else:
        os.makedirs(filesavepath)

    list0 = [f'ME{j}' for j in a]
    for i in list0:
        df = A[A["module"]==i]
        df = df.drop(df.columns[[1]], axis = 1)
        df.fillna(0, inplace=True)
        file = filesavepath + f"/{i}.csv"
        df.to_csv(file, index=False, header = True, mode='w')

    motif = pd.read_csv(loc + "/addinf.csv",header = 0,index_col=0).T
    data = pd.read_csv(filesavepath + "/ME2.csv", header=0, index_col=0)
    col = data.shape[1]
    
    for i in list0:
        b = pd.DataFrame(motif[i])
        MErep = pd.DataFrame(np.tile(b, (1, col)))
        MErep.index = motif.index
        data = pd.read_csv(filesavepath + f"/{i}.csv",header = 0,index_col=0)
        MErep.columns = data.columns
        bind = pd.concat([data, MErep])
        file2 = filesavepath + f"/{i}.csv"
        bind.to_csv(file2, index = True, header = True, mode='w')

def run_models():    
    if not os.path.exists(filesavepath):
        print(f'The folder {filesavepath} does not exist.')

    if not os.path.exists(outputpath):
        os.makedirs(outputpath)
        print(f'The folder {outputpath} does not exist. A new folder has been created in the {outputpath} you provided')
        
    #results_df = pd.DataFrame(columns=['module', 'm6A_site'])
    
    results_list = []
    
    for i in a:
        file_path = f"{filesavepath}/ME{i}.csv"
        X_input = pd.read_csv(file_path, index_col=0, header=0).T
        sc = StandardScaler()
        X_new = sc.fit_transform(X_input)
        folder_path = Path(loc + f"/model_save/ME{i}/")
        #folder_path = f'/Users/mac/pyProject/RBPtom6A/Stack_model/model_save2/ME{i}'
        pkl_files = [f for f in os.listdir(folder_path) if f.endswith('.pkl')]
        
        for pkl_file in pkl_files:
            file_path = os.path.join(folder_path, pkl_file)
            file_name, _ = os.path.splitext(pkl_file)
            
            with open(file_path, 'rb') as f:
                #best_models, stacked_model = joblib.load(f)
                xlf_lo = joblib.load(f)
                model_predictions = []
                predictions = xlf_lo.predict(X_new) # 模型使用
                model_predictions.extend(predictions)
            
           # for model_name, model in xlf_lo.items():
            #    predictions = model.predict(X_new)
            #    model_predictions.extend(predictions)
        
            prediction_dict = {'module': f'ME{i}', 'm6A_site': file_name}
            
            for col_name, prediction_value in zip(X_input.index, model_predictions):
                prediction_dict[col_name] = prediction_value
            results_list.append(prediction_dict)
            
    results_df = pd.DataFrame(results_list)
    
    file_name = os.path.splitext(os.path.basename(inputpath))[0]
   # results_df.to_csv(outputpath + f'/{file_name}_Scm6A_pred_results.csv', index=False)

    bed_df = pd.read_csv(loc +"/allpeaklist_bed.bed", sep='\t', header=None
                         ,usecols=[0, 1, 2, 3, 4, 5]
                         ,names=['chromosome', 'start', 'end','peak', 'gene', 'strand'])
    
    merged_df = results_df.merge(bed_df, left_on=results_df.columns[1], right_on='peak', how='left')
    merged_df = merged_df.drop("m6A_site", axis=1)
    last_six_columns = merged_df.columns[-6:]
    df1 = pd.concat([merged_df[last_six_columns], merged_df.drop(last_six_columns, axis=1)], axis=1)
    df1.to_csv(outputpath + f'/{file_name}_Scm6A_pred_results.csv', index=False)
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='manual to this script')
    parser.add_argument("-I","--input", required=True, type=str, default="0",help='Enter the absolute path to the folder where the input file')
    parser.add_argument("-O","--output", required=True, type=str, default="0",help='Enter the absolute path to the output file you want')
    args = parser.parse_args()
    inputpath = args.input
    outputpath = args.output
    preprocess()
    run_models()
    print('———————————— Successfully completed ————————————')
