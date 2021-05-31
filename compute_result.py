import pandas as pd
import os
import numpy as np
import joblib
import pickle
import os
import scipy.stats as st

#Issue for arguments in both the functions should be resolved
#Issue for converting normal gene expression file into GSVA should be resolved

# @app.route('/gexp_to_testData')  # ask where to put the .csv file and smiles list
# def gexp_smile_to_testSet(gene_expression, smiles_list):  #GSVA Gene expression dataframe and SMILES in form of a list
def gexp_smile_to_testSet():
    os.chdir('C:/Users/techuser/Downloads/Oncopretr-Git/Scripts')
    # df_gsva = pd.read_csv(gene_expression)
    df_gsva = pd.read_csv('temp_small.csv')
    smiles_list = ['COC1=C(C=C2C(=C1)N=CN=C2NC3=CC(=C(C=C3)F)Cl)OCCCN4CCOCC4', 'C1COCCN1C2=CC(=O)C=C(O2)C3=C4C(=CC=C3)SC5=CC=CC=C5S4','CC(C1=CC=C(C=C1)C(=O)NC2=C3C=CNC3=NC=C2)N.Cl','C1CCC(C1)C2=C(C=CC(=C2)C3=CNC4=NC=C(C=C34)C5=CC=CC=C5)C(=O)O','CC(C1=C(C=CC(=C1Cl)F)Cl)OC2=C(N=CC(=C2)C3=CN(N=C3)C4CCNCC4)N']

    df_gsva = df_gsva.dropna(axis = 1)
    if 'Unnamed: 0' in df_gsva.columns:
        df_gsva = df_gsva.drop(['Unnamed: 0'], axis = 1)
    
    smile_file = 'C:/Users/techuser/Downloads/Oncopretr-Git/SMILESVecProteinRepresentation/source/utils/smile_samples.txt'
    
    smiles = open(smile_file,'w')

    for element in smiles_list:
        smiles.write(element)
        smiles.write('\n')
    smiles.close()
    
    os.chdir('C:/Users/techuser/Downloads/Oncopretr-Git/SMILESVecProteinRepresentation/source')
    # %run getsmilesvec.py
    os.system('getsmilesvec.py')
    
    
    with open("smiles.vec", 'rb') as f:
        prots= pickle.load(f, encoding='bytes')
        
    for i,val in enumerate(smiles_list):
        prots[i].insert(0,val)
    
    col = ['SMILE']
    for x in range(1,101,1):
        col.append(f'X{x}')
        
    df_smileVec = pd.DataFrame(prots, columns = col)
    
    #temporary keys for cross join
    
    df_gsva['key'] = 1
    df_smileVec['key'] = 1
    
    df_test = df_gsva.merge(df_smileVec, on = 'key').drop('key', axis = 1)
    
    return df_test

# def compute_result(test_dataset):
def compute_result():
    test_dataset = gexp_smile_to_testSet()
     #making a list of 10 DataFrames for getting CI, Min, Max
    smiles = test_dataset['SMILE']
    main_test_dataset = test_dataset.drop(['SMILE'], axis = 1)
    os.chdir("C:/Users/techuser/Downloads/Oncopretr-Git/SVR Models")

    samples = []
    for i in range(len(main_test_dataset)):
        samples.append(f'Sample_{i}')

    df_results_list = []
    for i in range(10):
        model = joblib.load(f"model_svr_{i}.sav")
    #     print(f'Model Loaded {i}')
        z_score = model.predict(main_test_dataset)
    #     print(f'Predictions done {i}')
        df_result = pd.DataFrame(zip(samples, smiles, z_score), columns = ['Sample','SMILES','Z-Score'])
        df_results_list.append(df_result)

    print('Results List Saved')    

    # df = df_result.groupby(['SMILES']) #some changes are required at this level, try to add samples as a column too
    # df = df.agg(MIN_ZSCORE=('Z-Score', np.min))
    # df = df.reset_index()
    # df = df.sort_values('MIN_ZSCORE')

    #     return df_result.to_json()
   
    #computing results in form of a dataframe
    empty_array = np.zeros((len(samples), 6))

    ci = []
    minm = []
    maxm = []
    avgz = []
    sml = []

    for i,s in enumerate(samples):
        z = []

        for d in df_results_list:
            row = d.loc[d['Sample'] == s]
            z.append(row['Z-Score'].values[0])
            t = row['SMILES'].values[0]

        z = np.array(z)  
        ci.append(st.norm.interval(alpha=0.95, loc=np.mean(z), scale=st.sem(z)))
        minm.append(min(z))
        maxm.append(max(z))
        avgz.append(np.mean(z))  
        sml.append(t)

    df_final = pd.DataFrame(zip(samples, sml, ci, minm, maxm, avgz), columns = ['Sample','SMILE','CI','MAX','MIN','AVG'])
   
    print('Final List Computed')
    return df_final.to_json()

