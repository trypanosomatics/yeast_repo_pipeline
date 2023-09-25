#pipeline_functions
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

def read_dataframe_columns(file_name,sep,columns=0):
    '''
    Read a data table from a CSV file and return a dataframe 
    containing the columns specified as a parameter in list format
    '''
    data_raw = pd.read_csv(file_name,sep=sep)
    if columns==0:
        data_filtered = data_raw
    else: 
        data_filtered = data_raw[columns]
    return (data_filtered)


def data_filtering(df,column,criteria):
    '''
    Receives a dataframe and a filtering criteria, 
    returns the dataframe only with the parameters that match the criteria.    '''
    columns = df.columns
    for i in range(len(columns)):
        is_criteria = df[column] == criteria
        df_filtered = df[is_criteria]
    return (df_filtered)


def data_intersection(data,selection_data,column_df,column_selection_data,filter_name,data_type):
    '''
    Receives two dataframes, one with the data to analyze and another to select, 
    and the columns to use as filtering criteria returning the filtered df
    '''
    df = data[column_df].unique().tolist()
    df_filtering = selection_data[column_selection_data].unique().tolist()
    intersected_column = []
    not_intersected_column = []
    filter_name_column = []
    intersected_data_raw = []
    
    for i in range(len(df)):
        try:
            if str(df[i]) in df_filtering:
                intersected_column.append(df[i]) 
            else: 
                not_intersected_column.append(df[i]) 
                filter_name_column.append(filter_name) 
        except:
            not_intersected_column.append('N/A') 
            filter_name_column.append('N/A')
            pass
            
    intersected_data_raw = pd.DataFrame(columns=[data_type])
    intersected_data_raw[data_type] = intersected_column
    
    not_intersected_data_raw = pd.DataFrame(columns=[data_type, 'Filter Name'])
    not_intersected_data_raw[data_type] = not_intersected_column
    not_intersected_data_raw['Filter Name'] = filter_name_column
    
    pre_intersected_data = pd.merge(left=data,right=intersected_data_raw,how="inner",left_on=data_type, right_on=data_type)
    intersected_data = pd.merge(left=selection_data,right=pre_intersected_data,how="inner",left_on=column_selection_data, right_on=data_type)
   
    return (intersected_data, not_intersected_data_raw)

def data_not_in_intersection(data,selection_data,column_df,column_selection_data,filter_name):
    '''
    The function is designed to filter a DataFrame (data) by excluding rows 
    that have matching values in a specified column of the DataFrame 
    with a selection DataFrame (selection_data). 
    '''
    df_novel= pd.DataFrame()
    df_tested= pd.DataFrame()
    df_resultado = data.merge(selection_data, left_on=column_df, right_on=column_selection_data, how='left', indicator=True)
    df_resultado = df_resultado.drop(columns=['Unnamed: 0', 'smiles', 'mol', 'inchi', 'inchikey'])
    df_novel = df_resultado[df_resultado['_merge'] == 'left_only']
    df_tested = df_resultado[df_resultado['_merge'] == 'both']
    df_novel = df_novel.drop(columns=['_merge'])
    df_tested = df_tested.drop(columns=['_merge'])
    df_tested['Filter Name'] = filter_name
    return (df_novel, df_tested)

def resumen_contar_datos(df,column):
    '''
    Receives a dataframe and a filter criterion,
    returning a table with the count of each parameter
    '''
    columns = df.columns
    for i in range(len(columns)):
        n_datos = df[columns[i]].value_counts()
        print(n_datos)
    return (print(n_datos))
     
   
def smiles_to_inchikey(df,compound_column):
    '''
    Receives a dataframe and the column where the smiles of the compounds are located,
    returning the same df with a new column, with the compounds in inchikey format    '''
    pd.options.mode.chained_assignment = None
    mol_tested_list= []
    for element in df[compound_column]:
        mol_tested = Chem.MolFromSmiles(element)
        mol_tested_list.append(mol_tested)
    df['mol'] = mol_tested_list
    inchikey_list = []
    for element in df['mol']:
        try:
            inchikey = Chem.MolToInchiKey(element)
            inchikey_list.append(inchikey)
        except:
            inchikey_list.append('N/A')
            pass
    df['inchiKey'] = inchikey_list
    return(df)

def drug_likness(df,compound_column):
    pd.options.mode.chained_assignment = None
    mol_tested_list= []
    for element in df[compound_column]:
        mol_tested = Chem.MolFromSmiles(element)
        mol_tested_list.append(mol_tested)
    df['mol'] = mol_tested_list
    NumHDonors_list = []
    NumHAcceptors_list = []
    MW_list = []
    LogP_list = []
    for element in df['mol']:
        try:
            NumHDonors = Descriptors.NumHDonors(element)
            NumHDonors_list.append(NumHDonors)
        except:
            NumHDonors_list.append('N/A')
            pass
        try:
            NumHAcceptors = Descriptors.NumHAcceptors(element)
            NumHAcceptors_list.append(NumHAcceptors)
        except:
            NumHAcceptors_list.append('N/A')
            pass
        try:
            MW = Descriptors.rdMolDescriptors.CalcExactMolWt(element)
            MW_list.append(MW)
        except:
            MW_list.append('N/A')
            pass
        try:
            LogP = Descriptors.rdMolDescriptors.CalcCrippenDescriptors(element)[0]
            LogP_list.append(LogP)
        except:
            LogP_list.append('N/A')
            pass
    df['NumHDonors'] = NumHDonors_list
    df['NumHAcceptors'] = NumHAcceptors_list
    df['MW'] = MW_list
    df['logP'] = LogP_list
    #Cuento cuantas caracterísicas cumple cada compuesto
    countLipinski = lambda row: int(row['NumHDonors'] < 6) +  int(row['NumHAcceptors'] < 6) + int(row['MW'] < 500) + int(row['logP'] < 6)
    df['Lipinski'] = df.apply(countLipinski,axis=1)
    countRO3 = lambda row: int(row['NumHDonors']<= 3) +  int(row['NumHAcceptors']<= 3) + int(row['MW']< 300) + int(row['logP']<= 3)
    df['RO3'] = df.apply(countRO3,axis=1)
    #Selecciono aquellos compuestos que cumplen con 3 o más
    gdi_aviable_novel_druglike_lipinski = df[df['Lipinski']>2]
    gdi_aviable_novel_druglike_RO3 = df[df['RO3']>2]
    gdi_aviable_novel_druglike_not_lipinski = df[df['Lipinski']<3]
    gdi_aviable_novel_druglike_not_RO3 = df[df['RO3']<3]
    #Filtrado según drug like
    df_druglike = pd.concat([gdi_aviable_novel_druglike_lipinski,gdi_aviable_novel_druglike_RO3])
    df_not_druglike = pd.concat([gdi_aviable_novel_druglike_not_lipinski,gdi_aviable_novel_druglike_not_RO3])
    df_not_druglike['Filter Name'] = ["Not Druglike"] * len(df_not_druglike)
    return(df_druglike, df_not_druglike)

    