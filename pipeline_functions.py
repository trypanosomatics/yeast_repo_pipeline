#pipeline_functions
import pandas 
from rdkit import Chem
from rdkit.Chem import Descriptors

def leer_tabla_columnas(nombre_archivo,sep,columnas=0):
    '''
    Lee una tabla de datos en csv y devuelve un dataframe 
    que contiene las columnas indicadas como parametro en formato lista
    '''
    data_raw = pandas.read_csv(nombre_archivo,sep=sep)
    if columnas==0:
        data_filtered = data_raw
    else: 
        data_filtered = data_raw[columnas]
    return (data_filtered)


def filtrar_datos(df,columna,criterio):
    '''
    Recibe un dataframe y un criterio de filtrado, 
    devolviendo el df solo con los parametros que coinciden con el filtro
    '''
    columnas = df.columns
    for i in range(len(columnas)):
        is_criterio = df[columna] == criterio
        df_filtered = df[is_criterio]
    return (df_filtered)


def interseccion_datos(data,selection_data,columna_df,columna_selection_data,filter_name,data_type):
    '''
    Recibe dos dataframes, uno con los datos a analizar y otro para seleccionar, y las columnas con las que se va a realizar el filtrado
    devolviendo el df filtrado
    '''
    df = data[columna_df].unique().tolist()
    df_filtering = selection_data[columna_selection_data].unique().tolist()
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
            
    intersected_data_raw = pandas.DataFrame(columns=[data_type])
    intersected_data_raw[data_type] = intersected_column
    
    not_intersected_data_raw = pandas.DataFrame(columns=[data_type, 'Filter Name'])
    not_intersected_data_raw[data_type] = not_intersected_column
    not_intersected_data_raw['Filter Name'] = filter_name_column
    
    pre_intersected_data = pandas.merge(left=data,right=intersected_data_raw,how="inner",left_on=data_type, right_on=data_type)
    intersected_data = pandas.merge(left=selection_data,right=pre_intersected_data,how="inner",left_on=columna_selection_data, right_on=data_type)
   
    return (intersected_data, not_intersected_data_raw)


def resumen_contar_datos(df,columna):
    '''
    Recibe un dataframe y un criterio de filtrado, 
    devolviendo una tabla con el conteo de cada parametro
    '''
    columnas = df.columns
    for i in range(len(columnas)):
        n_datos = df[columnas[i]].value_counts()
        print(n_datos)
    return (print(n_datos))
     
   
def smiles_to_inchikey(df,compound_columna):
    '''
    Recibe un dataframe y la columna donde se encuentran los smiles de los compuestos, 
    devolviendo el mismo df con una columna nueva, con los compuestos en formato inchikey
    '''
    pandas.options.mode.chained_assignment = None
    gdi_aviable_novel_drug_like_analysis = df
    mol_tested_list= []
    for element in df[compound_columna]:
        mol_tested = Chem.MolFromSmiles(element)
        mol_tested_list.append(mol_tested)
    gdi_aviable_novel_drug_like_analysis['mol'] = mol_tested_list
    inchikey_list = []
    for element in gdi_aviable_novel_drug_like_analysis['mol']:
        try:
            inchikey = Chem.MolToInchiKey(element)
            inchikey_list.append(inchikey)
        except:
            inchikey_list.append('N/A')
            pass
    gdi_aviable_novel_drug_like_analysis['inchiKey'] = inchikey_list
    return(gdi_aviable_novel_drug_like_analysis)

def drug_likness(df,compound_columna):
    pandas.options.mode.chained_assignment = None
    gdi_aviable_novel_drug_like_analysis = df
    mol_tested_list= []
    for element in df[compound_columna]:
        mol_tested = Chem.MolFromSmiles(element)
        mol_tested_list.append(mol_tested)
    gdi_aviable_novel_drug_like_analysis['mol'] = mol_tested_list
    NumHDonors_list = []
    NumHAcceptors_list = []
    MW_list = []
    LogP_list = []
    for element in gdi_aviable_novel_drug_like_analysis['mol']:
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
    gdi_aviable_novel_drug_like_analysis['NumHDonors'] = NumHDonors_list
    gdi_aviable_novel_drug_like_analysis['NumHAcceptors'] = NumHAcceptors_list
    gdi_aviable_novel_drug_like_analysis['MW'] = MW_list
    gdi_aviable_novel_drug_like_analysis['logP'] = LogP_list
    #Cuento cuantas caracterísicas cumple cada compuesto
    countLipinski = lambda row: int(row['NumHDonors'] < 6) +  int(row['NumHAcceptors'] < 6) + int(row['MW'] < 500) + int(row['logP'] < 6)
    gdi_aviable_novel_drug_like_analysis['Lipinski'] = gdi_aviable_novel_drug_like_analysis.apply(countLipinski,axis=1)
    countRO3 = lambda row: int(row['NumHDonors']<= 3) +  int(row['NumHAcceptors']<= 3) + int(row['MW']< 300) + int(row['logP']<= 3)
    gdi_aviable_novel_drug_like_analysis['RO3'] = gdi_aviable_novel_drug_like_analysis.apply(countRO3,axis=1)
    #Selecciono aquellos compuestos que cumplen con 3 o más
    gdi_aviable_novel_druglike_lipinski = gdi_aviable_novel_drug_like_analysis[gdi_aviable_novel_drug_like_analysis['Lipinski']>2]
    gdi_aviable_novel_druglike_RO3 = gdi_aviable_novel_drug_like_analysis[gdi_aviable_novel_drug_like_analysis['RO3']>2]
    gdi_aviable_novel_druglike_not_lipinski = gdi_aviable_novel_drug_like_analysis[gdi_aviable_novel_drug_like_analysis['Lipinski']<3]
    gdi_aviable_novel_druglike_not_RO3 = gdi_aviable_novel_drug_like_analysis[gdi_aviable_novel_drug_like_analysis['RO3']<3]
    #Filtrado según drug like
    df_druglike = pandas.concat([gdi_aviable_novel_druglike_lipinski,gdi_aviable_novel_druglike_RO3])
    df_not_druglike = pandas.concat([gdi_aviable_novel_druglike_not_lipinski,gdi_aviable_novel_druglike_not_RO3])
    df_not_druglike['Filter Name'] = ["Not Druglike"] * len(df_not_druglike)
    return(df_druglike, df_not_druglike)

def generar_csv_filtrado(df,columna_df,criterio):
    datos_filtrados = filtrar_datos(df,columna_df,criterio)
    file_name = 'out_file.csv'
    datos_filtrados.to_csv(file_name)
    