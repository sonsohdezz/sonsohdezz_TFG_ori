# SCRIPT TFG SONSOLES HERNÁNDEZ PIÑEL


# __________________________________________________________________________________________________________________________________________


## MÓDULOS ##


# Importamos los módulos que se van a usar.
import sys
import time
import os
import csv
import re
import pandas as pd
import collections
from ete3 import NCBITaxa


# __________________________________________________________________________________________________________________________________________


## AYUDA Y COMPROBACIÓN DE ARGUMENTOS ##


# Función de ayuda.
def help():
    print("\nINFORMACIÓN")
    print("- Este programa sirve para crear una lista de los taxa encontrados y sus frecuencias a partir de un contig.")
    print("- El archivo necesario para ejecutar el programa es un archivo de anotaciones (.tsv) de los contigs, generado con eggNOG-mapper.")
    print("\nUso: python script.py <annotations_file>\n")
    time.sleep(5)
    sys.exit()


# Función de comprobación de argumentos.
def arg():

    # Se comprueba que el número de argumentos es mayor que 1 y/o si se requiere la ayuda.
    if len(sys.argv) > 1 and (sys.argv[1] == "-h" or sys.argv[1] == "-help"):
        help() 
        sys.exit() 

    # Se comprueba que se ha introducido al menos un argumentos y si no, se despliega la ayuda.
    if len(sys.argv) < 1:
        print("\nERROR: No se han introducido suficientes argumentos.")
        help() 
        sys.exit() 

    # Se comprueba que no se ha introducido más de un argumento y si no, se despliega la ayuda.
    if len(sys.argv) > 2:
        print("\nERROR: Se han introducido demasiados argumentos.")
        help() 
        sys.exit() 
    
    # Se comprueba que el argumento esté en el formato correcto y, si es el caso, se almacena en una variable.
    if sys.argv[1].endswith(".annotations.tsv"):
        input_file = sys.argv[1] 
        file_name = str(sys.argv[1])
        new_file_name = str(file_name[:-4])
    else:
        print("\nERROR: El archivo introducido no tiene un formato válido.")
        help()
        sys.exit()

    return input_file, new_file_name


# __________________________________________________________________________________________________________________________________________


## EXTRACCIÓN DE DATOS DEL ARCHIVO DE ANOTACIONES DE EGGNOG ##


# Función de generación de un diccionario con la columna 'query' y la columna 'eggNOG_OGs' a partir del archivo input.
def info(input_file):

    info_tmp = {}                                       # Creamos un diccionario temporal.

    # Leemos los datos del archivo de anotaciones como un dataframe y cambiamos el nombre de la columna de queries.
    data = pd.read_csv(input_file, header=4, delimiter="\t")
    data.rename(columns = {'#query':'query'}, inplace = True)

    # Recorremos los elementos del dataframe.
    for i in range(len(data)):
        contig_ID = (data.iloc[i]['query'])             # Almacenamos los datos de la columna de queries en una variable 'contig_ID'.
        taxas = (data.iloc[i]['eggNOG_OGs'])            # Almacenamos los datos de la columna de taxas en una variable 'taxas'.
        info_tmp[contig_ID] = taxas                     # Hacemos corresponder a los ID de los contigs los taxas en el diccionario temporal.

    info_dic = {}                                       # Creamos un nuevo diccionario, 'info'.                                   
    
    # Eliminamos los items del diccionario temporal con información extra y almacenamos la nueva información depurada en el diccionario nuevo.   
    for key, value in info_tmp.items():
        if not key.startswith('#'):
            info_dic[key] = value

    return info_dic


# Función de generación de una lista de taxas para cada contig en un nuevo diccionario.
def taxas(info_dic):

    taxas_dic = {}                                        # Creamos un diccionario.

    # Recorremos el diccionario con los IDs de los contigs y los taxa, y creamos una lista con los distintos taxa de cada contig.
    for query in info_dic:
        OGs = info_dic[query]
        taxas = str(OGs)
        list_of_taxas = taxas.split(",")
        
        list_of_new_taxas = []                            # Creamos otro diccionario.
        
        # Sacamos el nombre de los taxa de las líneas de información que hay en OGs, y eliminamos el resto de la información.
        for taxa in list_of_taxas:
            # ch = r"\@"
            # pattern  = ".*" + ch 
            # new_taxa = re.sub(pattern, '', taxa)
            # list_of_new_taxas.append(new_taxa)            # Almacenamos el nombre de la taxa depurado a una lista.

            new_taxa = taxa.split('@')[1]
            list_of_new_taxas.append(new_taxa)
            
        taxas_dic[query] = list(set(list_of_new_taxas))              # El nuevo diccionario hace corresponder {contig_ID:[lista con nombres de taxa]}

    return taxas_dic


# Función que genera un diccionario con la lista de taxas para cada ORF de cada contig.
def fusion(taxas_dic):

    # Creamos dos listas donde almacenaremos los contig_ID y taxas.
    list_of_new_queries = []
    list_of_taxas = []

    # Formateamos los nombres de las filas que corresponden a DISTINTOS ORFs de un MISMO contig, y les damos a todas el MISMO NOMBRE.
    for query in taxas_dic:
        taxas = taxas_dic[query]
        new_query = query[:-2]                      # Eliminamos los dos últimos carácteres del nombre (ej: FRAGMENTO_1_2 = FRAGMENTO_1).
        list_of_taxas.append(taxas)                 # Almacenamos la lista de taxas de cada ORF de un mismo contig a una misma lista.
        list_of_new_queries.append(new_query)       # Añadimos los distintos nombres nuevos para cada ORF a una misma lista.
        
    fusion_dic = {}                                 # Creamos un nuevo diccionario.

    # Creamos una función que devuelve las posisiciones de las filas con el mismo nombre (filas de los disintos ORFs de un MISMO contig).
    def list_duplicates_of(seq,item):
        start_at = -1
        locs = []
        while True:
            try:
                loc = seq.index(item,start_at+1)
            except ValueError:
                break
            else:
                locs.append(loc)
                start_at = loc
        return locs

    taxas_result_list = []                              # Creamos una nueva lista.    

    for query in list_of_new_queries:
        indexes = list_duplicates_of(list_of_new_queries, query)     # Almacenamos las posiciones de las filas con nombre de queries iguales.
        taxa_result_list = [list_of_taxas[i] for i in indexes]       # Extraemos la lista de taxas de esas filas y la añadimos a una misma lista.
        fusion_dic[query] = taxa_result_list                         # El diccionario tiene como values las listas de taxas de cada ORF de un contig.

    return fusion_dic


# Función que genera un diccionario con las estadísticas de cada contig.
def stat(fusion_dic):

    dic_stat = {}                          # Creamos un diccionario nuevo.
    ORFs = []                              # Creamos una lista donde almacenaremos el número de ORFs para cada contig.

    # Recorremos las keys del diccionario de la función anterior.
    for item in fusion_dic.keys():
    
        contig_ID = item                   # Variable con el ID del contig.
        taxa_lists = fusion_dic[item]      # Lista con las listas de taxa de cada ORF para cada contig (value del dic anterior).
        n_ORFs = len(taxa_lists)           # Número de ORFs = número de listas (value).

        # Añadimos el número de ORFs del contig a la lista general de ORFs.
        ORFs.append(n_ORFs)

        # Creamos una lista en la que vamos a recoger los conteos de cada taxa para los diferentes ORF del contig.
        count_list = []

        for nlist in taxa_lists:
            count_dic = {i:nlist.count(i) for i in nlist}
            count_list.append(count_dic)

        # Gracias al módulo counter(), sumamos los valores de cada diccionario de la lista creada para una misma key (=un mismo contig).
        counter = collections.Counter()
        for d in count_list: 
            counter.update(d)
        
        result = dict(counter)


        ## ________________________________________________________________

        ## PARTE NCBI PARA EXTRAER RANGOS TAXONÓMICOS

        # Almacenamos las taxas de cada contig en una variable independiente.
        taxas = result.keys()

        ncbi = NCBITaxa()
        #ncbi.update_taxonomy_database()

        rank_dic = {}

        # Gracias a una lista con los taxa, vamos a detectar el rango de cada taxa.
        for taxa in taxas:
            
            # name2taxid = ncbi.get_name_translator([taxa])           # Traducimos el taxa a su taxid.
            # taxid_tmp = str(name2taxid[taxa])                       # Extraemos el taxid del diccionario resultante.
            # taxid = taxid_tmp[1:-1]                                 # Eliminamos los corchetes que salen por defecto.
            
            # # En caso de que para el mismo nombre haya varios taxid, escogemos el primero, cuyo rango es mayor.
            # if taxid.find(',') >= 1:
            #     taxid = taxid.split(',', 1)[0]

            taxid = int(taxa.split('|')[0])
            taxid2rank = ncbi.get_rank([taxid])                     # Traducimos el taxid a su rango.
          
            if taxid not in taxid2rank:
                continue
            
            rank = str(taxid2rank[taxid])                             # Extraemos el rango del diccionario resultante.
            rank = rank.upper()                                     # Ponemos todo en mayúsculas.

            rank_dic[taxa] = rank                                     # Creamos el diccionario de rangos.

        ## ________________________________________________________________


        # Obtenemos los resultados estadísticos de los distintos taxa, pues sabemos cuantos ORFs hay para cada contig.
    
        # stat_dic = {key + " (" + rank_dic[key] + ") ": round(value/n_ORFs, 2) for key, value in result.items()}

        stat_dic = {rank_dic[key]: key + ": " + str(round(value/n_ORFs, 2)) for key, value in result.items() if key in rank_dic}
        
        dic_stat[item] = stat_dic

    # Mostramos los distintos datos estadísticos de cada contig de una sola vez.
    # for element in dic_stat.items():
    #     print(element)

    return dic_stat, ORFs

    
    

# __________________________________________________________________________________________________________________________________________


## EJECUCIÓN Y CREACIÓN DE DATAFRAME ##


# Ejecutamos las distintas funciones para obtener un DICCIONARIO con una lista de taxas para cada contigs, además de información estadística.
input_file, file_name = arg()
info = info(input_file)
taxas = taxas(info)
dic_wORFs = fusion(taxas)

# for i in dic_wORFs.items():
#     print(i)

dic_stat, ORFs = stat(dic_wORFs)

# for i in dic_stat.items():
#     print(i)


# Generamos un DATAFRAME a partir del diccionario final con los datos estadísticos de los taxa de cada contig.
df = pd.DataFrame.from_dict(dic_stat)

df_transposed = df.T                                                               # Invertimos columnas y filas para obtener los contigs en las filas.
df_transposed.insert(0, 'NUM_ORFs', ORFs, True)                                    # Añadimos una columna con el número de ORFs de cada contig.
df_transposed.index.name = 'CONTIG_ID'                                             # Le damos nombre a la primera columna.
df_transposed.rename(columns = {'SPECIES GROUP':'SPECIES'}, inplace = True)        # Cambiamos el nombre de la columna de especies.   

# Ordenamos las columnas con los 7 rangos deseados y guardamos la tabla en una variable.
def_table = df_transposed[['NUM_ORFs','SUPERKINGDOM','PHYLUM','CLASS','ORDER','FAMILY','GENUS', 'SPECIES']]


# Damos la opción de ver el dataframe en pantalla y/o de guardar en un documento .csv .
ans1 = input("> ¿Quiere visualizar el dataframe en el terminal? [y/n]\n")
ans1 = ans1.upper()
if ans1 == "Y":
    print(def_table)

ans2 = input("\n> ¿Quiere guardar el dataframe con los resultados en un archivo? [y/n]\n")
ans2 = ans2.upper()
if ans2 == "Y":
   def_table.to_csv('dataframe.csv')
   print("\nArchivo guardado como 'dataframe.csv'.")


# # # __________________________________________________________________________________________________________________________________________





