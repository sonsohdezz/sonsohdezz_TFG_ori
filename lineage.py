# SCRIPT LINAJES TFG


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
	for query in info_dic.keys():
		OGs = info_dic[query]
		taxas = str(OGs)
		list_of_taxas = taxas.split(",")
		
		# print(list_of_taxas)
		
		list_of_new_taxas = []

		# Sacamos el nombre de los taxa de las líneas de información que hay en OGs, y eliminamos el resto de la información.
		for taxa in list_of_taxas:
			new_taxa = taxa.split('@')[1]

			if new_taxa not in list_of_new_taxas:
				list_of_new_taxas.append(new_taxa)

		# print(list_of_new_taxas)

		taxas_dic[query] = list_of_new_taxas

	return taxas_dic


# Función que genera un diccionario con la lista de taxas Y SUS RANGOS para cada ORF de cada contig.
def ranks(taxas_dic):

	# Activamos el módulo del NCBI.
	ncbi = NCBITaxa()
	#ncbi.update_taxonomy_database()

	# Creamos el diccionario donde almacenaremos la información nueva.
	dic_ranks = {}

	# Sacamos los rangos de cada taxa en las listas de cada ORF (fila) del diccionario anterior.
	for item in taxas_dic.keys():

		# Creamos la lista donde almacenaremos los nuevos taxas con sus rangos.
		list_of_ranktaxas = []

		# Definimos el contig_ID para el diccionario.
		contig_ID = item
		# print(contig_ID)

		# Definimos la lista de taxas para poder recorrerla en un bucle.
		list_of_taxas = taxas_dic[item]
		# print(list_of_taxas)

		for taxa in list_of_taxas:

			taxid = int(taxa.split("|")[0])
			taxid2rank = ncbi.get_rank([taxid])

			if taxid not in taxid2rank:
				continue

			rank = str(taxid2rank[taxid]).upper()
			new_value = rank + ":" + taxa

			if new_value not in list_of_ranktaxas:
				list_of_ranktaxas.append(new_value)

		# print(contig_ID,list_of_ranktaxas)

		dic_ranks[contig_ID] = list_of_ranktaxas
	
	return dic_ranks

	
# Función que detecta los rangos de cada taxa en los diccionarios anteriores de cada ORF y saca la filogenia acorde.
def phylo(dic_ranks):

	# Creamos un diccionario nuevo.
	filo_dic = {}

	# Activamos el módulo del NCBI.
	ncbi = NCBITaxa()
	#ncbi.update_taxonomy_database()

	# Creamos una lista de referencia con la que comparar los rangos que extraigamos de dic_ranks.
	list_ref = ['SPECIES GROUP', 'GENUS', 'FAMILY', 'ORDER', 'CLASS', 'PHYLUM', 'SUPERKINGDOM']

	for i in dic_ranks.keys():

		ORF_id = i 							# Almacenamos la key (identificador del ORF) en una variable.
		list_of_taxas = dic_ranks[i]		# Almacenamos los values correspondientes a la key de turno en otra variable.

		# print(ORF_id)

		# Recorremos la lista de rangos de referencia.
		for k in list_ref:			

			# print(str(k), list_of_taxas)

			# Comprobamos si el rango de la lista de ref. se encuentra en la lista de taxas del ORF dado.
			match = [i for i in list_of_taxas if k in i]

			# Si no encontramos una coincidencia, seguimos recorriendo lista_ref.
			if len(match) == 0:
				
				# print(k, "is not found in", list_of_taxas, "\n")
				# print("Rango no encontrado.\n")
				continue
			
			# Si encontramos una coincidencia, realizamos una serie de operaciones.
			else:
				
				# print(k, "is found in", list_of_taxas)
				# print(match[0])
				# print("Rango encontrado.")
				
				taxid = match[0].split(":")[1].split("|")[0]				# Almacenamos en una variable el taxid de la value coincidente.
				
				lineage = ncbi.get_lineage(int(taxid))						# Extraemos el linaje completo del taxid de la value coincidente.	
				# print("Linaje creado.")
				# print(lineage)

				rank_dic = {}

				for j in lineage:
					taxid2rank = ncbi.get_rank([j])                     	# Traducimos cada taxid del linaje a su rango.
					rank = str(taxid2rank[j])                             	# Extraemos el rango del diccionario resultante.
					rank = rank.upper()
					rank_dic[j] = rank 										# Almacenamos en un diccionario los rangos para cada elemento del linaje.

				names = ncbi.get_taxid_translator(lineage)					# Convertimos los taxids del linaje a sus nombres respectivos.		
				
				# Almacenamos el linaje con sus rangos en un diccionario.
				lineage_dic = {rank_dic[taxid]: str(taxid) + "|" + str(names[taxid]) for taxid in lineage} 	
				
				# Creamos una entrada en el diccionario del tipo ORF_id:dic_rango_linaje.
				filo_dic[ORF_id] = lineage_dic													

				break

	return filo_dic


# Función que genera un diccionario con la lista de taxas para cada ORF de cada contig.
def fusion(filo_dic):

	# Creamos dos listas donde almacenaremos los ORFs (keys) y los linajes (values).
	contig_ids_list = []
	list_of_lineages = []

	# Formateamos los nombres de las filas que corresponden a DISTINTOS ORFs de un MISMO contig, y les damos a todas el MISMO NOMBRE.
	for ORF in filo_dic:
		lineage = filo_dic[ORF]
		contig_id = ORF[:-2]                      # Eliminamos los dos últimos carácteres del nombre (ej: FRAGMENTO_1_2 = FRAGMENTO_1).
		list_of_lineages.append(lineage)          # Almacenamos el linaje de cada ORF de un mismo contig a una misma lista.
		contig_ids_list.append(contig_id)      	 # Añadimos los distintos nombres nuevos para cada ORF a una misma lista.
		
	fusion_dic = {}                               # Creamos un nuevo diccionario.

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

	lineages_result_list = []                              # Creamos una nueva lista.    

	for contig in contig_ids_list:
		indexes = list_duplicates_of(contig_ids_list, contig)     					# Almacenamos las posiciones de las filas con nombre de queries iguales.
		lineages_result_list = [list_of_lineages[i] for i in indexes]       # Extraemos la lista de taxas de esas filas y la añadimos a una misma lista.
		fusion_dic[contig] = lineages_result_list                         		# El diccionario tiene como values las listas de taxas de cada ORF de un contig.

	return fusion_dic


# Función que genera un diccionario con las estadísticas de cada contig.
def format(fusion_dic):

	new_lineages_dic = {}

	# Recorremos las keys del diccionario de la función anterior.
	for item in fusion_dic.keys():
    
		contig_ID = item                   	# Variable con el ID del contig.
		lineages_lists = fusion_dic[item]    # Lista con las listas de linajes de cada ORF para cada contig (value del dic anterior).

		# print(contig_ID)
		# print(lineages_lists)
		# print(n_ORFs)

		# Creamos un diccionario de referencia que transforma los nombres de los rangos taxonómicos de interés en iniciales.
		ref_dic = {'SUPERKINGDOM':'SK', 'PHYLUM':'P', 'CLASS':'C', 'ORDER':'O', 'GENUS':'G', 'SPECIES GROUP':'SG'}	

		new_lineages_list = []

		# Para cada ORF de un mismo contig, transformamos el diccionario de linajes en una lista de rangos (indicando el rango y el nombre).
		for lineage in lineages_lists:
			
			lineage_list = []

			# Para CADA rango dentro del linaje de CADA ORF.
			for rank in lineage:
				
				# Si el rango se encuentra en los rangos del diccionario de referencia, se coge el valor y se transforma el rango en su inicial.
				if rank in ref_dic.keys():
					new_rank_name = ref_dic[rank]
					value = lineage[rank]
					new_value = value.split('|')[1] 
				else: 
					continue

				lineage_list.append(new_rank_name + '_' + new_value)			# Añadimos a la lista del linaje el rango y su valor en el formato nuevo.

		# 	print(lineage_list)
		# print("\n")

			new_lineages_list.append(lineage_list)
		
		new_lineages_dic[contig_ID] = new_lineages_list
	
	return new_lineages_dic


# Función que genera un diccionario en el cual se muestra el linaje más abundante para cada contig y su frecuencia.
def best_lineage(new_lineages_dic):

	best_lineage_dic = {}
	ORFs = []                              # Creamos una lista donde almacenaremos el número de ORFs para cada contig.

	# Recorremos el diccionario de la función anterior cuyas keys son los contig_IDs y las values son las listas de linajes de los ORFs del contig.
	for contig in new_lineages_dic:
		
		contig_ID = contig 											# Definimos el contig_ID.
		lists_of_lineages = new_lineages_dic[contig]			# En una variable, almacenamos las listas de linajes de los ORFs del contig.
		n_ORFs = len(lists_of_lineages)         				# Número de ORFs = número de listas (value).

		# Añadimos el número de ORFs del contig a la lista general de ORFs.
		ORFs.append(n_ORFs)

		lineage_count_dic = {tuple(i):lists_of_lineages.count(i) for i in lists_of_lineages}		# Contamos el número de apariciones de cada linaje.
		# print(lineage_count_dic)

		max_lineage = max(lineage_count_dic)
		max_value = max(lineage_count_dic.values())
		lineage_freq = round((max_value/n_ORFs), 2)

		best_lineage_dic[contig_ID] = (max_lineage, lineage_freq)

	return best_lineage_dic, ORFs


# _____________________________________________________________________________________________________________________________


## EJECUCIÓN Y CREACIÓN DEL DATAFRAME DE ORFs ##


# Ejecutamos las distintas funciones para obtener un DICCIONARIO con una lista de taxas para cada contigs, además de información estadística.
input_file, file_name = arg()

info = info(input_file)
# print(info)

taxas = taxas(info)
# print(taxas)

# dic_wORFs = fusion(taxas)
# for i in dic_wORFs.items():
# 	print(i)

dic_ranks = ranks(taxas)
# for i in dic_ranks.items():
# 	print(i)

filo_dic = phylo(dic_ranks)
# for i in filo_dic.items():
# 	print(i)


# Generamos un DATAFRAME a partir del diccionario final con los datos estadísticos de los taxa de cada contig.
# df = pd.DataFrame.from_dict(filo_dic)

# df_transposed = df.T                                                               # Invertimos columnas y filas para obtener los contigs en las filas.
# df_transposed.index.name = 'CONTIG_ID'                                             # Le damos nombre a la primera columna.
# # df_transposed.insert(0, 'NUM_ORFs', ORFs, True)                                    # Añadimos una columna con el número de ORFs de cada contig.

# # Ordenamos las columnas con los 7 rangos deseados y guardamos la tabla en una variable.
# def_table = df_transposed[['SUPERKINGDOM','PHYLUM','CLASS','ORDER','FAMILY','GENUS', 'SPECIES GROUP']]

# # Damos la opción de ver el dataframe en pantalla y/o de guardar en un documento .csv .
# ans1 = input("> ¿Quiere visualizar el dataframe en el terminal? [y/n]\n")
# ans1 = ans1.upper()
# if ans1 == "Y":
#     print(def_table)

# ans2 = input("\n> ¿Quiere guardar el dataframe con los resultados en un archivo? [y/n]\n")
# ans2 = ans2.upper()
# if ans2 == "Y":
#    def_table.to_csv('dataframe_lineage.csv')
#    print("\nArchivo guardado como 'dataframe_lineage.csv'.")


# _____________________________________________________________________________________________________________________________


## EJECUCIÓN Y CREACIÓN DEL DATAFRAME DE CONTIGs ##

fusion_dic = fusion(filo_dic)
# for i in fusion_dic.items():
# 	print(i)

new_lineages_dic = format(fusion_dic)
# for i in new_lineages_dic.items():
# 	print(i)

stat_dic, ORFs = best_lineage(new_lineages_dic)
# for i in stat_dic.items():
# 	print(i)


# Generamos un DATAFRAME a partir del diccionario final con los datos estadísticos de los taxa de cada contig.
df = pd.DataFrame.from_dict(stat_dic)

df_transposed = df.T                                                               # Invertimos columnas y filas para obtener los contigs en las filas.
df_transposed.index.name = 'CONTIG_ID'                                             # Le damos nombre a la primera columna.
df_transposed.insert(0, 'NUM_ORFs', ORFs, True)                                    # Añadimos una columna con el número de ORFs de cada contig.

df_transposed.rename(columns = {0:'TAXONOMY', 1:'FREQUENCY'}, inplace = True)

# print(df_transposed)

# # Ordenamos las columnas con los 7 rangos deseados y guardamos la tabla en una variable.
# def_table = df_transposed[['NUM_ORFs','BEST TAXONOMY','FREQ']]

# Damos la opción de ver el dataframe en pantalla y/o de guardar en un documento .csv .
ans1 = input("> ¿Quiere visualizar el dataframe en el terminal? [y/n]\n")
ans1 = ans1.upper()
if ans1 == "Y":
    print(df_transposed)

ans2 = input("\n> ¿Quiere guardar el dataframe con los resultados en un archivo? [y/n]\n")
ans2 = ans2.upper()
if ans2 == "Y":
   df_transposed.to_csv('dataframe_best_lineage.csv')
   print("\nArchivo guardado como 'dataframe_best_lineage.csv'.")


