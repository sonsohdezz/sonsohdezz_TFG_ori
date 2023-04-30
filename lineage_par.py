# SCRIPT LINAJES TFG PARAMETRADO


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
import getopt



## AYUDA Y COMPROBACIÓN DE ARGUMENTOS ##


# Función de ayuda.
def help():
	print("\nSINOPSIS\n")
	print("\tpython script.py -i <inputfile> \n")
	print("\nDESCRIPCIÓN\n")
	print("\t Este programa sirve para devolver la taxonomía de un genoma, predicha por eggNOG-mapper, pero ordenada y con datos estadísticos.")
	print("\t El archivo resultado consiste en una tabla con la taxonomía obtenida con eggNOG-mapper, ordenada por ORFs.")
	print("\t Las opciones de comando necesarios se describen a continuación:\n")
	print("\t [inputfile] \t Path al archivo que contiene el archivo de anotaciones .tsv obtenido tras correr eggNOG-mapper.")
	print("\nOPCIONES DE COMANDO\n")
	print("\t [-m] \t        Permite devolver como resultado una tabla adicional con la taxonomía ordenada por contigs.\n")

	time.sleep(2)
	sys.exit()


# Función de comprobación de argumentos.
def main(argv):

	inputfile = ''
	file_name = ''
	contig_table = ''

	try:
		opts, args = getopt.getopt(argv, "hi:m")
	except getopt.GetoptError:
		help()
		sys.exit()

	# Verificamos que se han introducido todos los argumentos necesarios.
	if len(opts) < 1 :
		help()
		sys.exit()

	# Para cada parámetro, comprobamos el argumento introducido y almacenamos la información si es correcta.
	for opt, arg in opts:

		if opt in ("-m"):
			contig_table = True

		elif opt in ("-h"):
			help()
			sys.exit()

		elif opt in ("-i"):

			if arg.endswith(".annotations.tsv"):
				inputfile = arg
				file_name = str(inputfile).split(".")[0]

		else:
			help()
			sys.exit()

	return inputfile, file_name, contig_table


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

	fusion_dic = {}                               # Creamos un nuevo diccionario.

	# Formateamos los nombres de las filas que corresponden a DISTINTOS ORFs de un MISMO contig, y les damos a todas el MISMO NOMBRE.
	for ORF in filo_dic:
		lineage = filo_dic[ORF]
		contig_id = "_".join(ORF.split("_")[:-1])                      # Eliminamos los dos últimos carácteres del nombre (ej: FRAGMENTO_1_2 = FRAGMENTO_1).
		if contig_id in fusion_dic:
			contig_list = fusion_dic[contig_id]
		else:
			contig_list = []
			fusion_dic[contig_id] = contig_list
		contig_list.append(lineage)

	return fusion_dic


# Función que genera un diccionario con las estadísticas de cada contig.
def format(fusion_dic):

	new_lineages_dic = {}

	# Recorremos las keys del diccionario de la función anterior.
	for contig_ID, lineages_lists in fusion_dic.items():
    
		# contig_ID = item                   	# Variable con el ID del contig.
		# lineages_lists = fusion_dic[item]    # Lista con las listas de linajes de cada ORF para cada contig (value del dic anterior).

		# print(contig_ID)
		# print(lineages_lists)
		# print(n_ORFs)

		# Creamos un diccionario de referencia que transforma los nombres de los rangos taxonómicos de interés en iniciales.
		ref_dic = {'SUPERKINGDOM':'SK', 'PHYLUM':'P', 'CLASS':'C', 'ORDER':'O', 'FAMILY':'F', 'GENUS':'G', 'SPECIES GROUP':'SG'}	

		new_lineages_list = []

		# Para cada ORF de un mismo contig, transformamos el diccionario de linajes en una lista de rangos (indicando el rango y el nombre).
		for lineage in lineages_lists:
			
			lineage_list = []

			# Para CADA rango dentro del linaje de CADA ORF.
			for rank in lineage:
				
				# Si el rango se encuentra en los rangos del diccionario de referencia, se coge el valor y se transforma el rango en su inicial.
				if rank in ref_dic:
					new_rank_name = ref_dic[rank]
					value = lineage[rank]
					new_value = value.split('|')[1] 
				else: 
					continue

				lineage_list.append(new_rank_name + '_' + new_value)			# Añadimos a la lista del linaje el rango y su valor en el formato nuevo.

		# print(lineage_list)
		# print("\n")

			new_lineages_list.append(lineage_list)
		
		new_lineages_dic[contig_ID] = new_lineages_list
	
	return new_lineages_dic


# Función que genera un diccionario en el cual se muestra el linaje más abundante para cada contig y su frecuencia.
def best_lineage(new_lineages_dic):

	best_lineage_dic = {}
	ORFs = []                              # Creamos una lista donde almacenaremos el número de ORFs para cada contig.

	# Recorremos el diccionario de la función anterior cuyas keys son los contig_IDs y las values son las listas de linajes de los ORFs del contig.
	for contig_ID, lists_of_lineages in new_lineages_dic.items():

		lineage_count_dic = {tuple(i):lists_of_lineages.count(i) for i in lists_of_lineages}		# Contamos el número de apariciones de cada linaje.
		
		# max_lineage = max(lineage_count_dic)
		max_value = max(lineage_count_dic.values())

		max_lineages = [lineage for lineage in lineage_count_dic if lineage_count_dic[lineage] == max_value]
		max_lineage = [lineage for lineage in max_lineages if len(lineage)==max([len(lin) for lin in max_lineages])]
		if len(max_lineage) > 1:
			continue

		n_ORFs = len(lists_of_lineages)         				# Número de ORFs = número de listas (value).
		lineage_freq = round((max_value/n_ORFs), 2)

		best_lineage_dic[contig_ID] = (max_lineage, lineage_freq, n_ORFs)

		# Añadimos el número de ORFs del contig a la lista general de ORFs.
		# ORFs.append(n_ORFs)

	return best_lineage_dic# , ORFs


# _____________________________________________________________________________________________________________________________


## EJECUCIÓN Y CREACIÓN DEL DATAFRAME DE ORFs ##

input_file, file_name, contig_table = main(sys.argv[1:])
info = info(input_file)
# print(info, file = sys.stderr)
taxas = taxas(info)
# print(taxas, file = sys.stderr)
dic_ranks = ranks(taxas)
# print(dic_ranks, file = sys.stderr)
filo_dic = phylo(dic_ranks)
# print(filo_dic, file = sys.stderr)


# Generamos un DATAFRAME DE ORFS a partir del diccionario final con los datos estadísticos de los taxa de cada ORF.
df = pd.DataFrame.from_dict(filo_dic)
df_transposed = df.T												# Invertimos columnas y filas para obtener los contigs en las filas.
df_transposed.index.name = 'CONTIG_ID'						# Le damos nombre a la primera columna.

# Ordenamos las columnas con los 7 rangos deseados y guardamos la tabla en un archivo .csv.
def_table = df_transposed[['SUPERKINGDOM','PHYLUM','CLASS','ORDER','FAMILY','GENUS', 'SPECIES GROUP']]
def_table.to_csv(file_name + '_lineage_ORFs.csv')

print("\n(1/2) Dataframe de ORFs creado y guardado como " + file_name + "_lineage_ORFs.csv")


# _____________________________________________________________________________________________________________________________


## EJECUCIÓN Y CREACIÓN DEL DATAFRAME DE CONTIGs ##


if contig_table == True:
	print(list(filo_dic.items())[15:19], file = sys.stderr)
	print()
	fusion_dic = fusion(filo_dic)
	print(list(fusion_dic.items())[5], file = sys.stderr)
	print()
	new_lineages_dic = format(fusion_dic)
	print(list(new_lineages_dic.items())[5], file = sys.stderr)
	print()
	stat_dic = best_lineage(new_lineages_dic)
	# print(stat_dic, file = sys.stderr)
	# print(ORFs, file = sys.stderr)

	# Generamos un DATAFRAME a partir del diccionario final con los datos estadísticos de los taxa de cada contig.
	df_2 = pd.DataFrame.from_dict(stat_dic)

	df_2_transposed = df_2.T											# Invertimos columnas y filas para obtener los contigs en las filas.
	df_2_transposed.index.name = 'CONTIG_ID'						# Le damos nombre a la primera columna.
	# df_2_transposed.insert(0, 'NUM_ORFs', ORFs, True)			# Añadimos una columna con el número de ORFs de cada contig.

	# Ordenamos las columnas con las 3 columnas deseados y guardamos la tabla en un archivo .csv.
	def_table_2 = df_2_transposed.rename(columns = {0:'TAXONOMY', 1:'FREQUENCY', 2:'NUM_ORFs'})
	def_table_2.to_csv(file_name + '_lineage_contigs.csv')

	print("(2/2) Dataframe de contigs creado y guardado como " + file_name + "_lineage_contigs.csv\n")
	print("\nFin de la ejecución.\n")
	exit()

else:

	print("\nFin de la ejecución.\n")
	exit()





