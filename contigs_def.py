# SCRIPT PARA GENERACIÓN DE CONTIGS PARAMETRADO



# Importamos los módulos que se van a usar.
import sys
import time
import os
import getopt
import shutil


# Función de ayuda.
def help():
	print("\nSINOPSIS\n")
	print("\tpython script.py -i <inputfile> -o <outputdir> -c <contigs_size> -w <scrollingwindow_size>\n")
	print("\nDESCRIPCIÓN\n")
	print("\t Este programa sirve para fragmentar un genoma de elección y simular contigs de un tamaño dado.")
	print("\t El archivo resultado consiste en un archivo multifasta con los distintos fragmentos (contigs) obtenidos con la fragmentación.") 
	print("\t Las opciones de comando necesarios se describen a continuación:\n")
	print("\t [inputfile]            \t Path al archivo que contiene el genoma a fragmentar.")
	print("\t [outputdir]            \t Path al directorio donde se desee almacenar el archivo resultado.")
	print("\t [contigs_size]         \t Tamaño deseado para los contigs simulados.")
	print("\t [scrollingwindow_size] \t Tamaño deseado para la ventana de desplazamiento (nivel de superposición de los contigs).\n")
	# print("\t El archivo resultado consiste en una tabla con la taxonomía obtenida con eggNOG-mapper, ordenada por ORFs.")
	# print("\nOPCIONES DE COMANDO\n")
	# print("-m\t Permite devolver como resultado una tabla adicional con la taxonomía ordenada por contigs.")

	time.sleep(5)
	sys.exit()


# Función de comprobación de argumentos.
def main(argv):

	inputfile = ''
	outputdir = ''
	fragment_size = ''
	window_size = ''

	try:
		opts, args = getopt.getopt(argv, "i:o:c:w:h")
	except getopt.GetoptError:
		help()
		sys.exit(2)

	# Verificamos que se han introducido todos los argumentos necesarios.
	if len(opts) != 4:
		help()
		sys.exit()

	# Para cada parámetro, comprobamos el argumento introducido y almacenamos la información si es correcta.
	for opt, arg in opts:

		if opt in ("-h"):
			help()
			sys.exit()

		elif opt in ("-i"):
			inputfile = arg

		elif opt in ("-o"):
			outputdir = arg

		elif opt in ("-c"):
			try:
				fragment_size = int(arg)
				try:
					x = 1/float(fragment_size)
				except ZeroDivisionError:
					print("\nERROR: El valor introducido [contigs_size] debe ser no nulo. Por favor, vuelva a intentarlo.\n")
					exit()
			except ValueError:
				print("\nERROR: El carácter introducido en [contigs_size] no es válido. Debe introducir un número. Por favor, vuelva a intentarlo.\n")
				exit()

		elif opt in ("-w"):
			try:
				window_size = int(arg)
				try:
					x = 1/float(window_size)
				except ZeroDivisionError:
					print("\nERROR: El valor introducido [scrollinwindow_size] debe ser no nulo. Por favor, vuelva a intentarlo.\n")
					exit()
			except ValueError:
				print("\nERROR: El carácter introducido en [scrollingwindow_size] no es válido. Debe introducir un número. Por favor, vuelva a intentarlo.\n")
				exit()

		else:
			help()
			sys.exit()


	return inputfile, outputdir, fragment_size, window_size


# Función que fragmenta el genoma.
def frag(inputfile, fragment_size, window_size):

	# Almacenamos el nombre del archivo genómico en una variable.
	path = str(inputfile)
	file_name = os.path.basename(path).split('/')[-1].split('.')[0]

	# Abrimos el archivo y lo leemos.
	with open(inputfile, 'r') as fin:
		data = fin.read().splitlines(True)
	with open('file_temp.txt', 'w') as fout:
		fout.writelines(data[1:])

	# Transformamos el texto del archivo en un string.
	text_file = open('file_temp.txt', "r")
	string = text_file.read()
	text_file.close()
	os.remove('file_temp.txt')

	# Eliminamos los saltos de línea del string y almacenamos el resultado en una variable 'genome'.
	genome = string.replace('\n', '').upper()
	print("\n(1/3) El genoma introducido tiene una longitud de " + str(len(genome)) + " pb.")

	# Recorremos el genoma en tandas del tamaño indicado de la ventana de desplazamiento y extraemos fragmentos del tamaño indicado.
	pos_ini = 0
	pos_end = pos_ini + fragment_size

	list_of_contigs = []
	fragment = ''

	while pos_end < len(genome):
        
		fragment = genome[pos_ini:pos_end]
		list_of_contigs.append(fragment)

		pos_ini += window_size + 1
		pos_end = pos_ini + fragment_size

	print("(2/3) El número de fragmentos generados con los parámetros introducidos es de:", len(list_of_contigs))
	return(list_of_contigs, file_name)


# Función que almacena los fragmentos generados en formato FASTA en un archivo multifasta en el directorio indicado.

def multifasta(fragments, file_name, outputdir):

	# Creamos un bucle que almacene en un archivo multifasta los diferentes fragmentos.
	with open('multifasta.fasta', 'w') as f:
		i = 0
		for fragment in fragments:
			i += 1
			f.writelines(">FRAGMENTO_" + str(i) + "\n")
			f.writelines(fragment + "\n\n")

	# Creamos una variable con el contenido del archivo multifasta.        
	fasta_file = open('multifasta.fasta')
	file_contents = fasta_file.read()

	# Renombramos el archivo con los fragmentos para que tenga el mismo nombre que el archivo con el genoma.
	new_file_name = '%s.fasta' % file_name
	os.rename('multifasta.fasta', new_file_name)

	# Definimos las rutas para mover el archivo multifasta creado al directorio destino.
	src_folder = os.getcwd()
	dst_folder = outputdir

	loc_1 = src_folder + "/" + new_file_name
	loc_2 = dst_folder + "/" + new_file_name

	# Verificamos si ya existe un archivo multifasta anterior en el directorio destino.
	if os.path.exists(loc_2):
		
		print("\nYa existe un archivo multifasta con esta información en el directorio destino.")
		print(">> ¿Desea reemplazarlo? [y/n]")
    
		ans = input().upper()
		if ans == 'Y':
			shutil.move(loc_1, loc_2)
			print("\n")
		else:
			new_file_name = new_file_name.split(".")[0] + "_new." + new_file_name.split(".")[1]
			new_loc_2 = dst_folder + "/" + new_file_name
			shutil.move(loc_1, new_loc_2)
			print("\n")

	else: 
		shutil.move(loc_1, loc_2)

	print("(3/3) Los diferentes contigs generados se encuentran en el archivo '%s'" % new_file_name)
	print('\nFin de la ejecución.\n')
	exit()


## ___________________________________________________________________________________________________________________________________


# EJECUCIÓN DE LAS DISTINTAS FUNCIONES


inputfile, outputdir, fragment_size, window_size = main(sys.argv[1:])
list_of_contigs, filename = frag(inputfile, fragment_size, window_size)
multifasta(list_of_contigs, filename, outputdir)


	

