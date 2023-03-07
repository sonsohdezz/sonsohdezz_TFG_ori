# SCRIPT PARA GENERACIÓN DE CONTIGS


# Importamos los módulos que se van a usar.
import sys
import time
import os


# Función de ayuda.
def help():
    print("\nINFORMACIÓN")
    print("- Este programa sirve para fragmentar un genoma de elección.")
    print("- El archivo necesario para ejecutar el programa es un archivo fasta con el genoma que se quiera fragmentar.")
    print("\nUso: python script.py <fasta_file>\n")
    time.sleep(5)
    sys.exit()


# Función de comprobación de argumentos.
def arg():

    # Se comprueba que el número de argumentos es mayor que 1 y si se requiere la ayuda.
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
    if sys.argv[1].endswith(".fna"):
        input_file = sys.argv[1] 
        file_name = str(sys.argv[1])
        new_file_name = str(file_name[:-4])
    else:
        print("\nERROR: El archivo introducido no tiene un formato válido.")
        help()
        sys.exit()

    return input_file, new_file_name


# Función de fragmentación del genoma contenido en el archivo input.
def contigs(input_file):

    # Eliminamos la primera fila del archivo con el genoma (que contiene la cabecera).
    with open(input_file, 'r') as fin:
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
    print("El genoma tiene una longitud de " + str(len(genome)) + " pb.")

    # Obtenemos el tamaño de fragmento que se desea.
    print("\n>> ¿Qué tamaño desea que tengan los fragmentos generados?")
    
    while True:
        try:
            fragment_size = int(input())
            break
        except ValueError:
            print("\nERROR: El carácter introducido no es válido. Debe introducir un número. Por favor, vuelva a intentarlo.\n")
            exit()

    if fragment_size == 0:
        print("\nERROR: El valor introducido debe ser no nulo. Por favor, vuelva a intentarlo.\n")
        exit()

    # Obtenemos el tamaño de la ventana de desplazamiento (grado de solapamiento o separación de los fragmentos) que se desea.
    print(">> ¿Qué tamaño desea que tenga la ventana de desplazamiento?")

    while True:
        try:
            window_size = int(input())
            break
        except ValueError:
            print("\nERROR: El carácter introducido no es válido. Debe introducir un número. Por favor, vuelva a intentarlo.\n")
            exit()
    
    if window_size == 0:
        print("\nERROR: El valor introducido debe ser no nulo. Por favor, vuelva a intentarlo.\n")
        exit()

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

    print("\nEl número de fragmentos generados es de:", len(list_of_contigs))
    return(list_of_contigs)


# __________________________________________________________________________________________________________________________________________


# Ejecutamos las distintas funciones para obtener una lista de fragmentos (contigs).
input_file, file_name = arg()
fragments = contigs(input_file)

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
print("Los diferentes fragmentos se encuentran en el archivo '%s'" % new_file_name)

# Damos la opción al usuario de visualizar los fragmentos en el terminal
print("\n>> ¿Desea visualizar los fragmentos generados? [Y/y o N/n]")

ans = input().upper()
if ans == 'Y':
    os.system('clear')
    print(file_contents)
else:
    print('\nFin de la ejecución.\n')
    exit()












