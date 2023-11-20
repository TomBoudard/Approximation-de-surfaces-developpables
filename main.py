import openmesh as om
import numpy as np

def read_off_file(file):
    with open(file, 'r') as file:
        first_line = file.readline().strip()
        print(first_line)
        if first_line != 'OFF':
            raise TypeError("Wrong file format, you need to use a .off file")
        line = file.readline().split()
        NB_vertices = line[0]
        NB_faces = line[1]
        NB_edges = line[2]
        print(NB_vertices)
        print(NB_faces)
        print(NB_edges)

# def write_off_file():

def main():
    print("---------- Début main\n")
    off_file = "Objects/test.off"
    #à la main:
    #read_off_file(off_file)
    # En utilisnt les librairies:
    mesh = om.read_trimesh(off_file) #créer un mesh à partir du fichier .off
    #om.write_mesh('test.off', mesh) #Ecrit un fichier off à partir d'un mesh
    # TODO implémenter algorithme
    print("---------- Fin main\n")

if __name__ == "__main__":
    main()  # Call the main function when this file is executed directly