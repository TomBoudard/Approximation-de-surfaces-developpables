import openmesh as om
import numpy as np
import math
import matplotlib.pyplot as plt


def cotan(x):
    return np.cos(x)/np.sin(x)


def dist(mesh,x1, x2):
    x1 = mesh.point(x1)
    x2 = mesh.point(x2)
    xDiff = x1[0] - x2[0]
    yDiff = x1[1] - x2[1]
    zDiff = x1[2] - x2[2]
    return np.sqrt(xDiff * xDiff + yDiff * yDiff + zDiff * zDiff)

def compute_angles(mesh, A, B, C):
    ''' Compute the angles of the given triangle'''
    # length of sides be a, b, c
    a = dist(mesh, B,C)
    b = dist(mesh, A,C)
    c = dist(mesh, A,B)

    # Cosinus law
    alpha = math.acos((b**2 + c**2 - a**2) / (2 * b * c))
    beta = math.acos((a**2 + c**2 - b**2) / (2 * a * c))
    gamma = math.acos((a**2 + b**2 - c**2) / (2 * a * b))

    return alpha,beta,gamma

def compute_area(mesh, face_vertices_dict, face_vertices_angles, face_id):
    ''' Compute the area of the face given'''
    A,B,C = face_vertices_dict[face_id]
    angle = face_vertices_angles[face_id][A.idx()]
    b = dist(mesh, A,C)
    c = dist(mesh, A,B)
    area = b * c * np.sin(angle) / 2
    return area

def map_curvature_to_color(value, max_val, alpha=1.0, colormap_name='viridis'):
    if value == None: #c'est un bord
        return[1, 1, 1, 1]
    # normalized_value = (value/max_val)
    # epsilon = 1e-5  # Small constant to avoid logarithm of zero
    # normalized_value = np.log(value + epsilon) / np.log(max_val + epsilon)
    power = 0.5  # Set the power factor for the normalization (adjust as needed)
    normalized_value = (value / max_val) ** power
    # Choose a colormap
    colormap = plt.get_cmap(colormap_name)
    # Map the normalized value to a color using the chosen colormap
    color = colormap(normalized_value)
    return color

def add_angles(mesh):
    '''This function aim to store for every face of our mesh, the vertices associated and the angle associated to
    each vertices in the face'''
    face_vertices_dict = {}
    face_vertices_angles = {}
    for face in mesh.faces(): #Iterate on the face of our mesh
        vertices = [v for v in mesh.fv(face)]  # Get vertex of the face
        face_vertices_dict[face.idx()] = vertices  # Store face ID as key and the verteices from the face as value

    for face in mesh.faces():
        vertices = face_vertices_dict[face.idx()] #get the 3 vertices our face
        a,b,c = compute_angles(mesh, vertices[0], vertices[1],vertices[2]) #comput the angles of our triangle
        vertices_angles = {vertices[0].idx(): a, vertices[1].idx() : b, vertices[2].idx(): c}
        face_vertices_angles[face.idx()] = vertices_angles #Store face ID as key and a list of the 3 [vertex ID, vertex's angle] as value
    return face_vertices_dict, face_vertices_angles

def get_opposite_vertices(face_vertices_dict, face_id1, face_id2, vertex):
    '''this function return the vertex id of the vertex being in face1 and face2 but not being vertex_id
    and also return the  two vertices ids of the two remaining vertices in the faces'''
    vertices_face1 = []
    vertices_face2 = []
    opp_vertices = []
    for i in range(3):
        vertex1 = face_vertices_dict[face_id1][i]
        vertex2 = face_vertices_dict[face_id2][i]
        if vertex1.idx() != vertex.idx():
            vertices_face1.append(vertex1.idx())
            if (vertex1.idx() in vertices_face2):
                opp_vertices.append(vertex1)
                vertices_face1.remove(vertex1.idx())
                vertices_face2.remove(vertex1.idx())
        if vertex2.idx() != vertex.idx():
            vertices_face2.append(vertex2.idx())
            if (vertex2.idx() in vertices_face1):
                opp_vertices.append(vertex2)
                vertices_face1.remove(vertex2.idx())
                vertices_face2.remove(vertex2.idx())
    opp_vertices.append(vertices_face1[0])
    opp_vertices.append(vertices_face2[0])
    return(opp_vertices)

def gradientColor(value):

    return (value, 0.5, 0.5)