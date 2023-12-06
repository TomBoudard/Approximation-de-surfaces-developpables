import openmesh as om
import numpy as np

def vector(mesh, vertexStart, vertexEnd):
    return (mesh.point(vertexEnd)[0] - mesh.point(vertexStart)[0],
            mesh.point(vertexEnd)[1] - mesh.point(vertexStart)[1],
            mesh.point(vertexEnd)[2] - mesh.point(vertexStart)[2])


def hsv_to_rgb(h, s, v):
    # Assurez-vous que les valeurs de H, S et V sont dans la plage correcte
    h = max(0, min(1, h))
    s = max(0, min(1, s))
    v = max(0, min(1, v))

    # Si la saturation est presque nulle, la couleur est une nuance de gris
    if s == 0:
        return int(v * 255), int(v * 255), int(v * 255)

    # Sinon, utilisez la formule de conversion HSV vers RGB
    i = int(h * 6.) # Sélectionne une région de l'espace couleur (0 à 5)
    f = (h * 6.) - i
    p = v * (1. - s)
    q = v * (1. - f * s)
    t = v * (1. - (1. - f) * s)

    # Mappez chaque composante à la plage des couleurs RVB (0 à 255)
    if i == 0:
        return int(v * 255), int(t * 255), int(p * 255)
    elif i == 1:
        return int(q * 255), int(v * 255), int(p * 255)
    elif i == 2:
        return int(p * 255), int(v * 255), int(t * 255)
    elif i == 3:
        return int(p * 255), int(q * 255), int(v * 255)
    elif i == 4:
        return int(t * 255), int(p * 255), int(v * 255)
    else:  # i == 5
        return int(v * 255), int(p * 255), int(q * 255)