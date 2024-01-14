import matplotlib.pyplot as plt
import numpy as np
import colorsys
from matplotlib.colors import Normalize

def main():
    filename = "../../output/test.off"

    file = open(filename, "r")

    # if file.readline() != "OFF":
    #     # raise SyntaxError(filename + " is not an OFF file")

    fileParams = file.readline().split()

    x = []
    y = []
    z = []
    verticesColor = []
    listIndex = []

    lineIndex = 1

    for line in file:
        lineContent = line.split()

        if lineIndex <= int(fileParams[0]):
            x.append(float(lineContent[0]))
            y.append(float(lineContent[1]))
            z.append(float(lineContent[2]))

            if len(lineContent) > 4:
                verticesColor.append(colorsys.rgb_to_hsv(float(lineContent[3])/256, float(lineContent[4])/256, float(lineContent[5])/256)[0])
                # verticesColor.append(lineIndex)

        elif lineIndex <= (int(fileParams[1]) + int(fileParams[0])):
            listIndex.append([int(index) for index in lineContent[1:]])
        else:
            break
        lineIndex += 1

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    norm = Normalize(min(verticesColor), max(verticesColor))
    mesh = ax.plot_trisurf(x, y, z, triangles=listIndex, cmap='plasma', shade=True, norm=norm)

    mesh.set_array(np.array(verticesColor))
    mesh.autoscale()

    fig.colorbar(mesh, label='Hue Value')

    plt.show()

    file.close()

    return 0


if __name__ == '__main__':
    main()