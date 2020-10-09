import trimesh
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from mpl_toolkits.mplot3d import Axes3D
import os
import math


def graphtri(vertexsorted, mesh):
    print(np.array(mesh["faces"]).shape)

    fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))

    ax.plot_trisurf(*vertexsorted, triangles=mesh["faces"])

    plt.show()


def graphscatter(matrix):
    fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))

    ax.scatter(*matrix)

    plt.show()


def writefile(file, returnminmax=True):

    # Takes mesh file, writes all triangles to a file - also writes bounding box and resolution to a file

    extension = file.split(".")[1]

    mesh = getattr(getattr(trimesh.exchange, extension), "load_" + extension)(open(file, "rb"))

    vertexarray = np.reshape(np.array(mesh["vertices"]), (len(mesh["vertices"]), 3))

    minvalues = np.min(vertexarray, axis=0)
    maxvalues = np.max(vertexarray, axis=0)

    vertexsorted = [[i[j] for i in vertexarray] for j in range(3)]

    trianglepoints = [[vertexarray[i] for i in j] for j in mesh["faces"]]

    # graphtri(vertexsorted, mesh)

    with open("triangles.txt", "a+") as f:
        for x, y, z in [minvalues, maxvalues]:
            f.write(format(x, ".6f") + ", " + format(y, ".6f") + ", " + format(z, ".6f") + ", ")
        f.write("\n")
        for triangle in trianglepoints:
            for point in triangle:
                for coord in point:
                    f.write(format(coord, ".6f") + ", ")
            f.write("\n")

    if returnminmax:
        return minvalues, maxvalues


def pyvoxel(resolution): # use binary to generate voxel file
    os.system('./voxelize {}'.format(int(resolution)))


def readfile(): # test that writefile() works
    trianglearray = []
    trianglematrix = [[], [], []]
    with open("triangles.txt", "r") as f:
        triangles = f.readlines()
        for triangle in triangles:
            pointarray = triangle.split(", ")
            pointarray = pointarray[:9]
            pointarray = [float(i) for i in pointarray]
            for i in range(3):
                for j in range(3):
                    trianglematrix[i].append(pointarray[(j * 3) + i])
            pointarray = np.reshape(pointarray, (3, 3))
            trianglearray.append(pointarray)

    graphscatter(trianglematrix)


def readvoxelfile(): # returns array of voxels in voxel form
    voxelarray = []
    with open("voxels.txt", "r") as f:
        voxels = f.readlines()
        for voxel in voxels:
            voxel = voxel[:-1]
            pointarray = voxel.split(" ")
            pointarray = [int(i) for i in pointarray]
            voxelarray.append(pointarray)

    voxelsorted = [[i[j] for i in voxelarray] for j in range(3)]

    graphscatter(voxelsorted)
    return voxelarray


def voxeltranslator(voxelarray, divisions, minvals, convfactor): # makes voxel array into units
    # voxtodist = lambda vox, maxvox, minval, maxval: ((1 - ((maxvox - vox)/maxvox))*(minval/maxval) + minval)
    # voxtodist = lambda vox, maxvox, minval, maxval: ((((maxval - minval)/maxvox)*(vox + 0.5)) + minval)
    return [[voxtodist(j[i], divisions[i], minvals[i], convfactor) for i in range(2)] for j in voxelarray]


def voxtodist(vox, division, minval, convfactor):
    comcoord = vox + 0.5
    unitconvert = division*comcoord + minval
    return unitconvert*convfactor


def calcvoxels(voxelarray, axiscoords, mass):
    dist = lambda dim0, dim1, axisdim0, axisdim1: math.sqrt(((dim0 - axisdim0)**2 + (dim1 - axisdim1)**2))
    rarray = [dist(x[0], x[1], axiscoords[0], axiscoords[1]) for x in voxelarray]
    # mois = np.multiply(rarray, (mass/len(voxelarray)))
    return np.sum(np.square(rarray)*(mass/len(voxelarray)))


def calcvolume(voxelarray, division):
    volumepervoxel = division[0]*division[1]*division[2]
    return volumepervoxel*len(voxelarray)


def resetfiles():
    os.system('rm triangles.txt')
    os.system('rm voxels.txt')


if __name__ == "__main__":

    resetfiles()

    args = {
        "file": "teapot.stl",
        "resolution": 50,
        "axis": "y",
        "mass": 0.66,
        "height": 0.152
    }

    axisdict = {
    "x": 0,
    "y": 1,
    "z": 2
    }

    axis = axisdict.get(args["axis"])
    assert (axis != None)
    axis = int(axis)

    resolution = args["resolution"]
    minvalues, maxvalues = writefile(args["file"])

    print("minvalues {}, maxvalues {}".format(minvalues, maxvalues))

    # axiscoordsall = (maxvalues + minvalues)/2

    height = maxvalues[1] - minvalues[1]

    unitmeter = float(args["height"])/height

    # axiscoords = np.delete(axiscoordsall, axis, axis=0)

    pyvoxel(resolution)
    voxelarray = readvoxelfile()
    dimvoxelarray = np.delete(voxelarray, axis, axis=1)

    divisions = (maxvalues - minvalues)/resolution

    transvox = voxeltranslator(dimvoxelarray, np.delete(divisions, axis, axis=0), np.delete(minvalues, axis, axis=0), unitmeter)

    axiscoords = np.mean(transvox, axis=0)

    # print("transvox {}".format(np.array(transvox)))
    print("transvox max {}".format(np.max(transvox, axis=0)))
    print("transvox min {}".format(np.min(transvox, axis=0)))

    print("axiscoords {}".format(axiscoords))
    print("axiscoords in units {}".format(axiscoords/unitmeter))

    print("volume {}".format(calcvolume(transvox, divisions*unitmeter)))

    moi = calcvoxels(transvox, axiscoords, args["mass"])

    print(moi)
