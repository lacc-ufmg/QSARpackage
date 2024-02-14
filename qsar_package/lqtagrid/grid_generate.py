#!/usr/bin/env
# coding: utf-8

from . import matrix_generate
from numpy import arange
import os
from pandas import DataFrame
from openbabel import pybel


class GridGenerate():

    def __init__(self, extension, coordinates=(), dimensions=(),
        atp="NH3+", directory="", step=1.0):

        #dataFile = open(files).read().splitlines()
        dataFile = os.listdir(directory)


        self.molecules = []

        if extension != "mol2":
            format_file = "g09" if extension == "log" else extension
            files = [pybel.readfile(format_file,os.path.join(directory,x))
                    for x in dataFile if x.endswith(extension)]
            for f in files:
                mol = list(f)[0]
                mol.write('mol2',mol.title[:-(len(extension))]+"mol2", overwrite=True)
                self.molecules.append(mol.title[:-(len(extension))]+"mol2")
        else:
            self.molecules = [os.path.join(directory,x) for x in dataFile if x.endswith('mol2')]

        self.molecules.sort()

        matrices = []

        minimos = [999999.0,999999.0,999999.0]
        maximos = [-999999.0,-999999.0,-999999.0]
        for i in range(len(self.molecules)):
            matrix = matrix_generate.MatrixGenerate(self.molecules[i])
            minimos[0] = min(minimos[0],matrix.minimos[0])
            minimos[1] = min(minimos[1],matrix.minimos[1])
            minimos[2] = min(minimos[2],matrix.minimos[2])
            maximos[0] = max(maximos[0],matrix.maximos[0])
            maximos[1] = max(maximos[1],matrix.maximos[1])
            maximos[2] = max(maximos[2],matrix.maximos[2])
            matrices.append(matrix)

        if coordinates != None and coordinates != ():
            x0, y0, z0 = coordinates
        else:
            x0 = int(minimos[0])-2
            y0 = int(minimos[1])-2
            z0 = int(minimos[2])-2

        if coordinates != None and dimensions != ():
            dim_x, dim_y, dim_z = dimensions
        else:
            dim_x = int(maximos[0]-minimos[0])+4
            dim_y = int(maximos[1]-minimos[1])+4
            dim_z = int(maximos[2]-minimos[2])+4

        if not step == 1:
            I = int((dim_x/step)+(1/step-1))
            J = int((dim_y/step)+(1/step-1))
            K = int((dim_z/step)+(1/step-1))
        else:
            I = dim_x + 1
            J = dim_y + 1
            K = dim_z + 1

        n = len(atp)
        coulomb = ""
        lj = ""

        self.cCoulomb = []
        self.cLJ = []
        # esse loop gera o cabeçalho do arquivo de saida
        for l in range(n):
            for i in arange(x0,dim_x+x0+step,step):
                for j in arange(y0,dim_y+y0+step,step):
                    for k in arange(z0,dim_z+z0+step,step):
                        coulomb += "%.2f_%.2f_%.2f_%s_C: \t" % (i, j,
                                                                k, atp[l])

                        lj += "%.2f_%.2f_%.2f_%s_LJ: \t" % (i, j,
                                                            k, atp[l])
                        self.cCoulomb.append("%.2f_%.2f_%.2f_%s_C:" % (i, j, k, atp[l]))
                        self.cLJ.append("%.2f_%.2f_%.2f_%s_LJ:" % (i, j, k, atp[l]))
        self.output = coulomb + lj

        self.coulombMatrix = []
        self.ljMatrix = []
        for i,matrix in enumerate(matrices):
            print("Processing molecule {}  of {}".format(i+1,len(matrices)))
            matrix.gridGenerate(dim_x, dim_y, dim_z, atp, x0, y0, z0, step)
            # valuesCoulomb = matrix.getMatrix("C")
            # valuesLj = matrix.getMatrix("L")
            textValuesCoulomb, textValuesLj, coulombMatrix, ljMatrix = matrix.getMatrix()
            self.output += "\n" + textValuesCoulomb + textValuesLj
            self.coulombMatrix.append(coulombMatrix)
            self.ljMatrix.append(ljMatrix)

    def saveGrid(self,output):
        arq = open(output+'.txt', "w")
        arq.write(self.output)
        arq.close()
        dfCoulomb = DataFrame(self.coulombMatrix, columns = self.cCoulomb,
            index = [os.path.basename(m[:-5]) for m in self.molecules])
        dfLj = DataFrame(self.ljMatrix, columns = self.cLJ,
            index = [os.path.basename(m[:-5]) for m in self.molecules])
        df = dfCoulomb.join(dfLj)
        df.to_csv(output+'.csv', sep =';')
