#!/usr/bin/env
# coding: utf-8

import math
import re
import os
from . import utils
#from . import utils
from numpy import arange
from numpy import zeros
from numpy import array
from numpy import reshape
from numpy import empty
from numpy import vstack
import pandas as pd
from rdkit import Chem
import openbabel as ob
import pybel


class MatrixGenerate():

    def __init__(self, filename):
        self.setX(filename)
        self.loadConstants(filename)
        self.loadAP()

    def setX(self, fileName):
        # filename_parts = fileName.split('.')
        # name = filename_parts[-2]
        print(fileName)
        molecule = pybel.readfile(fileName.split('.')[-1],fileName)

        self.X = empty((0,3))
        self.numberElements = 0
        self.cargas = []
        for i,c in enumerate(molecule):
            print(c)
            print(i)
            self.numberElements = len(c.atoms)
            filename_parts = (c.title).split('.')
            new_file = filename_parts[-2]+str(i)
            c.write("pdb",new_file+".pdb",overwrite=True)
            os.system('antechamber {} {} {} {} {} {} {} {} {} {}'.format(
                '-i',new_file+".pdb",
                '-fi',filename_parts[-1],
                '-o',new_file+".mol2",
                '-fo',
                'mol2',
                '-c',
                'bcc'))
            mol2 = pybel.readfile("mol2",new_file+".mol2")
            for m in mol2:
                for a in m:
                    self.X = vstack((self.X,a.coords))
                    self.cargas.append(a.partialcharge)
            os.remove(new_file+".mol2")
            os.remove(new_file+".pdb")
        self.minimos = self.X.min(axis=0)
        self.maximos = self.X.max(axis=0)
        self.m = len(self.X) # numero de linhas da matriz gerada

    def loadConstants(self, fileName):
        obConversion = ob.OBConversion() # objeto que converte entre formatos
        obConversion.SetInFormat(fileName.split('.')[-1])

        # abre molecula de referencia
        molecule = ob.OBMol()
        obConversion.ReadFile(molecule, fileName)
        ff = ob.OBForceField.FindForceField('GAFF')
        ff.Setup(molecule)
        ff.GetAtomTypes(molecule)
        self.typeConstants = []
        self.c6 = []
        self.c12 = []
        dfFF = pd.read_csv(os.path.join(os.path.dirname(__file__),"defaultsFiles/parametersVDWGAFF.csv"),index_col=0)
        for at in ob.OBMolAtomIter(molecule):
            tipo = at.GetData("FFAtomType").GetValue()
            self.typeConstants.append(tipo)
            s = dfFF.loc[tipo,"R"]/(2**(1/6))/5
            e = dfFF.loc[tipo,"E"]*4.184
            self.c6.append(4*e*s**6)
            self.c12.append(4*e*s**12)

        # molecule = Chem.MolFromMol2File(fileName,removeHs=False)
        # self.cargas = [a.GetProp('_TriposPartialCharge') for a in molecule.GetAtoms()]

    def loadAP(self):
        with open(os.path.dirname(__file__)+"/"+"defaultsFiles/AtomProva.atp") as f:
            input = f.readlines()
        currentLine = 1
        line = input[currentLine]
        self.ap = {}

        while len(input) > currentLine + 1:
            currentToken = 0
            currentLine += 1
            line = input[currentLine]
            tokens = line.split()
            #ap.insert(index, tokens[currentToken])
            group = tokens[currentToken]
            self.ap[group] = {}
            currentToken += 1
            #cargasap.insert(index, float(tokens[currentToken]))
            self.ap[group]['carga'] = float(tokens[currentToken])
            currentToken += 1
            #c6ap.insert(index, float(tokens[currentToken]))
            self.ap[group]['c6'] = float(tokens[currentToken])
            currentToken += 1
            #c12ap.insert(index, float(tokens[currentToken]))
            self.ap[group]['c12'] = float(tokens[currentToken])

    def gridGenerate(self, dimX, dimY, dimZ, atp, x0, y0, z0, step):
        self.DimX = dimX
        self.DimY = dimY
        self.DimZ = dimZ
        self.natp = len(atp)

        f = 138.935485
        nframes = self.m / self.numberElements
        #self.gridCoulomb = [[[[0 for x in range(self.natp)] for x in range(self.DimZ)]
        #                    for x in range(self.DimY)] for x in range(self.DimX)]

        #self.gridLJ = [[[[0 for x in range(self.natp)] for x in range(self.DimZ)]
        #                for x in range(self.DimY)] for x in range(self.DimX)]
        self.gridCoulomb = {}
        self.gridLJ = {}

        count = 0
        # esse loop roda sobre o número de sondas escolhidas
        for h in range(self.natp):
            #elem = self.search(self.ap, atp[h])
            #elem = self.ap.index(atp[h]) # encontra-se a posicao no vetor de elementos do elemento em questao
            # carrega-se as respectivas constantes
            q1 = self.ap[atp[h]]['carga'] #self.cargasap[elem]
            c6a = self.ap[atp[h]]['c6'] #self.c6ap[elem]
            c12a = self.ap[atp[h]]['c12'] #self.c12ap[elem]
            Vlj = 0
            Vc = 0
            npontos = 0
            #r1 = []
            r1 = [0.0,0.0,0.0]
            self.gridCoulomb[atp[h]] = {}
            self.gridLJ[atp[h]] = {}
            self.cutoffdistance = []
            # aqui começa o loop que gera as coordenadas cartesianas
            # acho que você pode gerar os pontos com o fecho convexo e substituir esses 3 loops por um loop
            # sobre os pontos gerados
            for i in arange(x0,self.DimX+x0+step,step):
                r1[0] = i
                self.gridCoulomb[atp[h]][i] = {}
                self.gridLJ[atp[h]][i] = {}
                for j in arange(y0,self.DimY+y0+step,step):
                    r1[1] = j
                    self.gridCoulomb[atp[h]][i][j] = {}
                    self.gridLJ[atp[h]][i][j] = {}
                    for k in arange(z0,self.DimZ+z0+step,step):
                        r1[2] = k
                        Vlj = 0
                        Vc = 0
                        npontos += 1
                        self.gridCoulomb[atp[h]][i][j][k] = {}
                        self.gridLJ[atp[h]][i][j][k] = {}
                        count += 1
                        # geradas as coordenadas cartesianas começa o loop sobre os átomos do PAC para
                        # calcular os descritores com base na distância entre a sonda e os átomos
                        for l in range(self.m):
                            r = utils.Distance(r1, self.X[l]) / 10
                            # self.cutoffdistance.append(r1)
                            index = l % self.numberElements
                            c6ij = math.sqrt(c6a * self.c6[index])
                            c12ij = math.sqrt(c12a * self.c12[index])

                            if r != 0:
                                Vlj = Vlj + (c12ij / (math.pow(r, 12))) - (c6ij / (math.pow(r, 6)))
                                Vc = Vc + f * float(q1) * float(self.cargas[l]) / r
                            else:
                                Vlj = float("inf")
                                Vc = float("inf")

                        self.gridCoulomb[atp[h]][i][j][k] = Vc / nframes
                        self.gridLJ[atp[h]][i][j][k] = Vlj / nframes
                        # self.gridLJ[atp[h]][i][j][k] = self.ljCut(Vlj / nframes,30)

    def getMatrix(self):
        textValuesCoulomb = ""
        textValuesLj = ""
        coulombMatrix = []
        ljMatrix = []
        count0 = 0
        count = 0
        for h in self.gridCoulomb:
            for i in self.gridCoulomb[h]:
                for j in self.gridCoulomb[h][i]:
                    for k in self.gridCoulomb[h][i][j]:
                        textValuesCoulomb += "%g\t" % (self.gridCoulomb[h][i][j][k])
                        textValuesLj += "%g\t" % (self.gridLJ[h][i][j][k])
                        coulombMatrix.append(self.gridCoulomb[h][i][j][k])
                        ljMatrix.append(self.gridLJ[h][i][j][k])
                        count += 1
        return textValuesCoulomb, textValuesLj, coulombMatrix, ljMatrix

    def ljCut(self,value,cut):
        calValue = value/4.18
        if calValue >= cut:
            calValue = cut + math.log10(calValue-(cut-1))
        calValue = calValue*4.18
        return calValue


    def activeConformation(self,matrix):
        # file = open(matrix,"r")
        # descritores = file.readline()
        # d = descritores.split()
        coord = []
        dfDesc = pd.read_csv(matrix,sep=';',header=None)
        regression_vector = array(list(dfDesc[1][:-1]))
        indT = list(dfDesc[1])[-1]
        d = list(dfDesc[0][:-1])
        for v in d:
            coord.append(v.split('_'))
        coordenadas = []
        probe = []
        types = []
        # file.close()
        for c in coord:
            coordenadas.append((float(c[0]),float(c[1]),float(c[2])))
            probe.append(c[3])
            types.append(c[4])
        nframes = (int(self.m/self.numberElements))
        V = zeros((nframes,len(coordenadas)))
        f = 138.935485
        soma = []
        confs = reshape(self.X,(nframes,self.numberElements,3))
        for i,c in enumerate(coordenadas):
            q1 = self.ap[probe[i]]['carga'] #self.cargasap[elem]
            c6a = self.ap[probe[i]]['c6'] #self.c6ap[elem]
            c12a = self.ap[probe[i]]['c12'] #self.c12ap[elem]
            for j,conf in enumerate(confs):
                for k,cAtom in enumerate(conf):
                    r = utils.Distance(list(c), cAtom) / 10
                    c6ij = math.sqrt(c6a * self.c6[k])
                    c12ij = math.sqrt(c12a * self.c12[k])

                    if r != 0:
                        if types[i] == "LJ":
                            V[j,i] += (c12ij / (math.pow(r, 12))) - (c6ij / (math.pow(r, 6)))
                        else:
                            V[j,i] += f * float(q1) * float(self.cargas[k]) / r
                    else:
                        V = float("inf")
                if types[i] == "LJ":
                    V[j,i] = self.ljCut(V[j,i],30)

        df = pd.DataFrame(V)
        df.to_csv('saida.csv', sep =';')
        y_conf = V.dot(regression_vector) + indT
        return y_conf
