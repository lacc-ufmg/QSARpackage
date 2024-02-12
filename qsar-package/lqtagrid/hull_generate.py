#!/usr/bin/env python3
# coding: utf-8

from ..lqtagrid import matrix_generate
import numpy as np
import pandas as pd
import os


class HullGenerate():

    def __init__(self, extension, atp, directory, delta_angle = 4, initial_distance = 2.5, delta_r = 1,
                 total_layers = 6):
        # self.n_pat = re.compile(r'(.*[0-9]+).+')
        # dataFile = open(files).read().splitlines()
        dataFile = os.listdir(directory)

        if delta_angle is None:
            delta_angle = 4
        if initial_distance is None:
            initial_distance = 2.5
        if delta_r is None:
            delta_r = 1
        if total_layers is None:
            total_layers = 6

        self.atp = atp
        self.delta_angle = delta_angle
        self.initial_distance = initial_distance
        self.delta_r = delta_r
        self.total_layers = total_layers

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

        self.cCoulomb = []
        self.cLJ = []

        self.coulombMatrix = []
        self.ljMatrix = []
        self.points = []
        for mols in self.molecules:
            print('Calculating ' + mols)
            matrix = matrix_generate.MatrixGenerate(mols)
            matrix.hullGenerate(atp, delta_angle, initial_distance,
                                total_layers, delta_r)

            points, coulombAtom, ljAtom = matrix.getHullList()
            self.points.append([','.join(map(str, p)) for p in points])
            self.coulombMatrix.append(coulombAtom)
            self.ljMatrix.append(ljAtom)

    def generateHeader(self):
        header = []
        step = self.delta_angle
        delta_r = self.delta_r
        r = self.initial_distance
        for p in self.atp:
            header += ["{}_{}_{}_{}".format(r+l*delta_r,0,0,p) for l in range(self.total_layers)]
            for theta in range(step, 180, step):
                for phi in range(0, 360, step):
                    for a in [r+l*delta_r for l in range(self.total_layers)]:
                        header.append("{}_{}_{}_{}".format(a,theta,phi,p))
        return header

    def saveGrid(self, output):
        np.savetxt(output + '.mol', self.molecules, fmt='%s')
        header = self.generateHeader()
        dfCoulomb = pd.DataFrame(self.coulombMatrix, index=self.molecules)
        dfCoulomb.columns = [a+"_C" for a in header]
        dfLj = pd.DataFrame(self.ljMatrix, index=self.molecules)
        dfLj.columns = [a+"_LJ" for a in header]
        dfPoints = pd.DataFrame(self.points, index=self.molecules)
        dfPoints.columns = header
        dfPoints.to_csv(output + '_Points.csv', sep=';')
        df = dfCoulomb.join(dfLj)
        df.index = self.molecules
        df.to_csv(output + '.csv', sep=';')
