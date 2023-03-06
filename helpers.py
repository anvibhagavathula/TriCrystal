import sys, csv, os
import os.path
import pandas
from operator import add
from crystals import Crystal, Atom, Element, distance_fractional, distance_cartesian
from shapely.geometry import Point, MultiPoint,  Polygon
from scipy.spatial import distance
from scipy.spatial.distance import cdist, pdist
from shapely.ops import nearest_points
from sklearn.neighbors import NearestNeighbors
from sklearn import neighbors
from sklearn.neighbors import KNeighborsRegressor
import numpy as np
import decimal
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import path
from datetime import date, datetime
import time


def seperate(my_crystal):
        """ Function isolates the z-coordinates (fractional) 
        from CRYSTAL and stores them in an array
        
        Args: 
        my_crystal: CRYSTAL structure 

        Returns: 
        Array with z-coordinates of bottom layer atoms"""

        nat = 0
        z_frac = []

        for atm in my_crystal:
                nat = nat + 1
                z_frac.append(atm.coords_fractional[2])

        z_sorted = sorted(z_frac)
        bottom_z = round(nat/2)
        answer = z_sorted[0:bottom_z]

        return answer



def lowest(my_crystal):
        """ Function finds the lowest z position in the crystal
        
        Args: 
        my_crystal: CRYSTAL structure 

        Returns: 
        Z coordinate in angstrom of lowest atom"""

        nat = 0
        z_cart = []

        for atm in my_crystal:
                nat = nat + 1
                z_cart.append(atm.coords_cartesian[2])

        return min(z_cart)



def highest(my_crystal):
        """ Function finds the highest z position in the crystal
        
        Args: 
        my_crystal: CRYSTAL structure 

        Returns: 
        Z coordinate in angstrom of highest atom"""

        nat = 0
        z_cart = []

        for atm in my_crystal:
                nat = nat + 1
                z_cart.append(atm.coords_cartesian[2])

        return max(z_cart)



def locate(test_atm,z_bot):
        """ Function checks if a given atom is from 
        the bottom or top layer
        
        Args: 
        test_atm: z-coordinate of a trial atom
        z_bot: array of z-coordinates of bottom layer atoms

        Returns: 
        Boolean indicating whether atom is from top
        or bottom layer """
        
        z = []

        test_atm = round(test_atm,2)

        for atm in z_bot:
                z.append(round(atm,2))

        if test_atm in z:
                answer = True
        else:
                answer = False

        return answer



def interlayer(my_crystal):
        """ Function finds the interlayer seperation 
        bettwen the layers
        
        Args: 
        my_crystal: CRYSTAL structure 

        Returns: 
        Interlayer seperation in bohr """

        atmz = []
        for atm in my_crystal:
            atmz.append(atm.coords_cartesian[2])

        z = sorted(atmz)
        z_bot = z[:len(z)//2]
        z_top = z[len(z)//2:]
        seperation = min(z_top) - max(z_bot)

        return seperation*0.529177249



def newcell(my_crystal,atoms,m,n):
        """ Function creates lattice vectors of a 
        rotated layer
        
        Args: 
        my_crystal: CRYSTAL structure 
        atoms: 
        m: int expansion parameter of desired supercell 
        n: int expansion parameter of desired supercell 

        Returns: 
        Supercell lattice vectors new_a; new_b """

        a1, a2, a3, alpha, beta, gamma = my_crystal.lattice_parameters

        old_a = np.array([0.5*a1*np.sqrt(3), -0.5*a2, 0])
        old_b = np.array([0, a2, 0])

        v1 = atoms[0] + (m+n)*old_a + (m+n)*old_b
        v2 = np.add(np.add(v1,(m+n)*old_a),n*old_b)
        v3 = np.add(np.add(v2,(m+n)*old_b),m*old_a)
        v4 = np.subtract(np.subtract(v3,(m+n)*old_a),n*old_b)

        new_a = v3 - v2
        new_b = v2 - v1

        return new_a, new_b, v1, v2, v3, v4



def rotcell(v1,v2,v3,v4,origin,R):
        """ Function gives the rotated coordinates 
        lattice vertices of newcell
        
        Args: 
        v1: vertex of supercell
        v2: vertex of supercell
        v3: vertex of supercell
        v4: vertex of supercell
        origin: origin atom 
        R: rotation matrix 

        Returns: 
        Rotated vertices"""
    
        v1 = v1 - origin
        v2 = v2 - origin
        v3 = v3 - origin
        v4 = v4 - origin
        vr1 = np.dot(v1,R)
        vr2 = np.dot(v2,R)
        vr3 = np.dot(v3,R)
        vr4 = np.dot(v4,R)

        return vr1,vr2,vr3,vr4



def poly(v1,v2,v3,v4):
        """ Function creates polygon / new unit 
        cell from vertices
        
        Args: 
        v1: vertex of supercell
        v2: vertex of supercell
        v3: vertex of supercell
        v4: vertex of supercell

        Returns: 
        Polygon vertices"""
        

        p1=v1[0], v1[1]
        p2=v4[0], v4[1]
        p3=v3[0], v3[1]
        p4=v2[0], v2[1]

        coords = [(p1), (p2), (p3), (p4)]
        ply = Polygon(coords)

        return ply, p1, p2, p3, p4



def inpoly(atz,ply):
        """ Function determines if an atom lies 
        within super unit cell
        
        Args: 
        atz: vertices of the supercell 
        ply: supercell polygon

        Returns: 
        Boolean indicating whether atom 
        is in supercell or not """

        atm = Point(atz[0], atz[1])
        check1 = atm.within(ply)
        check2 = atm.intersects(ply)
        if check1 == True:
                return check1
        if check2 == True:
                return check2



def central(ply):
        """ Function determines the centre of 
        unit cell
        
        Args: 
        ply: supercell polygon

        Returns: 
        Central coordinates of unit cell """

        center = ply.centroid

        return center



def swaped(layer):
        """ Function checks if rotation was 
        performed correctly
        
        Args: 
        layer: array of rotated atoms

        Returns: 
        Boolean indicating whether rotation
        was performed correctly """

        pos_count = 0
        neg_count = 0
        for atm in layer:
                if atm >= 0:
                        pos_count += 1
                else:
                        neg_count += 1
        if pos_count > neg_count:
                return True
        else:
                return False



def ntype(my_crystal):
        """ Function checks how many atomic species 
        are present in crystal
        
        Args: 
        my_crystal: CRYSTAL structure 

        Returns: 
        Integer number of atoms """
        
        nat = 0
        for atm in my_crystal.chemical_composition:
                nat+=1
        return nat



def bulk(my_crystal):
        """ Function creates arrays of top and 
        bottom atoms in the crystal
        
        Args: 
        my_crystal: CRYSTAL structure 

        Returns: 
        Arrays of top and bottom atom positions 
        and lists of chemical symbols for top and 
        bottom atoms """

        atoms_bot = []
        atoms_top = []
        atmxyz = []
        ele_bot = []
        ele_top = []
        z_bot = seperate(my_crystal)
        for atm in my_crystal:
                atm_frac = atm.coords_fractional
                atm_cart = atm.coords_cartesian
                atmxyz = atm_cart[0], atm_cart[1], atm_cart[2]
                if locate(atm_frac[2],z_bot) == True:
                        atoms_bot.append(atmxyz)
                        ELE = str(atm).lower()
                        ele_bot.append(ELE)
                else:
                        atoms_top.append(atmxyz)
                        ELE = str(atm).lower()
                        ele_top.append(ELE)
        return np.array(atoms_top), np.array(atoms_bot), ele_top, ele_bot



def scale_interlayer(vacuum_param): 
        """ Function scales the interlayer spacing between
        layers to ensure that they are ~3.6 angstroms apart
        scaled accordingly to the selected cell vacuum parameter 
        
        Args: 
        my_crystal: float of cell vacuum parameter 

        Returns: 
        Scaled value of interlayer spacing """

        spacing = (3.36/vacuum_param)*2
        spacing = np.round(spacing, 3)

        return spacing 