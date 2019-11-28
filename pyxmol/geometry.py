#!/usr/bin/python

import numpy as np
import math

def to_radians(degrees):
    return degrees/180*math.pi

def to_degrees(radians):
    return radians*180.0/math.pi

def distance(atom1, atom2):
    return np.linalg.norm(atom1.r-atom2.r)

def angle_rad(atom1, atom2, atom3):
    v1 = atom1.r - atom2.r
    v2 = atom3.r - atom2.r
    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)
    return math.acos(np.dot(v1.transpose(), v2))

def angle_deg(atom1, atom2, atom3):
    return to_degrees(angle_rad(atom1,atom2,atom3))


def dihedral_rad(a, b, c, d):
    ba = a.r - b.r
    bc = c.r - b.r
    cd = d.r - c.r
    abc = -np.cross(ba,bc)
    bcd = np.cross(bc,cd)
    return math.atan2(
        np.cross(abc,bcd).dot(bc) / np.linalg.norm(bc),
        abc.dot(bcd)
    )

def dihedral_deg(a, b, c, d):
    return to_degrees(dihedral_rad(a, b, c, d))

def geom_center(atoms):
    import copy
    c = np.array([0,0,0]).reshape((3,1))
    for a in atoms:
        c += a.r
    c /= len(atoms)
    return c


def translate(atoms, dr):
    for a in atoms:
        a.r += dr

def rotation_matrix(axis, radians):
    x,y,z = ( axis / np.linalg.norm(axis) ).tolist()
    s = math.sin(radians)
    c = 1.0 - math.cos(radians)
    c1 = 1.0 - c
    m = np.array([
        [x*x*c1 +  c , x*y*c1 - z*s, x*z*c1 + y*s],
        [y*x*c1 + z*s, y*y*c1 +  c , y*z*c1 - x*s],
        [z*x*c1 - y*s, z*y*c1 + x*s, z*z*c1 +  c ]
         ])
    return m


def rotate(atoms, m, fixed_point=None):
    if fixed_point is None:
        for a in atoms:
            a.r = m.dot(a.r)
    else:
        dr = fixed_point - m.dot(fixed_point)
        for a in atoms:
            a.r = m.dot(a.r) - dr
