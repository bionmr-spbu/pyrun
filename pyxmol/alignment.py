import numpy as np
import math

def get_alignment(atoms_ref, atoms_moved):
    """
        Returns tuple (dr, U) where
        dr --- translation vector
        U  --- rotation matrix
    :type atoms_ref: list of Atom
    :type atoms_moved: list of Atom
    """
    if len(atoms_moved) != len(atoms_ref):
        raise Exception("Atom sets must have the same length (%d != %d)" % (len(atoms_moved), len(atoms_ref)))

    X = np.matrix([np.reshape(a.r,(3,)) for a in atoms_ref])
    Y = np.matrix([np.reshape(a.r,(3,)) for a in atoms_moved])

    if len(atoms_moved) != len(atoms_ref):
        raise Exception("Atom sets must have the same number of valid atoms (%d != %d)" % (len(atoms_moved), len(atoms_ref)))

    xc = np.mean(X, 0)
    yc = np.mean(Y, 0)

    for x in X:
        x -= xc

    for y in Y:
        y -= yc

    C = X.T * Y

    d = 1. if np.linalg.det(C) > 0 else -1.
    P = np.matrix([
        [ 1.,  0., 0.],
        [ 0.,  1., 0.],
        [ 0.,  0., d]])

    V, s, W = np.linalg.svd(C, full_matrices=True)

    U = V * P * W
    dr = U*yc.reshape((3,1))

    return np.array(xc.reshape((3,1))-dr),U


def get_rmsd(atoms_ref, atoms_moved, alignment=None):
    """
    :type atoms_ref: list of Atom
    :type atoms_moved: list of Atom
    :type alignment: tuple
    :return: double
    """
    if len(atoms_moved) != len(atoms_ref):
        raise Exception("Atom sets must have the same length (%d != %d)" % (len(atoms_moved), len(atoms_ref)))

    rmsd = 0

    if alignment is None:
        for i in range(len(atoms_ref)):
            r1 = atoms_ref[i].r
            r2 = atoms_moved[i].r
            rmsd += np.linalg.norm(r1 - r2)**2
    else:
        dr, U = alignment
        for i in range(len(atoms_ref)):
            r1 = atoms_ref[i].r
            r2 = atoms_moved[i].r
            rmsd += np.linalg.norm(r1 - (dr + U*r2))**2

    return math.sqrt(rmsd/len(atoms_moved))

def align_with(atoms, alignment):
    dr, U = alignment
    for at in atoms:
        at.r = dr + U*at.r