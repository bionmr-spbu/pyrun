import numpy as np
from pdblib.num import *
from runmd.freeze import *
from pdblib.base import Atom


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

    X = np.matrix(getmat(atoms_ref))
    Y = np.matrix(getmat(atoms_moved))

    if len(atoms_moved) != len(atoms_ref):
        raise Exception("Atom sets must have the same number of valid atoms (%d != %d)" % (len(atoms_moved), len(atoms_ref)))


    xc = np.mean(X, 0)
    yc = np.mean(Y, 0)

    #
    # print "xc:", xc.shape
    # print "yc:", yc.shape
    #
    #
    #
    for x in X:
        x -= xc

    for y in Y:
        y -= yc


    # C --- covariance matrix
    # print "X shape = ", X.shape
    # print "Y shape = ", Y.shape

    C = X.T * Y

    d = 1. if np.linalg.det(C) > 0 else -1.
    P = np.matrix([
        [ 1.,  0., 0.],
        [ 0.,  1., 0.],
        [ 0.,  0., d]])

    V, s, W = np.linalg.svd(C, full_matrices=True)

    U = W.T * P * V.T
    dr = yc*U

    return np.array(xc-dr), U


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
            r1 = np.array([atoms_ref[i].r])
            r2 = np.array([atoms_moved[i].r])
            rmsd += np.linalg.norm(r1 - r2)**2
    else:
        dr, U = alignment
        for i in range(len(atoms_ref)):
            r1 = np.array([atoms_ref[i].r])
            r2 = np.array([atoms_moved[i].r])
            rmsd += np.linalg.norm(r1 - (dr + r2*U))**2

    return sqrt(rmsd/len(atoms_moved))


def align_with(atoms_moved, alignment):
    """
    Aligns atoms_moved atoms to alignment
    :param atoms_moved: list of moved atoms
    :type atoms_moved: list of Atom
    :param alignment: transformation applied to moved atoms
    :type alignment: tuple
    :return:
    """

    dr, U = alignment
    for at in atoms_moved:
        r2 = np.array([at.r])
        r = dr + r2*U
        at.r = (r[0, 0], r[0, 1], r[0, 2])