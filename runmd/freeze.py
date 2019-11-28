#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

from pdblib.num import *
from common.base import findcommon, range2list
from MD import AmberInputFile


def select_heavy(atoms):
    """
    :type atoms: list of Atom
    :return: list of Atom
    """
    heavy_atoms = [at for at in atoms if at.gettype() != "H"]
    return heavy_atoms


def select_all_secondary_structure_atoms(atoms, secondary_structure_string):
    """
    :type atoms: list of Atom
    :type secondary_structure_string: str
    :return: list of Atom
    """
    ss = range2list(secondary_structure_string)
    ss_atoms = [at for at in atoms if at.resi in ss]
    return ss_atoms


def select_heavy_secondary_structure_atoms(atoms, secondary_structure_string):
    """
    :type atoms: list of Atom
    :type secondary_structure_string: str
    :return: list of Atom
    """
    ss = set(range2list(secondary_structure_string))
    ss_atoms = [at for at in atoms if at.resi in ss and at.name in ["CA", "N", "C"] ]
    return ss_atoms


def select_backbone_atoms(atoms):
    """
    :type atoms: list of Atom
    :return: list of Atom
    """
    backbone_atoms = [at for at in atoms if at.name in ["CA", "N", "C"]]
    return backbone_atoms


def select_c_alpha_atoms(atoms):
    """
    :type atoms: list of Atom
    :return: list of Atom
    """
    backbone_atoms = [at for at in atoms if at.name == "CA"]
    return backbone_atoms


def freeze_atoms(parameters, atoms, spring_strength):
    """
    :type parameters: AmberInputFile
    :type spring_strength: float
    :type atoms: list of Atom
    """
    for i in range(len(atoms)):
        a = atoms[i]
        for j in range(i+1, len(atoms)):
            b = atoms[j]
            d = atdist(a, b)
            parameters.set_noe_restraint(a.atid, b.atid, 0.0, d, d, 999.0, spring_strength, spring_strength
                                         ,"%d:%s  <-->  %d:%s"%(a.resi, a.name, b.resi, b.name))


def freeze_heavy_atoms(parameters, atoms, spring_strength):
    freeze_atoms(parameters, select_heavy(atoms), spring_strength)


def freeze_all_secondary_structure_atoms(parameters, atoms, secondary_structure_string, spring_strength):
    freeze_atoms(parameters, select_all_secondary_structure_atoms(atoms, secondary_structure_string), spring_strength)


def freeze_heavy_secondary_structure_atoms(parameters, atoms, secondary_structure_string, spring_strength):
    freeze_atoms(parameters, select_heavy_secondary_structure_atoms(atoms, secondary_structure_string), spring_strength)


def freeze_backbone_atoms(parameters, atoms, spring_strength):
    freeze_atoms(parameters, select_backbone_atoms(atoms), spring_strength)


def freeze_c_alpha_atoms(parameters, atoms, spring_strength):
    freeze_atoms(parameters, select_c_alpha_atoms(atoms), spring_strength)
