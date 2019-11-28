from copy import deepcopy

from .alignment import get_alignment, get_rmsd, align_with
from pyxmol.base import *
from .predicate import *


def glue_chain(residues1, res_id_1, residues2, res_id_2, method='peptide_plain'):
    def select_within_residue(ats, resid):
        result = ats >> (aName.anyOf(["CA", "C", "N"]) & (rId==resid))
        assert(len(result)==3)
        return result

    def select_peptide_plain(ats, resid):
        result = ats >> ((aName.anyOf(["CA", "C"]) & (rId==resid-1)) |
                       (aName.anyOf(["CA", "N"]) & (rId == resid))
                       )
        assert(len(result)==4)
        return result

    if method=='peptide_plain':
        select = select_peptide_plain
    elif method=='within_residue':
        select = select_within_residue
    else:
        raise KeyError("no glue method `%s`" % method)

    head = residues1 >> (rId <= res_id_1)
    tail = residues2 >> (rId >= res_id_2-1)

    ats1 = select(head.asAtoms, res_id_1)
    ats2 = select(tail.asAtoms, res_id_2)

    assert(len(ats1) == len(ats2))
    assert(len(ats2) >= 3 )

    for i in range(len(ats1)):
        assert(ats1[i].aName == ats2[i].aName)

    al = get_alignment(ats1, ats2)

    print(get_rmsd(ats1,ats2,al))

    reses1 = ResidueSelection([ deepcopy(r) for r in head if r.id < res_id_1  ])
    reses2 = ResidueSelection([ deepcopy(r) for r in tail if r.id >= res_id_2 ])

    align_with(reses2.asAtoms, al)

    for i, r in enumerate(reses2):
        r.id = i+res_id_1

    last_atom_id = reses1[-1].asAtoms[-1].id+1
    for i, a in enumerate(reses2.asAtoms):
        a.id = i + last_atom_id

    return Chain(residues=reses1+reses2)
