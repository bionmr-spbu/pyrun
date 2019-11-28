import numpy
import copy

def bisect_contains(a, x):
    import bisect
    i = bisect.bisect_left(a, x)
    return i != len(a) and a[i] == x

def format_atom_name(aName):
    if len(aName)>3:
        return aName
    return " %-3s"%aName


def std_unique(a):
    i = 0
    for j in range(1,len(a)):
        if cmp(a[i],a[j]) != 0:
            i+=1
            a[i] = a[j]
        j += 1
    del a[i+1:]


def merge(left,right):
    result=[]
    i,j=0,0
    while i<len(left) and j<len(right):
        if left[i] < right[j]:
            result.append(left[i])
            i+=1
        else:
            result.append(right[j])
            j+=1
    while (i<len(left)):
        result.append(left[i])
        i+=1
    while (j<len(right)):
        result.append(right[j])
        j+=1
    return result


class Atom:
    def __init__(self,name="", id=None, r=None):
        self.name = name
        self.id = id
        self.r = r
        if (self.r is not None): self.r.shape=(3,1)
        self.residue = None

    def __repr__(self):
        return ("Atom("
                "name="+repr(self.name)+
                ", id="+repr(self.id)+
                ", r="+repr(self.r)+
                ")"
                )

    def __cmp__(self, other):
        return cmp((self.residue, self.id, id(self)), (other.residue, other.id, id(other)))

    def __deepcopy__(self, memo):
        if id(self) in memo: return memo[id(self)]
        a = Atom(name=self.name, id=self.id, r=self.r)
        memo[id(self)] = a
        return a

    def __copy__(self):
        return Atom(name=self.name, id=self.id, r=self.r)

    @property
    def aName(self):
        return self.name

    @property
    def aId(self):
        return self.id

    @property
    def rName(self):
        if self.residue:
            return self.residue.name
        return None

    @property
    def rId(self):
        if self.residue:
            return self.residue.id
        return None

    @property
    def cId(self):
        if self.residue:
            if self.residue.chain:
                return self.residue.chain.id
        return None

    def writeAsPdb(self, output):

        def default(value, d):
            return value if value is not None else d

        fmt = 'ATOM  %6s%4s%1s%-4s%1s%-5s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s  \n'
        output.write(fmt%(
            "%5d"%self.aId,
            format_atom_name(self.aName),
            '',  # location id
            default(self.rName,"XXX"),
            default(self.cId," "),
            "%4d"%default(self.rId, -1),
            self.r[0],
            self.r[1],
            self.r[2],
            1.0,   # occupancy
            0.0,   # temp factor
            '',    # elem
            '',    # charge
        ))


class Residue:
    def __init__(self, name="", id=None, atoms=None):
        self._atoms = []
        self._atom_name_map = {}
        self.name = name
        self.id = id
        self.chain = None
        self.add(atoms)

    @property
    def atoms(self):
        for a in self._atoms:
            yield a

    @property
    def asAtoms(self):
        return AtomSelection( self.atoms )

    def add(self, atoms):
        if atoms:
            for a in atoms:
                self._atom_name_map[a.name] = a
                a.residue = self
            atoms.sort()

            self._atoms = merge(self._atoms, atoms)
            std_unique(self._atoms)


    def remove(self, atoms):
        for x in atoms:
            del self._atom_name_map[x.name]
            self._atoms.remove(x)

    def __cmp__(self, other):
        return cmp((self.chain, self.id, id(self)), (other.chain, other.id, id(other)))

    def __deepcopy__(self, memo):
        if id(self) in memo: return memo[id(self)]
        r = Residue(name=self.name, id=self.id, atoms=copy.deepcopy(self._atoms, memo))
        memo[id(self)] = r
        return r

    def __copy__(self):
        return Residue(name=self.name, id=self.id, atoms=self._atoms)

    def __repr__(self):
        s  = "Residue("
        s += "id=" + repr(self.id)
        s += ", name=" + repr(self.id)
        s += ", atoms=[\n"
        for a in self.atoms:
            s += repr(a)+", \n"
        s += "]"
        s += ")"
        return s

    @property
    def rName(self):
        return self.name

    @property
    def rId(self):
        return self.id

    @property
    def cId(self):
        if self.chain:
            return self.chain.id
        return None

    def __getitem__(self, aname):
        """:rtype: Atom"""
        return self._atom_name_map[aname]

    def writeAsPdb(self, output):
        for a in self._atoms:
            a.writeAsPdb(output)


class Chain:
    def __init__(self, id='', residues=None):
        self._residue_id_map = {}
        self._residues = []
        self.id = id
        self.add(residues)

    @property
    def atoms(self):
        for r in self.residues:
            for a in r.atoms:
                yield a

    @property
    def residues(self):
        for r in self._residues:
            yield r

    @property
    def asAtoms(self):
        return AtomSelection( self.atoms )

    @property
    def asResidues(self):
        return ResidueSelection( self.residues )


    def add(self, residues):
        if residues:
            for r in residues:
                self._residue_id_map[r.id] = r
                r.chain = self

            self._residues = merge(self._residues, sorted(residues))
            std_unique(self._residues)

    def remove(self, residues):
        for x in residues:
            del self._residue_id_map[x.name]
            self._residues.remove(x)

    def __cmp__(self, other):
        if id(self)!=id(other) and self.id==other.id:
            c = cmp(len(self._residues),len(other._residues))
            if c==0 and len(self._residues)!=0:
                return cmp(self._residues[0].rId,other._residues[0].rId)
            return c
        return cmp( (self.id, id(self)), (other.id, id(other)) )

    def __deepcopy__(self, memo):
        if id(self) in memo: return memo[id(self)]
        c = Chain(id=self.id, residues=copy.deepcopy(self._residues, memo))
        memo[id(self)] = c
        return c

    def __copy__(self):
        return Chain(id=self.id, residues=self.residues)

    def __repr__(self):
        s  = "Chain("
        s += "id="+repr(self.id)
        s += ", residues=[\n"
        for r in self.residues:
            s += repr(r)+", \n"
        s += "]"
        s += ")"
        return s

    @property
    def cId(self):
        return self.id


    def __getitem__(self, rid):
        """:rtype: Residue"""
        return self._residue_id_map[rid]

    def __contains__(self, rid):
        return rid in self._residue_id_map

    def writeAsPdb(self, output):
        for r in self._residues:
            r.writeAsPdb(output)
        output.write("TER\n")

class AtomSelection:
    def __init__(self, atoms=None):
        self.elements = [a for a in atoms]
        self.elements.sort()
        std_unique(self.elements)
        self.hash = None

    def get_hash(self):
        if self.hash is None:
            self.hash = self.calc_hash()
        return self.hash

    def calc_hash(self):
        import hashlib
        m = hashlib.md5()
        for a in self.elements:
            m.update(str(a.rId))
            m.update(str(a.rName))
            m.update(str(a.aName))
        return m.hexdigest()

    @property
    def asResidues(self):
        """
        :rtype ResidueSelection
        """
        reses = [a.residue for a in self.elements]
        std_unique(reses)
        return ResidueSelection(reses)

    @property
    def asChains(self):
        """
        :rtype ChainSelection
        """
        return self.asResidues.asChains

    def __add__(self, other):
        """
        :rtype: AtomSelection
        """
        elems = self.elements[:]
        elems.extend(other)
        return AtomSelection(elems)

    def __contains__(self, item):
        return bisect_contains(self.elements, item)

    def __getitem__(self, item):
        if isinstance(item, int):
            return self.elements[item]
        else:
            return AtomSelection(self.elements[item])

    def __len__(self):
        return self.elements.__len__()

    def __rshift__(self, predicate):
        return AtomSelection([a for a in self.elements if predicate(a)])

    def writeAsPdb(self, output):
        for el in self.elements:
            el.writeAsPdb(output)

class ResidueSelection:
    def __init__(self, residues=None):
        self.elements = [a for a in residues]
        self.elements.sort()
        std_unique(self.elements)

    @property
    def asAtoms(self):
        result = []
        for r in self.elements:
            result += list(r.atoms)
        return AtomSelection(result)

    @property
    def asChains(self):
        chains = [r.chain for r in self.elements]
        std_unique(chains)
        return ChainSelection(chains)

    def __add__(self, other):
        elems = self.elements[:]
        elems.extend(other)
        return ResidueSelection(elems)

    def __contains__(self, item):
        return bisect_contains(self.elements, item)

    def __getitem__(self, item):
        if isinstance(item, int):
            return self.elements[item]
        else:
            return ResidueSelection(self.elements[item])

    def __len__(self):
        return self.elements.__len__()

    def __rshift__(self, predicate):
        return ResidueSelection([el for el in self.elements if predicate(el)])

    def writeAsPdb(self, output):
        for el in self.elements:
            el.writeAsPdb(output)

class ChainSelection:
    def __init__(self, chains=None):
        self.elements = [a for a in chains]
        self.elements.sort()
        std_unique(self.elements)

    @property
    def asResidues(self):
        result = []
        for c in self.elements:
            result += list(c.residues)
        return ResidueSelection(result)

    @property
    def asAtoms(self):
        return self.asResidues.asAtoms

    def __add__(self, other):
        elems = self.elements[:]
        elems.extend(other)
        return ChainSelection(elems)

    def __contains__(self, item):
        return bisect_contains(self.elements, item)

    def __getitem__(self, item):
        if isinstance(item, int):
            return self.elements[item]
        else:
            return ChainSelection(self.elements[item])

    def __len__(self):
        return self.elements.__len__()

    def __rshift__(self, predicate):
        return ResidueSelection([el for el in self.elements if predicate(el)])

    def writeAsPdb(self, output):
        for el in self.elements:
            el.writeAsPdb(output)

def pdbLineToTuple(line):
    """
    Pdb format:
    1-6:   "ATOM  ";
    7-11: ATOM ID
    13-16: Atom Name;
    17: Location indicator
    18-20: resname;
    22: chain identifier
    23-26: resnum;
    27: icode
    31-38: X,real(8.3)
    39-46: Y,real(8.3)
    47-54: Z,real(8.3)
    55-60: Occupancy
    61-66: TempFactor
    73-76: segID
    77-78: element
    79-80: Charge
    """
    assert(line.startswith("ATOM"))
    a = Atom(
        id=int(line[6:12]),
        name=line[12:16].strip(),
        r=numpy.array([
            float(line[30:38]),  # x
            float(line[38:46]),  # y
            float(line[46:54])   # z
        ])
    )
    resId = int(line[22:27])      # residue id
    resName = line[17:20].strip()   # residue name
    chainId = line[21]              # chain id
    return a, resId, resName, chainId

def readPdb(filename):
    models = []
    cids = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
    last_cid = 0
    last_chain_id = ""
    with open(filename) as f:
        chains = None
        chain = None
        residue = None
        residues = []
        atoms = []
        for line in f:
            if line.startswith("ATOM"):
                a, resId, resName, chainId = pdbLineToTuple(line)

                if chains is None:
                    chains = []
                    models.append(chains)

                if chain is None or last_chain_id != chainId:
                    if chain is not None: chain.add(residues)
                    residues = []
                    last_chain_id = chainId
                    if chainId == " ":
                        if (last_cid < len(cids)):
                            chainId = cids[last_cid]
                            last_cid += 1
                        else:
                            chainId = "z"
                    chain = Chain()
                    chain.id = chainId
                    chains.append(chain)

                if residue is None or residue.id != resId:
                    if residue is not None: residue.add(atoms)
                    atoms = []
                    residue = Residue()
                    residue.id = resId
                    residue.name = resName
                    residues.append(residue)

                atoms.append(a)

            if line.strip() in ("TER", "ENDMDL"):
                if residue is not None: residue.add(atoms)
                atoms = []
                if chain is not None: chain.add(residues)
                residues = []
                residue = None
                chain = None

        if residue is not None: residue.add(atoms)
        if chain is not None: chain.add(residues)

    return [ ChainSelection(cns) for cns in models ]

