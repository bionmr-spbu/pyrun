

class Predicate:
    def __init__(self, pred):
        self.pred = pred

    def __and__(self, other):
        return Predicate(lambda x: self.pred(x) and other(x))

    def __or__(self, other):
        return Predicate(lambda x: self.pred(x) or other(x))

    def __xor__(self, other):
        return Predicate(lambda x: self.pred(x) ^ other(x))

    def __invert__(self):
        return Predicate(lambda x: not self.pred(x))

    def __call__(self, *args, **kwargs):
        return self.pred(*args, **kwargs)


class _aNamePlaceholder:
    def anyOf(self, names):
        return Predicate(lambda x: x.aName in names)

    def __eq__(self, id):
        return Predicate(lambda x: x.aName == id)


class _rNamePlaceholder:
    def anyOf(self, names):
        return Predicate(lambda x: x.rName in names)

    def __eq__(self, id):
        return Predicate(lambda x: x.rName == id)


class _aIdPlaceholder:
    def anyOf(self, ids):
        return Predicate(lambda x: x.aId in ids)

    def __lt__(self, id):
        return Predicate(lambda x: x.aId < id)

    def __le__(self, id):
        return Predicate(lambda x: x.aId <= id)

    def __gt__(self, id):
        return Predicate(lambda x: x.aId > id)

    def __ge__(self, id):
        return Predicate(lambda x: x.aId >= id)

    def __eq__(self, id):
        return Predicate(lambda x: x.aId == id)

class _rIdPlaceholder:
    def anyOf(self, ids):
        return Predicate(lambda x: x.rId in ids)

    def __lt__(self, id):
        return Predicate(lambda x: x.rId < id)

    def __le__(self, id):
        return Predicate(lambda x: x.rId <= id)

    def __gt__(self, id):
        return Predicate(lambda x: x.rId > id)

    def __ge__(self, id):
        return Predicate(lambda x: x.rId >= id)

    def __eq__(self, id):
        return Predicate(lambda x: x.rId == id)


class _cIdPlaceholder:
    def anyOf(self, ids):
        return Predicate(lambda x: x.cId in ids)

    def __eq__(self, id):
        return Predicate(lambda x: x.cId == id)


rId = _rIdPlaceholder()
cId = _cIdPlaceholder()
aId = _aIdPlaceholder()
aName = _aNamePlaceholder()
rName = _rNamePlaceholder()

def within(d, points):
    def check_distance(x, d, points):
        from pyxmol.geometry import distance
        for p in points:
            if distance(x,p) < d:
                return True
        return False
    return Predicate(lambda a: check_distance(a,d,points))