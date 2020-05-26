# coding=utf-8
import c4d


class Node(object):

    def __init__(self, _key=None, _pos=None, _index=-1):
        self.key = _key
        self.i = _index
        self.pos = _pos
        self.nbl = None
        self.parent = None
        self.walkable = True
        self.g = 0
        self.h = 0
        self.f = 0

    def _setF(self, v):
        pass

    def _getF(self):
        return self.g + self.h

    def __add__(self, other):
        return self.pos + other.pos

    def __sub__(self, other):
        return self.pos - other.pos

    def __hash__(self):
        return hash(self.key)

    def __repr__(self):
        return "Node:" + str(self.i)

    f = property(_getF, _setF)


def mhdis(a, b):
    v = a - b
    return abs(v.x) + abs(v.y) + abs(v.z)


def setWalkable(l, b):
    for n in l:
        n.walkable = b


def dijkstra(nodes, sn, enl, maxlen, maxpath):
    # init
    openl = set([])
    found = set([])
    for n in nodes:
        n.g = float("inf")
        n.parent = None
        openl.add(n)
    if sn.key == "root":
        for n in sn.nbl:
            n.g = 0
            n.parent = sn
    else:
        sn.g = 0
    # find
    r = []
    encnt = len(enl)
    while len(openl) > 0:
        cn = min(openl, key=lambda n: n.g)
        openl.discard(cn)
        if cn.g == float("inf"):
            break
        if cn in enl:
            found.add(cn)
            path = genPath(sn, cn)
            if len(path) > maxlen:
                break
            r.append(path)
            if len(r) >= maxpath:
                return r

        if len(found) >= encnt:
            return r

        for nb in cn.nbl:
            if nb.walkable == False or nb in found:
                continue
            newg = cn.g + (cn - nb).GetLength()
            if newg < nb.g:
                nb.g = newg
                nb.parent = cn
    return r


def astar(nodes, sn, en, maxlen):
    for n in nodes:
        n.h = mhdis(n, en)
        n.g = float("inf")
        n.parent = None
    sn.g = 0
    openl = set([sn])
    closel = set([])
    while len(openl) > 0:
        cn = min(openl, key=lambda n: n.f)
        openl.discard(cn)
        closel.add(cn)
        if cn is en:
            path = genPath(sn, cn)
            if len(path) > maxlen:
                break
            return path

        for nb in cn.nbl:
            if nb.walkable == False or (nb in closel):
                continue
            newg = cn.g + (cn - nb).GetLength()
            if newg < nb.g or (nb not in openl):
                openl.add(nb)
                nb.g = newg
                nb.parent = cn
    return []


def genPath(sn, en):
    tmp = []
    cn = en
    while cn.parent != None:
        tmp.append(cn)
        cn = cn.parent
    if sn.key != "root":
        tmp.append(sn)
    return tmp
