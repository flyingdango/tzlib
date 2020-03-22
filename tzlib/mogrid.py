import c4d
from math import floor
from c4d.modules import mograph as mo


def pos2key(p):
    return str(round(p.x, 1)) + "_" + str(round(p.y, 1)) + "_" + str(round(p.z, 1))


def pos2index3d(p, gs):
    tmp = c4d.Vector(0)
    for i in xrange(3):
        if gs[i] > 0:
            tmp[i] = int(floor(p[i] / gs[i]))
    return tmp


def isFlatMatrix(m):
    return abs(m.v1.x + m.v2.y + m.v3.z - 3.0) < 0.01


def getGsByMarr(marr):
    sp = marr[0].off
    gs = c4d.Vector(999999)
    for j in xrange(1, len(marr)):
        p = marr[j].off
        for i in xrange(3):
            td = abs(sp[i] - p[i])
            if 0.01 < td < gs[i]:
                gs[i] = td
    for i in xrange(3):
        if gs[i] == 999999:
            gs[i] = 0
    return gs


def getGridBBX(nodes, attr_name="i3d"):
    minp, maxp = c4d.Vector(999999), c4d.Vector(-999999)
    for n in nodes:
        i = getattr(n, attr_name)
        for j in xrange(3):
            minp[j] = min(minp[j], i[j])
            maxp[j] = max(maxp[j], i[j])
    return minp, maxp, maxp - minp + c4d.Vector(1)


def gridIns_ez(pl, gs):
    a = pl[0]
    minp = None
    mind = float("inf")
    for i in xrange(1, len(pl)):
        b = pl[i]
        td = (a - b).GetLengthSquared()
        if td < mind:
            mind = td
            minp = b
    dv = minp - a
    dv = (abs(dv.x), abs(dv.y), abs(dv.z))
    if abs(max(dv) - min(gs.x, gs.y, gs.z)) > 1.0:
        return True
    return False


def getMarrFromClone(obj, minc=9):
    if not obj:
        return None
    md = mo.GeGetMoData(obj)
    if md is None:
        print "grid need a clone or matrix object"
        return None
    marr = md.GetArray(c4d.MODATA_MATRIX)
    farr = md.GetArray(c4d.MODATA_FLAGS)
    vl = map(lambda x: (x & (1 << 0)) and (not (x & (1 << 1))), farr)
    marr = [marr[i] for i in xrange(len(marr)) if vl[i]]
    if len(marr) < minc:
        print "no enough cube: ", len(marr)
        return None
    if not isFlatMatrix(marr[0]):
        print "not support clone rotation"
        return None
    return marr


def getMarrFromMd(obj, minc=9):
    md = mo.GeGetMoData(obj)
    if not md:
        return None
    marr = md.GetArray(c4d.MODATA_MATRIX)
    farr = md.GetArray(c4d.MODATA_FLAGS)
    vl = map(lambda x: (x & (1 << 0)) and (not (x & (1 << 1))), farr)
    marr = [marr[i] for i in xrange(len(marr)) if vl[i]]
    if len(marr) < minc:
        print "no enough cube: ", len(marr)
        return None
    if not isFlatMatrix(marr[0]):
        print "not support clone rotation"
        return None
    return marr


def isGridMode(gen):
    if not gen:
        return False
    if gen.GetType() not in (1018544, 1018545):
        print "only support clone/matrix"
        return False
    if gen[c4d.ID_MG_MOTIONGENERATOR_MODE] != 3:
        print "only support grid mode"
        return False
    return True


def getMarrPure(obj, minc=9):
    md = mo.GeGetMoData(obj)
    if not md:
        return None
    marr = md.GetArray(c4d.MODATA_MATRIX)
    if len(marr) < minc:
        print "no enough cube: ", len(marr)
        return None
    if not isFlatMatrix(marr[0]):
        print "not support clone rotation"
        return None
    return marr
