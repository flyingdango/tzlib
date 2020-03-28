# coding=utf-8
import c4d
import math
import utils as tzu
from c4d import utils
from c4d.utils import Neighbor
from collections import namedtuple
from itertools import combinations
from searcher import PointSearcher
"""
hasEdgeAB(a, b, tmp)
判断序号a,b是否在tmp里
tmp = [(0,3),(6,8)...]

pointPointsMap(el,ne=None)
建立边引索,结果的字典给定一个点序号,返回与之相连的点序号
ne为边是否为ngonEdge的标记列表
ne不给返回{0:(1,2,3),1:(0,2,5)..}
ne给,返回{0:{"i":(1,2,3),"ne":(True,True,False)}...}

getAllEdges(obj, nb=None, fl=None, nec=None)
返回{"el":edgeList,"ne":isneList}
el:物体所有边((1,2),(3,4)...)
ne:标记对应边是否为ngonEdge(False,False....)

selEdgePoint(obj, th=175)
选中边上多余的点,th角度阈值
自带的modelingCheck有时候无效,而且不好用脚本调用

el = ((0,1),(0,3)...)
结果可以根据任意一个点序,返回相连的所有点序
{0:(1,3,7,9),1:(2,5,6)...}

getCPolyEdge(f, i)
i:0-3,返回a-b,b-c,c-d,d-a

getVectorAvg(pil, pl)
根据面的所有点序获取面中心
pil：(0,3,5,6,7)

getCPolyNormal(f, pl)
获取cpoly的法线

deconstructSpline(obj, pa=None, th=5)
析构Spline
obj:spline对象,单/多分段都可以
pa:obj.GetAllPoints(),可以不给
th:如果>0,那么在这个范围内的点视为一点
返回{"p":pa,"e":el}
pa:点列表[c4d.Vector()...]
el:边列表((0,1),(1,2),(2,3)...)


"""


class Edge(object):

    def __init__(self, a, b):
        self.a = a
        self.b = b
        self.isNgonEdge = False
        self.isBoundary = False

    def __str__(self):
        return str(self.a) + "_" + str(self.b)

    def __repr__(self):
        return "Edge(" + str(self.a) + "," + str(self.b) + ")"

    def find(self, i):
        if i == self.a:
            return self.b
        if i == self.b:
            return self.a
        return -1

    def __eq__(self, other):
        return isEqEdgeEE(self, other)

    def __hash__(self):
        return hash((self.a, self.b))


class Ngon(object):

    def __init__(self, pil, eil, fil):
        self.pil = pil
        self.eil = eil
        self.fil = fil
        self.pc = len(pil)

    def __str__(self):
        return "pil:" + str(self.pil) + "fil:" + str(self.fil)

    def __repr__(self):
        return "Ngon" + str(self.pil)

    def __hash__(self):
        return hash(tuple(self.pil))


def hasEdgeAB(a, b, el):
    for i, e in enumerate(el):
        if isEqEdgeABE(a, b, e):
            return i
    return -1


def hasEdgeE(e, el):
    for i, te in enumerate(el):
        if isEqEdgeEE(e, te):
            return i
    return -1


def isEqEdgeEE(a, b):
    return a.a == b.a and a.b == b.b or a.a == b.b and a.b == b.a


def isEqEdgeABE(a, b, e):
    return a == e.a and b == e.b or a == e.b and b == e.a


def isEqEdgeABAB(a, b, c, d):
    return a == c and b == d or a == d and b == c


def pointPolyMap(fl):
    tmp = {}
    for i, f in enumerate(fl):
        tmp.setdefault(f.a, set([])).add(i)
        tmp.setdefault(f.b, set([])).add(i)
        tmp.setdefault(f.c, set([])).add(i)
        tmp.setdefault(f.d, set([])).add(i)
    for k in tmp:
        tmp[k] = tuple(tmp[k])
    return tmp


def pointNgonMap(nl):
    tmp = {}
    for i, n in enumerate(nl):
        for pid in n.pil:
            tmp.setdefault(pid, []).append(i)
    return tmp


def polyNgonMap(nl):
    tmp = {}
    for i, n in enumerate(nl):
        for j, fi in enumerate(n.fil):
            tmp[fi] = (i, j)
    return tmp


def pointEdgeMap(el):
    tmp = {}
    for i, e in enumerate(el):
        tmp.setdefault(e.a, []).append(i)
        tmp.setdefault(e.b, []).append(i)
    return tmp


def pointPointMap(fl):
    # include ngon link
    tmp = {}
    for f in fl:
        if f.IsTriangle():
            tmp.setdefault(f.a, set([])).add(f.b)
            tmp[f.a].add(f.c)
            tmp.setdefault(f.b, set([])).add(f.a)
            tmp[f.b].add(f.c)
            tmp.setdefault(f.c, set([])).add(f.a)
            tmp[f.c].add(f.b)
        else:
            tmp.setdefault(f.a, set([])).add(f.b)
            tmp[f.a].add(f.d)
            tmp.setdefault(f.b, set([])).add(f.a)
            tmp[f.b].add(f.c)
            tmp.setdefault(f.c, set([])).add(f.b)
            tmp[f.c].add(f.b)
            tmp.setdefault(f.d, set([])).add(f.a)
            tmp[f.d].add(f.c)
    for k in tmp:
        tmp[k] = tuple(tmp[k])
    return tmp


def edgeNgonMap(nl):
    tmp = {}
    for i, n in enumerate(nl):
        for ei in n.eil:
            tmp.setdefault(ei, set())
            tmp[ei].add(i)
    for k in tmp:
        tmp[k] = tuple(tmp[k])
    return tmp


def pointPointMap2(el):
    # not include ngon link
    tmp = {}
    for e in el:
        if e.isNgonEdge:
            continue
        tmp.setdefault(e.a, set([])).add(e.b)
        tmp.setdefault(e.b, set([])).add(e.a)
    for k in tmp:
        tmp[k] = tuple(tmp[k])
    return tmp


def pointNearMap(pl, th=1.0):
    tmp = {}
    ser = PointSearcher()
    ser.init(pl)
    for i, p in enumerate(pl):
        r = ser.search(p, r=th, n=20, HQ=True, sort=False)
        if not r:
            continue
        tmp[i] = set(r)
        tmp[i].discard(i)
    for k in tmp:
        tmp[k] = tuple(tmp[k])
    return tmp


def pointPointMap3(pl, el):
    # 破碎专用
    tmp = {}
    ppm = pointPointMap2(el)
    pnm = pointNearMap(pl)
    for k, v in ppm.iteritems():
        if k in tmp:
            continue
        npl = pnm.get(k)
        if not npl:
            tmp[k] = v
            continue
        npl = set(npl)
        npl.add(k)
        nbl = []
        for j in npl:
            nbl.extend(list(ppm[j]))
        nbl = tuple(set(nbl))
        for j in npl:
            tmp[j] = nbl
    return tmp


def getAllFaceNormal(pl, fl):
    tmp = []
    for f in fl:
        tmp.append(getCPolyNormal(f, pl))
    return tmp


def getAllPointNormal(pl, fnl, pfm):
    tmp = []
    for i, p in enumerate(pl):
        fil = pfm[i]
        if len(fil) == 0:
            tmp.append(c4d.Vector(0, 0, 1))
        tmp.append(getVectorAvg(fil, fnl).GetNormalized())
    return tmp


def getAllEdges(obj, nb=None, fl=None, nec=None):
    if not nb:
        nb = Neighbor()
        nb.Init(obj)
    if not fl:
        fl = obj.GetAllPolygons()
    if not nec:
        nec = obj.GetNgonEdgesCompact()
    el = []
    for i, f in enumerate(fl):
        pf = nb.GetPolyInfo(i)
        for ei in xrange(4):
            if pf["mark"][ei] or ei == 2 and f.c == f.d:
                continue
            a, b = getCPolyEdge(f, ei)
            te = Edge(a, b)
            te.isBoundary = True if pf["face"][ei] == c4d.NOTOK else False
            te.isNgonEdge = False if (nec[i] & (1 << ei) == 0) else True
            el.append(te)
    return el


def selEdgePoint(obj, pl=None, th=175):
    if not pl:
        pl = obj.GetAllPoints()
    el = getAllEdges(obj)
    esm = pointEdgeMap(el)

    sel = obj.GetPointS()
    sel.DeselectAll()
    for spi, d in esm.iteritems():
        t = [i for i in xrange(len(d)) if not el[d[i]].isNgonEdge]
        if len(t) == 2:
            e1, e2 = el[d[t[0]]], el[d[t[1]]]
            a, b, c = pl[spi], pl[e1.find(spi)], pl[e2.find(spi)]
            if utils.Deg(utils.VectorAngle(b - a, c - a)) > th:
                sel.Select(spi)
    c4d.EventAdd()


def getNgonTrans(obj):
    ngonCount, polyTrans = obj.GetPolygonTranslationMap()
    return obj.GetNGonTranslationMap(ngonCount, polyTrans)


def getAllNgons(obj, nb=None, fl=None, nec=None, nts=None):
    if not nb:
        nb = Neighbor()
        nb.Init(obj)
    if not fl:
        fl = obj.GetAllPolygons()
    if not nec:
        nec = obj.GetNgonEdgesCompact()
    if not nts:
        nts = getNgonTrans(obj)
    tmp = []
    el = []
    e_dict = {}
    for tfil in nts:
        # get all visiable edge
        tel = []
        teil = []
        for fi in tfil:
            f = fl[fi]
            #pf = nb.GetPolyInfo(fi)
            for ei in xrange(4):
                if ei == 2 and f.c == f.d:
                    continue
                if nec[fi] & (1 << ei) != 0:
                    continue
                a, b = getCPolyEdge(f, ei)
                #
                tek = str(a) + "_" + str(b)
                tek2 = str(b) + "_" + str(a)
                if tek in e_dict:
                    te, edgeIndex = e_dict[tek]
                elif tek2 in e_dict:
                    te, edgeIndex = e_dict[tek2]
                else:
                    te, edgeIndex = Edge(a, b), len(el)
                    el.append(te)
                    e_dict[tek] = (te, edgeIndex)
                teil.append(edgeIndex)
                tel.append(te)
                """
                edgeIndex = hasEdgeAB(a, b, el)
                if edgeIndex == -1:
                    te = Edge(a, b)
                    teil.append(len(el))
                    el.append(te)
                else:
                    te = el[edgeIndex]
                    teil.append(edgeIndex)
                tel.append(te)
                """
        tpil = [tel[0].a, tel[0].b]
        tsl = set([0])
        tsp = tel[0].b
        # get ngon point ring
        while True:
            for ei in xrange(1, len(tel)):
                if ei in tsl:
                    continue
                te = tel[ei]
                if te.a == tsp:
                    tpil.append(te.b)
                    tsp = te.b
                    tsl.add(ei)
                    break
                if te.b == tsp:
                    tpil.append(te.a)
                    tsp = te.a
                    tsl.add(ei)
                    break
            else:
                break
        tmp.append(Ngon(tuple(tpil[:-1]), tuple(teil), tuple(tfil)))

    return el, tmp


def getVectorAvg(pil, pl):
    sum = reduce(lambda x, y: x + y, map(lambda i: pl[i], pil))
    return sum / float(len(pil))


def getCPolyNormal(f, pl):
    a, b, c = pl[f.a], pl[f.b], pl[f.c]
    return ((c - b).Cross(a - b)).GetNormalized()


def getCPolyEdge(f, i):
    if i == 0:
        return f.a, f.b
    elif i == 1:
        return f.b, f.c
    elif i == 2:
        return f.c, f.d
    else:
        return f.d, f.a


def getThridPoint(f, e):
    ei = f.FindEdge(e.a, e.b)
    if ei == 0:
        return f.c
    if ei == 1:
        return f.a
    if ei == 3:
        return f.b
    return -1


def checkLineNear(ai, bi, pl, tmp, r2):
    ci, di = -1, -1
    for e in tmp:
        if tzu.isInDis2(pl[ai], pl[e.a], r2):
            ci = e.a
            break
        if tzu.isInDis2(pl[ai], pl[e.b], r2):
            ci = e.b
            break
    for e in tmp:
        if tzu.isInDis2(pl[bi], pl[e.a], r2):
            di = e.a
            break
        if tzu.isInDis2(pl[bi], pl[e.b], r2):
            di = e.b
            break
    if ci < 0 or di < 0:
        ai = ai if ci < 0 else ci
        bi = bi if di < 0 else di
        tmp.append(Edge(ai, bi))
        return
    if hasEdgeAB(ci, di, tmp) == -1:
        tmp.append(Edge(ci, di))


def searchNearEdge(i, el, pl, th2):
    s = el[i]
    sa, sb = pl[s.a], pl[s.b]
    ai, bi = -1, -1
    for j, e in enumerate(el):
        if j == i:
            continue
        ta, tb = pl[e.a], pl[e.b]
        dv = tb - ta
        if ai == -1 and utils.PointLineDistance(ta, dv, sa).GetLengthSquared() < th2:
            ai = j
            continue
        if bi == -1 and utils.PointLineDistance(ta, dv, sb).GetLengthSquared() < th2:
            bi = j
    return ai, bi


def findLineLink(el, pl, th):
    tmp = []
    th2 = th**2
    for i, e in enumerate(el):
        ai, bi = searchNearEdge(i, el, pl, th2)
        tmp.append((ai, bi))
    return tmp


def findLineJoint(el, pl, th=5):
    r2 = th**2
    spil = map(lambda x: x.a, el)
    epil = map(lambda x: x.b, el)
    oel = []
    nel = []
    for e in el:
        a, b = e.a, e.b
        seg = [a]
        p0, v = pl[a], pl[b] - pl[a]
        l = v.GetLength()
        step = l / 100.0
        n = v.GetNormalized()
        addset = set([])
        for k in xrange(100):
            tlen = step * k
            if tlen < r2 or (l - tlen) < r2:
                continue
            cp = p0 + n * (tlen)
            for spi in spil:
                if a == spi or b == spi or (spi in addset):
                    continue
                if tzu.isInDis2(cp, pl[spi], r2):
                    seg.append(spi)
                    addset.add(spi)
            for epi in epil:
                if a == epi or b == epi or (epi in addset):
                    continue
                if tzu.isInDis2(cp, pl[epi], r2):
                    seg.append(epi)
                    addset.add(epi)
        seg.append(b)
        if len(seg) == 2:
            continue
        oel.append(e)
        for i in xrange(len(seg) - 1):
            sp, ep = seg[i], seg[i + 1]
            nel.append(Edge(sp, ep))
    el = filter(lambda x: x not in oel, el)
    el.extend(nel)
    return el


def deconstructSpline(obj, pl=None, th=2):
    r2 = th**2
    tmp = []
    if not pl:
        pl = obj.GetAllPoints()
    segc = obj.GetSegmentCount()
    if th > 0:
        if segc >= 1:
            sindex = 0
            for i in xrange(segc):
                tsegd = obj.GetSegment(i)
                tcnt = tsegd["cnt"]
                if tcnt < 2:
                    continue
                for j in xrange(sindex, sindex + tcnt - 1):
                    ai, bi = j, j + 1
                    checkLineNear(ai, bi, pl, tmp, r2)
                if tsegd["closed"]:
                    ai, bi = sindex + tcnt - 1, sindex
                    checkLineNear(ai, bi, pl, tmp, r2)
                sindex += tcnt
        else:
            pc = len(pl)
            if pc < 2:
                return tmp
            for i in xrange(pc - 2):
                ai, bi = i, i + 1
                checkLineNear(ai, bi, pl, tmp, r2)
            if obj.IsClosed():
                ai, bi = pc - 1, 0
                checkLineNear(ai, bi, pl, tmp, r2)
    else:
        if segc >= 1:
            sindex = 0
            for i in xrange(segc):
                tsegd = obj.GetSegment(i)
                tcnt = tsegd["cnt"]
                if tcnt < 2:
                    continue
                for j in xrange(sindex, sindex + tcnt - 1):
                    ai, bi = j, j + 1
                    tmp.append((ai, bi))
                if tsegd["closed"]:
                    ai, bi = sindex + tcnt - 1, sindex
                    tmp.append((ai, bi))

                sindex += tcnt
        else:
            pc = len(pl)
            if pc < 2:
                return tmp
            for i in xrange(pc - 2):
                ai, bi = i, i + 1
                tmp.append((ai, bi))
            if obj.IsClosed():
                ai, bi = pc - 1, 0
                tmp.append((ai, bi))
    return tmp


#=====UV======
UVFace = namedtuple("UVFace", ["a", "b", "c", "d", "isTri"])


def isSameSide(A, B, C, P):
    AB = B - A
    AC = C - A
    AP = P - A
    v1 = AB.Cross(AC)
    v2 = AB.Cross(AP)
    # v1 and v2 should point to the same direction
    return v1.Dot(v2) >= 0


def isInTri(A, B, C, P):
    return isSameSide(A, B, C, P) and isSameSide(B, C, A, P) and isSameSide(C, A, B, P)


def cross2d(a, b, p):
    return (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x)


def toLeft(a, b, p):
    return cross2d(a, b, p) > 0


def offsetTri(a, b, c, th):
    m = (a + b + c) / 3.0
    na = a + (a - m).GetNormalized() * th
    nb = b + (b - m).GetNormalized() * th
    nc = c + (c - m).GetNormalized() * th
    return na, nb, nc


def isInTri2d(a, b, c, p, th=0):
    if th > 0:
        a, b, c = offsetTri(a, b, c, th)
    res = toLeft(a, b, p)
    if res != toLeft(b, c, p):
        return False
    if res != toLeft(c, a, p):
        return False
    if cross2d(a, b, c) == 0:  # ABC is in one line
        return False
    return True


def getTriUV(A, B, C, P):
    # v_AB,u_AC
    v0 = C - A
    v1 = B - A
    v2 = P - A
    dot00 = v0.Dot(v0)
    dot01 = v0.Dot(v1)
    dot02 = v0.Dot(v2)
    dot11 = v1.Dot(v1)
    dot12 = v1.Dot(v2)
    inverDeno = 1 / float(dot00 * dot11 - dot01 * dot01)
    u = (dot11 * dot02 - dot01 * dot12) * inverDeno
    v = (dot00 * dot12 - dot01 * dot02) * inverDeno
    return u, v


def getPosByUV(a, b, c, u, v, d=0, fn=None, an=None, bn=None, cn=None):
    if d == 0:
        return a + (b - a) * v + (c - a) * u

    if not an:
        r = a + (b - a) * v + (c - a) * u
        r = r + fn * d
    else:
        a = a + an * d
        b = b + bn * d
        c = c + cn * d
        r = a + (b - a) * v + (c - a) * u
    return r


def getAllUVFace(fl, uvtag, size, off):
    uvfl = []
    for i in xrange(uvtag.GetDataCount()):
        d = uvtag.GetSlow(i)
        p1 = d["a"] * size + off
        p1.y = p1.y * -1
        p2 = d["b"] * size + off
        p2.y = p2.y * -1
        p3 = d["c"] * size + off
        p3.y = p3.y * -1
        p4 = d["d"] * size + off
        p4.y = p4.y * -1
        uvfl.append(UVFace(a=p1, b=p2, c=p3, d=p4, isTri=fl[i].IsTriangle()))
    return uvfl


def getAllUVPoints(uvfl, fl):
    tmp = {}
    for i, uvf in enumerate(uvfl):
        f = fl[i]
        tmp.setdefault(f.a, {})[i] = uvf.a
        tmp.setdefault(f.b, {})[i] = uvf.b
        tmp.setdefault(f.c, {})[i] = uvf.c
        if not f.IsTriangle():
            tmp.setdefault(f.d, {})[i] = uvf.d

    return tmp


def findUVFace(p, uvfl, th):
    for i, f in enumerate(uvfl):
        if isInTri2d(f.a, f.b, f.c, p, th):
            return i
        if f.isTri == False and isInTri2d(f.a, f.d, f.c, p, th):
            return i
    return -1


#========spline connect========
class SplineSeg(object):

    def __init__(self, s):
        self.obj = s
        self.pc = s.GetPointCount()
        self.sp = s.GetPoint(0)
        self.ep = s.GetPoint(self.pc - 1)
        self.alive = True
        self.build = True

    def __repr__(self):
        return "SplineSeg{" + str(self.sp) + "_" + str(self.ep) + "}"


def joinSplineSeg(s, t, inva, invb, doc):
    s.alive = False
    t.alive = False
    s.build = False
    t.build = False
    p2i = s.pc
    a = tzu.getCurr(s.obj, doc)
    b = tzu.getCurr(t.obj, doc)
    if inva:
        tsp = s.ep
        tzu.invertSpline(a, doc)
    else:
        tsp = s.sp
    if invb:
        tep = t.sp
        tzu.invertSpline(b, doc)
    else:
        tep = t.ep
    # join
    b.InsertUnder(a)
    a = tzu.getJoin([a], doc)
    tzu.joinSplineSegment(a, doc)
    #del p2
    sel = a.GetPointS()
    sel.DeselectAll()
    sel.Select(p2i)
    tzu.delPoint(a, doc)

    tmp = SplineSeg(a)
    tmp.sp = tsp
    tmp.ep = tep
    return tmp


def connectSpline(splineList, th, doc):
    th2 = th**2
    sl = []
    for s in splineList:
        if s.GetSegmentCount() > 1:
            sl.extend(tzu.explodeObject(s))
        else:
            sl.append(s)
    sl = filter(lambda x: x.GetPointCount() > 1, sl)
    sl = map(lambda x: SplineSeg(x), sl)

    while len([s for s in sl if s.alive]) > 1:
        for i, s in enumerate(sl):
            if not s.alive:
                continue
            found = False
            for j, t in enumerate(sl):
                if i == j or t.alive == False:
                    continue
                if (s.sp - t.sp).GetLengthSquared() < th2:
                    sl.append(joinSplineSeg(s, t, True, False, doc))
                    found = True
                    break
                if (s.sp - t.ep).GetLengthSquared() < th2:
                    sl.append(joinSplineSeg(s, t, True, True, doc))
                    found = True
                    break
                if (s.ep - t.sp).GetLengthSquared() < th2:
                    sl.append(joinSplineSeg(s, t, False, False, doc))
                    found = True
                    break
                if (s.ep - t.ep).GetLengthSquared() < th2:
                    sl.append(joinSplineSeg(s, t, False, True, doc))
                    found = True
                    break
            if not found:
                s.alive = False
            break
    sl = filter(lambda x: x.build, sl)
    return map(lambda x: x.obj, sl)


def isSameRay(a, b, c):
    ba, bc = a - b, c - b
    bal, bcl = ba.GetLength(), bc.GetLength()
    # if bal < 0.1 or bcl < 0.1:
    #     return True
    return abs(ba.Dot(bc) / (bal * bcl)) > 0.9998


def isEdgeOverlap(a, b, c, d):
    ab = b - a
    cd = d - c
    ca = a - c
    cb = b - c
    da = a - d
    db = b - d
    l1, l2 = ca.GetLength(), cb.GetLength()
    l3, l4 = da.GetLength(), db.GetLength()
    if l1 < 0.1 or l4 < 0.1:
        cos = ab.Dot(cd) / (ab.GetLength() * cd.GetLength())
        if cos < -0.9998:
            return False
        elif abs(cos) > 0.9998:
            return True
        else:
            return False
    if l2 < 0.1 or l3 < 0.1:
        cos = ab.Dot(cd) / (ab.GetLength() * cd.GetLength())
        if cos > 0.9998:
            return False
        elif abs(cos) > 0.9998:
            return True
        else:
            return False
    c1 = ca.Dot(cb) / (l1 * l2)
    c2 = da.Dot(db) / (l3 * l4)
    return min(c1, c2) < 0.9998 and abs(c1) > 0.9998 and abs(c2) > 0.9998


def removeEdgeRayOverlap(el, pl, mode="first"):
    # short,long,both,first
    cb = combinations(list(range(len(el))), 2)
    removeList = set([])
    for i, j in cb:
        if i in removeList or j in removeList:
            continue
        s, t = el[i], el[j]
        a, b = pl[s.a], pl[s.b]
        c, d = pl[t.a], pl[t.b]
        if not isEdgeOverlap(a, b, c, d):
            continue

        if mode == "both":
            removeList.add(i)
            removeList.add(j)
        elif mode == "first":
            removeList.add(i)
        elif mode == "short":
            if (b - a).GetLengthSquared() < (d - c).GetLengthSquared():
                removeList.add(i)
            else:
                removeList.add(j)
        elif mode == "long":
            if (b - a).GetLengthSquared() > (d - c).GetLengthSquared():
                removeList.add(i)
            else:
                removeList.add(j)
    tmp = []
    for i, e in enumerate(el):
        if i in removeList:
            continue
        tmp.append(e)
    return tmp
