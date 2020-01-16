# coding=utf-8
import c4d
from c4d import utils
from c4d.utils.noise import Noise

"""
190721:
加入splinemask函数两个
1：用生成器塌陷
2：用utils里的函数
两个可靠性都一般

190921：
加入基于R21的splineOffset(测试)
加入getTagByName
"""


class DIR(object):
    up = c4d.Vector(0, 1, 0)
    down = c4d.Vector(0, -1, 0)
    left = c4d.Vector(-1, 0, 0)
    right = c4d.Vector(1, 0, 0)
    inn = c4d.Vector(0, 0, 1)
    out = c4d.Vector(0, 0, -1)
    zero = c4d.Vector(0)
    fx4 = ((-1, 0), (1, 0), (0, 1), (0, -1))
    fx4c = ((1, 1), (1, -1), (-1, 1), (-1, -1))


class RandomValue(object):

    def __init__(self, _d):
        self.d = _d
        if isinstance(self.d, list):
            self.t = "choice"
            self.get = self.gf_choice
        elif isinstance(self.d, tuple):
            self.t = "uniform"
            self.get = self.gf_uniform
        elif isinstance(self.d, (float, int)):
            self.t = "fixed"
            self.get = self.gf_fixed

    def gf_choice(self, rd):
        return rd.choice(self.d)

    def gf_uniform(self, rd):
        return rd.uniform(self.d[0], self.d[1])

    def gf_fixed(self, rd=None):
        return self.d

    def __repr__(self):
        return self.t + ":" + str(self.d)


def vector2matrix(zv):
    return utils.HPBToMatrix(utils.VectorToHPB(zv))


def clamp(v, minv, maxv):
    return max(min(v, maxv), minv)


def avg(l):
    sum = reduce(lambda x, y: x + y, l)
    return sum / float(len(l))


def getBlendList(cnt):
    tmp = []
    seg = 1.0 / float(cnt)
    for i in xrange(cnt):
        tmp.append(seg * i + 0.01)
    return tuple(tmp)


def safeDot(a, b):
    v = a.Dot(b)
    return v if v != 0 else 1


def checkType(o, tl):
    if not o:
        return False
    t = o.GetType()
    if isinstance(tl, tuple):
        return True if (t in tl) else False
    return True if t == tl else False


def isTP(obj):
    if not obj:
        return False
    i = obj.GetType()
    return i == 1001381


def isMG(obj):
    if not obj:
        return False
    i = obj.GetType()
    return i in (1018975, 1018791, 1018544, 1018545, 1019268, 1036557, 440000054)


def isPoly(obj):
    if not obj:
        return False
    i = obj.GetType()
    return i == 5100 or i == 1024542


def isSpline(obj):
    if not obj:
        return False
    i = obj.GetType()
    return i == 5101


def isInDis2(a, b, r2=16):
    return (a - b).GetLengthSquared() < r2


def centerAxis(obj):
    c = obj.GetMp()
    m = obj.GetMg()
    pa = obj.GetAllPoints()
    cg = m.Mul(c)
    pa = map(lambda x: x - cg, pa)
    obj.SetAllPoints(pa)
    m.off = m.off + cg
    obj.SetMg(m)
    obj.Message(c4d.MSG_UPDATE)


def setObjectColor(obj, c):
    obj[c4d.ID_BASEOBJECT_USECOLOR] = 2
    obj[c4d.ID_BASEOBJECT_COLOR] = c


def setVisibility(obj, v=2, r=2):
    obj[c4d.ID_BASEOBJECT_VISIBILITY_EDITOR] = v
    obj[c4d.ID_BASEOBJECT_VISIBILITY_RENDER] = r


def makeSelectionObject(tmp):
    s = c4d.BaseObject(5190)
    l = s[c4d.SELECTIONOBJECT_LIST]
    for o in tmp:
        l.InsertObject(o, 15)
    s[c4d.SELECTIONOBJECT_LIST] = l
    return s


def getTagByName(obj, name):
    tl = obj.GetTags()
    for t in tl:
        if t.GetName() == name:
            return t
    return None


def clearTags(obj, tn):
    if isinstance(tn, int):
        tn = (tn,)
    tags = obj.GetTags()
    for t in reversed(tags):
        if t.GetType() in tn:
            t.Remove()
    c4d.EventAdd()


def mergeSelTag(obj, typ, name):
    tmp = obj.GetTags()
    tags = []
    for t in tmp:
        if not t.CheckType(typ):
            continue
        tn = t.GetName()
        if tn.find(name) != -1:
            tags.append(t)
    if len(tags) <= 1:
        return
    first = tags[0]
    tsel = first.GetBaseSelect()
    for i in reversed(xrange(1, len(tags))):
        tt = tags[i]
        tsel.Merge(tt.GetBaseSelect())
        tt.Remove()
    first.SetName(name)
    c4d.EventAdd()


def getInExclusion(ie, doc):
    tmp = []
    for i in xrange(ie.GetObjectCount()):
        to = ie.ObjectFromIndex(doc, i)
        if to:
            tmp.append(to)
    return tmp


def showPoint(p, doc, col=None):
    if not col:
        col = c4d.Vector(1, 0, 0)
    n = c4d.BaseObject(5140)
    setObjectColor(n, col)
    n[c4d.NULLOBJECT_DISPLAY] = 0
    n[c4d.ID_BASEOBJECT_REL_POSITION] = p
    doc.InsertObject(n)
    c4d.EventAdd()
    return n


def showPoints(pl, doc, col=None):
    if not col:
        col = c4d.Vector(1, 0, 0)
    pn = c4d.BaseObject(5140)
    pn.SetName("Points:" + str(id(pl)))
    for p in pl:
        n = c4d.BaseObject(5140)
        setObjectColor(n, col)
        n[c4d.NULLOBJECT_DISPLAY] = 0
        n[c4d.ID_BASEOBJECT_REL_POSITION] = p
        n.InsertUnderLast(pn)
    doc.InsertObject(pn)
    c4d.EventAdd()
    return pn


def getBD(obj):
    r = obj.GetRad()
    c = obj.GetMp()
    m = obj.GetMg()
    return {"min": m.Mul((c - r)), "max": m.Mul(c + r), "r": r, "c": m.Mul(c), "m": m}
#==================================================


def CEdgeSelection(nb, mesh, faceCount, sel):
    tsel = sel.GetAll(faceCount * 4)
    tsel = [i for i in xrange(len(tsel)) if tsel[i] == 1]
    bs = c4d.BaseSelect()
    for i in tsel:
        ti = nb.GetPolyInfo(int(i / 4))["edge"][int(i % 4)]
        bs.Select(ti)
    return bs


def edge2Spline(obj, doc=None):
    utils.SendModelingCommand(
        command=c4d.MCOMMAND_EDGE_TO_SPLINE,
        list=[obj],
        mode=c4d.MODELINGCOMMANDMODE_EDGESELECTION,
        doc=doc)
    c = obj.GetDown()
    if not c or c.GetType() != 5101:
        return None
    return c


def getSplinePlane(s):
    # 返回 i012,v法线，f起始距离
    sp = s.GetPoint(0)
    b = s.GetRad()
    b = (abs(b.x), abs(b.y), abs(b.z))
    if b[2] < 0.1 and b[0] > 1 and b[1] > 1:
        return 0, DIR.inn, sp.z  # XY
    if b[0] < 0.1 and b[1] > 1 and b[2] > 1:
        return 1, DIR.right, sp.x  # YZ
    if b[1] < 0.1 and b[0] > 1 and b[2] > 1:
        return 2, DIR.up, sp.y  # XZ
    return -1, None, 0


def realSplineOffset(s, r, doc):
    plane, move, sp = getSplinePlane(s)
    if plane == -1:
        raise Exception("not flat plane")
    ta = (2, 0, 1)[plane]
    cs = CSTO(s, doc)
    if cs is None:
        raise Exception("error CSTO spline")
    e = c4d.BaseObject(5116)
    e[c4d.EXTRUDEOBJECT_MOVE] = move * 1000
    e[c4d.CAPSANDBEVELS_CAP_ENABLE_END] = False
    e[c4d.CAPSANDBEVELS_STARTBEVEL_TYPE] = 3
    e[c4d.CAPSANDBEVELS_STARTBEVEL_OFFSET] = -r
    e[c4d.CAPSANDBEVELS_STARTBEVEL_SEGMENTS] = 1
    e[c4d.CAPSANDBEVELS_AVOIDINTERSECTION] = True
    e[c4d.CAPSANDBEVELS_CAP_TYPE] = 2
    cs.InsertUnderLast(e)
    # return e
    # csto extrude
    c4d.EventAdd()
    doc.InsertObject(e)
    eo = CSTO(e, doc)
    e.Remove()
    if eo is None:
        raise Exception("error CSTO extrude")
    fc = eo.GetPolygonCount()
    if fc == 0:
        raise Exception("no polygon found")
    # sel edge
    st = getTagByName(eo, "ER1")
    if not st:
        raise Exception("no ER1 selectTag found")
    tsel = st.GetBaseSelect()
    if tsel.GetCount() < 3:
        return None
    tsel = tsel.GetAll(fc * 4)
    bs = c4d.BaseSelect()
    nb = utils.Neighbor()
    nb.Init(eo)
    for i, stat in enumerate(tsel):
        if stat == 0:
            continue
        tfi = int(i / 4)
        tei = int(i % 4)
        teil = nb.GetPolyInfo(tfi)["edge"]
        bs.Select(teil[tei])
    eo.SetSelectedEdges(nb, bs, c4d.EDGESELECTIONTYPE_SELECTION)
    # edge2spline
    rs = edge2Spline(eo, doc)
    if not rs:
        raise Exception("edge2spline failed")
    # rs = CSTO(rs, doc)这一步外部进行
    rs[c4d.SPLINEOBJECT_CLOSED] = True
    return rs


def splineMask(s1, s2, doc, mode=1, plane=0):
    # mode:0:并,1:1-2,2:2-1,3:交;4:交反,5:割
    # plane:0:xy,1:zy,2:xz
    sm = c4d.BaseObject(1019396)
    sm[c4d.MGSPLINEMASKOBJECT_MODE] = int(mode)
    sm[c4d.MGSPLINEMASKOBJECT_AXIS] = int(plane)
    sm[c4d.MGSPLINEMASKOBJECT_CREATECAP] = False
    t1 = s1.GetClone()
    t2 = s2.GetClone()
    t1.InsertUnderLast(sm)
    t2.InsertUnderLast(sm)
    # doc.InsertObject(sm)
    rs = tzu.getSpline(sm, doc)
    # sm.Remove()
    return rs


def splineMask2(s1, s2list, doc, mode=1, plane=0):
    if not isinstance(s2list, list):
        s2list = [s2list]
    BOOL_PLANE = [c4d.SPLINEBOOL_AXIS_XY, c4d.SPLINEBOOL_AXIS_ZY, c4d.SPLINEBOOL_AXIS_XZ]
    BOOL_MODE = [c4d.SPLINEBOOL_MODE_UNION, c4d.SPLINEBOOL_MODE_AMINUSB,
                 c4d.SPLINEBOOL_MODE_BMINUSA, c4d.SPLINEBOOL_MODE_AND,
                 c4d.SPLINEBOOL_MODE_OR, c4d.SPLINEBOOL_MODE_INTERSECTION]
    return utils.BooleanSplines(s1, s2list,
                                doc=doc,
                                bd=doc.GetActiveBaseDraw(),
                                projectionAxis=BOOL_PLANE[int(plane)],
                                booleanMode=BOOL_MODE[int(mode)])


def centerTool(obj, pos, doc=None):
    if not doc:
        doc = obj.GetDocument()
    settings = c4d.BaseContainer()
    settings[c4d.MDATA_CENTER_XAXIS] = pos.x
    settings[c4d.MDATA_CENTER_YAXIS] = pos.y
    settings[c4d.MDATA_CENTER_ZAXIS] = pos.z
    utils.SendModelingCommand(
        command=c4d.ID_MODELING_CENTER_TOOL,
        list=[obj],
        mode=c4d.MODELINGCOMMANDMODE_ALL,
        bc=settings,
        doc=doc)


def invertSpline(obj, doc=None):
    if not doc:
        doc = obj.GetDocument()
    utils.SendModelingCommand(
        command=c4d.MCOMMAND_SPLINE_REVERSE,
        list=[obj],
        mode=c4d.MODELINGCOMMANDMODE_ALL,
        doc=doc)


def splineChamfer(obj, r, flat, doc=None):
    if not doc:
        doc = obj.GetDocument()
    settings = c4d.BaseContainer()
    settings[c4d.MDATA_SPLINE_CHAMFERRADIUS] = r
    settings[c4d.MDATA_SPLINE_CHAMFERFLAT] = flat
    utils.SendModelingCommand(
        command=c4d.ID_MODELING_SPLINE_CHAMFER_TOOL,
        list=[obj],
        mode=c4d.MODELINGCOMMANDMODE_ALL,
        bc=settings,
        doc=doc)


def joinSplineSegment(obj, doc=None):
    if not doc:
        doc = obj.GetDocument()
    utils.SendModelingCommand(
        command=c4d.MCOMMAND_SPLINE_JOINSEGMENT,
        list=[obj],
        mode=c4d.MODELINGCOMMANDMODE_ALL,
        doc=doc)


def splineOffset(obj, offset, createNew=True, doc=None):
    if not doc:
        doc = obj.GetDocument()
    settings = c4d.BaseContainer()
    settings[c4d.MDATA_SPLINE_OUTLINE] = offset
    settings[c4d.MDATA_SPLINE_OUTLINESEPARATE] = createNew
    r = utils.SendModelingCommand(
        command=c4d.MCOMMAND_SPLINE_CREATEOUTLINE,
        list=[obj],
        mode=c4d.MODELINGCOMMANDMODE_ALL,
        bc=settings,
        doc=doc)
    if createNew:
        return r[0]


def splitPoly(obj, doc=None):
    if not doc:
        doc = obj.GetDocument()
    return utils.SendModelingCommand(
        command=c4d.MCOMMAND_SPLIT,
        list=[obj],
        mode=c4d.MODELINGCOMMANDMODE_POLYGONSELECTION,
        doc=doc
    )[0]


def disconnectPoly(obj, keepGroup=False, doc=None):
    if not doc:
        doc = obj.GetDocument()
    settings = c4d.BaseContainer()
    settings[c4d.MDATA_DISCONNECT_PRESERVEGROUPS] = keepGroup
    utils.SendModelingCommand(
        command=c4d.MCOMMAND_DISCONNECT,
        list=[obj],
        mode=c4d.MODELINGCOMMANDMODE_ALL,
        bc=settings,
        doc=doc)


def normalMove(obj, d, doc=None):
    if not doc:
        doc = obj.GetDocument()
    settings = c4d.BaseContainer()
    settings[c4d.MDATA_NORMALMOVE_VALUE] = d
    utils.SendModelingCommand(
        command=c4d.ID_MODELING_NORMALMOVE_TOOL,
        list=[obj],
        mode=c4d.MODELINGCOMMANDMODE_ALL,
        bc=settings,
        doc=doc
    )


def extrudePoly(obj, d, cap=True, doc=None):
    if not doc:
        doc = obj.GetDocument()
    settings = c4d.BaseContainer()
    settings[c4d.MDATA_EXTRUDE_OFFSET] = d
    settings[c4d.MDATA_EXTRUDE_PRESERVEGROUPS] = True
    settings[c4d.MDATA_EXTRUDE_CREATECAPS] = cap
    utils.SendModelingCommand(
        command=c4d.ID_MODELING_EXTRUDE_TOOL,
        list=[obj],
        mode=c4d.MODELINGCOMMANDMODE_POLYGONSELECTION,
        bc=settings,
        doc=doc
    )


def extrudeInner(obj, d, keepGroup=True, doc=None):
    if not doc:
        doc = obj.GetDocument()
    settings = c4d.BaseContainer()
    settings[c4d.MDATA_EXTRUDEINNER_OFFSET] = d
    settings[c4d.MDATA_EXTRUDEINNER_PRESERVEGROUPS] = keepGroup
    utils.SendModelingCommand(
        command=c4d.ID_MODELING_EXTRUDE_INNER_TOOL,
        list=[obj],
        mode=c4d.MODELINGCOMMANDMODE_POLYGONSELECTION,
        bc=settings,
        doc=doc
    )


def explodeObject(obj):
    # 这个不能传doc,否则结果必定为not alive,无解
    doc = obj.GetDocument()
    utils.SendModelingCommand(
        command=c4d.MCOMMAND_EXPLODESEGMENTS,
        list=[obj],
        mode=c4d.MODELINGCOMMANDMODE_ALL,
        doc=doc)
    ch = obj.GetChildren()
    return ch


def delPoly(obj, doc=None):
    if not doc:
        doc = obj.GetDocument()
    utils.SendModelingCommand(
        command=c4d.MCOMMAND_DELETE,
        list=[obj],
        mode=c4d.MODELINGCOMMANDMODE_POLYGONSELECTION,
        doc=doc)


def delPoint(obj, doc=None):
    if not doc:
        doc = obj.GetDocument()
    utils.SendModelingCommand(
        command=c4d.MCOMMAND_DELETE,
        list=[obj],
        mode=c4d.MODELINGCOMMANDMODE_POINTSELECTION,
        doc=doc)


def dissolvePoly(obj, doc=None):
    if not doc:
        doc = obj.GetDocument()
    utils.SendModelingCommand(
        command=c4d.MCOMMAND_MELT,
        list=[obj],
        mode=c4d.MODELINGCOMMANDMODE_POLYGONSELECTION,
        doc=doc)


def dissolvePoint(obj, doc=None):
    if not doc:
        doc = obj.GetDocument()
    utils.SendModelingCommand(
        command=c4d.MCOMMAND_MELT,
        list=[obj],
        mode=c4d.MODELINGCOMMANDMODE_POINTSELECTION,
        doc=doc)


def alignNormal(obj, doc=None):
    if not doc:
        doc = obj.GetDocument()
    utils.SendModelingCommand(
        command=c4d.MCOMMAND_ALIGNNORMALS,
        list=[obj],
        mode=c4d.MODELINGCOMMANDMODE_ALL,
        doc=doc)


def reverseNormal(obj, doc=None):
    if not doc:
        doc = obj.GetDocument()
    utils.SendModelingCommand(
        command=c4d.MCOMMAND_REVERSENORMALS,
        list=[obj],
        mode=c4d.MODELINGCOMMANDMODE_ALL,
        doc=doc)


def triangulate(obj, doc=None):
    if not doc:
        doc = obj.GetDocument()
    settings = c4d.BaseContainer()
    utils.SendModelingCommand(
        command=c4d.MCOMMAND_TRIANGULATE,
        list=[obj],
        mode=c4d.MODELINGCOMMANDMODE_ALL,
        doc=doc)


def unTriangulate(obj, r, ngon=False, doc=None):
    if not doc:
        doc = obj.GetDocument()
    settings = c4d.BaseContainer()
    settings[c4d.MDATA_UNTRIANGULATE_NGONS] = ngon
    settings[c4d.MDATA_UNTRIANGULATE_ANGLE_RAD] = r
    utils.SendModelingCommand(
        command=c4d.MCOMMAND_UNTRIANGULATE,
        list=[obj],
        mode=c4d.MODELINGCOMMANDMODE_ALL,
        bc=settings,
        doc=doc)


def optimize(obj, th, doc=None):
    if not doc:
        doc = obj.GetDocument()
    settings = c4d.BaseContainer()
    settings[c4d.MDATA_OPTIMIZE_TOLERANCE] = th
    settings[c4d.MDATA_OPTIMIZE_POINTS] = True
    settings[c4d.MDATA_OPTIMIZE_POLYGONS] = True
    settings[c4d.MDATA_OPTIMIZE_UNUSEDPOINTS] = True
    utils.SendModelingCommand(
        command=c4d.MCOMMAND_OPTIMIZE,
        list=[obj],
        mode=c4d.MODELINGCOMMANDMODE_ALL,
        bc=settings,
        doc=doc)


def getCurr(obj, doc=None):
    if not doc:
        doc = obj.GetDocument()
    obj = utils.SendModelingCommand(
        command=c4d.MCOMMAND_CURRENTSTATETOOBJECT,
        list=[obj],
        mode=c4d.MODELINGCOMMANDMODE_ALL,
        doc=doc
    )
    return None if not obj else obj[0]


def getJoin(l, doc=None):
    if not doc:
        doc = l[0].GetDocument()
    return utils.SendModelingCommand(
        command=c4d.MCOMMAND_JOIN,
        list=l,
        mode=c4d.MODELINGCOMMANDMODE_ALL,
        doc=doc
    )[0]


def CSTO(obj, doc):
    obj = utils.SendModelingCommand(
        command=c4d.MCOMMAND_CURRENTSTATETOOBJECT,
        list=[obj],
        mode=c4d.MODELINGCOMMANDMODE_ALL,
        doc=doc
    )
    r = utils.SendModelingCommand(
        command=c4d.MCOMMAND_JOIN,
        list=obj,
        mode=c4d.MODELINGCOMMANDMODE_ALL,
        doc=doc
    )
    if r:
        return r[0]
    return None


def getState(obj, doc=None, keepPSR=False):
    if not doc:
        doc = obj.GetDocument()
    m = obj.GetMg()
    cop = obj.GetClone()
    c = c4d.BaseObject(5140)

    if not keepPSR:
        cop.InsertUnder(c)
        cop.SetMg(m)
    else:
        cop.InsertUnder(c)
        c.SetMg(m)
        cop.SetMg(c4d.Matrix())
    obj = c

    obj = utils.SendModelingCommand(
        command=c4d.MCOMMAND_CURRENTSTATETOOBJECT,
        list=[obj],
        mode=c4d.MODELINGCOMMANDMODE_ALL,
        doc=doc
    )
    r = utils.SendModelingCommand(
        command=c4d.MCOMMAND_JOIN,
        list=obj,
        mode=c4d.MODELINGCOMMANDMODE_ALL,
        doc=doc
    )
    if r:
        return r[0]
    return None


def getAllState(l, tl, doc, keepPSR=False):
    tmp = []
    for o in l:
        s = getState(o, doc, keepPSR)
        if checkType(s, tl):
            tmp.append(s)
    return tmp


def getStateWithVTD(obj, tl, ti, doc, keepPSR=False):
    if not checkType(obj, tl):
        obj = getState(obj, doc, keepPSR)
    if not checkType(obj, tl):
        return None, None
    vtd = None
    vt = obj.GetTag(ti)
    if vt != None:
        vtd = vt.GetAllHighlevelData()
    return obj, vtd


def getAllStateWithVTD(l, tl, ti, doc, keepPSR=False):
    tmp = []
    tmp2 = []
    for o in l:
        to, td = getStateWithVTD(o, tl, ti, doc, keepPSR)
        tmp.append(to)
        tmp2.append(td)
    return tmp, tmp2
#==================================================


def getPoly(obj, doc=None, keepPSR=False):
    if not obj:
        return None
    if isPoly(obj):
        return obj
    if not doc:
        doc = obj.GetDocument()
    obj = getState(obj, doc, keepPSR)

    return obj if isPoly(obj) else None


def getSpline(obj, doc=None, keepPSR=False):
    if not obj:
        return None
    if isSpline(obj):
        return obj
    if not doc:
        doc = obj.GetDocument()
    obj = getState(obj, doc, keepPSR)

    return obj if isSpline(obj) else None


def clearPoints(obj):
    obj.ResizeObject(0)
    obj.Message(c4d.MSG_UPDATE)


def makeSpline(pl, close=False):
    pc = len(pl)
    s = c4d.SplineObject(pc, c4d.SPLINETYPE_LINEAR)
    s[c4d.SPLINEOBJECT_CLOSED] = close
    s[c4d.SPLINEOBJECT_TYPE] = 0
    s[c4d.SPLINEOBJECT_INTERPOLATION] = 0

    s.ResizeObject(pc, 1)
    s.SetSegment(0, pc, close)
    s.SetAllPoints(list(pl))
    s.Message(c4d.MSG_UPDATE)
    return s


def copySpline(f, t):
    pl = f.GetAllPoints()
    pc = len(pl)
    segc = f.GetSegmentCount()
    t.ResizeObject(pc, segc)
    if segc <= 1:
        t.SetSegment(0, pc, f.IsClosed())
    else:
        for i in xrange(segc):
            d = f.GetSegment(i)
            t.SetSegment(i, d["cnt"], d["closed"])
    t.SetAllPoints(list(pl))
    t.Message(c4d.MSG_UPDATE)


def makePoly(pl, fl=None):
    pc = len(pl)
    if not fl:
        if pc == 3:
            fl = ((0, 1, 2),)
        elif pc == 4:
            fl = ((0, 1, 2, 3),)
    if not fl:
        return None

    g = c4d.PolygonObject(pc, len(fl))
    g.SetAllPoints(pl)

    for i, f in enumerate(fl):
        if len(f) == 3:
            c = c4d.CPolygon(f[0], f[1], f[2])
            g.SetPolygon(i, c)
        elif len(f) == 4:
            c = c4d.CPolygon(f[0], f[1], f[2], f[3])
            g.SetPolygon(i, c)

    g.Message(c4d.MSG_UPDATE)
    return g


def copyPoly(f, t):
    pl = f.GetAllPoints()
    fl = f.GetAllPolygons()
    pc = len(pl)
    t.ResizeObject(pc, len(fl))
    t.SetAllPoints(pl)
    for i, c in enumerate(fl):
        t.SetPolygon(i, c)
    t.Message(c4d.MSG_UPDATE)

#==================================================


class CurlNoise(object):
    OFFSET = (c4d.Vector(2333), c4d.Vector(6666), c4d.Vector(-2333))

    def __init__(self, ss=5):
        self.sampleSize = ss
        self.offset = CurlNoise.OFFSET

    def __n(self, x, y, z, axis):
        pos = c4d.Vector(x, y, z)
        return Noise((pos + self.offset[axis]) * 0.001 * self.sampleSize)

    def sample(self, pos):
        d = 0.01
        v = c4d.Vector(0)
        x = pos.x
        y = pos.y
        z = pos.z
        v.x = self.__n(x, y + d, z, 2) - self.__n(x, y - d, z, 2) - \
            self.__n(x, y, z + d, 1) + self.__n(x, y, z - d, 1)
        v.y = self.__n(x, y, z + d, 0) - self.__n(x, y, z - d, 0) - \
            self.__n(x + d, y, z, 2) + self.__n(x - d, y, z, 2)
        v.z = self.__n(x + d, y, z, 1) - self.__n(x - d, y, z, 1) - \
            self.__n(x, y + d, z, 0) + self.__n(x, y - d, z, 0)
        return v / (d * 2) * 1000
