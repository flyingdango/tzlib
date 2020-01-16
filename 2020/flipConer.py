# -*- coding: UTF-8 -*-
"""
一个独特建模工具;用于反转角
可以做一些切割折叠效果
180917
初步完成,开始实践验证

by飞舞的团子
#"""

import c4d
import math
from c4d import utils


def debug(g):
    pn = c4d.BaseObject(5140)
    for i, p in enumerate(g.pl):
        n = c4d.BaseObject(5140)
        n.SetName(str(g.pil[i]))
        n[c4d.ID_BASEOBJECT_REL_POSITION] = p
        n[c4d.ID_BASEOBJECT_USECOLOR] = 1
        n[c4d.ID_BASEOBJECT_COLOR] = c4d.Vector(0, 1, 1)
        n[c4d.NULLOBJECT_DISPLAY] = 13
        n[c4d.NULLOBJECT_RADIUS] = 2
        n[c4d.NULLOBJECT_ORIENTATION] = 1
        n.InsertUnder(pn)
    doc.InsertObject(pn)
    c4d.EventAdd()


class FaceGroup(object):

    def __init__(self, n_=None):
        self.n = n_
        self.fl = []
        self.pil = set()
        self.pl = None
        self.pc = 0
        self.dl = None
        self.rm = None
        self.off = None
        self.anchor = None
        self.flag = None
        self.targetLength = 0

    def __repr__(self):
        return "G:" + str(len(self.fl))


def isPoly(obj):
    if not obj:
        return False
    return obj.CheckType(5100)


def disConnect(obj):
    settings = c4d.BaseContainer()
    settings[c4d.MDATA_DISCONNECT_PRESERVEGROUPS] = True
    utils.SendModelingCommand(
        command=c4d.MCOMMAND_DISCONNECT,
        list=[obj],
        mode=c4d.MODELINGCOMMANDMODE_POLYGONSELECTION,
        bc=settings,
        doc=doc)
    c4d.EventAdd()


def optimize(obj, th=0.1):
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
    c4d.EventAdd()


def getFaceS(obj):
    ps = obj.GetPolygonS()
    pl = obj.GetAllPoints()
    fl = obj.GetAllPolygons()
    fc = len(fl)
    sel = ps.GetAll(fc)
    sfl = [fl[i] for i in xrange(fc) if sel[i] == 1]
    return pl, fl, sfl


def toLeft(a, b, p):
    return (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x)


def disxy(a, b):
    v = b - a
    return math.hypot(v.x, v.y)


def toXY(v):
    return c4d.Vector(v.x, v.y, 0)


def isSameDir(v1, v2):
    if v1 == None:
        return True
    return v1.Dot(v2) > 0.998


def getFaceNormal(pl, f):
    a, b, c = pl[f.a], pl[f.b], pl[f.c]
    return ~((c - b).Cross(a - b))


def groupByNormal(pl, sfl):
    n1 = None
    n2 = None
    for f in sfl:
        n = getFaceNormal(pl, f)
        if isSameDir(n1, n):
            n1 = n
        elif isSameDir(n2, n):
            n2 = n
        else:
            print "groupByNormal_Error", f
            return None
    if n1 == None or n2 == None:
        print "groupByNormal_Error", n1, n2
        return None
    angle = utils.VectorAngle(n1, n2)
    if angle < utils.Rad(5) or angle > utils.Rad(175):
        print "angleError", utils.Deg(angle)
        return None
    axis = ~(n1.Cross(n2))
    g1 = FaceGroup(n1)
    g2 = FaceGroup(n2)
    return {"g1": g1, "g2": g2, "axis": axis, "angle": angle}


def getConerAnchor(obj, pl, fl, sfl, data):
    g1, g2, axis = data["g1"], data["g2"], data["axis"]
    angle = data["angle"]
    n1, n2 = g1.n, g2.n
    for f in sfl:
        n = getFaceNormal(pl, f)
        if isSameDir(n1, n):
            g1.fl.append(f)
            continue
        if isSameDir(n2, n):
            g2.fl.append(f)
    for g in (g1, g2):
        for f in g.fl:
            g.pil.update([f.a, f.b, f.c, f.d])
    ipi = g1.pil & g2.pil
    if len(ipi) == 0:
        print "no anchor found", g1.pil, g2.pil
        return None
    g2.pil -= ipi
    anchor = ipi.pop()
    data["anchor"] = anchor
    return True


def flip(obj, pl, data):
    g1, g2, axis = data["g1"], data["g2"], data["axis"]
    gl = [g1, g2]
    angle = math.pi - data["angle"]
    hpb = utils.VectorToHPB(axis)
    m = utils.HPBToMatrix(hpb)
    im = ~m
    anchor = im.MulV(pl[data["anchor"]])
    center = toXY(anchor)
    # to xy
    for g in gl:
        g.pil = list(g.pil)
        g.pc = len(g.pil)
        g.pl = map(lambda x: im.MulV(pl[x]), g.pil)
        g.dl = map(lambda x: disxy(anchor, x), g.pl)
        maxpi = max(xrange(g.pc), key=lambda x: g.dl[x])
        g.anchor = toXY(g.pl[maxpi])
        g.dl = map(lambda x: disxy(g.anchor, x), g.pl)
    rotAxis = c4d.Vector(0, 0, 1)
    side = toLeft(center, g1.anchor, g2.anchor)
    if side < 0:
        g1.rm = utils.RotAxisToMatrix(rotAxis, -angle)
        g2.rm = utils.RotAxisToMatrix(rotAxis, angle)
    else:
        g1.rm = utils.RotAxisToMatrix(rotAxis, angle)
        g2.rm = utils.RotAxisToMatrix(rotAxis, -angle)
    v1 = g1.anchor - center
    v2 = g2.anchor - center
    d1, d2 = v1.GetLength(), v2.GetLength()
    g1.targetLength = d2
    g2.targetLength = d1

    n1, n2 = ~v1, ~v2
    if d1 >= d2:
        g1.off = v1 + v2 - n1 * (d1 - d2)
        g2.off = c4d.Vector(0)
    else:
        g1.off = c4d.Vector(0)
        g2.off = v1 + v2 - n2 * (d2 - d1)

    # flip
    vis = set()
    for g in gl:
        for i, p in enumerate(g.pl):
            tpi, td = g.pil[i], g.dl[i]
            if tpi in vis:
                continue
            if td < 0.1:
                # print "static:",tpi,td
                continue
            if td < g.targetLength + 0.1:
                # print "rot:",tpi,td
                v = p - g.anchor
                v = g.rm.MulV(v)
                g.pl[i] = g.anchor + v
            else:
                # print "offset:",tpi,td
                g.pl[i] = p + g.off
            vis.add(tpi)
    # back,setPoints
    for g in gl:
        g.pl = map(lambda x: m.MulV(x), g.pl)
        for j, pid in enumerate(g.pil):
            obj.SetPoint(pid, g.pl[j])
    obj.Message(c4d.MSG_UPDATE)


def main():
    obj = doc.GetActiveObject()
    if not isPoly(obj):
        print "no polygon selected"
        return
    ps = obj.GetPolygonS()
    if ps.GetCount() < 2:
        print "no enough polygon selected"
        return
    pl, fl, sfl = getFaceS(obj)
    data = groupByNormal(pl, sfl)
    if data == None:
        return
    doc.StartUndo()
    doc.AddUndo(c4d.UNDOTYPE_CHANGE, obj)
    disConnect(obj)
    pl, fl, sfl = getFaceS(obj)
    d2 = getConerAnchor(obj, pl, fl, sfl, data)
    if d2 == None:
        optimize(obj)
        doc.EndUndo()
        return
    flip(obj, pl, data)
    optimize(obj)
    doc.EndUndo()


if __name__ == '__main__':
    main()
