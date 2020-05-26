# coding=utf-8
import c4d
import math
import utils as tzu
from c4d.utils import SplineHelp


class Voxel(object):

    def __init__(self, _i, _k):
        self.k = _k
        self.i = _i

    def __repr__(self):
        return self.k

    def __hash__(self):
        return hash(self.k)


class Shape(object):

    def __init__(self, _spline, _doc, _vd=1, _vs=5, _plane=0):
        self.doc = _doc
        self.vs = _vs
        self.ivs = 1.0 / float(_vs)
        self.voff = c4d.Vector(_vs) * .5
        self.plane = _plane

        self.spline = _spline
        if isinstance(_vd, dict):
            self.voxel_data = _vd
        elif _vd == 1:
            self.init_voxel_data()
        elif _vd == 2:
            self.init_voxel_data2()
        else:
            assert("unknow _vd")
        self.uuid = str(id(self.spline))
        self.context = {}

    def __repr__(self):
        return self.uuid

    def __hash__(self):
        return hash(self.uuid)

    def init_voxel_data(self):
        # 先获取边缘一圈的voxel，再交叉扫描填充，获取的voxel不会存在空洞
        vd = {}
        sh = SplineHelp()
        sh.InitSpline(self.spline)
        l = sh.GetSplineLength()
        # if len(vd) < 16 * self.vs:
        #     assert("spline too small")
        seg_cnt = int(math.ceil(l * self.ivs) * 3)
        step = 1.0 / float(seg_cnt)

        for i in xrange(0, seg_cnt + 5):
            t = step * i
            p = sh.GetPosition(offset=t, realoffset=True)
            # show_pos(p)
            vi, k = self.pos2vpos(p)
            if k in vd:
                continue
            v = Voxel(vi, k)
            vd[k] = v
        # fill voxel
        xl = {}
        yl = {}
        for k, v in vd.items():
            xk = int(v.i.x)
            yk = int(v.i.y)
            xl.setdefault(xk, set([]))
            yl.setdefault(yk, set([]))
            xl[xk].add(int(v.i.y))
            yl[yk].add(int(v.i.x))

        # 垂直扫描
        vset = set([])
        for xi, l in xl.items():
            syi = min(l)
            eyi = max(l)
            for yi in xrange(syi, eyi + 1):
                vset.add((xi, yi))
        # 水平扫描
        hset = set([])
        for yi, l in yl.items():
            sxi = min(l)
            exi = max(l)
            for xi in xrange(sxi, exi + 1):
                hset.add((xi, yi))
        # 交集求alpha_shape
        fset = vset & hset
        vd2 = {}
        for t in fset:
            vi = c4d.Vector(t[0], t[1], 0)
            tk = self.index2key(vi)
            vd2[tk] = Voxel(vi, tk)
        self.voxel_data = vd2

    def init_voxel_data2(self):
        # 备用方法，利用挤出获取voxel，能保留空洞
        s = self.spline.GetClone()
        s[c4d.SPLINEOBJECT_INTERPOLATION] = 4
        s[c4d.SPLINEOBJECT_MAXIMUMLENGTH] = self.vs * .5
        ex = c4d.BaseObject(5116)
        ex[c4d.EXTRUDEOBJECT_MOVE] = c4d.Vector(0)
        ex[c4d.CAP_END] = 0
        ex[c4d.CAP_TYPE] = 0
        ex[c4d.CAP_REGULAR] = True
        ex[c4d.CAP_REGULARWIDTH] = self.vs * .5
        s.InsertUnder(ex)
        co = tzu.getPoly(ex, self.doc)
        if not co:
            assert("extrude getPoly error")
        pl = co.GetAllPoints()
        vd = {}
        for p in pl:
            vi, k = self.pos2vpos(p)
            if k in vd:
                continue
            v = Voxel(vi, k)
            vd[k] = v
        self.voxel_data = vd

    def index2key(self, i):
        return str(int(i.x)) + "_" + str(int(i.y))

    # def checkbound(self, i):
    #     self.xmin = float("inf")
    #     self.xmax = -float("inf")
    #     self.ymin = float("inf")
    #     self.ymax = -float("inf")
    #     if i.x < self.xmin:
    #         self.xmin = i.x
    #     elif i.x > self.xmax:
    #         self.xmax = i.x
    #     if i.y < self.ymin:
    #         self.ymin = i.y
    #     elif i.y > self.ymax:
    #         self.ymax = i.y

    def pos2vpos(self, p):
        x = int(math.floor(p.x * self.ivs))
        y = int(math.floor(p.y * self.ivs))
        #z = int(math.floor(p.z * self.ivs))
        return c4d.Vector(x, y, 0), str(x) + "_" + str(y)

    def __sub__(self, other):
        # spline_bool
        rs = tzu.splineMask2(self.spline, other.spline, self.doc)
        # if not rs:
        #     assert("BooleanSplines Error")
        rvd = {}
        for k in self.voxel_data:
            if k not in other.voxel_data:
                rvd[k] = self.voxel_data[k]
        # return new
        return Shape(rs, self.doc, rvd)

    def getOverlap(self, other):
        tmp = {}
        for k, v in self.voxel_data.items():
            if k in other.voxel_data:
                tmp[k] = v
        return tmp

    def show(self, doc, col=None):
        if not col:
            col = c4d.Vector(0, 1, 1)
        ori = int(self.plane + 1)
        pn = c4d.BaseObject(5140)
        pn.SetName(self.uuid)
        for k, v in self.voxel_data.items():
            n = c4d.BaseObject(5140)
            tzu.setObjectColor(n, col)
            n[c4d.NULLOBJECT_DISPLAY] = 3
            n[c4d.NULLOBJECT_RADIUS] = self.vs * 0.68
            n[c4d.NULLOBJECT_ORIENTATION] = ori
            n[c4d.ID_BASEOBJECT_REL_POSITION] = v.i * self.vs + self.voff
            n.SetName(v.k)
            n.InsertUnder(pn)
        doc.InsertObject(pn)
        doc.InsertObject(self.spline)
        c4d.EventAdd()
