# coding=utf-8
import c4d
import math
from utils import clamp, isMG, isPoly, isTP, isSpline
from collections import namedtuple
from c4d.modules import mograph as mo


Sphere = namedtuple("Sphere", ["p", "r", "r2"])


def getD2(p1, p2):
    return sum([(p1[i] - p2[i])**2 for i in xrange(3)])


def SphereIntersect(a, b):
    return getD2(a.p, b.p) < a.r2 + b.r2 + 2 * a.r * b.r


def getMinPosOnTri(sp, tri0, tri1, tri2):
    edge0 = tri1 - tri0
    edge1 = tri2 - tri0
    v0 = tri0 - sp

    a = edge0.Dot(edge0)
    b = edge0.Dot(edge1)
    c = edge1.Dot(edge1)
    d = edge0.Dot(v0)
    e = edge1.Dot(v0)

    det = a * c - b * b
    s = b * e - c * d
    t = b * d - a * e

    if s + t < det:
        if s < 0:
            if t < 0:
                if d < 0:
                    s = clamp(-d / a, 0, 1)
                    t = 0
                else:
                    s = 0
                    t = clamp(-e / c, 0, 1)
            else:
                s = 0
                t = clamp(-e / c, 0, 1)
        elif t < 0:
            s = clamp(-d / a, 0, 1)
            t = 0
        else:
            invDet = 1 / det
            s *= invDet
            t *= invDet
    else:
        if s < 0:
            tmp0 = b + d
            tmp1 = c + e
            if tmp1 > tmp0:
                numer = tmp1 - tmp0
                denom = a - 2 * b + c
                s = clamp(numer / denom, 0, 1)
                t = 1 - s
            else:
                t = clamp(-e / c, 0, 1)
                s = 0
        elif t < 0:
            if a + d > b + e:
                numer = c + e - b - d
                denom = a - 2 * b + c
                s = clamp(numer / denom, 0, 1)
                t = 1 - s
            else:
                s = clamp(-e / c, 0, 1)
                t = 0
        else:
            numer = c + e - b - d
            denom = a - 2 * b + c
            s = clamp(numer / denom, 0, 1)
            t = 1 - s

    return tri0 + s * edge0 + t * edge1


class Searcher(object):

    def __init__(self):
        self.grid = {}
        self.result = None

    def pos2Index(self, p, s):
        return (int(math.floor(p[0] * s)), int(math.floor(p[1] * s)), int(math.floor(p[2] * s)))

    def index2Str(self, i):
        return str(i[0]) + "_" + str(i[1]) + "_" + str(i[2])

    def str2Index(self, s):
        return tuple(map(lambda a: int(a), s.split("_")))

    def vec2Tuple(self, p):
        return (int(p.x), int(p.y), int(p.z))

    def getNbKeys(self, kl):
        tmp = []
        for k in kl:
            ti = self.str2Index(k)
            tk = self.index2Str((ti[0] + 1, ti[1], ti[2]))
            tmp.append(tk)
            tk = self.index2Str((ti[0], ti[1] + 1, ti[2]))
            tmp.append(tk)
            tk = self.index2Str((ti[0], ti[1], ti[2] + 1))
            tmp.append(tk)
            tk = self.index2Str((ti[0] - 1, ti[1], ti[2]))
            tmp.append(tk)
            tk = self.index2Str((ti[0], ti[1] - 1, ti[2]))
            tmp.append(tk)
            tk = self.index2Str((ti[0], ti[1], ti[2] - 1))
            tmp.append(tk)
        return set(tmp)

    def find(self, p, HQ=False):
        ti = self.pos2Index(p, self.smul)
        tk = self.index2Str(ti)

        sk = tuple([tk])
        touch = False
        tmp = set([])

        for i in xrange(800):
            # 收集所有临近tile到tmp
            tmp = self.getNbKeys(sk)
            # 如果在grid中有tmp元素,说明找到最近一圈tile
            for k in tmp:
                if k in self.grid:
                    touch = True
                    break
            # 如果找到,移除空tile,扩一圈,跳出,否则继续查找
            if touch:
                tmp = set(filter(lambda k: k in self.grid, tmp))
                for j in xrange(int(HQ)):
                    tmp2 = self.getNbKeys(tmp)
                    tmp = tmp.union(tmp2)
                break
            else:
                if len(tmp) > 8000:
                    raise Exception("out of limit", len(tmp))
                sk = tuple(tmp)
                tmp = set([])

        il = []
        for k in tmp:
            if k in self.grid:
                il.extend(self.grid[k]["l"])

        return list(set(il))


class PointSearcher(Searcher):
    """
    PointSearcher最近点搜索by飞舞的团子
    ====接口说明====
    o = PointSearcher()
    创建

    o.init(obj,s=16)
    输入数据,支持vector列表,tp组对象,mg对象,polygon/spline
    s网格大小

    o.search(p,n=1,HQ=false,sort=True)
    搜索,p是搜索位置,n数量,HQ高质量,sort是否结果排序
    ====
    """

    def __init__(self):
        super(PointSearcher, self).__init__()
        self.data = {}

    def init(self, obj, s=16):
        self.gridSize = s
        self.smul = 1.0 / float(s)
        self.setData(obj)
        self.__buildGrid()

    def setData(self, obj):
        if not obj:
            raise Exception("no object found")
        self.data = {}
        # list
        if isinstance(obj, list):
            for i, v in enumerate(obj):
                self.data[i] = self.vec2Tuple(v)
            return

        # tp
        if isTP(obj):
            doc = obj.GetDocument()
            tp = doc.GetParticleSystem()
            il = obj.GetParticles()
            for i in il:
                self.data[i] = self.vec2Tuple(tp.Position(i))
            return
        # mg
        if isMG(obj):
            md = mo.GeGetMoData(obj)
            if md is None:
                raise Exception("no mg data")
            marr = md.GetArray(c4d.MODATA_MATRIX)
            for i, v in enumerata(marr):
                self.data[i] = self.vec2Tuple(v.off)
            return
        # poly/spline
        if isPoly(obj) or isSpline(obj):
            pa = obj.GetAllPoints()
            for i, v in enumerata(pa):
                self.data[i] = self.vec2Tuple(v)

        raise Exception("unknow object type", t)

    def __buildGrid(self):
        self.grid = {}
        for i, p in self.data.iteritems():
            ti = self.pos2Index(p, self.smul)
            tk = self.index2Str(ti)
            self.grid.setdefault(tk, {"i": ti, "l": []})["l"].append(i)

    def search(self, v, r=0, n=1, HQ=False, sort=True):
        if len(self.grid) <= 0:
            print "grid not init"
            return None

        p = self.vec2Tuple(v)
        il = self.find(p, HQ)
        self.result = None

        dis = {}
        for i in il:
            dis[i] = getD2(p, self.data[i])

        if r > 0:
            r2 = r**2
            il = filter(lambda x: dis[x] < r2, il)

        if len(il) <= 0:
            return None

        if sort:
            il.sort(key=lambda x: dis[x])
        self.result = il[0] if n == 1 else tuple(il[:n])
        return self.result


class PolySearcher(Searcher):
    """
    物体表面最近点搜索by飞舞的团子
    ====接口说明====
    o = PolySearcher()
    初始化

    o.init(obj,s=32,q=0.25)
    给模型,可以是polygon或者塌陷变形器
    s网格大小,q精度(越大越不精确,一般不用改)

    o.minpos(p,n=1,HQ=False,sort=True)
    获取给定点p到模型的最近点,n数量,HQ高质量,sort结果排序
    """

    def __init__(self):
        super(PolySearcher, self).__init__()
        self.polys = None
        self.points = None

    def sampleTri(self, i, a, b, c, mid=False):
        l = [(a, b, (a - b).GetLengthSquared()),
             (b, c, (b - c).GetLengthSquared()),
             (c, a, (c - a).GetLengthSquared())]

        l.sort(key=lambda x: x[2])

        if l[2][2] < self.s2:
            self.__pushPoint(i, a)
            if mid:
                self.__pushPoint(i, b)
                self.__pushPoint(i, c)
            return
        if l[1][2] < self.s2:
            self.__pushPoint(i, a)
            self.__pushPoint(i, (l[2][0] + l[2][1]) * .5)
            if mid:
                self.__pushPoint(i, b)
                self.__pushPoint(i, c)
            return

        d, e, f = (a + b) * .5, (b + c) * .5, (c + a) * .5
        self.sampleTri(i, e, f, d, True)
        self.sampleTri(i, a, d, f, False)
        self.sampleTri(i, b, e, d, False)
        self.sampleTri(i, c, f, e, False)

    def __pushPoint(self, i, v):
        ti = self.pos2Index(v, self.smul)
        tk = self.index2Str(ti)
        self.grid.setdefault(tk, {"l": set([]), "i": ti})["l"].add(i)

    def init(self, obj, s=32, Q=0.25):
        if not isPoly(obj):
            raise Exception("no object found")

        self.gridSize = s
        self.s2 = Q * s**2
        self.smul = 1.0 / float(s)
        m = obj.GetMg()
        self.polys = []
        self.points = obj.GetAllPoints()
        self.points = map(lambda x: m.Mul(x), self.points)
        tmp = obj.GetAllPolygons()
        for f in tmp:
            self.polys.append((f.a, f.b, f.c))
            if not f.IsTriangle():
                self.polys.append((f.a, f.c, f.d))

        if len(self.polys) == 0:
            raise Exception("no polygon found")

        for i, f in enumerate(self.polys):
            self.sampleTri(i, self.points[f[0]], self.points[f[1]], self.points[f[2]])

    def minpos(self, v, r=0, n=1, HQ=False, sort=True):
        if len(self.grid) <= 0:
            print "grid not init"
            return None

        p = self.vec2Tuple(v)
        il = self.find(p, HQ)
        self.result = None

        dis = {}

        for i in il:
            f = self.polys[i]
            thit = getMinPosOnTri(v, self.points[f[0]], self.points[f[1]], self.points[f[2]])
            dis[i] = {"d2": (v - thit).GetLengthSquared(), "hit": thit}

        if r > 0:
            r2 = r**2
            il = filter(lambda x: dis[x]["d2"] < r2, il)

        if len(il) <= 0:
            return None

        if sort:
            il.sort(key=lambda x: dis[x]["d2"])

        if n == 1:
            self.result = dis[il[0]]["hit"]
        else:
            il = il[:n]
            self.result = tuple(map(lambda x: dis[x]["hit"], il))
        return self.result
