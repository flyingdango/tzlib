# -*- coding: UTF-8 -*-
"""
相似选择
一个选择相似点/面的咸鱼脚本
by飞舞的团子

191216
基本成型

191217
优化对比方法
加入相似面选择
#"""
import c4d
from math import floor, ceil
from c4d import gui
from c4d import utils
from tzlib import utils as tzu
from tzlib import construct as cst
__TITLE__ = "SelectSimilar"
__VERSION__ = 1.0


class MyDialog(gui.GeDialog):

    def CreateLayout(self):
        self.SetTitle(__TITLE__ + "_v" + str(__VERSION__))
        self.GroupBegin(id=1000, flags=c4d.BFH_FIT, cols=2, title="Setting")
        self.GroupBorder(c4d.BORDER_GROUP_IN)
        self.GroupSpace(5, 5)
        self.GroupBorderSpace(10, 0, 10, 4)
        #
        self.AddStaticText(1002, flags=c4d.BFH_LEFT, initw=80, inith=10, name="lenTH:")
        self.AddEditNumberArrows(1003, flags=c4d.BFH_LEFT, initw=100, inith=10)
        self.SetFloat(1003, 10.0, min=1.0, max=100.0, step=0.1, format=c4d.FORMAT_FLOAT)
        self.AddStaticText(1004, flags=c4d.BFH_LEFT, initw=80, inith=10, name="angleTH:")
        self.AddEditNumberArrows(1005, flags=c4d.BFH_LEFT, initw=100, inith=10)
        self.SetFloat(1005, 10.0, min=1.0, max=30.0, step=0.1, format=c4d.FORMAT_FLOAT)
        #
        self.GroupEnd()
        self.AddButton(2020, flags=c4d.BFH_LEFT, initw=180, inith=10, name="GO")
        return True

    def Command(self, id, msg):
        if id == 2020:
            main()
        return True


def getPointInfo(spid, ppm, pl):
    sp = pl[spid]
    nbl = ppm[spid]
    #
    dl = []
    for j in nbl:
        td = int((pl[j] - sp).GetLength())
        dl.append(td)
    dl.sort()
    #
    al = []
    for i in nbl:
        tal = []
        tv = pl[i] - sp
        for j in nbl:
            if i == j:
                continue
            ta = int(utils.Deg(utils.VectorAngle(tv, pl[j] - sp)))
            tal.append(ta)
        tal.sort()
        al.append(tal)
    return len(nbl), dl, al


def getNgonInfo(n, pl):
    tpl = map(lambda x: pl[x], n.pil)
    sp = tzu.avg(tpl)
    #
    dl = []
    for j in tpl:
        td = int((j - sp).GetLength())
        dl.append(td)
    dl.sort()
    #
    al = []
    for i in tpl:
        tal = []
        tv = i - sp
        for j in tpl:
            if i == j:
                continue
            ta = int(utils.Deg(utils.VectorAngle(tv, j - sp)))
            tal.append(ta)
        tal.sort()
        al.append(tal)
    #
    return n.pc, dl, al


def cmplist(sl, tl, th):
    for i, v in enumerate(tl):
        if abs(sl[i] - v) > th:
            return False
    return True


def cmplist2(sll, tll, th):
    vis = set([])
    for i, sl in enumerate(sll):
        for j, tl in enumerate(tll):
            if j in vis:
                continue
            if not cmplist(sl, tl, th):
                continue
            vis.add(j)
            break
    return len(vis) == len(sll)


def sel_points():
    lth = dlg.GetFloat(1003)
    ath = dlg.GetFloat(1005)
    obj = doc.GetActiveObject()
    if not tzu.isPoly(obj):
        print "no polygon found"
        return
    sel = obj.GetPointS()
    sc = sel.GetCount()
    if sc == 0:
        print "no point selected"
        return
    if sc != 1:
        print "too many points selected:", sel.GetCount(), "/1"
        return

    pl, fl = obj.GetAllPoints(), obj.GetAllPolygons()
    pc, fc = len(pl), len(fl)
    ppm = cst.pointPointMap(fl)

    spid = sel.GetAll(fc).index(1)
    spf = getPointInfo(spid, ppm, pl)
    # print spf[2]
    sel.DeselectAll()
    for i, p in enumerate(pl):
        tpf = getPointInfo(i, ppm, pl)
        if tpf[0] != spf[0]:
            continue
        if lth != 1.0 and not cmplist(spf[1], tpf[1], lth):
            continue
        if ath != 1.0 and not cmplist2(spf[2], tpf[2], ath):
            # print tpf[2]
            continue
        sel.Select(i)

    c4d.EventAdd()


def getPolyIndexBySel(sel,n):
    d = sel.GetAll(n)
    tmp = []
    for i, v in enumerate(d):
        if v == 1:
            tmp.append(i)
    return tmp


def sel_polys():
    lth = dlg.GetFloat(1003)
    ath = dlg.GetFloat(1005)
    obj = doc.GetActiveObject()
    if not tzu.isPoly(obj):
        print "no polygon found"
        return

    sel = obj.GetPolygonS()
    if sel.GetCount() <= 0:
        print "no poly selected"
        return
    pl = obj.GetAllPoints()
    fl = obj.GetAllPolygons()
    fc = len(fl)
    el, nl = cst.getAllNgons(obj)
    pnm = cst.polyNgonMap(nl)
    pnm = {k:v[0] for k,v in pnm.items()}
    sfl = getPolyIndexBySel(sel,fc)
    spid = pnm[sfl[0]]
    for i in sfl:
        if pnm[i] != spid:
            print "too many poly selected", sfl
            return
    spf = getNgonInfo(nl[spid], pl)

    for n in nl:
        tpf = getNgonInfo(n, pl)
        if spf[0] != tpf[0]:
            continue
        if lth != 1.0 and not cmplist(spf[1], tpf[1], lth):
            continue
        if ath != 1.0 and not cmplist2(spf[2], tpf[2], ath):
            continue
        for i in n.fil:
            sel.Select(i)
    c4d.EventAdd()


def main():
    mode = doc.GetMode()
    if mode == c4d.Mpoints:
        sel_points()
    elif mode == c4d.Mpolygons:
        sel_polys()
    else:
        print "please use in point/poly mode"
        return

if __name__ == '__main__':
    __anchor__ = "\u98de\u821e\u7684\u56e2\u5b50"
    dlg = MyDialog()
    dlg.Open(c4d.DLG_TYPE_ASYNC, 450, 200, 180, 400)
