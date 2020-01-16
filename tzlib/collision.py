# coding=utf-8
"""
test doc
"""
import c4d
from base import *


def test():
    p = Plane()
    print p
    print "test ok"


def IsInDistance2(a, b, dis2):
    d = (b - a).GetLengthSquared()
    return d <= dis2


def IsPointInBox(p, box):
    return (p.x >= box.minp.x and p.x <= box.maxp.x) and \
        (p.y >= box.minp.y and p.y <= box.maxp.y) and \
        (p.z >= box.minp.z and p.z <= box.maxp.z)


def IsPointInSphere(p, sphere):

    return IsInDistance2(p, sphere.p, sphere.r2)


def IsPointInTri(p, tri):
    a = tri.p0
    b = tri.p1
    c = tri.p2

    v0 = a - p
    v1 = b - p
    v2 = c - p

    u = v1.Cross(v2)
    v = v2.Cross(v0)
    if (u.Dot(v) < 0):
        return False

    w = v0.Cross(v1)
    if (u.Dot(w) < 0):
        return False

    return True


def BoxPlaneIntersect(box, plane):
    center = box.p
    extents = box.maxp - center
    r = extents.x * abs(plane.normal.x) + extents.y * abs(plane.normal.y) + extents.z * abs(plane.normal.z)
    s = plane.normal.Dot(center) - plane.constant
    return abs(s) <= r


def BoxBoxIntersect(a, b):
    return (a.minp.x <= b.maxp.x and a.maxp.x >= b.minp.x) and \
        (a.minp.y <= b.maxp.y and a.maxp.y >= b.minp.y) and \
        (a.minp.z <= b.maxp.z and a.maxp.z >= b.minp.z)


def SphereSphereIntersect(a, b):
    d2 = (b.p - a.p).GetLengthSquared()
    return d2 <= (a.r2 + b.r2 + 2 * a.r * b.r)


def TriBoxIntersect(tri, box):
    a = tri.p0
    b = tri.p1
    c = tri.p2

    center = box.p
    extents = box.maxp - center
    v0 = a - center
    v1 = b - center
    v2 = c - center
    f0 = v1 - v0
    f1 = v2 - v1
    f2 = v0 - v2
    a00 = c4d.Vector(0, -f0.z, f0.y)
    a01 = c4d.Vector(0, -f1.z, f1.y)
    a02 = c4d.Vector(0, -f2.z, f2.y)
    a10 = c4d.Vector(f0.z, 0, -f0.x)
    a11 = c4d.Vector(f1.z, 0, -f1.x)
    a12 = c4d.Vector(f2.z, 0, -f2.x)
    a20 = c4d.Vector(-f0.y, f0.x, 0)
    a21 = c4d.Vector(-f1.y, f1.x, 0)
    a22 = c4d.Vector(-f2.y, f2.x, 0)

    # Test axis a00
    p0 = v0.Dot(a00)
    p1 = v1.Dot(a00)
    p2 = v2.Dot(a00)
    r = extents.y * abs(f0.z) + extents.z * abs(f0.y)
    if (max(-max(p0, p1, p2), min(p0, p1, p2)) > r):
        return False

    # Test axis a01
    p0 = v0.Dot(a01)
    p1 = v1.Dot(a01)
    p2 = v2.Dot(a01)
    r = extents.y * abs(f1.z) + extents.z * abs(f1.y)
    if (max(-max(p0, p1, p2), min(p0, p1, p2)) > r):
        return False

    # Test axis a02
    p0 = v0.Dot(a02)
    p1 = v1.Dot(a02)
    p2 = v2.Dot(a02)
    r = extents.y * abs(f2.z) + extents.z * abs(f2.y)
    if (max(-max(p0, p1, p2), min(p0, p1, p2)) > r):
        return False

    # Test axis a10
    p0 = v0.Dot(a10)
    p1 = v1.Dot(a10)
    p2 = v2.Dot(a10)
    r = extents.x * abs(f0.z) + extents.z * abs(f0.x)
    if max(-max(p0, p1, p2), min(p0, p1, p2)) > r:
        return False

    # Test axis a11
    p0 = v0.Dot(a11)
    p1 = v1.Dot(a11)
    p2 = v2.Dot(a11)
    r = extents.x * abs(f1.z) + extents.z * abs(f1.x)
    if max(-max(p0, p1, p2), min(p0, p1, p2)) > r:
        return False

    # Test axis a12
    p0 = v0.Dot(a12)
    p1 = v1.Dot(a12)
    p2 = v2.Dot(a12)
    r = extents.x * abs(f2.z) + extents.z * abs(f2.x)
    if max(-max(p0, p1, p2), min(p0, p1, p2)) > r:
        return False

    # Test axis a20
    p0 = v0.Dot(a20)
    p1 = v1.Dot(a20)
    p2 = v2.Dot(a20)
    r = extents.x * abs(f0.y) + extents.y * abs(f0.x)
    if max(-max(p0, p1, p2), min(p0, p1, p2)) > r:
        return False

    # Test axis a21
    p0 = v0.Dot(a21)
    p1 = v1.Dot(a21)
    p2 = v2.Dot(a21)
    r = extents.x * abs(f1.y) + extents.y * abs(f1.x)
    if max(-max(p0, p1, p2), min(p0, p1, p2)) > r:
        return False

    # Test axis a22
    p0 = v0.Dot(a22)
    p1 = v1.Dot(a22)
    p2 = v2.Dot(a22)
    r = extents.x * abs(f2.y) + extents.y * abs(f2.x)
    if max(-max(p0, p1, p2), min(p0, p1, p2)) > r:
        return False

    # Test the three axes corresponding to the face normals of AABB b
    if max(v0.x, v1.x, v2.x) < -extents.x or min(v0.x, v1.x, v2.x) > extents.x:
        return False
    if max(v0.y, v1.y, v2.y) < -extents.y or min(v0.y, v1.y, v2.y) > extents.y:
        return False
    if max(v0.z, v1.z, v2.z) < -extents.z or min(v0.z, v1.z, v2.z) > extents.z:
        return False

    # Test separating axis corresponding to triangle face normal
    plane = Plane()
    plane.normal = f1.Cross(f0).GetNormalized()
    plane.constant = plane.normal.Dot(a)
    return BoxPlaneIntersect(box, plane)


def SphereBoxIntersect(sphere, box):
    x = max(box.minp.x, min(sphere.p.x, box.maxp.x))
    y = max(box.minp.y, min(sphere.p.y, box.maxp.y))
    z = max(box.minp.z, min(sphere.p.z, box.maxp.z))
    return IsInDistance2(c4d.Vector(x, y, z), sphere.p, sphere.r2)


def lineLineIntersect(p1, p2, p3, p4):
    one = 0.998
    zero = 1.0 - one
    v1, v2 = p2 - p1, p4 - p3
    if v1.GetLengthSquared() < 1 or v2.GetLengthSquared() < 1:
        # 长度不足
        return False
    n1, n2 = ~v1, ~v2
    if abs(n1.Dot(n2)) > one:
        # 平行
        return False
    v3 = p3 - p1
    n3 = ~(v1.Cross(v3))
    if not abs(n3.Dot(n2)) < zero:
        # 不共面
        return False
    vs1 = v1.Cross(v2)
    vs2 = v3.Cross(v2)
    t = vs1.Dot(vs2) / vs1.GetLengthSquared()
    return p1 + v1 * t


def lineSegIntersect(p1, p2, p3, p4,ends = True):
    one = 0.998
    zero = 1.0 - one
    v1, v2 = p2 - p1, p4 - p3
    if v1.GetLengthSquared() < 1 or v2.GetLengthSquared() < 1:
        # 长度不足
        return False
    n1, n2 = ~v1, ~v2
    if abs(n1.Dot(n2)) > one:
        # 平行
        return False
    v3 = p3 - p1
    n3 = ~(v1.Cross(v3))
    if not abs(n3.Dot(n2)) < zero:
        # 不共面
        return False
    vs1 = v1.Cross(v2)
    vs2 = v3.Cross(v2)
    t = vs1.Dot(vs2) / vs1.GetLengthSquared()
    rp = p1 + v1 * t
    isinter = (~(p1 - rp)).Dot(~(p2 - rp))
    if ends and isinter == 0 or isinter < -one:
        return rp
    return False
