import c4d


class Plane(object):

    def __init__(self):
        self.normal = c4d.Vector(0, 1, 0)
        self.constant = 0


class Tri(object):

    def __init__(self, _p0, _p1, _p2):
        self.i = -1
        self.p0 = _p0
        self.p1 = _p1
        self.p2 = _p2

    def __str__(self):
        return "Tri(" + str(self.i) + ")"


class Sphere(object):

    def __init__(self, _p, _r):
        self.p = _p
        self.r = _r
        self.r2 = _r**2
        self.s = _r * 2

    def __str__(self):
        return "Sphere(" + str(self.r) + ")"


class Box(object):

    def __init__(self, _a, _b, initMode=0):
        # 0:pos,radius,1:min,max
        if initMode == 0:
            self.p = _a
            if str(type(_b)) == "<type 'c4d.Vector'>":
                self.r = _b
            else:
                self.r = c4d.Vector(_b)
            self.s = self.r * 2
            self.minp = _a - _b
            self.maxp = _a + _b
        else:
            self.minp = _a
            self.maxp = _b
            self.p = (_a + _b) * 0.5
            self.s = _b - _a
            self.r = self.s * 0.5

    def __str__(self):
        # return "Box(" + str(self.minp)+" "+str(self.maxp) + ")"
        return "Box()"
