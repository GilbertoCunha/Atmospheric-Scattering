import numpy as np


class Vector:

    def __init__(self, x, y, z,  coord):
        self.x = x
        self.y = y
        self.z = z
        self.coord = coord

    def __str__(self):
        return "[" + str(self.x) + ", " + str(self.y) + ", " + str(self.z) + ", " + str(self.coord) + "]"

    def __repr__(self):
        return "[" + str(self.x) + ", " + str(self.y) + ", " + str(self.z) + ", " + str(self.coord) + "]"

    def sph2cart(self):
        if self.coord == 'cart':
            return self
        r = self.x
        polar = self.y
        azimuth = self.z
        return Vector(r*np.sin(azimuth)*np.cos(polar), r*np.sin(azimuth)*np.sin(polar), r*np.cos(azimuth), 'cart')

    def cart2sph(self):
        if self.coord == 'sph':
            return self
        r = np.sqrt(self.x**2 + self.y**2 + self.z**2)
        polar = np.arctan2(self.y, self.x)
        azimuth = np.arctan2(np.sqrt(self.x**2 + self.y**2), self.z)
        return Vector(r, polar, azimuth, 'sph')

    def add(self, other):
        v1, v2 = self.sph2cart(), other.sph2cart()
        v3 = Vector(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z, 'cart')
        if self.coord == 'cart':
            return v3
        else:
            return v3.cart2sph()

    def dot(self, other):
        v1 = self.sph2cart()
        v2 = other.sph2cart()
        return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z

    def norm(self):
        return np.sqrt(self.dot(self))

    def normalize(self):
        if self.coord == 'cart':
            norm = self.norm()
            return Vector(self.x/norm, self.y/norm, self.z/norm, 'cart')
        else:
            return Vector(1, self.y, self.z, 'sph')

    def angle(self, other):
        v1 = self.sph2cart()
        v2 = other.sph2cart()
        aux = v1.dot(v2) / (v1.norm() * v2.norm())
        if aux > 1:
            aux = 1
        return np.arccos(aux)

    def __mull__(self, other):
        v = self.sph2cart()
        if self.coord == 'cart':
            return Vector(other*v.x, other*v.y, other*v.z, 'cart')
        else:
            return (Vector(other*v.x, other*v.y, other*v.z, 'cart')).cart2sph()

    def __rmul__(self, other):
        return self.__mull__(other)

    def dist(self, other):
        v1 = self.sph2cart()
        v2 = -1*other.sph2cart()
        return v1.add(v2).norm()


class Line:

    def __init__(self, point, vec):
        self.point = point
        self.vec = vec

    def __str__(self):
        return str(self.point) + "k * " + str(self.vec)


class Plane:

    def __init__(self, point, normal):
        self.point = point
        self.normal = normal

    def __str__(self):
        return "Point: " + str(self.point) + " | Normal: " + str(self.normal)

    def line_intersection(self, line):
        d = self.point.add(-1 * line.point).dot(self.normal) / line.vec.dot(self.normal)
        return line.point.add(d * line.vec)
