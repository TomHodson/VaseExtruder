from __future__ import division
from vectors import vec3, cross
from math import cos, sin, pi
from itertools import chain
def memodict(f):
    """ Memoization decorator for a function taking a single argument """
    class memodict(dict):
        def __missing__(self, key):
            ret = self[key] = f(key)
            return ret 
    return memodict().__getitem__

def memoize(f):
    """ Memoization decorator for a function taking one or more arguments. """
    class memodict(dict):
        def __getitem__(self, *key):
            return dict.__getitem__(self, key)

        def __missing__(self, key):
            ret = self[key] = f(*key)
            return ret

    return memodict().__getitem__

def flip(triangles):
	#flip the vertex order of a list of triangles
	return map(lambda x : x[::-1], triangles)

def circle(r, n, theta = 0):
	#return a list of n points describing a circle of radius r
	return tuple(vec3(r*cos(2*pi*i/n + theta), r*sin(2*pi*i/n + theta), 0) for i in range(n))

def naive_triangulisation(shape):
	return [[shape[0], shape[n+1], shape[n]] for n in range(1, len(shape) - 1)]

def polytotri(shape):
	z = shape[0].z
	from p2t import CDT, Point, Triangle
	polyline = [Point(p.x, p.y) for p in shape]
	cdt = CDT(polyline)
	triangles = cdt.triangulate()
	points = [[t.a, t.b, t.c] for t in triangles]
	return [[vec3(p.x, p.y, z) for p in t] for t in points]


def lerp(a,b,l):
	return a*(1-l) + b*l
def lerp_shapes(a,b,l):
	return [i*(1-l) + j*l for i,j in zip(a,b)]

#@memoize
def kochify(shape, iterations, max_iterations = None):
	#applys one iteration of the koch snowflake to an arbitray pointlist.
	if max_iterations == None: max_iterations = iterations
	if max_iterations <= 0: return tuple(shape)
	if 1 <= iterations: l = 1
	if 0 < iterations < 1: l = iterations
	if iterations <= 0: l = 0

	newshape = []
	for i in range(len(shape)):
		a = shape[i]
		b = shape[(i+1)%len(shape)]
		v = b - a
		newshape.append(a)
		newshape.append(lerp(a,b,1/3))
		newshape.append(lerp(a,b,1/2) + l * 1/3 * vec3(v.y, -v.x, 0))
		newshape.append(lerp(a,b,2/3))
	return tuple(kochify(newshape, iterations - 1, max_iterations - 1))

def translate(shape, vector):
	return tuple(point + vector for point in shape)

def join(a, b):
	#assert(len(a) == len(b))
	def half_triangles(a,b,shift): return [[ b[i], a[i], b[(i+shift) % len(a)] ] for i in range(len(a))]
	#return zip(b, a, chain(b[1:], [b[0],])) + zip(a, b, chain([a[-1],],a[:-1]))
	return half_triangles(a,b,+1) + half_triangles(b,a,-1)

def normal(triangle):
	a,b,c = triangle
	return cross(c-a, b-a)

def extrude(shape_func, samples = 50):
	shape_func = memodict(shape_func)
	triangles = []
	shapes = [shape_func(i/samples) for i in range(samples+1)]
	for i in range(samples):
		triangles += join(shapes[i], shapes[i+1])
	triangles += flip(polytotri(shapes[0]))# + polytotri( shapes[-1])
	return triangles

def triangles_to_stl(triangles):
	s = """
solid Model
{}
endsolid Model
"""
	vertex = \
"""
facet normal {n[0]:.5f} {n[1]:.5f} {n[2]:.5f}
	outer loop
  	vertex {t[0][0]:.5f} {t[0][1]:.5f} {t[0][2]:.5f}
  	vertex {t[1][0]:.5f} {t[1][1]:.5f} {t[1][2]:.5f}
  	vertex {t[2][0]:.5f} {t[2][1]:.5f} {t[2][2]:.5f}
	endloop
endfacet
"""
	vertices = "".join(vertex.format(n = normal(t), t = t) for t in triangles)
	return s.format(vertices)

def triangles_to_binary_stl(triangles):
	from struct import Struct, pack
	header = b"\x00" * 80 + pack("<I", len(triangles))
	out = header
	for t in triangles:
		#n = normal(t)
		n = vec3(0,0,0)
		data = list(n) + list(chain(*t)) + list([0,])
		s = Struct("<" + "f" * 12 + "H")
		out += s.pack(*data)
	return out

def rotating_koch(i):
	triangle = circle(sin(2*pi*i) + 1, 3, i*pi)
	koch = kochify(triangle, 3)
	return translate(koch, vec3(0,0,i*10))
def koch_to_circle(i):
	samples = 3
	iterations = 2
	height = 3
	radius = sin(0.8*pi*i)**2 + 0.2
	lerp_function = i**2
	spin = i* pi/6

	c = circle(radius, samples * 4**iterations, spin)
	koch = kochify(circle(radius, samples, spin), iterations)
	l = lerp_shapes(c, koch, lerp_function)
	return translate(l, [0,0,i*height])

def koch_circle_oscillations(i):
	samples = 3
	iterations = 2
	height = 5
	radius = (cos(2*pi*i)**2 + 1) / 2 if i < 0.5 else  cos(2*pi*i)
	lerp_function = sin(2*pi*i)**2
	spin = i* pi/2

	c = circle(radius, samples * 4**iterations, spin)
	koch = kochify(circle(radius, samples, spin), iterations)
	l = lerp_shapes(c, koch, lerp_function)
	return translate(l, [0,0,i*height])

def koch_growth(i):
	samples = 3
	iterations = 4
	height = 3
	radius = sin(pi*i)**2 + 0.4*2*(0.5-i)**2 if i < 0.5 else sin(pi*i) + 0.3*2*(0.5-i)**2
	lerp_function = i
	spin = i* pi/3

	c = circle(radius, samples * 4**iterations, spin)
	koch = kochify(circle(radius, samples, spin), sin(pi*i)*iterations, iterations)
	l = lerp_shapes(c, koch, lerp_function)
	return translate(l, [0,0,i*height])

def koch_plant_pot(i):
	samples = 8
	iterations = 2
	height = 100.0
	g=0.575; b=0.6; c=0.46; d=0.41
	radius = 75.0 * (g*sin(b*pi*i + c)**2 + d) / 2 / (g*sin(b*pi*1.0 + c)**2 + d)

	n = 4.0
	lerp_function = 0.6 * i**n * (1 - i)**n / 0.5**(2*n)
	spin = i* pi/2

	c = circle(radius, samples * 4**iterations, spin)
	koch = kochify(circle(radius, samples, spin), iterations)
	l = lerp_shapes(c, koch, lerp_function)
	return translate(l, [0,0,i*height])

def nice_swirl(i):
	#this one is really nice, it has small ribs going up the side
	samples = 50
	iterations = 2
	height = 100.0
	g=0.575; b=0.6; c=0.46; d=0.41
	radius = 74.0 * (g*sin(b*pi*i + c)**2 + d) / 2 / (g*sin(b*pi*1.0 + c)**2 + d)
	spin = i* pi/2

	n = 4.0
	lerp_function = 0.6 * i**n * (1 - i)**n / 0.5**(2*n)

	c = circle(radius, samples * 4**iterations, spin)
	koch = kochify(circle(radius, samples, spin), iterations*lerp_function, max_iterations = iterations)
	return translate(koch, [0,0,i*height])

def smooth_bulb(i):
	from math import sqrt
	samples = 100
	def top_radius(i):
			b=0.56; d=0.71
			c = pi/2.0*(1-b)
			return 0.5 * (sin(b*pi*i + c)**2 + d) / (sin(b*pi*0.5 + c)**2 + d)
	
	def bottom_radius(i):
		j = 0.36
		return sqrt(0.5**2 - (i - 0.5)**2/(1+j))
	
	radius = 74.0 / 2.0 * (top_radius(i) if i > 0.5 else bottom_radius(i)) / top_radius(1.0)
	spin = i* pi/2
	height = 74.0 / 2.0 * 1.0 / top_radius(1.0)

	c = circle(radius, samples, spin)
	return translate(c, [0,0,i*height])

def make_stl():
	with open("smooth_bulb.stl", 'w') as stl:
		surface = extrude(bulb, samples = 50)
		stl.write(triangles_to_binary_stl(surface))
if __name__ == '__main__':
	#import cProfile
	#cProfile.run("make_stl()")
	make_stl()