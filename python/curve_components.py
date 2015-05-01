from svg.path import Path, Line, Arc, CubicBezier, QuadraticBezier, parse_path
from point import Point

def FindIndex(dd, d):
	'''A helper method to find index i such that dd[i] < d < dd[i+1]'''

	i = 0
 	for j in xrange(0, len(dd)):
		if (d > dd[j]):
			i = j;

	return i;

class Component(object):
	'''
	A base class that all SVG components that are supported. Right now this is only lines and Bezier curves.
	'''

	def get_length(self):
		'''Get the length of this segment. For curves, this is the length around the curve.'''
		return 0

	def interpolate_points(self, fractions):
		'''Returns a list of points at the specified fractions along the line.
			For example if [0.5] was inputted, a point halfway along the line would be returned.'''
		return None

	def get_start_point(self):
		'''Get the first point in the component'''
		return None

	def get_end_point(self):
		'''Get the last point of the component'''
		return None

	def reverse(self):
		'''Reverses the start and end points of the component.'''
		return

class MyLine(Component):
	'''A wrapper around svg.path's Line class that adds some useful helper methods.'''

	def __init__(self, ln):
		'''
		:param ln: A Line object.
		'''

		self.ln = ln

	def get_length(self):
		return self.ln.length()

	def interpolate_points(self, fractions):
		QQ = []
		for fraction in fractions:
			p = self.ln.point(fraction)
			QQ.append( Point(p.real, p.imag) )

		return QQ

	def get_start_point(self):
		pt = self.ln.point(0)
		return Point(pt.real, pt.imag)

	def get_end_point(self):
		pt = self.ln.point(1)
		return Point(pt.real, pt.imag)

	def reverse(self):
		self.ln.start, self.ln.end = self.ln.end, self.ln.start

class MyBezCurve(Component):
	'''A wrapper around svg.path's CubicBezier class that adds some useful helper methods.'''

	def __init__(self, cb):
		'''
		:param cb: A CubicBezier object.
		'''

		self.cb = cb
		self.cached_arc_lengh = None

		self.precision = 1000

	def get_points(self):
		pts = []
		pts.append(Point(self.cb.start.real, self.cb.start.imag))
		pts.append(Point(self.cb.control1.real, self.cb.control1.imag))
		pts.append(Point(self.cb.control2.real, self.cb.control2.imag))
		pts.append(Point(self.cb.end.real, self.cb.end.imag))

		return pts

	def get_point_from_t(self, t):
		pos = self.cb.point(t);
		return Point(pos.real, pos.imag)

	def get_length(self):
		if self.cached_arc_lengh == None:
			n = self.precision
			s = 1.0/(n-1);

			PP = [None] * n; #An array of points
	
			self.cached_arc_lengh = 0

			for i in xrange(0, n):
				PP[i] = self.get_point_from_t(i * s);
				if (i > 0):
					self.cached_arc_lengh += PP[i].distance(PP[i-1]);

		return self.cached_arc_lengh


	def interpolate_points(self, fractions):
		'''Input is alist of values between 0 and 1 indicates distance along the arclength of the curve.
			This function returns a list of the points at each of these distances.'''

		n = self.precision

		'''n is the amount of precision desired. m is the number of points to generate'''
		#Construct polyline with large number of points
		s = 1.0/(n-1);
		t_Values = [0] * n;	#An array of doubles
		points_on_curve = [None] * n; #An array of points
		dist_to_prev = [0] * n;	#An array of doubles

		for i in xrange(0, n):
			t_Values[i] = i * s;
			points_on_curve[i] = self.get_point_from_t(t_Values[i]);
			if (i > 0):
				dist_to_prev[i] = points_on_curve[i].distance(points_on_curve[i-1]);

		#Get fractional arclengths along polyline
		arclength_to_point = [0] * n #An array of doubles
		arclength_to_point[0] = 0
		for i in xrange(1, n):
			arclength_to_point[i] = arclength_to_point[i-1] + dist_to_prev[i];

		for i in xrange(1, n):
			arclength_to_point[i] = arclength_to_point[i]/arclength_to_point[n-1];

		QQ = [] #An array of points

		for d in fractions:
			i = FindIndex(arclength_to_point, d);
			u = (d - arclength_to_point[i]) / (arclength_to_point[i+1] - arclength_to_point[i]);
			t = (i + u)*s;
			QQ.append(self.get_point_from_t(t))

		return QQ

	def get_start_point(self):
		return self.get_point_from_t(0)

	def get_end_point(self):
		return self.get_point_from_t(1)

	def reverse(self):
		self.cb.control1, self.cb.control2 = self.cb.control2, self.cb.control1
		self.cb.start, self.cb.end = self.cb.end, self.cb.start