import xml.etree.ElementTree as ET
import sys
from lxml import etree
import numpy
from Point import Point
import math
import svgwrite

from svg.path import Path, Line, Arc, CubicBezier, QuadraticBezier, parse_path

# Find index i such that dd[i] < d < dd[i+1]
def FindIndex(dd, d):
	'''dd is an array of doubles. d is a double '''

	i = 0
 	for j in xrange(0, len(dd)):
		if (d > dd[j]):
			i = j;

	return i;

class BezPoints:
	'''A wrapper around svg.path's CubicBezier class that adds some useful helper methods.'''

	def __init__(self, cb):
		''' Input is a CubicBezier object. '''

		self.cb = cb

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


	def interpolate_points(self, n=10000, m=10):
		'''n is the amount of precision desired. m is the number of points to generate'''
		#Construct polyline with large number of points
		s = 1.0/(n-1);
		tt = [0] * n;	#An array of doubles
		PP = [None] * n; #An array of points
		cc = [0] * n;	#An array of doubles

		for i in xrange(0, n):
			tt[i] = i*s;
			PP[i] = self.get_point_from_t(tt[i]);
			if (i > 0):
				cc[i] = PP[i].distance(PP[i-1]);

		#Get fractional arclengths along polyline
		dd = [0] * n #An array of doubles
		dd[0] = 0
		for i in xrange(1, n):
			dd[i] = dd[i-1] + cc[i];
		for i in xrange(1, n):
			dd[i] = dd[i]/dd[n-1];

		#Number of points to place on curve
		step = 1.0/(m-1)

		QQ = [None] * m; #An array of points

		for r in xrange(0, m):
			d = r*step;
			i = FindIndex(dd, d);
			u = (d - dd[i]) / (dd[i+1] - dd[i]);
			t = (i + u)*s;
			QQ[r] = self.get_point_from_t(t);

		return QQ


def get_path_from_svg(svg_file):
	print "Getting the path"

	tree = etree.parse(svg_file)

	all_paths = tree.findall("{*}path")
	print all_paths


	if len(all_paths) == 0:
		return None
	else:
		path_node = all_paths[0]
		path_string = path_node.get("d")

		path = parse_path(path_string)
		print(path)

		return path


def get_equidistant_points(path_points):
	li = Length_iterator(4, 0.01)
	li.initializeIterationOnCurve(path_points, 8)

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print "Input svg path"
		exit(1)

	path = get_path_from_svg(sys.argv[1])
	print("INFO")
	print(dir(path))
	#get_equidistant_points(path, 20)

	print("")

	#Accumulate all of the cubic bezier curves in the svg.
	bez_curves = []
	for p in path:
		if type(p) is CubicBezier:
			print "CUBEY: " + str(p) + " -!--! " + str(type(p.start))
			bez_curves.append(BezPoints(p))
		else:
			print "NOPE"

	#Paramaterize the curve
	for curve in bez_curves:
		interped_points = curve.interpolate_points()

		for q in interped_points:
			print str(q.x) + "\t" + str(q.y)
		print "\n"

		#write output svg 
		output_svg = svgwrite.Drawing('output.svg')
		for i in xrange(0, len(interped_points)-1):
			p1 = interped_points[i]
			p2 = interped_points[i+1]
			output_svg.add(output_svg.line((p1.x, p1.y), (p2.x, p2.y), stroke=svgwrite.rgb(10, 10, 16, '%')))
		output_svg.save()
	

