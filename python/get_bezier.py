import xml.etree.ElementTree as ET
import sys
from lxml import etree
import numpy
from point import Point
import math
import svgwrite

from svg.path import Path, Line, Arc, CubicBezier, QuadraticBezier, parse_path
from curve_components import *

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
	components = []
	for p in path:
		if type(p) is CubicBezier:
			print "CUBEY: " + str(p) + " -!--! " + str(type(p.start))
			components.append(MyBezCurve(p))
		elif type(p) is Line:
			components.append(MyLine(p))
			print "LINE: " + str(p)

	#The desired length of each segment of the output.
	seg_length = 20;
	#Iterate over the components in the svg, finding points where necesary.

	#Paramaterize the curve
	for curve in components:
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
	

