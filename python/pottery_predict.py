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
	'''Given a file path, reads in the SVG at that path.'''

	print "Reading in the file: " + str(svg_file)

	try:
		tree = etree.parse(svg_file)
	except Error:
		print "Unable to parse the file. Are you sure it is an SVG?" 
		exit(1)

	#Find all paths in the SVG
	all_paths = tree.findall("{*}path")

	if len(all_paths) == 0:
		return None
	if len(all_paths) > 1:
		print "WARNING: The input SVG contains multiple paths. Only the first one will be considered."

	path_node = all_paths[0]
	path_string = path_node.get("d")

	#The path node in the SVG contains the string storing the path info, but this is not easy to
	#work with directly. Therefore parse this string into a path datastructure.
	path = parse_path(path_string)

	return path

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print "USAGE: python pottery_predict.py <svg file>"
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

	#The desired distance between points along the curve
	seg_length = 5;

	#Iterate over the components in the svg, finding equidistant points along the curve.
	total_length = sum(c.get_length() for c in components)	

	interped_points = []

	cur_length = 0
	i = 0
	leftovers = 0
	for comp in components:
		cur_comp_length = comp.get_length()
		d = leftovers
		fractions = []
		while d < cur_comp_length:
			fractions.append(d / float(cur_comp_length))
			d += seg_length
		leftovers = d - cur_comp_length

		interped_points += comp.interpolate_points(fractions)


	#write an output svg with a circle marker at each chosen point position
	output_svg = svgwrite.Drawing('output.svg')
	for p in interped_points:
		output_svg.add(output_svg.circle(center=(p.x, p.y), r=1.2, fill=svgwrite.rgb(100, 0, 0, '%')))
	output_svg.save()


