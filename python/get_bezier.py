import xml.etree.ElementTree as ET
import sys
from lxml import etree
import numpy
from param import Length_iterator

from svg.path import Path, Line, Arc, CubicBezier, QuadraticBezier, parse_path

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
	//li = Length_iterator(4, 0.01)
	//li.initializeIterationOnCurve(path_points, 8)

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print "Input svg path"
		exit(1)

	path = get_path_from_svg(sys.argv[1])
	print("INFO")
	print(dir(path))
	#get_equidistant_points(path, 20)

	print("")

	for p in path:
		if type(p) is CubicBezier:
			print "CUBEY: " + str(p) + " -!--! " + str(type(p.start))
			points.append(p.start)
			points.append(p.control1)
			points.append(p.control2)
			points.append(p.end)
		else:
			print "NOPE"

	print "\nIn the end there were points"
	print points
	get_equidistant_points(points)

