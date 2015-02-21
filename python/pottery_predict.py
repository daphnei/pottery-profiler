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

def get_points_along_path(components, seg_length = 5):
	''' Inputs are
		components list of components forming the path
		seg_length: The desired distance between points along the curve
	'''

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

	return interped_points

def split_curve(components):
	#Find the component with the topmost start/end point
	index_of_highest = 0
	index_of_lowest = 0
	highest_y = sys.maxint #The y value closest to 0 (the top of the page)
	lowest_y = -sys.maxint #The largest y (furthest from the top of the page)
	for i in xrange(0, len(components)):
		start_point = components[i].get_start_point()
		end_point = components[i].get_end_point()

		if min(start_point.y, end_point.y) < highest_y:
			index_of_highest = i
			highest_y = min(start_point.y, end_point.y)

		if max(start_point.y, end_point.y) > lowest_y:
			index_of_lowest = i
			lowest_y = max(start_point.y, end_point.y)

	left_profile = [] #The inner side of the profile.
	right_profile = [] #The outside side of the profile.

	#Figure out whether this component should be considered part of the
	#inner profile or outer profile, and figure out which direction
	#the curve was digitized in (clockwise or counterclockwise)
	top_comp = components[index_of_highest]
	start = top_comp.get_start_point()
	end = top_comp.get_end_point()
	if (start.x < end.x):
		clockwise = True
		if (start.y > end.y):
			left_profile.append(top_comp)
		else:
			right_profile.append(top_comp)
	else:
		clockwise = False
		if (start.y > end.y):
			right_profile.append(top_comp)
		else:
			left_profile.append(top_comp)

	#Find all the components that belong to the right curve.
	i = index_of_highest;
	while i != index_of_lowest:
		i = i + 1 if clockwise else i - 1
		if i == len(components):
			i = 0
		elif i == -1:
			i = len(components) - 1

		right_profile.append(components[i])

	#Find all the components that belong to the left curve.
	i = index_of_highest;
	while i != index_of_lowest:
		i = i - 1 if clockwise else i + 1
		if i == len(components):
			i = 0
		elif i == -1:
			i = len(components) - 1

		left_profile.append(components[i])

	return left_profile, right_profile


if __name__ == "__main__":
	if len(sys.argv) != 2:
		print "USAGE: python pottery_predict.py <svg file>"
		exit(1)

	path = get_path_from_svg(sys.argv[1])

	#Accumulate all of the lines and cubic bezier curves in the svg.
	components = []
	for p in path:
		print p
		if type(p) is CubicBezier:
			#print "CUBIC BEZ: " + str(p) + " -!--! " + str(type(p.start))
			components.append(MyBezCurve(p))
		elif type(p) is Line:
			components.append(MyLine(p))
			#print "LINE: " + str(p)

	#Split the components into those that are part of the left edge
	#and those that are part of the right edge.
	inner_profile, outer_profile = split_curve(components)

	interped_points_inner = get_points_along_path(inner_profile)
	interped_points_outer = get_points_along_path(outer_profile)


	#write an output svg with a circle marker at each chosen point position
	output_svg = svgwrite.Drawing('output.svg')
	for p in interped_points_inner:
		output_svg.add(output_svg.circle(center=(p.x, p.y), r=1.2, fill=svgwrite.rgb(100, 0, 0, '%')))
	for p in interped_points_outer:
		output_svg.add(output_svg.circle(center=(p.x, p.y), r=1.2, fill=svgwrite.rgb(0, 0, 100, '%')))
	output_svg.save()


