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
	except:
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

	#One path consists of multiple components. These components can either be Bezier curves
	#or straight lines. We want to choose a point every "seg_length" units along the entire
	#path, which means we can't necessarily start the points at the very beginning of each
	#component.
	leftovers = 0
	old_comp_length = 0
	for comp in components:
		#Gets the arc length of the component
		comp_length = comp.get_length()

		print "\nWorking on component with length " + str(comp_length)

		if comp_length != 0:
			d = leftovers
			goal_lengths = []
			while d <= comp_length:

				goal_lengths.append(d)
				d += seg_length
			print "The comp_length is " + str(comp_length) + " and d went up to " + str(d)

			leftovers = d - comp_length
			print "The leftovers are " + str(leftovers)

			fractions = map(lambda x: x / comp_length, goal_lengths)
			print fractions
			interped_points += comp.interpolate_points(fractions)

	return interped_points

def find_components_on_side(components, clockwise, index_of_highest, index_of_lowest):
	components_on_side = []

	i = index_of_highest;
	while True:
		i = i + 1 if clockwise else i - 1
		if i == len(components):
			i = 0
		elif i == -1:
			i = len(components) - 1

		if i == index_of_lowest:
			break

		components_on_side.append(components[i])

	return components_on_side

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

	#figure out which direction the curve was digitized in (clockwise or counterclockwise)
	#Also, decide whether the topmost component should be in the left or the right profile.
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

	print "This path is " + str("CLOCKWISE" if clockwise else "COUNTER-CLOCKWISE")

	#Find all the components that belong to the right curve.
	right_profile += find_components_on_side(components, clockwise, index_of_highest, index_of_lowest)
	left_profile += find_components_on_side(components, not clockwise, index_of_highest, index_of_lowest)

	#Now decide if the bottommost component should be in the left or the right profile.
	bottom_comp = components[index_of_lowest]
	start = bottom_comp.get_start_point()
	end = bottom_comp.get_end_point()
	if (start.x < end.x):
		if (start.y > end.y):
			right_profile.append(bottom_comp)
		else:
			left_profile.append(bottom_comp)
	else:
		if (start.y > end.y):
			left_profile.append(bottom_comp)
		else:
			right_profile.append(bottom_comp)

	#If the curve was digitized in a clockwise direction than the order of the paths
	#in the left profile needs to be reversed so that the 0th interpolated point
	#is at the very top of the profile, not the very bottom.
	if clockwise:
		for c in left_profile:
			c.reverse()

	#Likewise, the paths in the right profile need to be reverse if the digitization order
	#was counter-clockwise.
	if not clockwise:
		for c in right_profile:
			c.reverse()

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
		output_svg.add(output_svg.circle(center=(p.x, p.y), r=0.7, fill=svgwrite.rgb(100, 0, 0, '%')))
	for p in interped_points_outer:
		output_svg.add(output_svg.circle(center=(p.x, p.y), r=0.7, fill=svgwrite.rgb(0, 0, 100, '%')))
	output_svg.save()
