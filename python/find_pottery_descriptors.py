import xml.etree.ElementTree as ET
import sys
from lxml import etree
import numpy
import svgwrite
import os
import pickle

from svg.path import Path, Line, Arc, CubicBezier, QuadraticBezier, parse_path
from settings import *
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
		print "WARNING: The SVG is empty."
		return []
	if len(all_paths) > 1:
		print "WARNING: The SVG contains multiple paths. Only the first one will be considered."

	path_node = all_paths[0]
	path_string = path_node.get("d")

	#The path node in the SVG contains the string storing the path info, but this is not easy to
	#work with directly. Therefore parse this string into a path datastructure.
	path = parse_path(path_string)

	#we only care about lines and Bezier curves. Create a list of just these two types, and also
	#wrap each of them in an object that provides some helpful utility methods.
	components = []
	for p in path:
		if type(p) is CubicBezier:
			components.append(MyBezCurve(p))
		elif type(p) is Line:
			components.append(MyLine(p))

	return components

def get_points_along_path(components):
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

		if comp_length != 0:
			d = leftovers
			goal_lengths = []
			while d <= comp_length:
				goal_lengths.append(d)
				d += SEG_LENGTH

			leftovers = d - comp_length

			fractions = map(lambda x: x / comp_length, goal_lengths)
			interped_points += comp.interpolate_points(fractions)

	return interped_points

def split_profile_points(points):
	points_y = list(p.y for p in points)

	#Find the index of the topmost point
	top_index = points_y.index(min(points_y))

	#Find the index of the bottommost point.
	bottom_index = points_y.index(max(points_y))

	if (bottom_index < top_index):
		left_profile = points[top_index:] + points[0:bottom_index+1]
		right_profile = points[bottom_index:top_index+1]
	elif (top_index < bottom_index):
		left_profile = points[top_index:bottom_index+1]
		right_profile = points[bottom_index:] + points[0:top_index+1]
	else:
		print "TODO: deal with this"

	#Ensure that the 0th element of each profile is the element at the very top
	left_profile.reverse()

	return left_profile, right_profile

def draw_points_to_output_file(left_profile_points, right_profile_points):
	'''	write an output svg with a circle marker at each chosen point position '''

	output_svg = svgwrite.Drawing('output.svg')
	for p in left_profile_points:
		output_svg.add(output_svg.circle(center=(p.x, p.y), r=0.7, fill=svgwrite.rgb(100, 0, 0, '%')))
	for p in right_profile_points:
		output_svg.add(output_svg.circle(center=(p.x, p.y), r=0.7, fill=svgwrite.rgb(0, 0, 100, '%')))
	output_svg.save()

def compute_fft_points(points):
	#We only care about the x values since
	#points_x = p.x for p in points

	coords = numpy.array(points).transpose()

	# Position the curve so that the top most point is at the origin (0,0)
	x = list(p.x - coords[0].x for p in coords)
	y = list(p.y - coords[0].y for p in coords)

	# z = x + jy
	z = [complex(x[i],y[i]) for i in range(len(x))]

	# Fourier descriptors
	z_f = numpy.fft.fft(numpy.asarray(z))

	# scale invariance 
	z_f1 = [m/abs(z_f[1]) for m in z_f]

	return z_f1

def compute_curvature(points):
	'''
	:param points:
	:return: curvatures for each inner point. A smaller radius of curvaturs means that the the section is more bendy.
			A larger value means that the line section is closer to linear.

			Math for this is based on:
			http://www.intmath.com/applications-differentiation/8-radius-curvature.php
	'''

	curvature = [0] * len(points)

	#It is only possible to compute the curvature at a point if there are points
	#both to the left and to the right of it. Therefore, it is not possible to compute curvature
	#for the first and last points. Leave these at 0
	for i in xrange(1, len(points) - 1):
		a = points[i-1]
		b = points[i]
		c = points[i+1]

		try:
			#find the slope between points a and b
			m1 = ((a.y - b.y) / (float) (a.x - b.x))
			#find the slope between points b and c
			m2 = ((b.y - c.y) / (float) (b.x - c.x))
			#find the slope between points a and c
			m3 = ((a.y - c.y) / (float) (a.x - c.x))

			#Find the center of the circle passing through these three points.
			center_x = ((m1 * m2 * (a.x - c.x)) + (m2 * (a.x + b.x)) - (m1 * (b.x + c.x))) / (2.0 * (m2 - m1))
			center_y = ((-1 / m1) * (center_x - ((a.x + b.x) / 2.0))) + (a.y + b.y) / 2.0

			#Find the distance between the center of the circle and any of the three points. This is the radius
			#of curvature.
			curvature[i] = a.distance(Point(center_x, center_y))

			#Curvature should also encode the direction of the curve. The direction can be gotten from the slope
			#between points a and c
			curvature[i] *= (1 if m3 >= 0 else -1)
		except ZeroDivisionError:
			#The curve around this point is linear. This means that the curvature radius is technically infinity,
			#but I'll instead just use a very large number.
			curvature[i] = sys.maxint;
	return curvature

def save_fft_for_all_svgs(dir):
	fft_data = {}

	''' Goes through every .svg in the input directory, and saves its fourier transform'''
	for filename in os.listdir(dir):
		if filename.endswith(".svg"):
			path = get_path_from_svg(dir + filename)

			#If the svg didn't contain any path, then skip doing any calculation on it.
			if path == []:
				data = {}
				data[LEFT_FFT_KEY] = []
				data[RIGHT_FFT_KEY] = []
				data[LEFT_CURVATURE_KEY] = []
				data[RIGHT_CURVATURE_KEY] = []

				fft_data[filename] = data
			else:
				points = get_points_along_path(path)
				left_profile_points, right_profile_points = split_profile_points(points)

				data = {}
				#Calculature the fourier transforms for each profile
				data[LEFT_FFT_KEY] = compute_fft_points(left_profile_points)
				data[RIGHT_FFT_KEY] = compute_fft_points(right_profile_points)

				#calculate the curvature along each profile.
				data[LEFT_CURVATURE_KEY] = compute_curvature(left_profile_points)
				data[RIGHT_CURVATURE_KEY] = compute_curvature(left_profile_points)

				fft_data[filename] = data

	#pickle_string = jsonpickle.encode(fft_data)
	#with open("fft_data.json", "w") as pickle_file:
	#	pickle_file.write(pickle_string)

	pickle.dump(fft_data,  open(DESC_OUTPUT_FILE, "wb"))

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print "USAGE: python find_pottery_descriptors.py <svg file>"
		exit(1)

	save_fft_for_all_svgs(sys.argv[1])

	path = get_path_from_svg("/Users/daphne/Documents/School/CSC494/pottery-profiler/TestPottery/test_curve_07.svg")
	points = get_points_along_path(path)
	left_profile_points, right_profile_points = split_profile_points(points)
	draw_points_to_output_file(left_profile_points, right_profile_points)

	print "The pottery descriptors have been written to the pickle: " + DESC_OUTPUT_FILE
