import xml.etree.ElementTree as ET
import sys
from lxml import etree
import numpy
from point import Point
import math
import svgwrite
import jsonpickle
import os
import pickle
import cmath

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

		if comp_length != 0:
			d = leftovers
			goal_lengths = []
			while d <= comp_length:
				goal_lengths.append(d)
				d += seg_length

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
		right_profile = points[top_index:] + points[0:bottom_index+1]
		left_profile = points[bottom_index:top_index+1]
	elif (top_index < bottom_index):
		right_profile = points[top_index:bottom_index+1]
		left_profile = points[bottom_index:] + points[0:top_index+1]
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
	output_svg.save(left_profile_points, right_profile_points)

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

def save_fft_for_all_svgs(dir):
	fft_data = {}

	''' Goes through every .svg in the input directory, and saves its fourier transform'''
	for filename in os.listdir(dir):
		if filename.endswith(".svg"):
			path = get_path_from_svg(dir + filename)

			#If the svg didn't contain any path, then skip doing any calculation on it.
			if path == []:
				data = {}
				data["left_fft"] = []
				data["right_fft"] = []

				fft_data[filename] = data
			else:
				points = get_points_along_path(path)
				left_profile_points, right_profile_points = split_profile_points(points)

				data = {}
				data["left_fft"] = compute_fft_points(left_profile_points)
				data["right_fft"] = compute_fft_points(right_profile_points)

				fft_data[filename] = data

	#pickle_string = jsonpickle.encode(fft_data)
	#with open("fft_data.json", "w") as pickle_file:
	#	pickle_file.write(pickle_string)

	pickle.dump(fft_data,  open("fft_data.json", "wb"))

def check_one_svg(target_name, profile_side="right"):
	'''Compares one svg to all the other svgs in the directory. Prints out match scores.
	   I tried to base this as much as possible on what was done in the old crane script.
	'''
	print "Comparing " + target_name

	profile_side += "_fft"

	# with open("fft_data.json", "r") as pickle_file:
	# 	all_data_string = pickle_file.read()
	# 	all_data = jsonpickle.decode(all_data_string)
	all_data = pickle.load( open( "fft_data.json", "rb" ) )

	#Find the one we are comparing against.
	target_obj = all_data[target_name]

	target_descriptors = numpy.asarray(target_obj[profile_side])
	phase = numpy.angle(target_descriptors)
	target_len = len(target_obj[profile_side])

	fd_dist = {} #store the different in magnitudes
	phase_dist = {} #store the difference in phases

	for shape_name, shape_data in all_data.items():
		# gets euclid distance for 2 sets of points. Only consider the first n descriptors where n is the minimum
		# length of either list
		descriptors = numpy.asarray(shape_data[profile_side])
		n = min(descriptors.size, target_descriptors.size)

		fd_dist[shape_name] = (cmath.sqrt(sum((descriptors[:n] - target_descriptors[:n])**2))).real
		phase_dist[shape_name] = ((cmath.sqrt(sum((numpy.angle(descriptors[:n])-phase[:n])**2))).real)

	#Sort the items by how close they are to the target one.
	sorted_fd_dist = sorted(fd_dist.items(), key=lambda fd_dist: fd_dist[1])

	print("{:25}\t{:10}\t{:3}".format("name", "FD", "phase"))
	for i in xrange(len(sorted_fd_dist)):
		print("{:25}\t{:.6f}\t{:.3f}".format(sorted_fd_dist[i][0], sorted_fd_dist[i][1], phase_dist[sorted_fd_dist[i][0]]))


if __name__ == "__main__":
	if len(sys.argv) != 2:
		print "USAGE: python pottery_predict.py <svg file>"
		exit(1)

	# save_fft_for_all_svgs(sys.argv[1])

	print("Done saving")
	check_one_svg("AS_4A_2012_1.svg");
	# path = get_path_from_svg(sys.argv[1])

	# points = get_points_along_path(path)
	# left_profile_points, right_profile_points = split_profile_points(points)

	# l = compute_fft_points(left_profile_points)
	# r = compute_fft_points(right_profile_points)

	# draw_points_to_output_file(l, r)
