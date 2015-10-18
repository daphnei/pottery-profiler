import xml.etree.ElementTree as ET
import sys
from lxml import etree
import numpy
import svgwrite
import os
import pickle
import math

from svg.path import Path, Line, Arc, CubicBezier, QuadraticBezier, parse_path
from settings import *
from curve_components import *


def get_path_from_svg(svg_file):
	'''Given a file path, reads in the SVG at that path.'''

	print "Reading in the file: " + str(svg_file)

	try:
		tree = etree.parse(svg_file)
	except:
		print "Unable to parse the file. Are you sure it exists, and it is an SVG?"
		exit(1)

	# Find all paths in the SVG
	all_paths = tree.findall("{*}path")

	if len(all_paths) == 0:
		print "WARNING: The SVG is empty."
		return []
	if len(all_paths) > 1:
		print "WARNING: The SVG contains multiple paths. Only the first one will be considered."

	path_node = all_paths[0]
	path_string = path_node.get("d")

	# The path node in the SVG contains the string storing the path info, but this is not easy to
	# work with directly. Therefore parse this string into a path datastructure.
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

	# Iterate over the components in the svg, finding equidistant points along the curve.
	total_length = sum(c.get_length() for c in components)

	interped_points = []

	# Find the component that contains that contains the height y value, and start processing
	# that component first. doing this only helps to ensure that the top most curve has a
	#SCRAP THIS FOR NOW

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
	points_x = list(p.x for p in points)
	# Find the index of the topmost point
	top_index = points_y.index(min(points_y))

	# Find the index of the bottommost point.
	bottom_index = points_y.index(max(points_y))
	starts_on_left = points_x[0] < points_x[len(points_x) - 1]

	if (bottom_index < top_index):
		pr1 = points[top_index:] + points[0:bottom_index + 1]
		pr2 = points[bottom_index:top_index + 1]
	elif (top_index < bottom_index):
		pr1 = points[top_index:bottom_index + 1]
		pr2 = points[bottom_index:] + points[0:top_index + 1]

	else:
		print "TODO: deal with paththis"

	# Ensure that the 0th element of each profile is the element at the very top
	pr2.reverse()

	#Figure out which profile is the left and which is the right by checking the
	#second element of each list.
	if pr1[1].x < pr2[1].x:
		left_profile = pr1
		right_profile = pr2
	else:
		left_profile = pr2
		right_profile = pr1

	return left_profile, right_profile


def draw_points_to_output_file(left_profile_points, right_profile_points, output_file_name='output.svg'):
	'''	write an output svg with a circle marker at each chosen point position '''

	output_svg = svgwrite.Drawing(output_file_name)
	for p in left_profile_points:
		output_svg.add(output_svg.circle(center=(p.x, p.y), r=0.7, fill=svgwrite.rgb(100, 0, 0, '%')))
	for p in right_profile_points:
		output_svg.add(output_svg.circle(center=(p.x, p.y), r=0.7, fill=svgwrite.rgb(0, 0, 100, '%')))
	output_svg.save()


def compute_fft_points(points):
	'''
	Compute a fourier descriptor for each of the input points.
	'''

	coords = numpy.array(points).transpose()

	# Position the curve so that the top most point is at the origin (0,0)
	x = list(p.x - coords[0].x for p in coords)
	y = list(p.y - coords[0].y for p in coords)

	# z = x + jy
	z = [complex(x[i], y[i]) for i in range(len(x))]

	# Fourier descriptors
	z_f = numpy.fft.fft(numpy.asarray(z))

	# scale invariance 
	z_f1 = list([m / abs(z_f[1]) for m in z_f])

	return z_f1


def dfourier(x, T):
	#TRANSLATED OVER FROM THE MATLAB CODE

	l = 2 * math.pi / T;
	t = len(x);
	f = numpy.fft.fft(x);
	a1 = range(0, int(numpy.fix((t - 1) / 2) + 1))
	a2 = range(-int(numpy.fix(t / 2)), -1 + 1)
	a = numpy.concatenate((a1, a2), 0)
	b = numpy.array([0j] * t)
	for k in xrange(0, t):
		b[k] = (1j * l) * (f[k] * a[k])

	d = numpy.real(numpy.fft.ifft(b))
	return d


def compute_first_derivative(x):
	# TRANSLATED OVER FROM THE MATLAB CODE
	# This function should calculate the first derivative of a function x using
	# Fourier transform.
	# x doesn't have to be periodic.

	t = len(x)
	tt = numpy.array(list(range(1, t + 1)))

	if x[len(x) - 1] == x[0]:
		mean_change = numpy.zeros([1, t])
		dp = numpy.zeros([1, t])
		p = dp
	else:
		q = (x[1] - x[0]) - (x[len(x) - 1] - x[len(x) - 2])
		a = float(q) / (2 * (1 - t))
		b = (x[0] - x[len(x) - 1]) / float(1 - t) - a * (t + 1)
		c1 = x[0] - a - b
		c2 = x[len(x) - 1] - a * t * t - b * t
		c = (c1 + c2) / 2.0
		p = a * tt * tt + b * tt + c
		mean_change = p
		dp = 2 * a * tt + b

	X = x - mean_change
	tt = int(numpy.fix(t / 2))
	XX = numpy.concatenate((X[tt: len(X) - 2], X, X[1:tt]), 0)
	DXX = dfourier(XX, len(XX))
	DX = numpy.array(DXX[numpy.fix((t + 1) / 2): numpy.fix((t + 1) / 2) + len(X)])
	d1 = dfourier(X[0:len(X) - 1], t - 1)
	d2 = dfourier(X[1:], t - 1)
	DX = numpy.concatenate(([d1[0]], d2), 0)
	dx = DX + dp

	return dx


def compute_tangent(points):
	#TRANSLATED OVER FROM THE MATLAB CODE
	#I don't really understand exactly what this is doing, but I think it is supposed to get
	#the tangent along the curve. I copied this over from curve_tangent.m in the Matlab code.

	x = list(p.x for p in points)
	y = list(p.y for p in points)

	x = x - numpy.mean(x)
	y = y - numpy.mean(y)

	X = numpy.array(x)
	Y = numpy.array(y)

	t1 = numpy.angle(x[0] + (y[0] * 1j), False)
	M = numpy.array([[numpy.cos(t1), numpy.sin(t1)], [-numpy.sin(t1), numpy.cos(t1)]]);
	for k in xrange(len(x)):
		Q = numpy.dot(M, numpy.array([[x[k]], [y[k]]]))
		X[k] = Q[0]
		Y[k] = Q[1]

	dx = compute_first_derivative(X);
	dy = compute_first_derivative(Y);

	d = numpy.divide(dy, dx)
	teta = list(math.atan(x) for x in d)

	I = numpy.multiply(teta[:len(teta) - 1], teta[1:len(teta)])
	I = numpy.where(I < -1.2)[0]
	if len(I) > 0:
		k = len(I) - 1
		while k >= 0:
			II = I[k] + 1;
			if teta[II] < 0:
				teta[II:] = list((teta[x] + math.pi) for x in range(II, len(teta)))
			else:
				teta[II:] = list((teta[x] - math.pi) for x in range(II, len(teta)))

			k -= 1

	circle_step = (2 * math.pi) / len(x)
	circle = (math.pi / 2) + numpy.arange(0, (2 * math.pi) * (1 - 1.0 / len(x)) + 0.001, circle_step)
	teta = teta - circle
	if numpy.mean(teta) < -math.pi / 2:
		teta = teta + math.pi;

	return teta


def compute_tangent_alt(points):
	"""
	Computes the slope of the curve at each point. For vertical slopes, NaN is returned.
	:param points: List of points forming the curve.
	:return:
	"""
	tangent = [0] * len(points)

	for i in xrange(1, len(points) - 1):
		p1 = points[i - 1]
		p2 = points[i + 1]

		if p2.x - p1.x == 0:
			slope = float('NaN')
		else:
			slope = (p2.y - p1.y) / float(p2.x - p1.x)

		tangent[i] = slope

	return tangent


def compute_thickness(target_profile, other_profile, tangents_target, tangents_other):
	'''

	An attempt to compute the thickness of the sherd at each point in the target_profile.
	It doesn't work very well.

	TODO: Improve this function by making use of the medial axis.

	:param target_profile: For each point in this profile, the width will be calculated.
	:param other_profile: The opposite profile.
	:param tangents_target: The slope of the tangent to the curve at each of the points in target_profile.
	:param tangents_other: The slope of the tangent to the curve at each of the points in other_profile.
	:return:
	'''
	thicknesses = [0] * len(target_profile)

	output_svg = svgwrite.Drawing("temp.svg")

	for i in xrange(1, len(target_profile) - 1):
		try:
			# Find the intersection between the line normal to the tangent, and each of the line segments making up the
			#opposite profile.

			p1 = target_profile[i - 1]
			p2 = target_profile[i + 1]
			p_on_line = Point((p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2.0)

			#The slope and y-intercept of the line normal to the tangent.
			s1 = - 1.0 / tangents_target[i]
			b1 = target_profile[i].y - s1 * target_profile[i].x

			# TO REMOVE
			s0 = tangents_target[i]
			b0 = target_profile[i].y - s0 * target_profile[i].x

			min_dist = sys.maxint
			best_point = None

			#TODO: Do something smarter here. If the database grows too large, this will have crappy runtime.
			#Maybe restrict to only look at a subset of points that are nearby?
			for j in xrange(1, len(other_profile) - 1):
				p1 = other_profile[j - 1]
				p2 = other_profile[j + 1]
				p_on_line = Point((p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2.0)

				s2 = tangents_other[j]
				b2 = p_on_line.y - s2 * p_on_line.x

				#determine the x-y position of the intercept between the two lines
				x = (b2 - b1) / (s1 - s2)
				y = s1 * x + b1

				#Is this point closer than any others found previously?
				intercept = Point(x, y)
				distance = intercept.distance(target_profile[i])
				if distance < min_dist and intercept.is_between(p1, p2):
					min_dist = distance
					best_point = other_profile[j]

			thicknesses[i] = min_dist


			#Visualize the results to see how good they are.
			#TODO: Remove this.
			#output_svg.add(output_svg.line((0, b0), (target_profile[i].x, target_profile[i].y), stroke=svgwrite.rgb(80, 0, 160, '%')))
			#output_svg.add(output_svg.line((0, b1), (target_profile[i].x, target_profile[i].y), stroke=svgwrite.rgb(10, 10, 16, '%')))
			if (best_point != None):
				output_svg.add(output_svg.line((best_point.x, best_point.y), (target_profile[i].x, target_profile[i].y),
											   stroke=svgwrite.rgb(10, 10, 16, '%')))
		except ZeroDivisionError:
			thicknesses[i] = float('NaN')

	#Save the visualization.
	output_svg.save()
	return thicknesses


def compute_curvature(points):
	'''
	:param points:
	:return: curvatures for each inner point. A smaller radius of curvaturs means that the the section is more bendy.
			A larger value means that the line section is closer to linear.

			Math for this is based on:
			http://www.intmath.com/applications-differentiation/8-radius-curvature.php
	'''

	curvature = [0] * len(points)

	# It is only possible to compute the curvature at a point if there are points
	# both to the left and to the right of it. Therefore, it is not possible to compute curvature
	# for the first and last points. Leave these at 0
	for i in xrange(1, len(points) - 1):
		a = points[i - 1]
		b = points[i]
		c = points[i + 1]

		try:
			# find the slope between points a and b
			m1 = ((a.y - b.y) / (float)(a.x - b.x))
			# find the slope between points b and c
			m2 = ((b.y - c.y) / (float)(b.x - c.x))
			# find the slope between points a and c
			m3 = ((a.y - c.y) / (float)(a.x - c.x))

			# Find the center of the circle passing through these three points.
			center_x = ((m1 * m2 * (a.x - c.x)) + (m2 * (a.x + b.x)) - (m1 * (b.x + c.x))) / (2.0 * (m2 - m1))
			center_y = ((-1 / m1) * (center_x - ((a.x + b.x) / 2.0))) + (a.y + b.y) / 2.0

			# Find the distance between the center of the circle and any of the three points. This is the radius
			# of curvature.
			curvature[i] = a.distance(Point(center_x, center_y))

		# Curvature should also encode the direction of the curve. The direction can be gotten from the slope
		# between points a and c.
		# Now that we also have a descriptor for the slope of the tangent, this becomes less important.
		# It can be commented out.
		# curvature[i] *= (1 if m3 >= 0 else -1)

		except ZeroDivisionError:
			# The curve around this point is linear. This means that the curvature radius is 0,
			curvature[i] = 0;
	return curvature

def save_descriptors_for_al_svgs(dir):
	""" Goes through every .svg in the input directory, paramaterizes its curve into a series of points,
	    and then computed a bunch of different metrics on those points.

	    The results are stored as a dictionary (keyed by file name) of dictionaries (keyed by Metric name), pointing
	    to a list of values, one for each generated point. Since computation is pretty slow, this dictionary is written
	    out to a pick then can be imported by other scripts.

	    When applicable, the various metrics are computed separately for the left and right profiles of the pottery sherd.

	    :param dir: The directory in which to search for SVG files.

	"""

	all_data = {}

	for filename in os.listdir(dir):
		if filename.endswith(".svg"):
			path = get_path_from_svg(dir + filename)

			data = {}
			data[Metric.LEFT_FFT_KEY] = []
			data[Metric.RIGHT_FFT_KEY] = []
			data[Metric.LEFT_CURVATURE_KEY] = []
			data[Metric.RIGHT_CURVATURE_KEY] = []
			data[Metric.LEFT_TANGENT_KEY] = []
			data[Metric.RIGHT_TANGENT_KEY] = []
			data[Metric.LEFT_TANGENT_ALT_KEY] = []
			data[Metric.RIGHT_TANGENT_ALT_KEY] = []
			data[Metric.THICKNESS_KEY] = []

			all_data[filename] = data

			# If the svg didn't contain any path, then skip doing any calculation on it.
			if path:
				points = get_points_along_path(path)
				left_profile_points, right_profile_points = split_profile_points(points)

				draw_points_to_output_file(left_profile_points, right_profile_points, "out/output_" + filename);

				# Calculature the fourier transforms for each profile
				data[Metric.LEFT_FFT_KEY] = compute_fft_points(left_profile_points)
				data[Metric.RIGHT_FFT_KEY] = compute_fft_points(right_profile_points)

				# calculate the curvature along each profile.
				data[Metric.LEFT_CURVATURE_KEY] = compute_curvature(left_profile_points)
				data[Metric.RIGHT_CURVATURE_KEY] = compute_curvature(right_profile_points)

				data[Metric.LEFT_TANGENT_ALT_KEY] = compute_tangent_alt(left_profile_points)
				data[Metric.RIGHT_TANGENT_ALT_KEY] = compute_tangent_alt(right_profile_points)

				#calculate the direction of the tangent to the curbe along each profile.
				data[Metric.LEFT_TANGENT_KEY] = compute_tangent(left_profile_points)
				data[Metric.RIGHT_TANGENT_KEY] = compute_tangent(right_profile_points)

				#also keep trick of the point location.
				data[Metric.X_KEY] = list(p.x for p in points)
				data[Metric.Y_KEY] = list(p.y for p in points)

				data[Metric.THICKNESS_KEY] = compute_thickness(left_profile_points, right_profile_points, \
															   data[Metric.LEFT_TANGENT_ALT_KEY],
															   data[Metric.RIGHT_TANGENT_ALT_KEY])

	pickle.dump(all_data, open(DESC_OUTPUT_FILE, "wb"))


if __name__ == "__main__":
	if len(sys.argv) != 2:
		print "USAGE: python find_pottery_descriptors.py <directory>"
		exit(1)

	save_descriptors_for_al_svgs(sys.argv[1])

	# path = get_path_from_svg("/Users/daphne/Documents/School/CSC494/pottery-profiler/Pottery/AS_99B_2012_1.svg")
	# points = get_points_along_path(path)
	# left_profile_points, right_profile_points = split_profile_points(points)
	# draw_points_to_output_file(left_profile_points, right_profile_points)
	#
	# left_tan = compute_tangent_alt(left_profile_points)
	# right_tan = compute_tangent_alt(right_profile_points)
	# thickness = compute_thickness(left_profile_points, right_profile_points, left_tan, right_tan)
	# print thickness
	# print right_tan

	print "The pottery descriptors have been written to the pickle: " + DESC_OUTPUT_FILE
