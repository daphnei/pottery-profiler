__author__ = 'daphne'
import sys
import numpy
from point import Point
import math
import pickle
import cmath
from settings import *

def check_one_svg(target_name, desc_type=RIGHT_CURVATURE_KEY):
	'''Compares one svg to all the other svgs in the directory. Prints out match scores.
	   I tried to base this as much as possible on what was done in the old crane script.
	'''
	print "Comparing " + target_name

	# with open("fft_data.json", "r") as pickle_file:
	# 	all_data_string = pickle_file.read()
	# 	all_data = jsonpickle.decode(all_data_string)
	all_data = pickle.load( open( "sherd_data.pickle", "rb" ) )

	#Find the one we are comparing against.
	target_obj = all_data[target_name]

	target_descriptors = numpy.asarray(target_obj[desc_type])

	#for each sherd in the database, store its average distance from the target
	dists = {}

	for shape_name, shape_data in all_data.items():
		# gets euclid distance for 2 sets of points. Only consider the first n descriptors where n is the minimum
		# length of either list
		tomatch_descriptors = numpy.asarray(shape_data[desc_type])
		n = min(tomatch_descriptors.size, target_descriptors.size)

		dists[shape_name] = numpy.linalg.norm(target_descriptors[:n] - tomatch_descriptors[:n])

	#Sort the items by how close they are to the target one.
	sorted_dists = sorted(dists.items(), key=lambda dists: dists[1])

	print("{:25}\t{:10}".format("name", "FD"))
	for i in xrange(len(sorted_dists)):
		print("{:25}\t{:.6f}".format(sorted_dists[i][0], sorted_dists[i][1]))

if __name__ == "__main__":
	check_one_svg(sys.argv[1]);

	# path = get_path_from_svg("/Users/daphne/Documents/School/CSC494/pottery-profiler/Pottery/" + sys.argv[1])
	# points = get_points_along_path(path)
	# left_profile_points, right_profile_points = split_profile_points(points)
	# draw_points_to_output_file(left_profile_points, right_profile_points)
