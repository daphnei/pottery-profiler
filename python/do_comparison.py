__author__ = 'daphne'
import sys
import numpy
from point import Point
import math
import pickle
import cmath

def check_one_svg(target_name, profile_side="right"):
	'''Compares one svg to all the other svgs in the directory. Prints out match scores.
	   I tried to base this as much as possible on what was done in the old crane script.
	'''
	print "Comparing " + target_name

	profile_side += "_curvature"

	# with open("fft_data.json", "r") as pickle_file:
	# 	all_data_string = pickle_file.read()
	# 	all_data = jsonpickle.decode(all_data_string)
	all_data = pickle.load( open( "sherd_data.pickle", "rb" ) )

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
	check_one_svg(sys.argv[1]);

	# path = get_path_from_svg("/Users/daphne/Documents/School/CSC494/pottery-profiler/Pottery/" + sys.argv[1])
	# points = get_points_along_path(path)
	# left_profile_points, right_profile_points = split_profile_points(points)
	# draw_points_to_output_file(left_profile_points, right_profile_points)
