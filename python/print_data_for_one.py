'''
This script takes the name of one of the sherds stored in the database,
and prints out all of the descriptors for that sherd.
'''

__author__ = 'daphne'

import pickle
import sys
from settings import *

if __name__ == "__main__":
	'''
		Prints the metrics computed for the specified pottery sherd.
	'''

	if len(sys.argv) != 2:
		print "USAGE: python print_data_for_one.py <name>"
		exit(1)

	all_data = pickle.load(open(DESC_OUTPUT_FILE, "rb" ) )
	one_data = all_data[sys.argv[1]]

	num_points_left = len(one_data[Metric.LEFT_FFT_KEY])
	num_points_right = len(one_data[Metric.RIGHT_FFT_KEY])
	max_either_side = max(num_points_left, num_points_right)
	min_either_side = min(num_points_left, num_points_right)

	c1 = one_data[Metric.LEFT_FFT_KEY]
	c2 = one_data[Metric.RIGHT_FFT_KEY]
	c3 = one_data[Metric.LEFT_CURVATURE_KEY]
	c4 = one_data[Metric.RIGHT_CURVATURE_KEY]
	print ("Lengths: " + str(len(c1)) + ", " + str(len(c2)) + ", " + str(len(c3)) + ", " + str(len(c4)))


	# print "Data for " + sys.argv[1]
	# print("\t{:10}\t\t{:10}\t\t{:10}\t{:10}".format("FD (L)", "FD (R)", "Curve (L)", "Curve (R)"))
	# for i in xrange(max_either_side):
	# 	if i < min_either_side:
	# 		print(str(i) + "\t{:.4f}\t{:16.4f}\t{:8.4f}\t{:.4f}".format(c1[i], c2[i], c3[i], c4[i]))
	# 	elif i < num_points_left:
	# 		print(str(i) + "\t{:.4f}{:29.4f}".format(c1[i], c3[i]))
	# 	else:
	# 		print(str(i) + "\t\t\t\t{:20.4f}\t\t\t\t{:.4f}".format(c2[i], c4[i]))

	print one_data[Metric.X_KEY]
	print one_data[Metric.Y_KEY]
