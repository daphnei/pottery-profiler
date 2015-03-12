'''
This script takes the name of one of the sherds stored in the database,
and prints out all of the descriptors for that sherd.
'''

__author__ = 'daphne'

import pickle
import sys
from settings import *

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print "USAGE: python print_data_for_one.py <name>"
		exit(1)

	all_data = pickle.load( open(DESC_OUTPUT_FILE, "rb" ) )
	one_data = all_data[sys.argv[1]]

	num_points = len(one_data[LEFT_FFT_KEY])

	c1 = one_data[LEFT_FFT_KEY]
	c2 = one_data[RIGHT_FFT_KEY]
	c3 = one_data[LEFT_CURVATURE_KEY]
	c4 = one_data[RIGHT_CURVATURE_KEY]

	print "Data for " + sys.argv[1]
	print("{:10}\t\t{:10}\t\t{:10}\t{:10}".format("FD (L)", "FD (R)", "Curve (L)", "Curve (R)"))
	for i in xrange(num_points):
		print("{:.4f}\t{:16.4f}\t{:8.4f}\t{:.4f}".format(c1[i], c2[i], c3[i], c4[i]))
