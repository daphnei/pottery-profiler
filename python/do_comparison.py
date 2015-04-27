import scipy

__author__ = 'daphne'
import sys
import numpy
from point import Point
import math
import pickle
import cmath
from settings import *
import numpy as np
from scipy.sparse import *
from sklearn.cluster import KMeans
from sklearn.preprocessing import Imputer
import matplotlib.pyplot as plt
from sklearn.datasets import load_digits
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
from scipy.cluster.vq import kmeans, vq

def get_max_desc_length(all_data, desc_type):
	inner_dicts = all_data.values()
	desc_lists = map(lambda x: x[desc_type], inner_dicts)
	return max(len(sublist) for sublist in desc_lists)

def kmeans_attempt():
	#coo_matrix

	all_data = pickle.load( open( "sherd_data.pickle", "rb" ) )


	#get the max length of any curve descriptors and the max length of any fft descriptors
	max_length = get_max_desc_length(all_data, RIGHT_FFT_KEY)

	names = []

	data_for_kmeans = np.zeros((len(all_data), 3 * max_length))  #Instead should be initing a Numpy array?
	for i, (shape_name, shape_data), in enumerate(all_data.items()):
		curve_descs = shape_data[RIGHT_CURVATURE_KEY]
		fft_descs = shape_data[RIGHT_FFT_KEY]
		fft_descs_real = list(x.real for x in fft_descs)
		fft_descs_imag = list(x.imag for x in fft_descs)

		n = len(curve_descs)

		if n != 0:
			#Pad the lists up to their maximum length.
			curve_descs = curve_descs + ((max_length - len(curve_descs)) * [np.mean(curve_descs)])
			fft_descs_real = fft_descs_real + ((max_length - len(fft_descs)) * [np.mean(fft_descs_real)])
			fft_descs_imag = fft_descs_imag+ ((max_length - len(fft_descs)) * [np.mean(fft_descs_imag)])

			#Append all the lists together.
			combined_data = np.array(curve_descs + fft_descs_real + fft_descs_imag) #This should be a row of the matrix

			# print
			# print combined_data

			#Do something to add the row to the matrix
			data_for_kmeans[i, :] = combined_data

			names.append(shape_name)

	print ("Starting PCA")

	#Reduce the dimensionality by some amount
	reduced_data = data_for_kmeans#PCA(n_components=20).fit_transform(data_for_kmeans)

	codebook, distortion = kmeans(reduced_data, 4)

	clusters, n = vq(reduced_data, codebook)

	l = list(str(index) + ',\t ' + name for name, index in zip(names, clusters))
	l.sort()
	for label in l:
		print label

def do_comp_with_all_metrics(target_name, weights=None):
	print "Comparing: " + target_name

	all_data = pickle.load( open( "sherd_data.pickle", "rb" ) )

	#Find the one we are comparing against.
	target_obj = all_data[target_name]

	#for each sherd in the database, store its average distance from the target
	dists = {}

	#If no weights were provided, default to equal weights.
	if weights == None:
		weights = [1.0 / len(ALL_KEYS)] * len(ALL_KEYS)

	for metric_index in xrange(len(ALL_KEYS)):
		metric = ALL_KEYS[metric_index]

		target_descriptors = numpy.asarray(target_obj[metric])

		for shape_name, shape_data in all_data.items():
			tomatch_descriptors = numpy.asarray(shape_data[metric])

			if shape_name not in dists:
				dists[shape_name] = 0

			dists[shape_name] += compare(target_descriptors, tomatch_descriptors) * weights[metric_index]



	#Print out the results.
	print_comparison_result(dists)



def do_comp_with_one_metric(target_name, metric, normalize_by_length=False):
	'''Compares one svg to all the other svgs in the directory. Prints out match scores.

		Inputs
		target_name: The name of the svg that should be compared with all other entries in the database.
		metric: Any one of the keys found in settings.py
		normalize_by_length: whether or not match scores should be divided by the number of datapoints used in the match


	'''
	print "Comparing " + target_name

	all_data = pickle.load( open( "sherd_data.pickle", "rb" ) )

	#Find the one we are comparing against.
	target_obj = all_data[target_name]

	target_descriptors = numpy.asarray(target_obj[metric])

	#for each sherd in the database, store its average distance from the target
	dists = {}

	for shape_name, shape_data in all_data.items():
		#We are checking how well the target_descriptors match against the tomatch_descriptors
		tomatch_descriptors = numpy.asarray(shape_data[metric])

		dists[shape_name] = compare(target_descriptors, tomatch_descriptors)

	#Print out the results.
	print_comparison_result(dists)

def print_comparison_result(dists):
	#Sort the items by how close they are to the target one.
	sorted_dists = sorted(dists.items(), key=lambda dists: dists[1])

	print("{:25}\t{:10}".format("name", "matchValue"))
	for i in xrange(len(sorted_dists)):
		print("{:25}\t{:.6f}".format(sorted_dists[i][0], sorted_dists[i][1]))

def compare(target_descriptors, tomatch_descriptors, normalize_by_length = False):

	#I attempted to do comparisons within a moving window in order to help resolve
	#the issue of small phase changes causing completely different matches. However,
	#setting the window size to 2-5 didn't end up drastically improving comparisons.
	window_size = 1

	#Check as many points as we have in both sherds.
	n = min(tomatch_descriptors.size, target_descriptors.size)

	distance = -1

	if n != 0:
		means_targ = []
		means_tomatch = []
		#For each 5 points find the average.
		for i in xrange(n - window_size):
			subset1 = target_descriptors[i:i + window_size]
			subset2 = tomatch_descriptors[i:i + window_size]

			#Find the average of each subset
			means_targ.append(numpy.mean(subset1))
			means_tomatch.append(numpy.mean(subset2))

		distance = (numpy.linalg.norm(numpy.asarray(means_targ) - numpy.asarray(means_tomatch)))

	#In order to not give preference to shorter profiles, normalize by the number of points being operated on.
	#This is more useful for some metrics (such as curvature) than for others (such as fourier descriptors)
	if normalize_by_length and n > 0:
		distance = distance / n

	return distance

if __name__ == "__main__":
	if len(sys.argv) == 3:
		print "Comparing " + sys.argv[1] + " using the " + sys.argv[2] + " metric:"
		do_comp_with_one_metric(sys.argv[1], sys.argv[2]);
	elif len(sys.argv) == 2:
		print "Comparing " + sys.argv[1] + " using all metrics:"
		do_comp_with_all_metrics(sys.argv[1])
	else:
		print "USAGE python do_comparison.py SVG_NAME [METRIC]"

	#kmeans_attempt()
