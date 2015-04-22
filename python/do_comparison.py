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
	reduced_data = PCA(n_components=2).fit_transform(data_for_kmeans)

	kmeans = KMeans(init='k-means++', n_clusters=4, n_init=10)
	kmeans.fit(reduced_data)

	# labeler = KMeans(n_clusters=4)
	# # labeler.fit(reduced_data)
	#
	# print labeler
	# labels = labeler.predict(data_for_kmeans)
	#
	# for i in xrange(len(names)):
	# 	print names[i] + ": " + str(labels[i])


	# Plot the decision boundary. For that, we will assign a color to each
	x_min, x_max = reduced_data[:, 0].min() + 1, reduced_data[:, 0].max() - 1
	y_min, y_max = reduced_data[:, 1].min() + 1, reduced_data[:, 1].max() - 1
	xx, yy = np.meshgrid(np.arange(0, 1, 0.02), np.arange(-1000, 1000, 5))

	# Obtain labels for each point in mesh. Use last trained model.
	Z = kmeans.predict(np.c_[xx.ravel(), yy.ravel()])

	# Put the result into a color plot
	Z = Z.reshape(xx.shape)
	plt.figure(1)
	plt.clf()
	plt.imshow(Z, interpolation='nearest',
			   extent=(xx.min(), xx.max(), yy.min(), yy.max()),
			   cmap=plt.cm.Paired,
			   aspect='auto', origin='lower')

	plt.plot(reduced_data[:, 0], reduced_data[:, 1], 'k.', markersize=2)
	# Plot the centroids as a white X
	centroids = kmeans.cluster_centers_
	plt.scatter(centroids[:, 0], centroids[:, 1],
				marker='x', s=169, linewidths=3,
				color='w', zorder=10)
	plt.title('K-means clustering on the digits dataset (PCA-reduced data)\n'
			  'Centroids are marked with white cross')
	plt.xlim(x_min, x_max)
	plt.ylim(y_min, y_max)
	plt.xticks(())
	plt.yticks(())
	plt.show()

def check_one_svg(target_name, desc_type=RIGHT_CURVATURE_KEY):
	'''Compares one svg to all the other svgs in the directory. Prints out match scores.
	   I tried to base this as much as possible on what was done in the old crane script.
	'''
	print "Comparing " + target_name

	all_data = pickle.load( open( "sherd_data.pickle", "rb" ) )

	#Find the one we are comparing against.
	target_obj = all_data[target_name]

	target_descriptors = numpy.asarray(target_obj[desc_type])

	#for each sherd in the database, store its average distance from the target
	dists = {}

	window_size = 1

	for shape_name, shape_data in all_data.items():
		# gets euclid distance for 2 sets of points. Only consider the first n descriptors where n is the minimum
		# length of either list
		tomatch_descriptors = numpy.asarray(shape_data[desc_type])
		n = min(tomatch_descriptors.size, target_descriptors.size)

		if n != 0:
			#dists[shape_name] = numpy.linalg.norm(target_descriptors[:n] - tomatch_descriptors[:n])

			means_targ = []
			means_tomatch = []
			#For each 5 points find the average.
			for i in xrange(n - window_size):
				subset1 = target_descriptors[i:i + window_size]
				subset2 = tomatch_descriptors[i:i + window_size]

				#Find the average of each subset
				means_targ.append(numpy.mean(subset1))
				means_tomatch.append(numpy.mean(subset2))

			# means_targ = target_descriptors[:n]
			# means_tomatch = tomatch_descriptors[:n]

			dists[shape_name] = (numpy.linalg.norm(numpy.asarray(means_targ) - numpy.asarray(means_tomatch)))


		#In order to not give preference to shorter profiles, normalize by the number of points being operated on.
		#if n > 0:
		#	dists[shape_name] = dists[shape_name] / n

	#Sort the items by how close they are to the target one.
	sorted_dists = sorted(dists.items(), key=lambda dists: dists[1])

	print("{:25}\t{:10}".format("name", "matchValue"))
	for i in xrange(len(sorted_dists)):
		print("{:25}\t{:.6f}".format(sorted_dists[i][0], sorted_dists[i][1]))

if __name__ == "__main__":
	check_one_svg(sys.argv[1], RIGHT_TANGENT_KEY);

	#kmeans_attempt() sys.argv[1]

	# path = get_path_from_svg("/Users/daphne/Documents/School/CSC494/pottery-profiler/Pottery/" + sys.argv[1])
	# points = get_points_along_path(path)
	# left_profile_points, right_profile_points = split_profile_points(points)
	# draw_points_to_output_file(left_profile_points, right_profile_points)
