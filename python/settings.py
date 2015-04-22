__author__ = 'daphne'

SEG_LENGTH = 3

LEFT_FFT_KEY = "left_fft"               # A fourier transform of the points on the left profile.
RIGHT_FFT_KEY = "right_fft"				# A fourier transform of the points on the right profile.
LEFT_CURVATURE_KEY = "left_curvature"	# Radius of curvature along the left profile.
RIGHT_CURVATURE_KEY = "right_curvature" # Radius of curvature along the right profile.
LEFT_TANGENT_KEY = "left_tangent"
RIGHT_TANGENT_KEY = "right_tangent"
X_KEY = "x"								# Points generated along the right profile.
Y_KEY = "y"								# Points generated along the right profile.

DESC_OUTPUT_FILE = "sherd_data.pickle"