__author__ = 'daphne'

SEG_LENGTH = 3

#Note, when a new key is added, it should also be appended to the list bel
LEFT_FFT_KEY = "left_fft"               # A fourier transform of the points on the left profile.
RIGHT_FFT_KEY = "right_fft"				# A fourier transform of the points on the right profile.
LEFT_CURVATURE_KEY = "left_curvature"	# Radius of curvature along the left profile.
RIGHT_CURVATURE_KEY = "right_curvature" # Radius of curvature along the right profile.
LEFT_TANGENT_KEY = "left_tangent"
RIGHT_TANGENT_KEY = "right_tangent"
LEFT_TANGENT_ALT_KEY = "left_tangent_alt"
RIGHT_TANGENT_ALT_KEY = "right_tangent_alt"
THICKNESS_KEY = "thickness"

ALL_KEYS = [LEFT_FFT_KEY, RIGHT_FFT_KEY, LEFT_CURVATURE_KEY, RIGHT_CURVATURE_KEY, LEFT_TANGENT_KEY, RIGHT_TANGENT_KEY,
			LEFT_TANGENT_ALT_KEY, RIGHT_TANGENT_ALT_KEY]

X_KEY = "x"								# Points generated along the right profile.
Y_KEY = "y"								# Points generated along the right profile.

DESC_OUTPUT_FILE = "sherd_data.pickle"