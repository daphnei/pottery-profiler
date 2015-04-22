import numpy
import math

def dfourier(x,T):
	l = 2 * math.pi/T;
	t = len(x);
	f = numpy.fft.fft(x);
	a1 = range(0, int(numpy.fix((t-1)/2)+1))
	a2 = range( -int(numpy.fix(t/2)), -1 + 1)
	a = numpy.concatenate((a1, a2), 0)
	b = numpy.array([0j] * t)
	for k in xrange (0,t):
		b[k] = (1j*l) * (f[k] * a[k])

	d = numpy.real(numpy.fft.ifft(b))
	return d

print dfourier([1, 4, 6, 3, 19, 20], 6)