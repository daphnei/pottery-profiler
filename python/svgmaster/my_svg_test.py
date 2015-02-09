import sys, math
import svg

def Distance(p1, p2):
	return math.sqrt((p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2)

# Find index i such that dd[i] < d < dd[i+1]
def FindIndex(dd, d):
	'''dd is an array of doubles. d is a double '''

	i = 0
 	for j in xrange(0, len(dd)):
		if (d > dd[j]):
			i = j;

	return i;


if __name__ == "__main__":
	f = svg.parse(sys.argv[1])

	#print f
	#print f.items
	g = f.items[0]
	svg_components = g.items[0]
	print svg_components

	#for it in svg_components.items:
		#if type(it) is svg.geometry.Bezier:
	#	print dir(it)

	bez_0 = svg_components.items[1]

	#Construct polyline with large number of points
	n = 1000;
	s = 1.0/(n-1);
	tt = [0] * n;	#An array of doubles
	PP = [None] * n; #An array of points
	cc = [0] * n;	#An array of doubles

	for i in xrange(0, n):
		tt[i] = i*s;
		PP[i] = bez_0._bezierN(tt[i]);
		if (i > 0):
			cc[i] = Distance(PP[i], PP[i-1]);

	#Get fractional arclengths along polyline
	dd = [0] * n #An array of doubles
	dd[0] = 0
	for i in xrange(1, n):
		dd[i] = dd[i-1] + cc[i];
	for i in xrange(1, n):
		dd[i] = dd[i]/dd[n-1];

	#Number of points to place on curve
	m = 10; # m is an int
	step = 1.0/(m-1)

	QQ = [None] * m; #An array of points

	for r in xrange(0, m):
		d = r*step;
		i = FindIndex(dd, d);
		u = (d - dd[i]) / (dd[i+1] - dd[i]);
		t = (i + u)*s;
		QQ[r] = bez_0._bezierN(t);

	for q in QQ:
		print str(q.x) + "\t" + str(q.y)
	
