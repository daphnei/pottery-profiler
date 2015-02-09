import sys
import math

class Length_iterator:
	def __init__(self, reclimit, err):
		self.Side = ["LEFT", "RIGHT"]
		self.curLeafCtrlPolyLengths = [0] * 3

		self.limit = reclimit;
		self.minTincrement = float(1) / (1 << self.limit);
		self.ERR = err;
		self.recCurveStack = [[0 for x in range(8)] for x in range(reclimit+1)];
		self.sides = self.Side[0]; #Fix this line?
		#if any methods are called without first initializing this object on
		#a curve, we want it to fail ASAP.
		self.nextT = sys.float_info.max
		self.lenAtNextT = sys.float_info.max
		self.lenAtLastSplit = sys.float_info.min
		self.recLevel = sys.float_info.min
		self.lastSegLen = sys.float_info.max
		self.done = True;

		self.cachedHaveLowAcceleration = -1

		self.nextRoots = [float(0)] * 4;
		self.flatLeafCoefCache = [0, 0, -1, 0];

	def initializeIterationOnCurve(self, pts, type_id):
		#self.recCurveStack[0:type_id] = pts[0:type_id]
		self.recCurveStack[0] = pts

		self.curveType = type_id;
		self.recLevel = 0;
		self.lastT = 0;
		self.lenAtLastT = 0;
		self.nextT = 0;
		self.lenAtNextT = 0;
		self.goLeft(); # initializes nextT and lenAtNextT properly
		self.lenAtLastSplit = 0;
		if (recLevel > 0):
			self.sides[0] = self.Side[0];
			self.done = False;
		else:
			#the root of the tree is a leaf so we're done.
			self.sides[0] = self.Side[1];
			self.done = True;

		self.lastSegLen = 0;

	def haveLowAcceleration(self, err):
		if cachedHaveLowAcceleration == -1:
			len1 = curLeafCtrlPolyLengths[0];
			len2 = curLeafCtrlPolyLengths[1];
			# the test below is equivalent to !within(len1/len2, 1, err).
			# It is using a multiplication instead of a division, so it
			# should be a bit faster.
			if not within(len1, len2, err*len2):
				self.cachedHaveLowAcceleration = 0;
				return False;
		
			if self.curveType == 8:
				len3 = self.curLeafCtrlPolyLengths[2];
				# if len1 is close to 2 and 2 is close to 3, that probably
				# means 1 is close to 3 so the second part of this test might
				# not be needed, but it doesn't hurt to include it.
				if not (within(len2, len3, err*len3) and within(len1, len3, err*len3)):
					self.cachedHaveLowAcceleration = 0;
					return False;

			self.cachedHaveLowAcceleration = 1;
			return True;

		return (self.cachedHaveLowAcceleration == 1);
	

	# returns the t value where the remaining curve should be split in
	# order for the left subdivided curve to have length len. If len
	#is >= than the length of the uniterated curve, it returns 1.
	def next(self, len1):
		targetLength = self.lenAtLastSplit + len1;
		while lenAtNextT < targetLength:
			if self.done:
				self.lastSegLen = self.lenAtNextT - self.lenAtLastSplit;
				return 1;
			
			self.goToNextLeaf();
		
		self.lenAtLastSplit = targetLength;
		leaflen = self.lenAtNextT - self.lenAtLastT;
		t = (targetLength - self.lenAtLastT) / leaflen;

		# cubicRootsInAB is a fairly expensive call, so we just don't do it
		# if the acceleration in this section of the curve is small enough.
		if not haveLowAcceleration(float(0.05)):
			# We flatten the current leaf along the x axis, so that we're
			# left with a, b, c which define a 1D Bezier curve. We then
			# solve this to get the parameter of the original leaf that
			# gives us the desired length.

			if self.flatLeafCoefCache[2] < 0:
				x = 0 + self.curLeafCtrlPolyLengths[0],
				y = x + self.curLeafCtrlPolyLengths[1];

				if self.curveType == 8:
					z = y + self.curLeafCtrlPolyLengths[2];
					self.flatLeafCoefCache[0] = 3*(x - y) + z;
					self.flatLeafCoefCache[1] = 3*(y - 2*x);
					self.flatLeafCoefCache[2] = 3*x;
					self.flatLeafCoefCache[3] = -z;
				elif curveType == 6:
					self.flatLeafCoefCache[0] = float(0);
					self.flatLeafCoefCache[1] = y - 2*x;
					self.flatLeafCoefCache[2] = 2*x;
					self.flatLeafCoefCache[3] = -y;
		   
			a = self.flatLeafCoefCache[0];
			b = self.flatLeafCoefCache[1];
			c = self.flatLeafCoefCache[2];
			d = t * self.flatLeafCoefCache[3];

			# we use cubicRootsInAB here, because we want only roots in 0, 1,
			# and our quadratic root finder doesn't filter, so it's just a
			# matter of convenience.
			n = cubicRootsInAB(a, b, c, d, self.nextRoots, 0, 0, 1);
			if n == 1 and not math.isnan(self.nextRoots[0]):
				t = self.nextRoots[0];
			
		
		# t is relative to the current leaf, so we must make it a valid parameter
		# of the original curve.
		t = t * (self.nextT - self.lastT) + self.lastT;
		if (t >= 1):
			t = 1;
			self.done = True;
		
		# even if done = true, if we're here, that means targetLength
		# is equal to, or very, very close to the total length of the
		# curve, so lastSegLen won't be too high. In cases where len
		# overshoots the curve, this method will exit in the while
		# loop, and lastSegLen will still be set to the right value.
		self.lastSegLen = len1;
		return t;
	
	def lastSegLen(self):
		return self.lastSegLen

	# go to the next leaf (in an inorder traversal) in the recursion tree
	# preconditions: must be on a leaf, and that leaf must not be the root.
	def goToNextLeaf(self):
		# We must go to the first ancestor node that has an unvisited
		# right child.
		self.recLevel -= 1;
		while self.sides[recLevel] == self.Side[1]:
			if self.recLevel == 0:
				self.done = True;
				return;
			self.recLevel -= 1;

		self.sides[ self.recLevel ] = self.Side[1];
		
		#System.arraycopy(recCurveStack[recLevel], 0, recCurveStack[recLevel+1], 0, curveType);
		self.recCurveStack[recLevel] = self.recCurveStack[self.recLevel+1]

		self.recLevel += 1;
		self.goLeft();
	
	# go to the leftmost node from the current node. Return its length.
	def goLeft(self):
		len1 = self.onLeaf();
		if len1 >= 0:
			self.lastT = nextT;
			self.lenAtLastT = self.lenAtNextT;
			self.nextT += (1 << (self.limit - self.recLevel)) * self.minTincrement;
			self.lenAtNextT += len1;
			# invalidate caches
			self.flatLeafCoefCache[2] = -1;
			self.cachedHaveLowAcceleration = -1;
		else:
			subdivide(self.recCurveStack[self.recLevel], 0,
							  self.recCurveStack[self.recLevel+1], 0,
							  self.recCurveStack[self.recLevel], 0, self.curveType);
			self.sides[recLevel] = self.Side[0];
			self.recLevel += 1;
			self.goLeft();

	# this is a bit of a hack. It returns -1 if we're not on a leaf, and
	# the length of the leaf if we are on a leaf.
	def onLeaf(self):
		curve = self.recCurveStack[self.recLevel];
		polyLen = 0;

		x0 = curve[0]
		y0 = curve[1];
		#for (int i = 2; i < curveType; i += 2) {
		i = 2
		while i < self.curveType:
			x1 = curve[i]
			y1 = curve[i+1];
			
			len1 = linelen(x0, y0, x1, y1);
			polyLen += len1;
			self.curLeafCtrlPolyLengths[i/2 - 1] = len1;
			x0 = x1;
			y0 = y1;

			i += 2
		

		lineLen = linelen(self.curve[0], self.curve[1], self.curve[self.curveType-2], self.curve[self.curveType-1]);
		if olyLen - lineLen < self.ERR or self.recLevel == self.limit:
			return (polyLen + lineLen)/2;
		
		return -1
	

	def subdivide(src, srcoff, left, leftoff, right, rightoff, type):
		if type == 6:
			subdivideQuad(src, srcoff, left, leftoff, right, rightoff);
		elif type == 8:
			subdivideCubic(src, srcoff, left, leftoff, right, rightoff);
		else:
			print("Unsupported curve type")
			exit(1)

	#  // Most of these are copied from classes in java.awt.geom because we need
	# // float versions of these functions, and Line2D, CubicCurve2D,
	# // QuadCurve2D don't provide them.
	# /**
	#  * Subdivides the cubic curve specified by the coordinates
	#  * stored in the <code>src</code> array at indices <code>srcoff</code>
	#  * through (<code>srcoff</code>&nbsp;+&nbsp;7) and stores the
	#  * resulting two subdivided curves into the two result arrays at the
	#  * corresponding indices.
	#  * Either or both of the <code>left</code> and <code>right</code>
	#  * arrays may be <code>null</code> or a reference to the same array
	#  * as the <code>src</code> array.
	#  * Note that the last point in the first subdivided curve is the
	#  * same as the first point in the second subdivided curve. Thus,
	#  * it is possible to pass the same array for <code>left</code>
	#  * and <code>right</code> and to use offsets, such as <code>rightoff</code>
	#  * equals (<code>leftoff</code> + 6), in order
	#  * to avoid allocating extra storage for this common point.
	#  * @param src the array holding the coordinates for the source curve
	#  * @param srcoff the offset into the array of the beginning of the
	#  * the 6 source coordinates
	#  * @param left the array for storing the coordinates for the first
	#  * half of the subdivided curve
	#  * @param leftoff the offset into the array of the beginning of the
	#  * the 6 left coordinates
	#  * @param right the array for storing the coordinates for the second
	#  * half of the subdivided curve
	#  * @param rightoff the offset into the array of the beginning of the
	#  * the 6 right coordinates
	#  * @since 1.7
	#  */
	def subdivideCubic(src, srcoff, left, leftoff, right, rightoff):
		x1 = src[srcoff + 0];
		y1 = src[srcoff + 1];
		ctrlx1 = src[srcoff + 2];
		ctrly1 = src[srcoff + 3];
		ctrlx2 = src[srcoff + 4];
		ctrly2 = src[srcoff + 5];
		x2 = src[srcoff + 6];
		y2 = src[srcoff + 7];
		if (left != null):
			left[leftoff + 0] = x1;
			left[leftoff + 1] = y1;
		
		if (right != null):
			right[rightoff + 6] = x2;
			right[rightoff + 7] = y2;
		
		x1 = (x1 + ctrlx1) / float(2.0);
		y1 = (y1 + ctrly1) / float(2.0);
		x2 = (x2 + ctrlx2) / float(2.0);
		y2 = (y2 + ctrly2) / float(2.0);
		centerx = (ctrlx1 + ctrlx2) / float(2.0);
		centery = (ctrly1 + ctrly2) / float(2.0);
		ctrlx1 = (x1 + centerx) / float(2.0);
		ctrly1 = (y1 + centery) / float(2.0);
		ctrlx2 = (x2 + centerx) / float(2.0);
		ctrly2 = (y2 + centery) / float(2.0);
		centerx = (ctrlx1 + ctrlx2) / float(2.0);
		centery = (ctrly1 + ctrly2) / float(2.0);
		if (left != null):
			left[leftoff + 2] = x1;
			left[leftoff + 3] = y1;
			left[leftoff + 4] = ctrlx1;
			left[leftoff + 5] = ctrly1;
			left[leftoff + 6] = centerx;
			left[leftoff + 7] = centery;
		
		if (right != null):
			right[rightoff + 0] = centerx;
			right[rightoff + 1] = centery;
			right[rightoff + 2] = ctrlx2;
			right[rightoff + 3] = ctrly2;
			right[rightoff + 4] = x2;
			right[rightoff + 5] = y2;
		
	def subdivideCubicAt(t, src, srcoff, left, leftoff, right, rightoff):
		x1 = src[srcoff + 0];
		y1 = src[srcoff + 1];
		ctrlx1 = src[srcoff + 2];
		ctrly1 = src[srcoff + 3];
		ctrlx2 = src[srcoff + 4];
		ctrly2 = src[srcoff + 5];
		x2 = src[srcoff + 6];
		y2 = src[srcoff + 7];
		if (left != None):
			left[leftoff + 0] = x1;
			left[leftoff + 1] = y1;
		
		if (right != None):
			right[rightoff + 6] = x2;
			right[rightoff + 7] = y2;
		
		x1 = x1 + t * (ctrlx1 - x1);
		y1 = y1 + t * (ctrly1 - y1);
		x2 = ctrlx2 + t * (x2 - ctrlx2);
		y2 = ctrly2 + t * (y2 - ctrly2);
		centerx = ctrlx1 + t * (ctrlx2 - ctrlx1);
		centery = ctrly1 + t * (ctrly2 - ctrly1);
		ctrlx1 = x1 + t * (centerx - x1);
		ctrly1 = y1 + t * (centery - y1);
		ctrlx2 = centerx + t * (x2 - centerx);
		ctrly2 = centery + t * (y2 - centery);
		centerx = ctrlx1 + t * (ctrlx2 - ctrlx1);
		centery = ctrly1 + t * (ctrly2 - ctrly1);
		if (left != None):
			left[leftoff + 2] = x1;
			left[leftoff + 3] = y1;
			left[leftoff + 4] = ctrlx1;
			left[leftoff + 5] = ctrly1;
			left[leftoff + 6] = centerx;
			left[leftoff + 7] = centery;
		
		if (right != None):
			right[rightoff + 0] = centerx;
			right[rightoff + 1] = centery;
			right[rightoff + 2] = ctrlx2;
			right[rightoff + 3] = ctrly2;
			right[rightoff + 4] = x2;
			right[rightoff + 5] = y2;
	   
   

	def subdivideQuad(src, srcoff, left, leftoff, right, rightoff):
		x1 = src[srcoff + 0];
		y1 = src[srcoff + 1];
		ctrlx = src[srcoff + 2];
		ctrly = src[srcoff + 3];
		x2 = src[srcoff + 4];
		y2 = src[srcoff + 5];
		if (left != None):
			left[leftoff + 0] = x1;
			left[leftoff + 1] = y1;
		
		if (right != None):
			right[rightoff + 4] = x2;
			right[rightoff + 5] = y2;
		
		x1 = (x1 + ctrlx) / float(2.0);
		y1 = (y1 + ctrly) / float(2.0);
		x2 = (x2 + ctrlx) / float(2.0);
		y2 = (y2 + ctrly) / float(2.0);
		ctrlx = (x1 + x2) / float(2.0);
		ctrly = (y1 + y2) / float(2.0);
		if (left != None):
			left[leftoff + 2] = x1;
			left[leftoff + 3] = y1;
			left[leftoff + 4] = ctrlx;
			left[leftoff + 5] = ctrly;
		
		if (right != None):
			right[rightoff + 0] = ctrlx;
			right[rightoff + 1] = ctrly;
			right[rightoff + 2] = x2;
			right[rightoff + 3] = y2;
		
	

	def subdivideQuadAt(t, src, srcoff, left, leftoff, right, rightoff):
		x1 = src[srcoff + 0];
		y1 = src[srcoff + 1];
		ctrlx = src[srcoff + 2];
		ctrly = src[srcoff + 3];
		x2 = src[srcoff + 4];
		y2 = src[srcoff + 5];
		if (left != None):
			left[leftoff + 0] = x1;
			left[leftoff + 1] = y1;
		
		if (right != None):
			right[rightoff + 4] = x2;
			right[rightoff + 5] = y2;
		
		x1 = x1 + t * (ctrlx - x1);
		y1 = y1 + t * (ctrly - y1);
		x2 = ctrlx + t * (x2 - ctrlx);
		y2 = ctrly + t * (y2 - ctrly);
		ctrlx = x1 + t * (x2 - x1);
		ctrly = y1 + t * (y2 - y1);
		if (left != None):
			left[leftoff + 2] = x1;
			left[leftoff + 3] = y1;
			left[leftoff + 4] = ctrlx;
			left[leftoff + 5] = ctrly;
		
		if (right != None):
			right[rightoff + 0] = ctrlx;
			right[rightoff + 1] = ctrly;
			right[rightoff + 2] = x2;
			right[rightoff + 3] = y2;
	   

	def subdivideAt(t, src, srcoff, left, leftoff, right, rightoff, size):
		if (size == 8):
			subdivideCubicAt(t, src, srcoff, left, leftoff, right, rightoff);
		elif (size== 6):
			subdivideQuadAt(t, src, srcoff, left, leftoff, right, rightoff);
		
	


