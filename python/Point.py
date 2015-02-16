import math

class Point:
    def __init__(self, x=None, y=None):
        '''A Point is defined either by a tuple/list of length 2 or
           by 2 coordinates
        >>> Point(1,2)
        (1.000,2.000)
        >>> Point((1,2))
        (1.000,2.000)
        >>> Point([1,2])
        (1.000,2.000)
        >>> Point('1', '2')
        (1.000,2.000)
        >>> Point(('1', None))
        (1.000,0.000)
        '''
        if (isinstance(x, tuple) or isinstance(x, list)) and len(x) == 2:
            x,y = x

        # Handle empty parameter(s) which should be interpreted as 0
        if x is None: x = 0
        if y is None: y = 0

        try:
            self.x = float(x)
            self.y = float(y)
        except:
            raise TypeError("A Point is defined by 2 numbers or a tuple")

    def __add__(self, other):
        '''Add 2 points by adding coordinates.
        Try to convert other to Point if necessary
        >>> Point(1,2) + Point(3,2)
        (4.000,4.000)
        >>> Point(1,2) + (3,2)
        (4.000,4.000)'''
        if not isinstance(other, Point):
            try: other = Point(other)
            except: return NotImplemented
        return Point(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        '''Substract two Points.
        >>> Point(1,2) - Point(3,2)
        (-2.000,0.000)
        '''
        if not isinstance(other, Point):
            try: other = Point(other)
            except: return NotImplemented
        return Point(self.x - other.x, self.y - other.y)

    def __mul__(self, other):
        '''Multiply a Point with a constant.
        >>> 2 * Point(1,2)
        (2.000,4.000)
        >>> Point(1,2) * Point(1,2) #doctest:+IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        TypeError:
        '''
        if not isinstance(other, numbers.Real):
            return NotImplemented
        return Point(self.x * other, self.y * other)
    def __rmul__(self, other):
        return self.__mul__(other)

    def __eq__(self, other):
        '''Test equality
        >>> Point(1,2) == (1,2)
        True
        >>> Point(1,2) == Point(2,1)
        False
        '''
        if not isinstance(other, Point):
            try: other = Point(other)
            except: return NotImplemented
        return (self.x == other.x) and (self.y == other.y)

    def __repr__(self):
        return '(' + format(self.x,'.3f') + ',' + format( self.y,'.3f') + ')'

    def __str__(self):
        return self.__repr__();

    def coord(self):
        '''Return the point tuple (x,y)'''
        return (self.x, self.y)

    def length(self):
        '''Vector length, Pythagoras theorem'''
        return math.sqrt(self.x ** 2 + self.y ** 2)

    def rot(self, angle):
        '''Rotate vector [Origin,self] '''
        if not isinstance(angle, Angle):
            try: angle = Angle(angle)
            except: return NotImplemented
        x = self.x * angle.cos - self.y * angle.sin
        y = self.x * angle.sin + self.y * angle.cos
        return Point(x,y)

    def distance(self, other):
        return math.sqrt((self.x - other.x) ** 2 + (self.y - other.y) ** 2)

