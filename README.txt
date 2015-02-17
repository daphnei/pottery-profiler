# pottery-profiler
An attempt to classify pottery sherds based on SVGs of their profiles. 

USAGE: python pottery_predict.py <svg file>
The interpolated points will be outputted to file "output.svg" which has a circle at each point location. 

NOTE: pottery_predict relies on a bunch of libraries being installed. Before running it, you will need to do
sudo pip install svgwrite
sudo pip install svg.path
sudo pip install numpy
