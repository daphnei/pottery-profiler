import pickle

def qwert():
	print "Hello"

if __name__ == "__main__":
	a = qwert
	a()

	f = open("p.pickle", "w")

	pickle.dump(a, f, 1)

	f.close()