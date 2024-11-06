from multiprocessing import Pool
import numpy as np

def square(a):
	return a*a, a
	
if __name__ == "__main__":
	a = np.random.rand(100)
	with Pool(processes=40) as pool:
		results = pool.map(square, a)
	results = np.array(results)
	a = results[:,0]
	b = results[:,1]
	print(a, b)
