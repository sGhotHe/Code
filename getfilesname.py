import os
import sys

def file_name(file_dir):
	fn = []
	try:
		for root, dirs, files in os.walk(file_dir):
			r = root
			for f in files:
				fn.append(f)
	except:
		print('No such directory. Please check and try again.')
		sys.exit()
	return r, fn
		
if __name__ == '__main__':
	file_dir = "data/test"
	root, files = file_name(file_dir)
	print(files)
