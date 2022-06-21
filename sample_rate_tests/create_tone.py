#!env python3
import sys
import math

# Number of samples should be a multiple of 1024 to match bladeRF
# buffer size constraints
n_samples = 1024

# IQ values are in the range [-2048, 2047]. Clamp to 1800 just to 
# avoid saturating
scale = 1800


if (len(sys.argv) < 2):
	print('Usage: ' + sys.argv[0] + ': <output file> [n_samples]\n')
	sys.exit(1)

if (len(sys.argv) > 2):
	try: 
		n_samples = int(sys.argv[2])
	except ValueError:
		print('Invalid value for n_samples: ' + sys.argv[2] + '\n')
		sys.exit(1)

	if n_samples < 1024 or n_samples % 1024 != 0:
		print('n_samples must be a multiple of 1024\n')
		sys.exit(1)

with open(sys.argv[1], 'w') as out_file:
	for n in range(0, n_samples):
		theta = n * (2 * math.pi) / n_samples 
		i = int(scale * math.cos(theta))
		q = int(scale * math.sin(theta))

		out_file.write(str(i) + ', ' + str(q) + '\n')