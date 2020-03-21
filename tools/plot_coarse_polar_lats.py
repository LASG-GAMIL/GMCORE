#!/usr/bin/env python3

import matplotlib
matplotlib.use('agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np

num_lat = 180

lat = np.linspace(-np.pi/2, np.pi/2, num_lat)
wgt = np.exp(-100 * (np.abs(lat[1:len(lat)-1]) - np.pi/2)**2)

dlat = (1 + 2 * wgt) / np.sum(1 + 2 * wgt) * np.pi

with PdfPages('coarse_polar_lats.pdf') as pdf:
	fig = plt.figure(figsize=(10, 5))
	ax = fig.add_subplot(1, 1, 1)
	plt.plot(lat[1:len(lat)-1], np.degrees(dlat))
	pdf.savefig()
	plt.close()
