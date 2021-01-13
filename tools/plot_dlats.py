#!/usr/bin/env python3

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np

num_lat = 1800

lat = np.linspace(-np.pi/2, np.pi/2, num_lat)
wgt = np.exp(-50 * (np.abs(lat[1:len(lat)-1]) - np.pi/2)**4)

dlat = (1 + 2 * wgt) / np.sum(1 + 2 * wgt) * np.pi

fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(1, 1, 1)
plt.plot(np.degrees(lat[1:len(lat)-1]), np.degrees(dlat))
plt.show()
plt.close()
