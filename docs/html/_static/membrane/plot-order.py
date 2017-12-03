#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

plt.style.use('ggplot')
mpl.rc('lines', linewidth=2, markeredgewidth=1, markersize=10)
mpl.rc('axes', labelsize='x-large', titlesize="xx-large")
mpl.rc('xtick', labelsize="x-large")
mpl.rc('ytick', labelsize="x-large")

if len(sys.argv) != 2:
    print("USAGE: {} order.dat".format(sys.argv[0]))
    sys.exit(0)

dat_file = sys.argv[1]
png_file = "{}.png".format(os.path.splitext(dat_file)[0])

beads = []
order_values = []

with open(dat_file) as fp:
    for lino, line in enumerate(fp):
        if lino == 0:
            line = line.split()

            beads = [val.strip() for val in line[1:]]

        if lino > 2:
            if line.startswith("-"):
                break

            line = line.split()

            order_values.append([float(val) for val in line[1:]])

order_values = np.array(order_values)
x = np.arange(len(beads))
y = order_values.mean(axis=0)
y_err = order_values.std(axis=0)

plt.errorbar(x, y, yerr=y_err)
plt.gca().xaxis.set_ticklabels(beads)
plt.xlabel("Martini bond")
plt.ylabel("P_2 order")
plt.title("Order profile")
plt.tight_layout()
plt.savefig(png_file)

print("INFO: {} successfully converted to {}".format(dat_file, png_file))
