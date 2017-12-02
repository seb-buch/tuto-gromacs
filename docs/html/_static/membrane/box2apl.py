#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, division
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


if len(sys.argv) != 3:
    print("USAGE: box2apl.py FILE.xvg NLIPIDPERLEAFLET")
    sys.exit(0)

xvg_file = sys.argv[1]
png_file = "{}_APL.png".format(os.path.splitext(xvg_file)[0])

title = "Area per lipids"
xaxis_title = "Time (ns)"
yaxis_title = "Area (nm^2)"
x = []
y = []
labels = []

try:
    with open(xvg_file) as fp:
        for line in fp:
            line = line.strip()

            if line.startswith("#"):
                continue

            if line.startswith("@"):
                line = line.split()
                if len(line) == 2:
                    continue

                if line[1] == "title":
                    title = " ".join(line[2:]).strip("\"")
                elif line[1] == "xaxis":
                    xaxis_title = " ".join(line[3:]).strip("\"")
                elif line[1] == "yaxis":
                    yaxis_title = " ".join(line[3:]).strip("\"")
                elif line[1].startswith("s") and line[2] == "legend":
                    labels.append(" ".join(line[3:]).strip("\""))
            else:
                line = [float(val) for val in line.split()]
                x.append(line[0])
                y.append(line[1])

except IOError:
    print("ERROR: Could not read '{}'".format(xvg_file))
    sys.exit(1)
except (ValueError, TypeError):
    print("ERROR: '{}' does not seem to be a valid .xvg file".format(xvg_file))
    sys.exit(1)

x = np.array(x) * 0.001
y = np.array(y)**2 / float(sys.argv[2])

if len(y) > 500:
    marker = ""
else:
    marker = "+"

plt.plot(x, y, marker=marker)
plt.xlabel(xaxis_title)
plt.ylabel(yaxis_title)
plt.title(title)
plt.tight_layout()

plt.savefig("{}".format(png_file))
print("INFO: Area per lipid saved to '{}'".format(png_file))