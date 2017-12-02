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


if len(sys.argv) == 1:
    print("USAGE: xvg2png.py file1.xvg [file2.xvg...]")
    sys.exit(0)


for fname in sys.argv[1:]:
    title = "Title"
    xaxis_title = "X axis"
    yaxis_title = "Y axis"
    x = []
    ys = []
    labels = []

    new_fname = "{}.png".format(os.path.splitext(fname)[0])

    try:
        with open(fname) as fp:
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

                    for i, val in enumerate(line[1:]):
                        try:
                            ys[i].append(val)
                        except IndexError:
                            ys.append([val, ])
    except IOError:
        print("ERROR: Could not read '{}' -> skipped".format(fname))
        continue
    except (ValueError, TypeError):
        print("ERROR: '{}' does not seem to be a valid .xvg file -> skipped".format(fname))
        continue

    if len(labels) != 0:
        if len(ys) != len(labels):
            print("ERROR: Incomplete legend in '{}' -> skipped".format(fname))
            continue

    x = np.array(x)
    if len(x) > 500:
        marker = ""
    else:
        marker = "+"

    for i, vals in enumerate(ys):
        if len(labels) == 0:
            label = "value #%i" % (i + 1)
        else:
            label = labels[i]

        plt.plot(x, vals, label=label, marker=marker)

    if len(labels) != 0:
        plt.legend(loc="best")

    plt.xlabel(xaxis_title)
    plt.ylabel(yaxis_title)
    plt.title(title)
    plt.tight_layout()

    plt.savefig("{}".format(new_fname))
    print("INFO: '{}' succesfully converted to '{}'".format(fname, new_fname))