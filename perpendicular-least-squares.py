__author__ = 'Jacob Bieker'
import os, sys
import numpy
from multiprocessing import Pool


def line_solve():
    # TODO: Find the Least Squares for a line
    return 0


def plane_solve():
    # TODO: Find the least Squares for a plane
    return 0


def read_clusters():
    # TODO: Read in the files containing the cluster points
    return 0


def determine_type():
    # TODO: Determine the type (i.e. Plane or Line) to solve for
    return 0


def initial_guess(largest_cluster, type):
    if type == "line":
        return 0
    elif type == "plane":
        return 0
    # TODO: Get an inital guess from the perpendicular least squares method from the largest cluster
    return


def check_guess(cluster, type, guess):
    if type == "line":
        return 0
    elif type == "plane":
        return 0
    # TODO: check guess with another cluster


def bootstrap_cluster(cluster):
    # TODO: bootstrap a cluster to get a different distribution to check with check_guess
    return 0


def determine_uncertainty(solutions):
    # TODO: Take a list or dict of solutions and determine the uncertainty in them
    return 0

if __name__ == "__main__":
    arguments = str(sys.argv)
    print(arguments)

