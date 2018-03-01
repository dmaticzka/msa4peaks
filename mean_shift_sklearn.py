#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
from itertools import cycle
import numpy as np
from sklearn.cluster import MeanShift

information = """
Calculates the mean_shift/peaks for the given dataset.

Input:
* bed6 file containing dataset
* bandwidth (an integer value) to calculate the cluster radius/range
* Output file name (optional)
* Output_Graph file name (optional)

Output:
* bed6 file with new_centroids/peaks for the given data set.
 default name "output.bed".
* graph file, default name "graph.png".

Example:
- read file input.bed and write new peaks/centroids into output.bed
- mean_shift.py input.bed 200 -o data_output.bed -g output-graph.png

"""


class Mean_Shift():

    """
    applying mean shift algorithm to an array of integers
    >>> a = Mean_Shift(3)
    >>> x = np.array([[1, 2, 3, 4, 5, 9, 10, 11, 12, 13, 20, 21, 22, 23, 24]])
    >>> a.peak(x)
    array([[11.],
           [ 3.],
           [22.]])
    """

    def __init__(self, bandwidth):
        self.bandwidth = bandwidth

    def peak(self, chroms_start):
        # ####################################################################
        # Computes clustering with MeanShift
        # ####################################################################

        ms = MeanShift(bandwidth=self.bandwidth, bin_seeding=True)
        # Performs clustering on input(ndarray, shape (n_samples, n_features))
        ms.fit(chroms_start.transpose())
        # retrieving clusters centers (array, [n_clusters, n_features])
        cluster_centers = ms.cluster_centers_

        return cluster_centers

    def output_graph(self, cluster_centers, chroms_start, graph):
        # ####################################################################
        # Plots result and saves into a file
        # ####################################################################

        # creates a new figure
        plt.figure()
        # clears the current figure
        plt.clf()
        # cycle through colors to assign different color
        # to consecutive centroids
        colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')

        for k, color in zip(range(len(cluster_centers)), colors):
            cluster_center = cluster_centers[k]
            plt.plot(chroms_start, chroms_start, color + '.',
                     cluster_center[0], cluster_center[0], 'o',
                     markerfacecolor=color,
                     markeredgecolor='k', markersize=16)
        plt.title('Estimated number of clusters: %d' % len(cluster_centers))
        # plt.show()
        plt.savefig(graph)

    def output_bed(self, centroids, output):
        # ####################################################################
        # writes output to a bed file
        # ####################################################################

        # BED file fields
        chromosome = 'chrX'
        name = 'X'
        score = '255'
        strand = '+'
        centroids_sorted = np.sort(centroids, axis=None, kind='quicksort')

        with open(output, "w") as wr:
            for cent in centroids_sorted:
                bed_data = ("%s\t%s\t%s\t%s\t%s\t%s"
                            % (chromosome, int(cent), int(cent + 1),
                               name, score, strand))
                wr.write(bed_data + '\n')


def main():
    parser = argparse.ArgumentParser(description=information,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument("data_file", type=str, help='Enter the bed file name')
    parser.add_argument("bandwidth", type=int, help='Enter the bandwidth')
    parser.add_argument("-o", "--output_file", default='output.bed',
                        help='Enter the output bed file name')
    parser.add_argument("-g", "--graph_file", default="graph.png",
                        help='Enter the output graph file name')

    args = parser.parse_args()

    ms = Mean_Shift(args.bandwidth)
    # retrieves chromosomes start positions from input bed file
    # converts them into an ndarray.
    chromosomes_start = np.array([np.loadtxt(args.data_file, usecols=[1])])
    # retrieves centroids/peaks using mean shift algorithm
    centroids = ms.peak(chromosomes_start)
    print("number of estimated clusters : %d" % len(centroids))
    ms.output_graph(centroids, chromosomes_start, args.graph_file)
    ms.output_bed(centroids, args.output_file)


if __name__ == "__main__":
    main()
