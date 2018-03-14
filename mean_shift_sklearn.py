#!/usr/bin/env python

import argparse
<<<<<<< HEAD
=======
import matplotlib.pyplot as plt
from itertools import cycle
>>>>>>> 9208fec6bb243aa8f0e060afcd3d13ae54c49468
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
<<<<<<< HEAD

Example:
- reads file input.bed and writes output.bedwith new peaks/centroids
- mean_shift.py input.bed -b 3 -o data_output.bed 
=======
* graph file, default name "graph.png".

Example:
- reads file input.bed and writes output.bedwith new peaks/centroids
- mean_shift.py input.bed 200 -o data_output.bed -g output-graph.png
>>>>>>> 9208fec6bb243aa8f0e060afcd3d13ae54c49468

"""


class Mean_Shift():

    """
    applying mean shift algorithm to an array of integers
    >>> a = Mean_Shift(3)
    >>> x = np.array([[1, 2, 3, 4, 5, 9, 10, 11, 12, 13, 20, 21, 22, 23, 24]])
    >>> a.peak(x)
    array([[ 3.],
           [11.],
           [22.]])
    """

    def __init__(self, bandwidth):
        self.bandwidth = bandwidth

    def peak(self, chroms_start):
        # ####################################################################
        # Computes clustering with MeanShift
        # ####################################################################

<<<<<<< HEAD
        ms = MeanShift(bandwidth=self.bandwidth, bin_seeding=True,
                       min_bin_freq=1, cluster_all=True)
=======
        ms = MeanShift(bandwidth=self.bandwidth, bin_seeding=True)
>>>>>>> 9208fec6bb243aa8f0e060afcd3d13ae54c49468
        # Performs clustering on input(ndarray, shape (n_samples, n_features))
        ms.fit(chroms_start.transpose())
        # retrieving clusters centers (array, [n_clusters, n_features])
        cluster_centers = ms.cluster_centers_

        return cluster_centers

<<<<<<< HEAD
    def output_bed(self, centroids, output, chromosome, strand):
=======
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
>>>>>>> 9208fec6bb243aa8f0e060afcd3d13ae54c49468
        # ####################################################################
        # writes output to a bed file
        # ####################################################################

        # BED file fields
<<<<<<< HEAD
        name = 'X'
        score = '255'

        with open(output, "a") as wr:
            for cent in centroids:
=======
        chromosome = 'chrX'
        name = 'X'
        score = '255'
        strand = '+'
        centroids_sorted = np.sort(centroids, axis=None, kind='quicksort')

        with open(output, "w") as wr:
            for cent in centroids_sorted:
>>>>>>> 9208fec6bb243aa8f0e060afcd3d13ae54c49468
                bed_data = ("%s\t%s\t%s\t%s\t%s\t%s"
                            % (chromosome, int(cent), int(cent + 1),
                               name, score, strand))
                wr.write(bed_data + '\n')

<<<<<<< HEAD
    def segregate_data(self, bed_file):
        # ####################################################################
        # writes output to a bed file
        # ####################################################################

        dict_chr = dict()

        for j in bed_file:
            chrm_name_strand = (j[0], j[2])
            chrm_pos = j[1]

            if chrm_name_strand in dict_chr:
                # append to an existing key
                dict_chr[chrm_name_strand].append(chrm_pos)
            else:
                # create a new array in this slot
                dict_chr[chrm_name_strand] = [chrm_pos]

        return dict_chr

    def clear_output_file(self, output_file):
        # ####################################################################
        # delete contents of output_file if already exist
        # ####################################################################

        with open(output_file, 'w'):
            pass

    def computation(self, chr_details, output_file):
        # ####################################################################
        # computes peaks and write it to an output file
        # ####################################################################

        for d in chr_details:
            # converts retrieved chromosomes start positions list to int
            chr_details[d] = list(map(int, chr_details[d]))
            # converts them into an ndarray.
            chromosomes_start = np.array([chr_details[d]])
            # retrieves centroids/peaks using mean shift algorithm
            centroids = self.peak(chromosomes_start)
            chromosomes, strand = d[0], d[1]
            print("number of estimated clusters in (%s,%s) are %d :"
                  % (chromosomes, strand, len(centroids)))
            self.output_bed(centroids, output_file, chromosomes, strand)

=======
>>>>>>> 9208fec6bb243aa8f0e060afcd3d13ae54c49468

def main():
    parser = argparse.ArgumentParser(description=information,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument("data_file", type=str, help='Enter the bed file name')
<<<<<<< HEAD
    parser.add_argument("-b", "--bandwidth", type=int, default=5,
                        help='Enter the bandwidth')
    parser.add_argument("-o", "--output_file", default='output.bed',
                        help='Enter the output bed file name')

    args = parser.parse_args()
    output_file = args.output_file
    ms = Mean_Shift(args.bandwidth)
    print("Bandwidth : ", args.bandwidth)
    # delete the contents of output file if already exist
    ms.clear_output_file(output_file)
    # retrieves chromosome name, start and strand from input bed file
    bed_file = np.loadtxt(args.data_file, dtype=str, usecols=(0, 1, 5))
    # segregates data based on chromosome and strand for mean_shift computation
    chr_details = ms.segregate_data(bed_file)
    # computes the peaks and writes to a file
    ms.computation(chr_details, output_file)
=======
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
>>>>>>> 9208fec6bb243aa8f0e060afcd3d13ae54c49468


if __name__ == "__main__":
    main()
