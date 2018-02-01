#! /usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np

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

    def __init__(self, bandwidth):
        self.bandwidth = bandwidth

    def peak(self, data):
        # initializing an empty dictionary for centroids
        centroids = {}
        # setting the initial centroids
        for i in range(len(data)):
            # location and data at the centroids; key and value in a dictionary
            centroids[i] = data[i]

        while True:
            # new found centroids are being saved here 'new_centroids'
            new_centroids = []
            # here we are going to cycle through our known centrods
            for i in centroids:
                # Features list will have all the features/datapoints
                # in our bandwidth
                Features = []
                # storing the value at each location
                # i.e. data we have at that location
                centroid = centroids[i]
                # Now we are going to iterate through the data
                # and decide whether the features/data are within
                # that bandwidth or not
                for featuresset in data:
                    # calculates what will be in one cluster
                    if abs(featuresset - centroid) < self.bandwidth:
                        Features.append(featuresset)
                # provides the mean as new_centroid
                new_centroid = np.mean(Features)
                new_centroids.append(new_centroid)
            # calculating unique centroids we have got
            unique_cent = sorted(set(new_centroids))
            # copying the centroids dictionary without taking the attributes
            prev_centroids = dict(centroids)
            # new centroids going to define as new dictionary
            # and defining when convergence is reached
            centroids = {}
            for i in range(len(unique_cent)):
                centroids[i] = np.array(unique_cent[i])
            convergence = True

            for i in centroids:
                if not np.array_equal(centroids[i], prev_centroids[i]):
                    convergence = False
                # break out of the for loop
                if not convergence:
                    break
            # break out of the while loop
            if convergence:
                break
        # finally the centroids are reset
        return centroids


def main():
    parser = argparse.ArgumentParser(description=information,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument("dataFile", type=str, help='Enter the data file name')
    parser.add_argument("Bandwidth", type=int, help='Enter the bandwidth')
    parser.add_argument("-o", "--output", default='output.bed',
                        help='Enter the output data file name')
    parser.add_argument("-g", "--graph", default="graph.png",
                        help='Enter the output graph file name')

    Pars = parser.parse_args()

    clf = Mean_Shift(Pars.Bandwidth)

    X = []
    end_pos = []
    chromosome = 'chrX'
    name = 'X'
    score = '255'
    strand = '+'

    with open(Pars.dataFile, "r") as f:
        for line in f:
            line = line.strip().split()
            if line[0] == "chrX":
                X.append(int(line[1]))
    X = np.array(X)
    centroids = clf.peak(X)
    print("No. of Clusters : ", len(centroids))

    plt.plot(X, np.zeros_like(X), 'ys')
    with open(Pars.output, "w") as wr:
        for c in centroids:
            end_pos.append((centroids[c] + 1))
            cl = ("%s\t%s\t%s\t%s\t%s\t%s" % (chromosome, int(centroids[c]),
                                              int(end_pos[c]),
                                              name, score, strand))
            wr.write(cl + '\n')
            plt.plot(centroids[c], np.zeros_like(centroids[c]), 'X')
        plt.savefig(Pars.graph)


if __name__ == "__main__":
    main()
