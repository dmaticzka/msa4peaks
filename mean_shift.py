#! /usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
<<<<<<< HEAD
import doctest

information = """
Calculates the mean_shift/peaks for the given dataset.

Input:
* bed6 file containing dataset
* bandwidth (an integer value) to calculate the cluster radius/range
* Output file name (optional)
* Output_Graph file name (optional)

Output:
* bed6 file with new_centroids/peaks for the given data set. default name "output.bed".
* graph file, default name "graph.png".

Example:
- read file input.bed and write new peaks/centroids into output.bed
- mean_shift.py input.bed 200 -o data_output.bed -g output-graph.png
=======
information = """
Calculates the mean_shift/peaks for the given dataset.
Input:
* bed6 file containing dataset
* bandwidth (an integer value) to calculate the cluster radius/range
Output:
* bed6 file with new_centroids/peaks for the given data set
Example:
- read file input.bed and write new peaks/centroids into output.bed
- mean_shift.py input.bed 200
>>>>>>> bdda4a3ba7e5203d191dc38c64cde3903f779c9d
"""


class Mean_Shift():
<<<<<<< HEAD
    """ calcualtes the mean shift
    >>> a = Mean_Shift(3)
    >>> x = [1, 2, 3, 4, 5, 9, 10, 11, 12, 13, 20, 21, 22, 23, 24]
    >>> a.peak(x)
    {0: array(3.0), 1: array(11.0), 2: array(22.0)}
    """

=======
>>>>>>> bdda4a3ba7e5203d191dc38c64cde3903f779c9d
    def __init__(self, bandwidth):
        self.bandwidth = bandwidth

    def peak(self, data):
        # initializing an empty dictionary for centroids
        centroids = {}
<<<<<<< HEAD
=======

>>>>>>> bdda4a3ba7e5203d191dc38c64cde3903f779c9d
        # setting the initial centroids
        for i in range(len(data)):
            # location and data at the centroids; key and value in a dictionary
            centroids[i] = data[i]

        while True:
            # new found centroids are being saved here 'new_centroids'
            new_centroids = []
            # here we are going to cycle through our known centrods
            for i in centroids:
                # Features list will have all the features/datapoints in our bandwidth
                Features = []
                # storing the value at each location i.e. data we have at that location
                centroid = centroids[i]
<<<<<<< HEAD
=======

>>>>>>> bdda4a3ba7e5203d191dc38c64cde3903f779c9d
                # Now we are going to iterate through the data
                # and decide whether the features/data are within that bandwidth or not
                for featuresset in data:
                    # print("second loop related to dictionary containing featureset", featuresset)
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
<<<<<<< HEAD
=======

>>>>>>> bdda4a3ba7e5203d191dc38c64cde3903f779c9d
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
                if convergence==False:
                    break
            # break out of the while loop
            if convergence==True:
                break
        # finally the centroids are reset
<<<<<<< HEAD
        # self.centroids = centroids
        return centroids
=======
        self.centroids = centroids
>>>>>>> bdda4a3ba7e5203d191dc38c64cde3903f779c9d

def main():
    parser = argparse.ArgumentParser(description=information,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
<<<<<<< HEAD
    parser.add_argument("dataFile", type=str, help='Enter the data file name')
    parser.add_argument("Bandwidth", type=int, help='Enter the bandwidth')
    parser.add_argument("-o", "--output", default='output.bed', help='Enter the output data file name')
    parser.add_argument("-g", "--graph", default="graph.png", help='Enter the output graph file name')

    Pars = parser.parse_args()

    clf = Mean_Shift(Pars.Bandwidth)

=======
    parser.add_argument("dataFile", type=str)
    parser.add_argument("Bandwidth", type=int)
    Pars = parser.parse_args()

    clf = Mean_Shift(Pars.Bandwidth)
>>>>>>> bdda4a3ba7e5203d191dc38c64cde3903f779c9d
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
<<<<<<< HEAD

    X = np.array(X)
    centroids = clf.peak(X)
    print("No. of Clusters : ", len(centroids))
    plt.plot(X, np.zeros_like(X), 'ys')

    with open(Pars.output, "w") as wr:
=======
    X = np.array(X)
    clf.peak(X)
    centroids = clf.centroids
    print("No. of Clusters : ", len(centroids))

    with open("output.bed", "w") as wr:
>>>>>>> bdda4a3ba7e5203d191dc38c64cde3903f779c9d
        for c in centroids:
            end_pos.append((centroids[c] + 1))
            cl = ("%s\t%s\t%s\t%s\t%s\t%s" % (chromosome, int(centroids[c]), int(end_pos[c]), name, score, strand))
            wr.write(cl + '\n')
<<<<<<< HEAD
            plt.plot(centroids[c], np.zeros_like(centroids[c]), 'X')
        plt.savefig(Pars.graph)

if __name__ == "__main__":
    doctest.testmod()
=======
    # plotting the output
    plt.plot(X, np.zeros_like(X), 'ys')
    for c in centroids:
        plt.plot(centroids[c], np.zeros_like(centroids[c]), 'X')
    plt.show()

if __name__ == "__main__":
>>>>>>> bdda4a3ba7e5203d191dc38c64cde3903f779c9d
    main()
