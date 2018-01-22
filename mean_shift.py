#! /usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
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
"""

# X = np.array([1, 2, 3, 4, 5, 9, 10, 11, 12, 13, 20, 21, 22, 23, 24])
"""
X = []
with open('test_2.bed', "r") as f:
    for line in f:
        line = line.strip().split()
        if line[0]=="chrX":
            X.append(int(line[1]))
X = np.array(X)
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
                # Features list will have all the features/datapoints in our bandwidth
                Features = []
                # storing the value at each location i.e. data we have at that location
                centroid = centroids[i]

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
        self.centroids = centroids

def main():
    parser = argparse.ArgumentParser(description=information,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("dataFile", type=str)
    parser.add_argument("Bandwidth", type=int)
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
    clf.peak(X)
    centroids = clf.centroids
    print("No. of Clusters : ", len(centroids))

    with open("output.bed", "w") as wr:
        for c in centroids:
            end_pos.append((centroids[c] + 1))
            cl = ("%s\t%s\t%s\t%s\t%s\t%s" % (chromosome, int(centroids[c]), int(end_pos[c]), name, score, strand))
            wr.write(cl + '\n')
    # plotting the output
    plt.plot(X, np.zeros_like(X), 'ys')
    for c in centroids:
        plt.plot(centroids[c], np.zeros_like(centroids[c]), 'X')
    plt.show()

if __name__ == "__main__":
    main()
