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
    calcualtes the mean shift
    >>> a = Mean_Shift(3)
    >>> x = np.array([[1, 2, 3, 4, 5, 9, 10, 11, 12, 13, 20, 21, 22, 23, 24]])
    >>> a.peak(x)
    (array([[11.],
           [ 3.],
           [22.]]), 3)
    """
    def __init__(self, bandwidth):
        self.bandwidth = bandwidth
        self.labels = []

    def peak(self, dt):
        # #############################################################################
        # Compute clustering with MeanShift
        ms = MeanShift(bandwidth=self.bandwidth, bin_seeding=True)
        ms.fit(dt.transpose())
        self.labels = ms.labels_
        cluster_centers = ms.cluster_centers_
        labels_unique = np.unique(self.labels)
        n_clusters_ = len(labels_unique)

        return cluster_centers, n_clusters_

    def plot_graph(self, cluster_centers, n_clusters_, X, graph):
        # #############################################################################
        # Plot result
        plt.figure(1)
        plt.clf()
        colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')

        for k, col in zip(range(n_clusters_), colors):
            self.labels == k
            cluster_center = cluster_centers[k]
            plt.plot(X, X, col + '.', cluster_center[0],
                     cluster_center[0], 'o',
                     markerfacecolor=col,
                     markeredgecolor='k', markersize=14)
        plt.title('Estimated number of clusters: %d' % n_clusters_)
        # plt.show()
        plt.savefig(graph)

    def output_file(self, centroids, output):

        chromosome = 'chrX'
        name = 'X'
        score = '255'
        strand = '+'
        centroid = np.sort(centroids, axis=None, kind='quicksort')

        with open(output, "w") as wr:
            for c in centroid:
                cl = ("%s\t%s\t%s\t%s\t%s\t%s" % (chromosome,
                                                  int(c),
                                                  int(c + 1),
                                                  name, score, strand))
                wr.write(cl + '\n')


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
    dt = np.array([np.loadtxt(Pars.dataFile, usecols=[1])])
    cntr, no_clust = clf.peak(dt)
    print("number of estimated clusters : %d" % no_clust)
    clf.plot_graph(cntr, no_clust, dt, Pars.graph)
    clf.output_file(cntr, Pars.output)

if __name__ == "__main__":
    main()
