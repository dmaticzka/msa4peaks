#!/usr/bin/env python

import argparse
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

Example:
- reads file input.bed and writes output.bedwith new peaks/centroids
- mean_shift.py input.bed -b 3 -o data_output.bed

"""


class Mean_Shift():

    """
    bandwidth to perform clustering
    >>> ms = Mean_Shift(3)
    >>> ms.bandwidth
    3
    """

    def __init__(self, bandwidth):
        self.bandwidth = bandwidth

    def group_chrm(self, bed_file):
        # ####################################################################
        # group chromosomes with respect to name and strand
        # ####################################################################

        """
        groups bed_file data(ndarray) w.r.t chromosome and strand for centroids
        >>> ms = Mean_Shift(3) # bandwidth:  3
        >>> bd = [['chrX', '1', '+'], ['chrX', '2', '+'],
        ... ['chrX', '3', '+'], ['chrX', '4', '+'],['chrX', '5', '+'],
        ... ['chrX', '9', '-'], ['chrX', '10', '-'], ['chrX', '11', '-'],
        ... ['chrX', '12', '-'], ['chrX', '13', '-'], ['chrY', '20', '+'],
        ... ['chrY', '21', '+'], ['chrY', '22', '+'], ['chrY', '23', '+'],
        ... ['chrY', '24', '+']]
        >>> sorted(ms.group_chrm(bd).items())
        [(('chrX', '+'), ['1', '2', '3', '4', '5']),\
 (('chrX', '-'), ['9', '10', '11', '12', '13']),\
 (('chrY', '+'), ['20', '21', '22', '23', '24'])]
        """
        chrm_details = dict()

        for b in bed_file:
            chrm_name_strand = (b[0], b[2])
            chrm_pos = b[1]

            if chrm_name_strand in chrm_details:
                # append to an existing key
                chrm_details[chrm_name_strand].append(chrm_pos)
            else:
                # create a new array in this slot
                chrm_details[chrm_name_strand] = [chrm_pos]

        return chrm_details

    def centroids_chrm(self, chroms_start):
        # ####################################################################
        # Determines chromosomes centres with MeanShift
        # ####################################################################

        """
        applying mean shift algorithm to compute centres w.r.t bandwidth(3)
        >>> ms = Mean_Shift(3)
        >>> x = np.array([[1, 2, 3, 4, 5, 9, 10, 11, 12, 13, 20, 21, 22,\
         23, 24]])
        >>> data = ms.centroids_chrm(x)
        >>> sorted(data, key=lambda tup: tup[0])
        [array([3.]), array([11.]), array([22.])]
        """
        ms = MeanShift(bandwidth=self.bandwidth, bin_seeding=True,
                       min_bin_freq=1, cluster_all=True)
        # Performs clustering on input(ndarray, shape (n_samples, n_features))
        ms.fit(chroms_start.transpose())
        # retrieving clusters centers (array, [n_clusters, n_features])
        centroids = ms.cluster_centers_

        return centroids

    def clear_output_file(self, output_file):
        # ####################################################################
        # delete contents of output_file already exists
        # ####################################################################

        with open(output_file, 'w'):
            pass

    def output_bed(self, centroids, output, chromosome, strand):
        # ####################################################################
        # writes output to a bed file format
        # ####################################################################
        """
        writes to a bed_file
        >>> Mean_Shift(3).output_bed('1', "test.bed", "chrX", "+") is None
        True
        """

        # BED file fields
        name = 'X'
        score = '255'
        # print(centroids)
        with open(output, "a") as wr:
            for chrm_start in centroids:
                bed_data = ("%s\t%s\t%s\t%s\t%s\t%s" %
                            (chromosome, int(chrm_start), int(chrm_start) + 1,
                             name, score, strand))
                wr.write(bed_data + '\n')

    def call_centroid_and_output(self, chrm_details, output_file):
        # ####################################################################
        # calls to centroids_chrm method to determine centroids
        # calls to output_bed method to write to an output file
        # ####################################################################

        """
        determined centroid(s) and chromosomes information is written to file
        >>> import collections
        >>> chrm =  {('chrX', '+'): ['1', '2', '3', '4', '5'],
        ... ('chrX', '-'): ['9', '10', '11', '12', '13'],
        ... ('chrY', '+'): ['20', '21', '22', '23', '24']}
        >>> dt = collections.OrderedDict(sorted(chrm.items()))
        >>> Mean_Shift(3).call_centroid_and_output(dt, "test.txt")
        number of estimated clusters in (chrX,+) are 1
        number of estimated clusters in (chrX,-) are 1
        number of estimated clusters in (chrY,+) are 1
        """

        for d in chrm_details:
            # converts retrieved chromosomes start positions list to int
            chrm_details[d] = list(map(int, chrm_details[d]))
            # converts them into an ndarray.
            chromosomes_start = np.array([chrm_details[d]])
            # retrieves centroids/peaks using mean shift algorithm
            centroids = self.centroids_chrm(chromosomes_start)
            chromosomes, strand = d[0], d[1]
            print("number of estimated clusters in (%s,%s) are %d"
                  % (chromosomes, strand, len(centroids)))
            self.output_bed(centroids, output_file, chromosomes, strand)


def main():
    parser = argparse.ArgumentParser(description=information,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument("data_file", type=str, help='Enter the bed file name')
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
    # groups data based on chromosome and strand
    chrm_details = ms.group_chrm(bed_file)
    # calls method for centroids computation and writes output to a file
    ms.call_centroid_and_output(chrm_details, output_file)


if __name__ == "__main__":
    main()
