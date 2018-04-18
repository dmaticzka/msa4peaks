#!/usr/bin/env python

import argparse
import numpy as np
from sklearn.cluster import MeanShift
import multiprocessing
from itertools import zip_longest
from pybedtools import BedTool

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
        >>> bd = [['chrX', '1', '+', '1'], ['chrX', '2', '+', '1'],
        ... ['chrX', '3', '+', '1'], ['chrX', '4', '+', '1'],
        ... ['chrX', '5', '+', '1'], ['chrX', '9', '-', '9'],
        ... ['chrX', '10', '-', '9'], ['chrX', '11', '-', '9'],
        ... ['chrX', '12', '-', '9'], ['chrX', '13', '-', '9'],
        ... ['chrY', '20', '+', '20'], ['chrY', '21', '+', '20'],
        ... ['chrY', '22', '+', '20'], ['chrY', '23', '+', '20'],
        ... ['chrY', '24', '+', '20']]
        >>> sorted(ms.group_chrm(bd).items())
        [(('chrX', '+', '1'), ['1', '2', '3', '4', '5']),\
 (('chrX', '-', '9'), ['9', '10', '11', '12', '13']),\
 (('chrY', '+', '20'), ['20', '21', '22', '23', '24'])]
        """

        chrm_details = dict()

        for b in bed_file:
            chrm_name_strand = (b[0], b[2], b[3])
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
        >>> x = np.array([1, 2, 3, 4, 5, 9, 10, 11, 12, 13, 20, 21, 22,\
         23, 24])
        >>> data = ms.centroids_chrm(x)
        >>> sorted(data, key=lambda tup: tup[0])
        [array([3.]), array([11.]), array([22.])]
        """
        ms = MeanShift(bandwidth=self.bandwidth, bin_seeding=False,
                       min_bin_freq=1, cluster_all=True)

        chroms_start = np.array([chroms_start])

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

    def output_bed(self, centroids, output, chromosome):
        # ####################################################################
        # writes output to a bed file format
        # ####################################################################
        """
        writes to a bed_file
        >>> Mean_Shift(3).output_bed('1', "test.bed",
        ... [("chrX", '+', '1')]) is None
        True
        """

        # BED file fields
        name = 'X'
        score = '255'

        with open(output, "a") as wr:
            for key, chrm_start in zip(chromosome, centroids):
                chromosome, strand, _ = key

                for k in range(len(chrm_start)):
                    bed_data = ("%s\t%s\t%s\t%s\t%s\t%s" %
                                (chromosome, int(chrm_start[k]),
                                 int(chrm_start[k]) + 1, name,
                                 score, strand))
                    wr.write(bed_data + '\n')

    def call_centroid_and_output(self, chrm_details, output_file):
        # ####################################################################
        # calls to centroids_chrm method to determine centroids
        # calls to output_bed method to write to an output file
        # ####################################################################
        """
        determined centroid(s) and chromosomes information is written to file
        >>> import collections
        >>> chrm =  {('chrX', '+', '1'): ['1', '2', '3', '4', '5'],
        ... ('chrX', '-', '9'): ['9', '10', '11', '12', '13'],
        ... ('chrY', '+', '20'): ['20', '21', '22', '23', '24']}
        >>> dt = collections.OrderedDict(sorted(chrm.items()))
        >>> Mean_Shift(3).call_centroid_and_output(dt, "test.txt")
        number of estimated clusters in ('chrX', '+', '1') are 1
        number of estimated clusters in ('chrX', '-', '9') are 1
        number of estimated clusters in ('chrY', '+', '20') are 1
        """

        chrom_keys = []
        chrom_values = []

        for v in chrm_details:
            # converts retrieved chromosomes start positions list to int
            chrm_details[v] = list(map(int, chrm_details[v]))

        for k in chrm_details.keys():
            chrom_keys.append(k)

        for i in chrm_details.values():
            chrom_values.append(i)

        pool = multiprocessing.Pool(2)
        centroids = pool.map(self.centroids_chrm, chrom_values)
        for c in range(len(centroids)):
            print("number of estimated clusters in %s are %s" %
                  (chrom_keys[c], len(centroids[c])))
        self.output_bed(centroids, output_file, chrom_keys)


if __name__ == "__main__":
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

    # segregating bed file for independent data chunks
    input_file = BedTool(args.data_file)
    merged_data = BedTool(input_file.sort().
                          merge(s=True, c=(2, 6), d=500, o=('max', 'first')))
    intersected_data = merged_data.\
        intersect(s=True, loj=True, a=input_file, b=merged_data)

    # retrieves chromosome name, start and strand from input bed file
    bed_data = []
    for data in intersected_data:
        bed_col = str(data).strip().split()
        bed_data.append([bed_col[0], bed_col[1], bed_col[5], bed_col[7]])

    # groups data based on chromosome and strand
    chrm_details = ms.group_chrm(bed_data)

    detail_1, detail_2 = zip(*zip_longest(*[iter(chrm_details.items())] * 2))
    chrm_details1 = dict(item for item in detail_1 if item is not None)
    chrm_details2 = dict(item for item in detail_2 if item is not None)

    # calls method for centroids computation and writes output to a file
    # ms.call_centroid_and_output(chrm_details, output_file)
    proc_1 = multiprocessing.Process(target=ms.call_centroid_and_output,
                                     args=(chrm_details1, output_file))
    proc_2 = multiprocessing.Process(target=ms.call_centroid_and_output,
                                     args=(chrm_details2, output_file))
    proc_1.start()
    proc_2.start()

    proc_1.join()
    proc_2.join()

    # sorting and saving as bedfile
    BedTool(output_file).sort().saveas(args.output_file)
