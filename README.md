# GradHC: Highly Reliable Gradual Hash-based Clustering for DNA Storage Systems

## Background
An implementation of the Gradual Hash-based clustering (GradHC) algorithm for DNA storage systems. The primary strength of GradHC lies in its capability to cluster
with excellent accuracy various types of designs, including varying strand lengths, cluster sizes (including extremely small clusters), and different error ranges.

## Input
The tool receives an 'evyat.txt' text file, in the following format: 
```
<original data string #1> 
*****************************
<erroneous copies of the above string>

(two blank line)
<original data string #2> 
etc...
```
For example:
```
CGAGGTTGTGGGCTTCTGTATATGCTCCAATGAAAGTGCCAACTTTTCATAGTCTTCCCAGTGGATTCAATGACGACATCGCACACATACCGCAGTGCGGAAGGCCTAG
*****************************
CGAGGGTTATGGGCTTCTGTATATGCCTCCAATGAACAGTGCCAACTTTTCATAGTCTTCCCAGTGGATATGACGACATCGCACACATACCGCAGTGCGGAAGGCCTAG
CGAGGTTGTGGGGCTTCTGTATATGCTCCAATGAAACCAACTTTTACATAGTCTTCCCAGTGGATTAAATGACGACATCGCACACATACCGCAGTGGCGGATAGGCCTAG
CGAGGGTTGTGTGGCCTTCTGTATATGCCAATGAAAGTGCCAACTTTTCATAGTCTTCCCAGTGGATTCAATGACGACATCGCACACATACCGCAGTGCGGAAGGCCTAG


TCATCAGTGTTAAAATCTTGTGTAGGCAGACGCTTCCTGGAAAACCCGTCCTGGGTATACACAACGGTATGTACACTCTAAGAATTGGTTGCCACTGCGCACTTCTAGG
*****************************
TCATCAGTAGCTAAATCTTGTGTAGGCAGACGCTTCCTGGAAAAACCCGTCCTGGGTATACATCAACGGTATGTACACTTTACGAATTAGTTGCCACTGCGCACTTCTAGG
TCATCAGGTGTTAAAATCTTGTGGCAGGCAGACGCTTCCTGGAAAACCCGTCCTGGGTATACACAACCGGTATGTACACTCTAAGATATTGGTTGCCACTGCGCACTTCTAGG
TCATCATTGTTAAAATCTTGTAGTAGGCAGACGCTTCCTGGAAAACCCGTCCTGGGTATACACAACGGTTATGTACACTCTAAGAATATGGTTGCCACATGCGCACTTCTAGG
```
The file consists of noisy copies generated from the original data. In other words, the sequences *after* the synthesis.
For our purposes, as we aim to test the clustering result, the file is expeted to include the parition to clusters beforehand, so it can be later compared to our result, using different metrics.

## Output
After the clustering successfully finishes, the result can exported as an evyat.txt file, if the 'export' argument of the GradHCBasedCluster's init is set to True.

## How-To
```
usage: GradHC_clustering.py [-h] -e EVYAT
GradHC_clustering.py: error: the following arguments are required: -e/--evyat
```
The log is printed to the standard output as a default, can be turned off manually (set 'enable' argument to False in the 'info' function). Pay attention in case you redirect the output to a file, since longer runs may result with a heavier log file.

## References
[1] C. Rashtchian et al. “Clustering billions of reads for DNA data storage,” Advances in Neural Information Processing Systems, vol. 30, 2017.\
[2] P. L. Antkowiak, J. Lietard, M. Z. Darestani et al. ”Low cost DNA data storage using photolithographic synthesis and advanced information reconstruction and error correction,” Nature Communications, vol. 11, 2020.

