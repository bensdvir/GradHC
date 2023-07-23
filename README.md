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
After the clustering successfully finishes, the following output is printed to the standard output in the following structure:
```
Total time: 90.47252297401428
Total Clusters: 5016, Singles: 290
Metric Accrcy:
0.6: 0.9760765550239234
0.7: 0.9726874003189793
0.8: 0.9712918660287081
0.9: 0.9667065390749602
0.95: 0.9593301435406698
0.99: 0.9545454545454546
1.0: 0.9545454545454546
Metric FalsePos:
Total num. of strands: 75009
(FP) False Positives: 0
(TN) True Negatives: 75009
(FN) False Negatives: 3081
(TP) True Positives: 71928
(TS) Threat Score / (CSI) Critical Success Index: 0.958924929008519
Metric NMI:
0.9294847141008635
Metric RandIndex:
0.5754018335777495
```
Few metrics are used in order to evaulate the result.
* Accuracy: defined in [1]
* FalsePos: refer to [Sensitivity and specificity](https://en.wikipedia.org/wiki/Sensitivity_and_specificity)
* NMI: refer to [Normalized Mutual Info](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.normalized_mutual_info_score.html#sklearn.metrics.normalized_mutual_info_score)
* RandIndex: refer to [Adjusted Rand index](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.adjusted_rand_score.html)
* Purity: can be added by inserting 'Purity" to the array in self.results_str(), but take into account that it requires a long run time. Uses implementation from [jhumigas](https://gist.github.com/jhumigas/010473a456462106a3720ca953b2c4e2).

## Export
The result can exported as an evyat.txt file, if the 'export' argument of the LSHBasedCluster's init is set to True.

## How-To
```
usage: lsh_based_clustering.py [-h] -e EVYAT
lsh_based_clustering.py: error: the following arguments are required: -e/--evyat
```
The log is printed to the standard output as a default, can be turned off manually (set 'enable' argument to False in the 'info' function). Pay attention in case you redirect the output to a file, since longer runs may result with a heavier log file.

## References
[1] C. Rashtchian et al. “Clustering billions of reads for DNA data storage,” Advances in Neural Information Processing Systems, vol. 30, 2017.\
[2] P. L. Antkowiak, J. Lietard, M. Z. Darestani et al. ”Low cost DNA data storage using photolithographic synthesis and advanced information reconstruction and error correction,” Nature Communications, vol. 11, 2020.

## Note - DNAsimulator
(Refers to [DNASimulator](https://github.com/gadihh/DNASimulator))\
For possible future usages, if the tool is to be assimilated in the DNAsimulator, the needed files are attached in 'DNAsimulator_files' file. Replace them with exisitng ones. \
Changes to the original DNAsimulator code:
* Modified:
	- \DNASimulator\dnaSimulator\app.py
	- \DNASimulator\dnaSimulator\dnaSimulator_ui2.py
* Added:
	- \DNASimulator\dnaSimulator\lsh_based_clustering.py

Notes:
- Multiprocessing is not supported in the Windows version
- Progress bar is currently unreponsive.
