# The Method QuMinS for Labeling and Summarization
    
<p>
    The algorithm QuMinS focuses on two distinct data mining tasks - the tasks of labeling and summarizing large sets of moderate-to-high dimensionality data, such as features extracted from Gigabytes of complex data elements. Specifically, QuMinS is a fast and scalable solution to two problems (a) low-labor labeling - given a large collection of data objects, very few of which are labeled with keywords, find the most suitable labels for the remaining ones, and (b) mining and attention routing - in the same setting, find clusters, the top-NO outlier objects, and the top-NR representative objects. The algorithm is fast and it scales linearly with the data size, besides working even with tiny initial label sets.
</p>
<p>
    <b>Remark:</b> A first implementation of QuMinS was initially named as the method QMAS (after Querying, Mining And Summarizing Multi-modal Databases) in an earlier <a href="http://doi.ieeecomputersociety.org/10.1109/ICDM.2010.150"><b>Conference Publication</b></a> of this work. Later, it was renamed to QuMinS for clarity, since several improvements on the initial implementation were included into a <a href="http://dx.doi.org/10.1016/j.ins.2013.11.013"><b>Journal Publication</b></a>.
</p>
<p>
    <a href="http://doi.ieeecomputersociety.org/10.1109/ICDM.2010.150"><b>Conference Publication</b></a><br/>
    <a href="http://dx.doi.org/10.1016/j.ins.2013.11.013"><b>Journal Publication</b></a><br/>
</p>

<h3>QMAS - basic information</h3>

Please run 'make demo' to see a demo.

Input:
  feature.input -- the input image features adopted from MrCC code. Each line corresponds to an image tile.
  location.input -- each line contains the location info for an image tile, it is a triple of (index, x, y). Two tiles a, b are adjacent if a.index = b.index (they belong to the same satellite snapshot), and |a.x-b.x| + |a.y-b.y| = 1.
  label.input -- each line contains a pair of (image-id, label-id). All such id's should be positive integers
  query.input -- each line contains an image-id which we want to label
Output:
  cluster.output -- the output clustering results, each line contains a pair of (image-id, cluster-id). Outliers are not included.
  label.output -- the output labeling results, each line contains a triple of (image-id, label-id, label-score).
  label.detail -- the output labeling details, each line contains image-id, followed by scores of each label.
Source:
  gcap.m -- the MATLAB source which compiles MrCC source, execute the clustering algorithm (soft clustering), and performs RWR.
