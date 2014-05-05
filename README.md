VP Tree
=======


Introduction
------------
The VP tree is particularly useful in dividing data in a [non-standard metric space](https://en.wikipedia.org/wiki/Metric_space#Examples_of_metric_spaces) into a
[BSP tree](https://en.wikipedia.org/wiki/Binary_space_partitioning).
Tree construction executes in O(n&nbsp;log(n)) time, and search is under certain circumstances and in the limit, O(log(n))
expected time. This makes it suitable when distance computations are expensive.


Construction
------------

 It is most efficient to simply build the tree by repeatedly partitioning the data. We build the tree from the top down from an array of items. For each node, we first choose a point at random, and then partition the list into two sets: The left children contain the points farther away than the median, and the right contains the points that are closer than the median. Then we recursively repeat this until we have run out of points. 


[![Bitdeli Badge](https://d2weczhvl823v0.cloudfront.net/priyankt68/vptree/trend.png)](https://bitdeli.com/free "Bitdeli Badge")

