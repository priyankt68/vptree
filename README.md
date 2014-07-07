vptree
======

different implementations of vptree


Introduction
------------

The way a VP tree stores data can be represented by a circle. First, understand that each node of this tree contains 
an input point and a radius. All the left children of a given node are the points inside the circle and all the right
children of a given node are outside of the circle. The tree itself does not need to know any other information about 
what is being stored. All it needs is the distance function that satisfies the properties of the metric space.
The left children are all located inside the circle and the right children are located outside the circle.
