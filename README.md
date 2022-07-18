# Delaunay-Triangulation
Delaunay Triangulation with Chew's Algorithm 

The input file needs to cointain a polygon with vertices in counterclockwise order currently it does not support the update on a given delaunay triangulation. 
From what I have understood from the documentation of triagulating the cavity that has been created by the constraint edge that we are given to include in the triangulation the program needs to support a DCEL as an input.
So, the cavity can be computed and the algorithm can triangulate the cavity which should be a polygon and the order is found as counterclockwise by traversing the DCEL.
I am not so sure about how the cavity can be computed in a reasonable time, without comparing every edge of the previous triangulation. Given we are dealing with a single constraint at the moment it should be reasonable to consider all the edges in the DCEL for collision with the given constraint. 
