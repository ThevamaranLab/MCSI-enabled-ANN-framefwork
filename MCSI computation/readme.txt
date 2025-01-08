A simple MATLAB script for computing Hu's geometric moment invariant-based multi-component shape invariants (MCSI).


The skeleton of the code was taken from Mathworks.inc and the code was further improved by Tejas.K
Formulas for computing Hu's geometric moment invariant-based MCSI's were incorporated and further modified by Bhanu.M
August 2023

If you make any modifications to the codes - including using the skeleton to create your own version - you must rename the file. Do not distribute modified versions with the same file name.


Step 1: Have the set of images in the same path as the code. Images set should have the multi-component image, individual component images. (counting the number of components can be done separately or incorporated to the code)
Images are recommended to be grey scale. Works with RGB as well, but analysis fails with images having too many colors/ color gradient.
Step 2: Make sure centerOFMass.m and central_moments.m are in the same path.
Step 3: Counts of each component can be incorporated in to the count_vec (line 106) in order or as mentioned earlier can be coded together with this script.
 

You can also include the formulas for higher order MCSIs. 