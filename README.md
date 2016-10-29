# SVPGen_and_Rad

This repository contains the codes I used to generate the results for my Master's thesis available [here][1], and a similar journal publication available [here][2]. They produce geometries such as the one seen below:

![Sample Porous Block](test/porous_block.png?raw=true "Sample Porous Block")

The python script is for use with [YADE DEM][3] software to generate periodic domains of spherical-void-phase Representative Elemental Volumes (REVs). It exports text files containing lists of sphere locations and their associated radii.

The porous block macro file is actually a [Solidworks][4] `.swp` macro file. Just copy-paste the contents into a new Solidworks macro.

Further details on usage and issues are given in the comments of each of the scripts, and the test folder contains example files which can be used to test the solidworks macro. It should produce the model shown in the above image.

I have applied MIT Licencses to both of the scripts.

[1]: http://ir.lib.uwo.ca/etd/2506/
[2]: http://www.sciencedirect.com/science/article/pii/S0017931014009028
[3]: https://yade-dem.org/doc/index.html
[4]: http://www.solidworks.com/