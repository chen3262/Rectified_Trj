# Rectified_Trj

<img src ="https://github.com/chen3262/Rectified_Trj/blob/master/pic.png" width="750">

This is a C++ code to calculat velocity profiles from [GROMACS](http://www.gromacs.org/) trajectory files (.trr) trajectory files (.trr). This code used advanced algorithm to obtained velocity profiles by means of Rectified Trajectory, detailed in our paper [J. Phys. Chem. B, 2014, 118 (28), pp 8170–8178](http://pubs.acs.org/doi/abs/10.1021/jp5012523). [OpenMp](http://www.openmp.org) are implented to fulfill parallel computing in multiple-processors computers. This code requires the [XTC Libray](http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library)

## Requirements
[xdrfile](http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library)>=1.1.4 is required.

## Compilation

After compiling the xdrfile-1.1.4 library on the local machine, cd into this repository. Then:

```bash
icpc -DCPLUSPLUS -I /export/singer_group3/Sihan/Tools/xdrfile-1.1.4/_install/include/xdrfile -L /export/singer_group3/Sihan/Tools/xdrfile-1.1.4/_install/lib test.cpp -lxdrfile
```

## Testing

To test your build, do:

```bash
./a.out TRR GRO molinfo output.dens output.vx log.txt
```

## License

Copyright (C) 2017 Si-Han Chen chen.3262@osu.edu
