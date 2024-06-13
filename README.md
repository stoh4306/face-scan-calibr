## Library for stereo calibration and 3D reconstruction

The repository contains the header and the source files to create a static library for stereo calibration and 3D reconstruction.

* OS & buld environment : Windows10, Visual stuido 2017
* dependencies : [OpenCV 4.1.0](https://github.com/opencv/opencv/tree/4.1.0), [Eigen 3.3.5](http://bitbucket.org/eigen/eigen/get/3.3.5.tar.bz2)

---

## Test the library

The folder "test_calibr" has the main source file and the data of stereo check board images to compute the calibration of stereo camera using the calibration library.
The execution binary creates the following files:

- corner points on the check images
- camera pose and intrinsic camera parameters
- reconstructed corner points in 3D by triangulation
