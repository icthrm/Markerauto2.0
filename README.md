# Overview

Here we present a fast, robust and fully automatic software, Markerauto2, for efficient and effective cryo-electron tomography (cryo-ET) tilt series alignment. Firstly, Markerauto2 implement the fiducial marker detection based on an fast accelerated wavelet template-based detection algorithm. After detecting fiducial markers, Markerauto2 utilizes a novel two-stage algorithm to establish correspondences in electron tomography, enhancing both speed and accuracy in this crucial step.
Then, Markerauto2 employs a high-angle fiducial marker supplementation strategy designed to identify undetected markers in high-angle projections and generate more complete tracks. Finally, Markerauto2 utilizes a group-weighted parameter optimization algorithm for precise calibration of projection parameters. In comparison to other alignment software and method, Markerauto2 stands out through faster, more robust alignment while it is also fully automated.

# Install

The sections below explain how to download and install Markerauto2 on your computer.

## Prerequisites

Note that Markerauto2 depends on and uses several external programs and libraries.

- **Linux or Unix-like operating systems**
- **GCC**
  ```bash
  sudo apt install gcc
- **CMake**
    ```bash
  sudo apt install cmake
- **Opencv4.5.5 (C++ version)**
    ```bash
  sudo apt-get install build-essential cmake git libgtk2.0-dev pkg-config libavcodec-dev libavformat-dev libswscale-dev
    Download the Opencv4.5.0 installation package (https://opencv.org/releases/).
    unzip opencv-4.5.0.zip
    cd opencv-4.5.0
    mkdir build
    cd build 
    cmake -DCMAKE_BUILD_TYPE=Release -DOPENCV_GENERATE_PKGCONFIG=ON -DCMAKE_INSTALL_PREFIX=/usr/local ..
    make -j2
    sudo make install
- **Ceres-Solver2.1 (C++ version)**
    ```bash
    sudo apt-get install libeigen3-dev libatlas-base-dev libsuitesparse-dev
    git clone https://ceres-solver.googlesource.com/ceres-solver
    tar zxf ceres-solver-2.1.0.tar.gz
    mkdir ceres-bin
    cd ceres-bin
    cmake ../ceres-solver-2.1.0
    make -j3
    sudo make install
## Installtion of Markerauto2

We store the public release versions of Markerauto2 on GitHub, a site that provides code-development with version control and issue tracking through the use of git. We will not describe the use of git in general, as you will not need more than very basic features. Below we outline the few commands needed on a UNIX-system, please refer to general git descriptions and tutorials to suit your system. To get the code, you clone or download the repository. We recommend cloning, because it allows you very easily update the code when new versions are released. To do so, use the shell command-line:

- **Download**
  ```bash
  git clone https://github.com/icthrm/Markerauto2.git
- **Compilation**
    ```bash
    cd ./Markerauto2/
    mkdir build
    cd ./build
    cmake ..
    make -j2
# Explanation of parameters and useful examples

## Parameter explanation of Markerauto2

Markerauto2: command line arguments

- **==== Required options ====**

  - **--input (-i):** Specify the input file (The input file extension is .mrc or .st).
    
    Example: `--input inputfile.st` or `-i inputfile.st`

  - **--output (-o):** Specify the output file.
    
    Example: `--output outputfile.st` or `-o outputfile.st`

  - **--initangle (-a):** Specify the initial angle file.
    
    Example: `--initangle initangle.rawtlt` or `-a initangle.rawtlt`

  - **--newangle (-n):** Specify the new angle filename.
    
    Example: `--newangle newanglefile.tlt` or `-n newanglefile.tlt`

  - **--diameter (-d):** Specify the marker diameter (pixel).
    
    Example: `--diameter n` or `-d n` (n is the marker diameter. If -d -1, the software will automatically initialize the diameter)

- **==== Optional options ====**

  - **--help (-h):** Print help content, no extra argument needed.

  - **--rotationangle (-r):** Image rotation.
    
    Example: `--rotationangle 0` or `-r 0`

  - **--verbose (-v):** Display verbose information. A number greater than 1 enables this mode. If you do not provide this argument, it will run as '-v 0'.
    
    Example: `--verbose 0` or `-v 1`

  - **-t:** Enable fast mode. In the fast mode, a two-stage algorithm for the fiducial marker correspondence establishment is used. Otherwise, the software uses the global algorithm to match the fiducial markers.


## Explanation of the content in the output folder

- **'xtiltangle.txt'**

    - Record the rotation angle of each micrograph $\alpha/\pi*180$.

- **'points.mot'**

    - The final coordinates of 3D points in space.

- **'cameras.params'**

    - Save camera parameters for different images: $s, \alpha, \beta, \gamma, t_0, t_1$.

- **'BBa\_new.tlt'**

    - Optimized new angle file.

- **'BBa\_fin.xf'**

    - The file to be applied to the linear transformation of the image, calculated from the optimized camera parameters. 
      The file is processed with Autom or IMOD to perform the image transformation to obtain the aligned image sequence (see subsection \ref{Final}).

    - Six columns of data are stored. The first four columns are the values obtained by multiplying the rotation matrix by the parameter $s$ and the focal length. The last two columns are the values obtained by multiplying the translation matrix by the parameter $s$, the focal length, and the scale.

- **'fids.txt'**

    - Output the width, height, number of photos, and ratio of the micrograph.
    - Output the $i$-th micrograph (frame).
    - Output the number of fiducial markers in the $i$-th micrograph.
    - The number of fiducial fiducial markers in each micrograph and output its coordinates $(x,y)$ in the micrograph.

- **'proadd.txt'**

    - Output high-angle fiducial markers coordinates $(x, y)$ supplemented by reprojection.

- **'transmx folder'**

    - Stores affine matrices between pairs of micrographs. The number of affine matrices may not be unique.

- **'matches folder'**

    - Store the coordinates of the corresponding fiducial markers in the two matched micrographs.


## Examples

When using Markerauto2 for alignment, the necessary parameters provided by the user include the initial micrograph, the tilt angle of the micrograph, the linear transformation *.xf file storage path to align the images, and the optimized tilt angle storage path. Other parameters are optional and users can choose according to their needs. The following are serveral useful examples of input from the command line:

- **The software will automatically initialize the diameter, no fast mode, no log output:**
    ```bash
    ./markerauto2 -i BBa.st -a BBa.rawtlt -n BBa_new.tlt -o BBa_fin.xf -d -1 
- **User gives the initial diameter value, no fast mode, no log output:**
    ```bash
    ./markerauto2 -i BBa.st -a BBa.rawtlt -n BBa_new.tlt -o BBa_fin.xf -d 8
- **The software will automatically initialize the diameter, with log output:**
    ```bash
    ./markerauto2 -i BBa.st -a BBa.rawtlt -n BBa_new.tlt -o BBa_fin.xf -d -1 -v
- **The software will automatically initialize the diameter, with fast mode:**
    ```bash
    ./markerauto2 -i BBa.st -a BBa.rawtlt -n BBa_new.tlt -o BBa_fin.xf -d -1 -t
## Generate final stack

After users use Markerauto2 to obtain optimized projection parameters, they can linearly transform the image sequence through the projection parameters to achieve image alignment. Files that perform linear transformations on images are stored in *.xf files. Users can use Autom or IMOD to process the *.xf file to realize the linear transformation with the following commands:

- **For Autom:**
  ```bash
  mrcstack -i BBa.st -o BBa_fin.mrc -x BBa_fin.xf
  
- **For IMOD:**
  ```bash
  newsatck -i BBa.st -o BBa_fin.mrc -x BBa_fin.xf

