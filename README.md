
<h1 align="center">
  <br>
  <a href="http://dual.sphysics.org/"><img src="http://design.sphysics.org/img/logo_dualsphysics.png" alt="DualSPHysics" width="300"></a>
  <br>
  DualSPHysics
  <br>
</h1>

<h4 align="center"><a href="http://www.dual.sphysics.org" target="_blank">DualSPHysics</a> is based on the Smoothed Particle Hydrodynamics model named <a href="http://www.sphysics.org" target="_blank">SPHysics</a>.</h4>

<h4 align="center">The code is developed to study free-surface flow phenomena where Eulerian methods can be difficult to apply, such as waves or impact of dam-breaks on off-shore structures. DualSPHysics is a set of C++, <a href="https://developer.nvidia.com/cuda-zone" target="_blank">CUDA</a> and Java codes designed to deal with real-life engineering problems.</h4>

# Instructions for regular users

If you only want a copy of DualSPHysics to create and run cases in your system, you probably want <a href="http://www.dual.sphysics.org/index.php/downloads/" target="_blank">the full DualSPHysics package</a> from the official website. There you will find documentation and packages of different versions for different Operating Systems.

It is possible that you want the latest version in this repository that is not yet uploaded in the official web page. In this case check the [Building the project](#building-the-project) section to build an executable.

Have in mind that DualSPHysics needs a case already created to execute the SPH solver, so you need to use GenCase, which is included the main package on the <a href="http://www.dual.sphysics.org/index.php/downloads/" target="_blank">DualSPHysics webpage</a>.

If you need help check out the wiki for this project.

# Instructions for developers

If you are a developer and want to use this code check the following guides.

Take into account that for pre- and post-processing you will need to use GenCase and the different post-processing tools included in the main DualSPHysics package, <a href="http://www.dual.sphysics.org/index.php/downloads/" target="_blank"> here</a>. If you compile your own DualSPHyiscs version just overwrite the one in the package.

## Developing a modified version to fit your own needs.

You can fork this repository and change or add anything you want. Keep in mind that your changes will not be taken into account into the main versions. If your objective is to implement your changes/improvements to the main code, check the next section.

## Developing a modified or improved version to contribute to the project.

We appreciate your efforts! But please, if you are trying to develop/implement a functionality to be added to the main repository, be sure to follow the steps described in the [CONTRIBUTING.md](CONTRIBUTING.md) file. 

# Building the project

## Microsoft Windows

This application is being developed in Visual Studio Community 2015 since it is free and compatible with CUDA 9.2 (<a href="https://www.visualstudio.com/vs/older-downloads/" target="_blank">download web</a>). The repository contains project files.

Make sure that you install the CUDA SDK beforehand if you want to compile the GPU version, and configure the Visual Studio project to point to the CUDA libraries directory to compile (now prepared for CUDA 9.2).

You can also use the [Makefile](src/source/Makefile). It is possible that you'll need to edit it. Check the GNU/Linux guide on how to compile if you're using the makefile, as it is mostly the same, just installing the things manually by yourself.

## GNU/Linux

### Using Makefile

You can build the project in GNU/Linux using the [Makefile](src/source/Makefile) included in the source folder. Follow these steps (for the GPU version):

1. Clone this repository into your system `git clone https://github.com/DualSPHysics/DualSPHysics.git`
2. In a terminal, go to the folder `cd DualSPHysics/src/source/`
3. Edit the `Makefile` file with a text editor and then:
   * Set the `DIRTOOLKIT` variable with the path to CUDA in your system e.g. `DIRTOOLKIT=/opt/cuda`
   * Make sure that your `G++` compiler version is compatible with the CUDA version installed in your system (e.g. CUDA 9.1 only supports G++ versions up to 6.x). If you want to use an specific version, you should modify the variable `CC`, for example: `CC=/usr/local/bin/g++-6`
   * If you are using G++ version 5 of greater make sure to set `USE_GCC5=YES`
4. Execute `make clean` to make sure the environment is clean and ready to compile.
5. Execute `make`

After compiling you should see a message like `--- Compiled Release GPU/CPU version ---`. Go to `bin/linux/` to check that `DualSPHyiscs_linux64` or `DualSPHyiscsCPU_linux64` is there and build correctly.

**For the CPU version**: If you want to compile de CPU version just ignore CUDA and use the makefile `Makefile_cpu`. To specify a different file to `make`, use the `-f` parameter: `make -f Makefile_cpu`

### Using CMake

Alternatively you can use [CMake](https://cmake.org) to compile DualSPHysics following these steps:

1. Clone this repository into your system `git clone https://github.com/DualSPHysics/DualSPHysics.git`
2. In a terminal, go to the folder `cd DualSPHysics/src/source/`
3. Create a temporal directory were the build files will be placed `mkdir build`
4. Make sure that your `GCC` compiler version is compatible with the CUDA version installed in your system (e.g. CUDA 9.1 only supports G++ versions up to 6.x). If you want to use an specific version, you should modify the environment variable `CC` and `CXX`, for example: `export CC=/usr/local/bin/gcc-6; export CXX=/usr/local/bin/g++6`
4. Go to the build folder and execute cmake `cd build; cmake ..`
5. Execute `make`

If CUDA is not installed in your system, only the CPU version will be compiled.

If you want the binaries be placed into `bin/linux` directory, just type `make install`

# Graphical user interface (GUI)

<h1 align="center">
  <a href="http://design.sphysics.org/"><img src="https://i.imgur.com/D5dabSD.png" alt="DualSPHysics" width="1000"></a>
</h1>

Please check [DesignSPHysics website](http://design.sphysics.org/).

# Advanced visualisation with Blender

<h1 align="center">
  <a href="http://visual.sphysics.org/"><img src="https://i.imgur.com/I1psKiT.png" alt="DualSPHysics" width="1000"></a>
</h1>

Please check [VisualSPHysics website](http://visual.sphysics.org/).

# Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) if you want to make changes to the code.

# Authors and people involved

* **Dr Jose M. Dominguez Alonso** - *Main Developer*
* **Dr Alejandro J.C. Crespo**
* **Prof. Moncho Gomez Gesteira**
* **Prof. Benedict D. Rogers**
* **Dr Georgios Fourtakas**
* **Prof. Peter Stansby**
* **Dr Renato Vacondio**
* **Dr Corrado Altomare**
* **Dr Angelo Tafuni**
* **Orlando Garcia Feal**
* **Ivan Martinez Estevez**


## Former people involved
* **Dr Jose Gonzalez Cao**
* **Dr Ricardo Canelas**
* **Dr Athanasios Mokos**
* **Dr Stephen Longshaw**
* **Dr Anxo Barreiro**

See also the list of [contributors](https://github.com/dualsphysics/DualSPHysics/contributors) who participated in this project.

## License

This project is licensed under the LGPL License - see the [LICENSE](LICENSE) file for details.
