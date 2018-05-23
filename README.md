
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

This application is being developed in Visual Studio Community 2015 since it is free and compatible with CUDA 9.1 (<a href="https://www.visualstudio.com/vs/older-downloads/" target="_blank">download web</a>). The repository contains project files.

Make sure that you install the CUDA SDK beforehand if you want to compile the GPU version, and configure the Visual Studio project to point to the CUDA libraries directory to compile (now prepared for CUDA 9.1).

You can also use the [Makefile](Source/Makefile). It is possible that you'll need to edit it. Check the GNU/Linux guide on how to compile if you're using the makefile, as it is mostly the same, just installing the things manually by yourself.

## GNU/Linux

You can build the project in GNU/Linux using the [Makefile](Source/Makefile) included in the source folder. Follow these steps (for the GPU version):

1. Clone this repository into your system
2. Ensure you have GCC version 4.X installed. Usually there are packages in your distro like `gcc49` that provides the `g++-4.9` executable.
3. In a terminal, go to the folder `src/source/`
4. Execute `make clean` to make sure the environment is clean and ready to compile
5. Execute `make CC=g++-4.9 CPP=g++-4.9 CXX=g++-4.9 LD=g++-4.9 -f ./Makefile`. Be sure to replace `g++-4.9` for the executable name you have in your system (previously installed in step 2)

After compiling you should see a message like `--- Compiled Release GPU/CPU version ---`. Go to `bin/linux/` to check that `DualSPHyiscs_linux64` or `DualSPHyiscsCPU_linux64` is there and build correctly.

**For the GPU version**: Install the CUDA package and edit the makefile to point the CUDA libs directories to the paths on your system. Also, check the CUDA version and adapt the Makefile according to that.

**For the CPU version**: If you want to compile de CPU version just ignore CUDA and use the makefile `Makefile_cpu`

# Graphical user interface (GUI)

Please check [DesignSPHysics website](http://design.sphysics.org/).

# Advanced visualisation with Blender

Please check [VisualSPHysics website](http://visual.sphysics.org/).

# Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) if you want to make changes to the code.

# Authors and people involved

* **Dr Jose M. Dominguez** - *Main Developer*
* **Dr Alejandro Crespo**
* **Prof. Moncho Gomez Gesteira**
* **Dr Anxo Barreiro**
* **Orlando G. Feal**
* **Dr Benedict Rogers**
* **Prof. Peter Stansby**
* **Dr Georgios Fourtakas**
* **Dr Athanasios Mokos**
* **Dr Renato Vacondio**
* **Dr Ricardo Canelas**
* **Dr Stephen Longshaw**
* **Dr Corrado Altomare**

See also the list of [contributors](https://github.com/dualsphysics/DualSPHysics/contributors) who participated in this project.

## License

This project is licensed under the LGPL License - see the [LICENSE](LICENSE) file for details.
