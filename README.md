
<h1 align="center">
  <br>
  <a href="http://dual.sphysics.org/"><img src="http://design.sphysics.org/img/logo_dualsphysics.png" alt="DualSPHysics" width="300"></a>
  <br>
  DualSPHysics
  <br>
</h1>

<h4 align="center">A combined <a href="https://developer.nvidia.com/cuda-zone" target="_blank">CUDA</a> and <a href="http://www.openmp.org/" target="_blank">OpenMP</a> implementation of the Smoothed Particle Hydrodynamics method based on the advanced <a href="https://wiki.manchester.ac.uk/sphysics/index.php/Main_Page" target="_blank">SPHysics</a> code.</h4>

<p align="center">
<img src="http://design.sphysics.org/img/dualsphysics_demonstration.gif" alt="DualSPHysics" width="400">
</p>

# Instructions for regular users

If you only want a copy of DualSPHysics to make and run cases in your system, you probably want <a href="http://www.dual.sphysics.org/index.php/downloads/" target="_blank">the full DualSPHysics package</a> from the official webpage. There you will find documentation and a different set of packages for different versions of the application and different Operating Systems.

It is possible that you want the latest version in this repository that is not yet uploaded in the official web page. In this case check the [Building the project](#building-the-project) section to build an executable.

Have in mind that DualSPHysics needs a generated case to execute the SPH solver, so you want to check out GenCase too, included the main package on the <a href="http://www.dual.sphysics.org/index.php/downloads/" target="_blank">DualSPHysics webpage</a>.

If you need help check out the wiki for this project, or you can also check the previously mentioned download page to get a PDF guide.

# Instructions for developers

If you are a developer and want to use this code check the following guides.

Take into account that for generating cases and post-processing them you will need to use GenCase and the different post-processing tools included in the main DualSPHyisics package, <a href="http://www.dual.sphysics.org/index.php/downloads/" target="_blank">downloadable here</a>. If you compile your own DualSPHyiscs version just overwrite (or don use) the one in the package, but yours.

## Developing a modified version to fit your own needs.

You can fork this repository and change or add anything you want in it. Keep in mind that your changes will not be taken into account into the main versions. If your objective is to add your changes/improvements to the main code, check the next section.

## Developing a modified or improved version to contribute to the project.

We appreciate your efforts! But please, if you are trying to develop a feature to add to the main repository, be sure to follow the steps on the [CONTRIBUTING.md](CONTRIBUTING.md) file, to make the process easier for both parts and optimize the process.

# Building the project

## On Microsoft Windows

This application is being developed in Visual Studio 2013. The repository contains project files for that IDE and it is recomended to use them to make the process easier.

If you prefer it, you can also use the [Makefile](Source/Makefile). It is possible that you'll need to change some things in it to make it work.

## On GNU/Linux

You can build the project in GNU/Linux using the [Makefile](Source/Makefile) included in the source folder. It's meant to be used in this OS but may need a bit of tinkering to make it work.

# Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) if you want to make changes to the code that would be later merged with the main project.

# Authors and people involved

* **Dr Jose M. Dominguez** - *Main Developer*
* **Dr Alejandro Crespo**
* **Prof. Moncho Gomez Gesteira**
* **Dr Anxo Barreiro**
* **Dr Benedict Rogers**
* **Dr Georgios Fourtakas**
* **Dr Athanasios Mokos**
* **Dr Renato Vacondio**
* **Dr Ricardo Canelas**
* **Dr Stephen Longshaw**
* **Dr Corrado Altomare**

See also the list of [contributors](https://github.com/dualsphysics/DualSPHysics/contributors) who participated in this project.

## License

This project is licensed under the LGPL License - see the [LICENSE](LICENSE) file for details.
