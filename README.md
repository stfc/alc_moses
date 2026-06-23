## About the code
ALC_MOSES is the acronym for [**A**da **L**ovelace **C**entre]_[**MO**delling and **S**ettings for **E**lectromical **S**imulations]. It is a software that generates atomistic models of electrified interfaces automatically, together with input files for Density Functional Theory (DFT) simulations of such interfaces, including the possibility of grand-canonical DFT (GC-DFT). The implemented options also allow the setting of implicit models, which are crucial to account for the electrostatic screening of the solvent (and electrolyte) and reduce the computational cost of MD simulations. In addition, the software can also generate script files for HPC submission. It is important to mention this tool does not execute simulations but set the required input files.

## Disclaimer
The ALC does not fully guarantee the code is free of errors and assumes no legal responsibility for any incorrect outcome or loss of data.

## Contributors
**Original author:** Ivan Scivetti (SCD, STFC)  

## Structure of files and folders
ALC_MOSES contains the following set of files and folders (in italic-bold):

* [***CI-tests***](./CI-tests): contains the tests files (in .tar format) needed for testing. The user should execute the available scripts of the [***tools***](./tools) folder to run the test automatically and verify the code has been installed properly (see the [build_code.md](./build_code.md) file for instructions).
* [***examples***](./examples): folder with example cases to help familiarising with the code (see [manual.pdf](./manual.pdf)).  
* [***scripts***](./scripts): contains scripts for data processing.
* [***source***](./source): contains the source code. Files have the *.F90* extension
* [***tools***](./tools): shell files for building, compiling and testing the code automatically.
* [.gitignore](./.gitignore): instructs Git which file to ignore.  
* [CMakeLists.txt](./CMakeLists.txt): sets the framework for code building and testing with CMake. This file must ONLY be modified to add test cases.  
* README.md: this file.
* [build_code.md](./build_code.md): steps to build, compile and run tests using the CMake platform.  
* [manual.pdf](./manual.pdf): ALC_MOSES manual.  

## Dependencies
The user must have access to the following software (locally):

* GNU-Fortran (11.2.0) or Intel-Fortran (ifx 2025.1.0)
* Gnuplot (5.2)
* Cmake (3.16)
* Make (4.2.1)
* Git (2.34.1)

Information in parenthesis indicates the minimum version tested during the development of the code. The specification for the minimum versions is not fully rigorous but indicative, as there could be combinations of other minimum versions that still work. 

## Getting started

### Obtaining the code
The user can clone the code locally by executing the following command with the SSH protocol
```sh
$ git clone git@github.com:stfc/alc_moses.git
```
Instead, if the user wants to use the HTTPS protocol it must execute
```sh
$ git clone https://github.com/stfc/alc_moses.git
```
Both ways generate the ***alc_moses*** folder as the root directory. Alternatively, the code can be downloaded from any of the availab
le assets.

### Building and testing the code with CMake
Details can be found in file [build_code.md](./build_code.md)

### Using the code
The user is referred to the [manual](./manual.pdf) of the code for a detailed description of the implemented functionalities. The examples cases within the [***examples***](./examples) folder complement sec. 3 of the manual.

## Acknowledgements
* ALC for funding.  
* Gilberto Teobaldi (SCD,STFC) for support.  
* Ziwei Chai for help in setting up GC-DFT simulation files using CP2K.
* Manuel dos Santos Dias, Brad Ayers and Jacek Dziedzic for support with the ONETEP code.
* Ciaran O'Brien for discussions on fundamental aspects and applications of the GC-DFT formulation.
* Lesley Mansfield for management support.  
* The Innovation Department of the STFC for assistance with the licensing process. 
