# foamPlasmaToolkit

**foamPlasmaToolkit** is a collection of solvers, models, boundary conditions, and utilities
developed for plasma applications in **OpenFOAM**. It aims to extend the capabilities of OpenFOAM to
simulate plasma-related phenomena, including plasma–fluid interactions, plasma chemistry, and 
coupled electromagnetic effects.

##	Description

This toolkit provides an evolving framework for implementing plasma physics features within OpenFOAM.  
It is designed for research and development purposes, focusing on modularity and compatibility with 
the OpenFOAM v2412 (OpenCFD Ltd.).  

##	Status

The project is currently in its **initial development phase**.  
Functionality, structure, and documentation are expected to change as new components are added.

##	License

This project is distributed under the terms of the **GNU General Public License v3.0 (GPLv3)**.  
For full license details, see the [LICENSE](LICENSE) file.

© 2025 Rention Pasolari

This software is not part of OpenFOAM, but it is developed using the OpenFOAM framework and linked 
against OpenFOAM libraries (v2412).

##	Requirements

- **OpenFOAM v2412 (OpenCFD Ltd.)**  

This toolkit has been **tested and verified to work on**:
- **Ubuntu 24.04 LTS**
- **OpenFOAM v2412 (OpenCFD Ltd. release)**

### Optional Dependencies

- **PETSc** (for users who want PETSc-based solvers)

> PETSc is not required for basic installation.
> See [`docs/petsc.md`](docs/petsc.md) for setup instructions.

## Installation / Setup

### 1. Install OpenFOAM v2412 (required)

This toolkit is designed **for OpenFOAM v2412 (OpenCFD Ltd.)**. It has not been tested in other 
versions.

If OpenFOAM is not already installed, follow the official installation instructions from OpenCFD:

➡ https://develop.openfoam.com/Development/openfoam/-/wikis/precompiled

Once installed, ensure OpenFOAM is correctly loaded in your terminal:

```bash
source /opt/openfoam2412/etc/bashrc
foamVersion   # should print: OpenFOAM-v2412
```

### 2. Download foamPlasmaToolkit

Download the toolkit by cloning the GitHub repository:

```bash
git clone https://github.com/rpasolari/foamPlasmaToolkit.git
```

Enter the toolkit directory:

```bash
cd foamPlasmaToolkit
```

Activate the toolkit environment (so OpenFOAM can detect its extensions):
```bash
source etc/bashrc
```

:warning: **IMPORTANT:** To avoid having to run this manually every time, add it to 
your `~/.bashrc` as shown below.

Open the .bashrc file in the user's directory in an editor, e.g. by typing in a terminal window 
(note the dot):

```bash
gedit ~/.bashrc
```

Then add the following line at the end of the file (change the path to match where you installed the
toolkit, e.g. /home/USERNAME):

```bash
source /ABSOLUTE_PATH_TO_FOAMPLASMATOOLKIT/etc/bashrc
```
To verify that the environment was successfully loaded, open a **new terminal** and run:

```bash
echo $FOAM_PLASMA
```

### 3. Build foamPlasmaToolkit
From inside the foamPlasmaToolkit directory, compile everything:

```bash
./Allwmake
```
If the compilation finishes without errors, the installation is complete.

##	Contributions

All development and maintenance are currently handled by:
**Rention Pasolari <r.pasolari@gmail.com>**  
