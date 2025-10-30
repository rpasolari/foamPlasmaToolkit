# PETSc Support
The `foamPlasmaToolkit` optionally supports **PETSc** to enable advanced linear solvers and
preconditioners, suitable for stiff plasma systems and large-scale parallel simulations.

PETSc is **not required** for basic use.  
Use it only if you need advanced iterative solvers or HPC acceleration.

## References

This documentation and functionality rely on the following projects:

- **PETSc** - Portable, Extensible Toolkit for Scientific Computation  
  https://petsc.org/

- **petsc4foam** - PETSc interface for OpenFOAM  
  https://gitlab.com/petsc/petsc4foam

## Tested Configuration

This toolkit has been tested with:

- **PETSc 3.24.4**
- **petsc4foam**
- **OpenFOAM v2412 (OpenCFD Ltd.)**
- **Ubuntu 24.04 LTS**

> Other versions may work but are not guaranteed.

## 1. Installing PETSc

### A. Go to ThirdParty directory and download PETSc 3.24.4

```bash
foam
cd ThirdParty
```

#### Option I - Clone PETSc 3.24.4 using Git

```bash
git clone -b release https://gitlab.com/petsc/petsc.git petsc-3.24.4
cd petsc-3.24.4
git checkout v3.24.4    # pin exactly to v3.24.4
```

#### Option II - Download PETSc 3.24.4

```bash
wget https://gitlab.com/petsc/petsc/-/archive/v3.24.4/petsc-v3.24.4.tar.gz
tar -xvf petsc-v3.24.4.tar.gz
mv petsc-v3.24.4 petsc-3.24.4
cd petsc-3.24.4
```

### B. Clean old builds

```bash
make distclean >/dev/null 2>&1 || true
rm -rf DPInt32-cuda
```

### C. Configure PETSc

#### I. Configure without CUDA

```bash
./configure \
  --with-cuda=0 \
  --with-shared-libraries=1 \
  --with-precision=double \
  --with-debugging=0 \
  --with-cc=mpicc \
  --with-cxx=mpicxx \
  --with-fc=mpif90 \
  --download-fblaslapack \
  --download-metis \
  --download-parmetis \
  --download-superlu_dist \
  --with-hypre=1 \
  --with-petsc-arch=DPInt32 \
  COPTFLAGS="-O3" \
  CXXOPTFLAGS="-O3"
```

#### II. Configure with CUDA

:warning: **IMPORTANT:** Before configuring, **CUDA Toolkit** must be installed.
> See [`docs/cuda.md`](docs/cuda.md) for setup instructions.

```bash
./configure \
  --with-cuda=1 \
  --with-cudac=nvcc \
  --with-cuda-dir=/usr/local/cuda-12.8 \
  --with-shared-libraries=1 \
  --with-precision=double \
  --with-debugging=0 \
  --with-cc=mpicc \
  --with-cxx=mpicxx \
  --with-fc=mpif90 \
  --download-fblaslapack \
  --download-metis \
  --download-parmetis \
  --download-superlu_dist \
  --with-hypre=1 \
  --with-petsc-arch=DPInt32-cuda \
  COPTFLAGS="-O3" \
  CXXOPTFLAGS="-O3" \
  CUDAOPTFLAGS="-O3"
```

### D. Build and check

#### I. Without CUDA

```bash
make PETSC_DIR=$PWD PETSC_ARCH=DPInt32 all -j$(nproc)
make PETSC_DIR=$PWD PETSC_ARCH=DPInt32 check
```

#### II. With CUDA

```bash
make PETSC_DIR=$PWD PETSC_ARCH=DPInt32-cuda all -j$(nproc)
make PETSC_DIR=$PWD PETSC_ARCH=DPInt32-cuda check
```

#### E. Add PETSc to Environment

Open the .bashrc file in the user's directory in an editor, e.g. by typing in a terminal window 
(note the dot):

```bash
gedit ~/.bashrc
```

Then add the following lines:

#### I. Without CUDA

```bash
export PETSC_DIR=/opt/openfoam2412/ThirdParty/petsc-3.24
export PETSC_ARCH=DPInt32
export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
```

#### II. With CUDA

```bash
export PETSC_DIR=/opt/openfoam2412/ThirdParty/petsc-3.24
export PETSC_ARCH=DPInt32-cuda
export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
```

Reload:

```bash
source ~/.bashrc
```

## 2. Installing petsc4foam

### A. Go to modules directory

```bash
foam
cd modules
```

### B. Remove external-solver and clone from git

```bash
rm -rf external-solver
git clone https://gitlab.com/petsc/petsc4foam.git external-solver
```

### C. Make petsc4foam

```bash
./Allwclean
./Allwmake
```

### D. Confirm installation

```bash
foamHasLibrary -verbose petscFoam   # should print: Can load "petscFoam"
```

## 3. Use petsc4foam

:warning: **IMPORTANT:** Add libpetscFoam.so to the optional keyword entry libs of the
controlDict file:

```bash
libs
(
   "libpetscFoam.so"
);
```

For usage examples, see the `foamPlasmaToolkit` [`tutorials`](tutorials) or the official `petsc4foam` tutorials
at <https://gitlab.com/petsc/petsc4foam>.

