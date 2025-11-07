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

- **PETSc 3.24.0**
- **petsc4foam**
- **OpenFOAM v2412 (OpenCFD Ltd.)**
- **Ubuntu 24.04 LTS**

> Other versions may work but are not guaranteed.

## 1. Installing PETSc

### A. Go to ThirdParty directory (*OF Directory*) and download PETSc

```bash
foam
cd ThirdParty
```

```bash
git clone -b release https://gitlab.com/petsc/petsc.git petsc-3.24.0
cd petsc-3.24.0
git checkout v3.24.0    # pin exactly to v3.24.0
```
### B. Clean old builds

```bash
make distclean >/dev/null 2>&1 || true
rm -rf DPInt32
rm -rf DPInt32-cuda
```

### C. Configure PETSc

#### I. Configure *WITHOUT* CUDA

```bash
./configure \
  --with-shared-libraries=1 \
  --with-precision=double \
  --with-debugging=0 \
  --with-cc=mpicc \
  --with-cxx=mpicxx \
  --with-fc=mpif90 \
  --with-mpi=1 \
  --download-fblaslapack \
  --download-metis \
  --download-parmetis \
  --download-superlu_dist \
  --download-hypre \
  --with-petsc-arch=DPInt32 \
  COPTFLAGS="-O3" \
  CXXOPTFLAGS="-O3" \
  CUDAOPTFLAGS="-O3"
```

#### II. Configure *WITH* CUDA

:warning: **IMPORTANT:** Before configuring, **CUDA Toolkit** must be installed.
> See [`docs/cuda.md`](docs/cuda.md) for setup instructions.

:information_source: **GPU Compute Capability Required**

The flag `--with-cuda-arch=<value>` must match the **SM architecture** of your GPU.  
See NVIDIA's Compute Capability list here:  
https://developer.nvidia.com/cuda-gpus  

Example:
- RTX 5090 / 5080 / 5070 / 5060 / 5050 → `--with-cuda-arch=120` (Blackwell)
- RTX 4090 / 4080 / 4070 / 4060 / 4050 → `--with-cuda-arch=89` (Ada Lovelace)

```bash
./configure \
  --with-cuda=1 \
  --with-cudac=nvcc \
  --with-cuda-dir=/usr/local/cuda-12.8 \
  --with-cuda-arch=120 \
  --with-shared-libraries=1 \
  --with-precision=double \
  --with-debugging=0 \
  --with-cc=mpicc \
  --with-cxx=mpicxx \
  --with-fc=mpif90 \
  --with-mpi=1 \
  --download-fblaslapack \
  --download-metis \
  --download-parmetis \
  --download-superlu_dist \
  --download-hypre \
  --with-petsc-arch=DPInt32-cuda \
  COPTFLAGS="-O3" \
  CXXOPTFLAGS="-O3" \
  CUDAOPTFLAGS="-O3"
```

#### D. Add PETSc to Environment

Open the .bashrc file in the user's directory in an editor, e.g. by typing in a terminal window 
(note the dot):

```bash
gedit ~/.bashrc
```

Then add the following lines:

#### I. *WITHOUT* CUDA

```bash
export PETSC_DIR=/opt/openfoam2412/ThirdParty/petsc-3.24.0
export PETSC_ARCH=DPInt32
export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
```

#### II. *WITH* CUDA

```bash
export PETSC_DIR=/opt/openfoam2412/ThirdParty/petsc-3.24.0
export PETSC_ARCH=DPInt32-cuda
export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
```

:warning: **IMPORTANT: WSL2 GPU Notice**  

When using **OpenFOAM + PETSc with CUDA under WSL2**, you probably need to disable CUDA-aware
shared-memory MPI.
WSL2 does **not fully support CUDA IPC / pinned memory**, which can cause crashes.

Add this to your terminal or `~/.bashrc`:

```bash
export OMPI_MCA_btl=^smcuda
export OMPI_MCA_pml=ob1
export OMPI_MCA_opal_cuda_support=true
```

Reload:

```bash
source ~/.bashrc
```

### E. Build and check

#### I. *WITHOUT* CUDA

```bash
make PETSC_DIR=$PWD PETSC_ARCH=DPInt32 all -j$(nproc)
make PETSC_DIR=$PWD PETSC_ARCH=DPInt32 check
```

#### II. *WITH* CUDA

```bash
make PETSC_DIR=$PWD PETSC_ARCH=DPInt32-cuda all -j$(nproc)
make PETSC_DIR=$PWD PETSC_ARCH=DPInt32-cuda check
```

## 2. Installing petsc4foam
> **Note:** *`petsc4Foam` is automatically compiled during the build of the `foamPlasmaToolkit`
**as long as PETSc is already installed on your system**. When you install the toolkit after
installing PETSc, it will detect PETSc and build `petsc4Foam` automatically.*

> *If the `foamPlasmaToolkit` was installed **before** PETSc, then `petsc4Foam` will not have been
built. In that case, after installing PETSc manually (following the steps above), you need to
compile `petsc4Foam` yourself by following the steps below.*

⚠️ **IMPORTANT:** The `petsc4Foam` version included in the `foamPlasmaToolkit` is 
**not the same** as the version available from the official upstream repository.  
The toolkit version contains additional modifications that enable **implicit (monolithic) 
multi-region coupling**, allowing multi-processing simulations to solve tightly coupled regions in a
single PETSc system.  

Therefore, even if you have installed the official petsc4Foam, you must use the version bundled with
the `foamPlasmaToolkit` when performing multi-region implicit simulations (see below).

### A. Go to petsc4Foam ThirdParty directory (foamPlasmaToolkit)

```bash
cd $FOAM_PLASMA_THIRD_PARTY/petsc4Foam
```


### B. Build petsc4foam

```bash
./Allwmake -prefix=openfoam
```

### C. Ensure that petsc4Foam is installed
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

For usage examples, see the `foamPlasmaToolkit` [`tutorials`](tutorials) or the official
`petsc4foam` tutorials at <https://gitlab.com/petsc/petsc4foam>.

