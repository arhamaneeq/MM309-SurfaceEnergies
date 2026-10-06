# Contributing

This document covers the development environment, LAMMPS installation, project workflow, and basic contribution guidelines.

## Development Environment

I recommend working with Linux or WSL2 (Windows Subsystem for Linux).

If you're on Windows, install WSL2 with Ubuntu. Microsoft provides the installation instructions here:

https://learn.microsoft.com/en-us/windows/wsl/install

Once Ubuntu is installed, make sure your packages are up to date:

```sh
sudo apt update
sudo apt upgrade
```

You'll also need the basic development tools:

```sh
sudo apt install git cmake build-essential
```

### Installing LAMMPS without GPU Support

We will use LAMMPS for the Molecular Dynamics simulations in this project.

Clone the LAMMPS repository:

```sh
git clone https://github.com/lammps/lammps.git
```

Then enter the repository

```sh
cd lammps
```

I recommend building LAMMPS using CMake

Create a separate Build directory:

```sh
mkdir build
cd build
```

Configure the build:

```sh
cmake ../cmake
```

Compile:

```sh
cmake --build . -j$(nproc)
```

The resulting lmp executable will be located in `lammps/build/lmp`

You can verify the installation using

```sh
./lmp -h
```

### Building LAMMPS with NVIDIA GPU Support

If you have an NVIDIA GPU, you can build LAMMPS with KOKKOS to accelerate supported calculations on the GPU.

First, make sure your NVIDIA drivers and CUDA toolkit are installed and working:

```sh
nvidia-smi
```

If this command does not work, fix your NVIDIA/CUDA installation before attempting to build LAMMPS.

From the LAMMPS repository:

```sh
mkdir -p build
cd build
```

Configure LAMMPS with KOKKOS:

```sh
cmake ../cmake \
    -D PKG_KOKKOS=on \
    -D Kokkos_ENABLE_CUDA=on \
    -D Kokkos_ENABLE_OPENMP=on
```

Then compile

```sh
cmake --build . -j$(nproc)
```

Verify the build:
```sh
./lmp -h
```

You should see KOKKOS listed amongst the installed packages.

> [!NOTE]
> The exact KOKKOS/CUDA configuration may depend on your GPU architecture and CUDA installation. If the configuration fails, check the LAMMPS and KOKKOS documentation before modifying the build configuration.

## Cloning This Project

Once LAMMPS is installed, clone this repository:

```sh
git clone https://github.com/arhamaneeq/MM309-SurfaceEnergies.git
cd MM309-SurfaceEnergies
```

You're now set up.