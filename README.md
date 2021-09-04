# simple-md
Basic molecular dynamics implementation

## Introduction

Classical molecular dynamics simulations assume a set of 
particles interacting due to the forcefields the particles
exert on each other. This repository contains a simple 
demonstration program for molecular dynamics written in
[SYCL](https://www.khronos.org/sycl/).

## Equations of Motion

In classical molecular dynamics, particles follow 
[Newton's laws of motion](http://en.wikipedia.org/wiki/Newton%27s_laws_of_motion)

  Force = Mass x Acceleration

The great challenge in using this model is in finding
an appropriate force function.  A commonly used 
energy function is the [Lennard Jones interatomic potential
function](http://en.wikipedia.org/wiki/Lennard-Jones_potential)
which describes the energy between two interacting particles
as a function only of the distance between the particles:

  E = 4 \epsilon * [  (\sigma/r)^12 -  ( \sigma /r ) ^6 ]
     
where \epsilon and \sigma are model parameters that can be fitted
for different particles.  This empirical function expresses the
observation that particles which get too close to each other
will repel each other, particles which are very far apart
have no influence on each other, and particles an intermediate
distance away, attract each other. Thus for each particle,
equations of motion are

    acceleration of particle i 
  = [sum over all other particles F(r_ij) ]/[mass of particle i]

wher r_ij is the distance between particle i and particle j and
F(r_ij) is the force between particle i and particle j.
These equations of motion can be integrated in time using the [velocity
verlet algorithm](https://en.wikipedia.org/wiki/Verlet_integration)

   v_i^(n+1/2) = v_i^n + 0.5 x [delta t] x F_i^n

   x_i^(n+1) =  x_i^n ] + [delta t] x v_i^(n+1/2) 

   v_i^(n+1) = v_i^(n+1/2) + [delta t] x F_i^(n+1)

## Pseudo code

A pseudo code for this algorithm is

0. For each particle, calculate pairwise distances between it and other particles,
then update the force arrays
1. For each particle, calculate intermediate velocities using a half step
2. For each particle, calculate new position using intermediate velocities
3. For each particle, calculate the new pairwise distances between it and 
other particles, then update the force arrays
4. For each particle calculate final velocities using a half step
5. Go to step 1 and repeat for further timesteps

The main level of parallelism available, is over the particles.

## Serial C implementation

See the file [md.c](md.c)

## Sycl implementation

See the file [md_sycl.cpp](md_sycl.cpp)
`
### Running on intel devcloud

#### Build hipSYCL
```bash
cd $HOME
git clone https://github.com/VCCA2021HPC/simple-md
cd simple-md
qsub -l nodes=1:gen9:ppn=2 -d . get_hipSYCL.sh
```

#### Build programs with hipSYCL and execute

Once the build has finished, compile and run the programs
```bash
cd $HOME/simple-md
$HOME/hipSYCL-install/bin/syclcc -O3 md.c -o gcc_md
echo "time ./gcc_md" > job.sh
qsub -l nodes=1:gen9:ppn=2 -d . gcc_job.sh
$HOME/hipSYCL-install/bin/syclcc -O3 md_sycl.c -o gcc_md_sycl
echo "time ./gcc_md_sycl" > gcc_job_sycl.sh
qsub -l nodes=1:gen9:ppn=2 -d . gcc_job_sycl.sh
```

#### Build programs with DPCPP and execute

This is at present very slow when run in parallel
```bash
cd $HOME/simple-md
dpcpp md.c -o dpcpp_md
echo "time ./dpcpp_md" > dpcpp_job.sh
qsub -l nodes=1:gen9:ppn=2 -d . dpcpp_job.sh
dpcpp -lOpenCL -lsycl md_sycl.cpp -o dpcpp_md_sycl
echo "time ./dpcpp_md_sycl" > job_sycl.sh
qsub -l nodes=1:gen9:ppn=2 -d . dpcpp_job_sycl.sh
```

## Discussion

The current algorithm requires O(p^2 x n) running time, where p is the number
of particles and n is the number of timesteps. This can be quite 
prohibitive when we want to simulate many particles. For many forcefields, one can
either ignore the contributions from particles that are far away and keep
track of a list of particles that are nearby, or one can decompose the force
into a nearby and a far away component. The far away component can usually be
approximated using a fast algorithm to get an computational cost of O(p x n) 

## References

- [James Vance et al, "Molecular Dynamics Collaborative Project"](https://github.com/jnvance/ljmd-c)
- [Furio Ercolessi, "A Molecular Dynamics Primer"](https://www.glennklockwood.com/materials-science/molecular-dynamics/ercolessi-1997.pdf)
- [Michael P. Allen, "Introduction to Molecular Dynamics Simulation"](https://udel.edu/~arthij/MD.pdf)
- [Sina Kazemi and Peter Güntert, "Molecular Dynamics Simuation Tutoral"](http://www.bpc.uni-frankfurt.de/guentert/wiki/images/9/96/180618_TutorialMD.pdf)
- [Eugene Klyshko "Molecular Dynamics Simulations in Python"](https://klyshko.github.io/teaching/2019-03-01-teaching)
- ["Molecular Modeling Practical"](http://www.cgmartini.nl/~mdcourse/pepmd/index.html)
- [Eijkhout, "Molecular Dynamics"](https://pages.tacc.utexas.edu/~eijkhout/istc/html/md.html)
- [Yi Yao, "Molecular Dynamics Simulations by Gromacs, LAMMPS, CPMD, CP2K, CP.X and Qbox from Scratch"](https://yaoyi92.github.io/molecular-dynamics-simulations-by-gromacs-lammps-cpmd-cp2k-cpx-qbox-from-scratch.html)
- [Schroeder "Interactive Molecular Dynamics", website](https://physics.weber.edu/schroeder/md/)
- [Schroeder "Interactive Molecular Dynamics", paper](https://physics.weber.edu/schroeder/md/InteractiveMD.pdf)
- [Allister Beharry, "oneMD: A Molecular Dynamics Simulator in Modern C++ and SYCL"](https://www.codeproject.com/Articles/5295109/oneMD-A-Molecular-Dynamics-Simulator-in-Modern-Cpl)
- [Wes Barnett, "Lennard Jones Molecular Dynamics"](https://github.com/wesbarnett/lennardjones)
- [Wes Barnett, "Gromacs Tutorials"](https://group.miletic.net/en/tutorials/gromacs/)
- [Jeff Hammond, "A comparative analysis of Kokkos and Sycl"](https://www.iwocl.org/wp-content/uploads/iwocl-2019-dhpcc-jeff-hammond-a-comparitive-analysis-of-kokkos-and-sycl.pdf)
- [Brian Homerding, "Investigation of the performance of SYCL kernels across various architectures"](https://p3hpcforum2020.alcf.anl.gov/wp-content/uploads/sites/8/2020/09/P3HPC_Homerding_Day-1.pdf)
- [Beau Johnston, Josh Milthorpe, Jeffrey S. Vetter, "Evaluating the Performance and Portability of Contemporary SYCL Implementations"](https://www.researchgate.net/publication/345990610_Evaluating_the_Performance_and_Portability_of_Contemporary_SYCL_Implementations)
- [Tom Deakin and Tom Lin, "Experiences with SYCL for the hydrodynamics mini-app Cloverleaf"](http://uob-hpc.github.io/2020/01/06/cloverleaf-sycl.html)
- [Zheming Jin, "A Case Study with the HACCmk Kernel in SYCL"](https://publications.anl.gov/anlpubs/2019/12/157540.pdf)
- [Ramesh Peri, "Analyzing the Performance of Reduction Operations in DPC++"](https://software.intel.com/content/www/cn/zh/develop/articles/analyzing-performance-reduction-operations-dpc.html#gs.9kg8aw)
- [Zheming Jin, "Evaluating the Performance of Integer Sum Reduction in SYCL on GPUs"](https://oaciss.uoregon.edu/icpp21/views/includes/files/ppss_pap101s3-file2.pdf)
