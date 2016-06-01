# CloverLeaf_OpenACC - feature/heterogeneous

This is the OpenACC version of CloverLeaf version 1.3. 

## Branch Notes

### Explicit mesh decomposition
This branch lets you explicitly control the mesh decomposition ratios between MPI ranks.

Simply state:
ratio 0.75 0.2 0.05
in clover.in to get 70% of the mesh allocated to MPI rank 0, 20% to rank 1 and 5% to rank 2.
This is especially useful for running OpenACC heterogeneously - 1 MPI rank on the GPU and one on the host.
Internally a 1D processor decomposition is used to evenly balance the cells.


Note: When using the explicit mode the responsibility is on the user to decompose the problem correctly - performance will be poor with a bad configuration.
A message is printed to STDOUT to outline the resulting mesh decomposition, for the user to check.

### Explicit PGI targets

To enable support for the PGI multicore target architecture we have enabled three PGI compiler targets:

* `make COMPILER=PGI` - no OpenACC support
* `make COMPILER=PGI_GPU` - for GPU target (specifically nvidia,cc35)
* `make COMPILER=PGI_HOST` - for multicore target (running on the host)

To support these different options the generated binary is renamed to clover_leaf_COMPILER.



## Release Notes

### Version 1.3

CloverLeaf 1.3 contains a number of optimisations over previous releases.
These include a number of loop fusion optimisations and the use of scalar variables over work arrays.
Overall this improves cache efficiency.

This OpenACC version is based off of the CAPS Kernels version.

Due to limitations in data movement for Fortran data types, this version does not currently support tiling.




## Performance

Expected performance is give below.

If you do not see this performance, or you see variability, then is it recommended that you check MPI task placement and OpenMP thread affinities, because it is essential these are pinned and placed optimally to obtain best performance.

Note that performance can depend on compiler (brand and release), memory speed, system settings (e.g. turbo, huge pages), system load etc. 

### Performance Table

| Test Problem  | Time                         | Time                        | Time                        |
| ------------- |:----------------------------:|:---------------------------:|:---------------------------:|
| Hardware      |  E5-2670 0 @ 2.60GHz Core    | E5-2670 0 @ 2.60GHz Node    | E5-2698 v3 @ 2.30GHz Node   |
| Options       |  make COMPILER=INTEL         | make COMPILER=INTEL         | make COMPILER=CRAY          |
| Options       |  mpirun -np 1                | mpirun -np 16               | aprun -n4 -N4 -d8           |
| 2             | 20.0                         | 2.5                         | 0.9                         |
| 3             | 960.0                        | 100.0                       |                             |
| 4             | 460.0                        | 40.0                        | 23.44                       |
| 5             | 13000.0                      | 1700.0                      |                             |

### Weak Scaling - Test 4

| Node Count | Time         |
| ---------- |:------------:|
| 1          |   40.0       |
| 2          |              |
| 4          |              |
| 8          |              |
| 16         |              |


### Strong Scaling - Test 5

| Node Count | Time          | Speed Up |
| ---------- |:-------------:|:--------:|
| 1          |   1700        |  1.0     |
| 2          |               |          |
| 4          |               |          |
| 8          |               |          |
| 16         |               |          |


