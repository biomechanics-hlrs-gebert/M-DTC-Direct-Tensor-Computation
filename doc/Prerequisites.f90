!> \page prerequisites Prerequisites 
!>
!> struct-process depends on the following libraries:
!> <ol>
!>   <li>CMAKE<br>
!>   &nbsp;&nbsp;&nbsp;To build METIS.<br>
!>   &nbsp;&nbsp;&nbsp;On recent Linux operating systems like Fedora 33
!>   a suitable cmake version should be available
!>   via the package manager.</li>
!>   <li>BLAS<br>
!>   &nbsp;&nbsp;&nbsp;To build PETSc.<br>
!>   &nbsp;&nbsp;&nbsp;For testing purposes BLAS can be
!>   downloaded along with LAPCK from http://www.netlib.org/lapack.<br>
!>   &nbsp;&nbsp;&nbsp;For production runs we strongly recomend to use a tuned version
!>   like Intel-MKL or AOCL.</li>
!>   <li>LAPACK<br>
!>   &nbsp;&nbsp;&nbsp;To build PETSc. <br>
!>   &nbsp;&nbsp;&nbsp;For testing purposes LAPACK can be
!>   downloaded from http://www.netlib.org/lapak. <BR>
!>   &nbsp;&nbsp;&nbsp;For production runs we
!>   strongly recomend to use a tuned version like Intel-MKL or AOCL.
!>       <ul>
!>         <li>Build with the netlib.org version with:
!>         <pre>
!>           tar -xf lapack-3.9.0.tar.gz</li>
!>           mkdir lapack-3.9.0-build</li>
!>           cd lapack-3.9.0-build</li>  
!>           cmake ../lapack-3.9.0 -DCMAKE_INSTALL_PREFIX=<install_path></li>
!>           make install</li>  
!>         </pre>
!>         </li>
!>       </ul>
!>   </li>
!>   <li>MPI<br>
!>   &nbsp;&nbsp;&nbsp;To build PETSc.
!>     <ul>
!>       <li>mpi-compilers mpif90 and mpicc</li>
!>       <li>mpirun</li>
!>     </ul>
!>   </li>
!>   <li>METIS
!>     <ol>
!>       <li>In <metis_src_dir>/include/metis.h change
!>         <pre>
!>           #define IDXTYPEWIDTH 32 to #define IDXTYPEWIDTH 64</li>
!>           #define REALTYPEWIDTH 32 to #define REALTYPEWIDTH 64</li>
!>         </pre>
!>       </li>
!>       <li>Build with
!>         <pre>
!>           cd <path_to_metis_source>/build
!>           cmake -DGKLIB_PATH=../GKlib \
!>                 -DCMAKE_INSTALL_PREFIX=<metis_install_path> ../
!>           make install
!>         </pre>
!>       </li>
!>     </ol>
!>   </li>
!>   <li>PETSc
!>     <ul>
!>       <li> Get PETSc with
!>       <pre>git clone -b release https://gitlab.com/petsc/petsc.git petsc</pre>
!>       </li>
!>       <li> Configure with:
!>       <pre>./configure -prefix=<petsc_install_prefix>                      \
!>            --with-fortran-datatypes --with-fortran-interfaces=1            \
!>            --with-x=0  --with-64-bit-indices=1                             \
!>            CC_LINKER_FLAGS="-O3 -Wall" CFLAGS="-O3 -Wall" LDFLAGS="-O3"    \
!>            CXXFLAGS="-O3 -Wall" CXX_LINKER_FLAGS="-O3 -Wall"               \
!>            FFLAGS="-O3 -Wall" FC_LINKER_FLAGS="-O3 -Wall"                  \
!>            --with-precision=double --with-fortran-datatypes                \
!>            --with-shared-libraries=0 </pre>
!>       </li>
!>       <li> Build according to the configure instructions.
!>     </ul>
!>   </li>
!> </ol>
!>
