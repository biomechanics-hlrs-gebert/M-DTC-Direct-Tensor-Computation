  !> \page regression Regression
  !>
  !> The cases described below are to be used for struct-process regression and
  !> correctness checks. Base data and reference results can be found in the
  !> regression directory.
  !>
  !> <ol>
  !>   <li><b><u> solid_cube </u></b><br>
  !>   A puredat datafiled representing a cube shaped solid with 21³ voxels.
  !>   Two input files are provided along with the corresponding reference
  !>   results.
  !>
  !>   The cases can be executed from the struct-process base directory.
  !>   <pre>
  !>   mpirun -n 3 struct_process_x86_64 regression/solid_cube/solid_cube_HEX08.input
  !>   diff -r regression/solid_cube/solid_cube_Hex08 regression/solid_cube_reg/solid_cube_Hex08
  !>
  !>   mpirun -n 3 struct_process_x86_64 regression/solid_cube/solid_cube_HEX20.input
  !>   diff -r regression/solid_cube/solid_cube_Hex20 regression/solid_cube_reg/solid_cube_Hex20
  !>   </pre>
  !>   </li>
  !> </ol>
  !>
