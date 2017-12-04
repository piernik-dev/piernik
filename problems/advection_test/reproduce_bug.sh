./setup advection_test -o BUG_advection -d MAGNETIC,IONIZED
cd runs/advection_test_BUG_advection
./piernik -n '&OUTPUT_CONTROL dt_res = 1.0 /'
mpirun -np 5 ./piernik -n '&END_CONTROL tend = 2. / &OUTPUT_CONTROL vars(:) = "magx", "magy"/'

# The problem is that not all magnetic field gets properly initialized from the restart file for certain thread counts.
# if occurs also when the second call is
# mpirun -np 3 ./piernik -n '&END_CONTROL tend = 2. / &MPI_BLOCKS dd_rect_quality = 0./ &OUTPUT_CONTROL vars(:) = "magx", "magy"/'
