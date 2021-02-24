# Configure
`python configure.py --prob=streaming_instability -p --eos=isothermal --nghost=3 -hdf5 --hdf5_path /usr/lib/x86_64-linux-gnu/hdf5/openmpi -h5double -mpi`

# Run
`mpiexec -n 4 $ATHENA_DUST -i athinput.si > log`
