# NAS Workflow
## Login
### Pleiades
```bash
$ ./scripts/nas
```
### Lou
```bash
$ ./scripts/lou
```

## File Transfers
### Intra-enclave (between two NAS hosts) with [Shift Transfer Tool](https://www.nas.nasa.gov/hecc/support/kb/shift-transfer-tool-overview_300.html)
```bash
lfeX:~> shiftc --hosts=8 /nobackup/$USERNAME/dir ~/
```
For large directories (> 1 GB), archive first:
```bash
lfeX:~> shiftc --hosts=8 --create-tar /nobackup/$USERNAME/../athdf/ ./athdf.tar
```


### [Secure Unattended Proxy (SUP)](https://www.nas.nasa.gov/hecc/support/kb/entry/145)
1. Start SUP on local (remote) host:
```bash
$ eval `sup -s bash -u $USERNAME -ols=--color=always`
```
2. Authorize front end (FE) host:
```bash
pfeXX:~> touch ~/.meshrc
```
3. Authorize writable directories:
```bash
pfeXX:~> echo /nobackup/$USERNAME >> ~/.meshrc
```
4. Execute local commands:
```bash
$ sup shiftc -r lfe:/nobackup/ /mnt/c/Users/xzele/Downloads/sup/
```


## [Archiving to Lou](https://www.nas.nasa.gov/hecc/support/kb/using-shift-for-transfers-and-tar-operations-between-two-nas-hosts_513.html)
- Check quota status:
```bash
lfs quota -h -u $USERNAME /nobackupp12
```
- Snapshotting:
```bash
lfeX:~$ mkdir snapshots/$(date +"%Y-%m-%d")
lfeX:~$ cd /nobackup/$USERNAME
lfeX:/nobackup/$USERNAME$ shiftc --hosts=8 --create-tar --index-tar github lfe:~/snapshots/$(date +"%Y-%m-%d")/github.tar
```
- Archiving:
```bash
lfeX:~$ mkdir archives/.../dir
lfeX:~$ cd /nobackup/$USERNAME
lfeX:/nobackup/$USERNAME$ shiftc --hosts=8 --create-tar --index-tar github/.../dir lfe:~/archives/.../dir/dir.$(date +"%Y-%m-%d").tar
```


## [PBS](https://www.nas.nasa.gov/hecc/support/kb/running-jobs-with-pbs-121/)
### Interactive Jobs
Using [VNC Xterm](https://www.nas.nasa.gov/hecc/support/kb/vnc-a-faster-alternative-to-x11_257.html): 
1. Start VNC server on PFE:
```bash
pfeXX:~> vncserver -localhost
```
2. Copy-paste this line into PFE prompt: `~C`
3. In FE prompt: `-L 2222X:localhost:592X` (or any available local port)
4. Locally (same FE port above): 
```bash
$ vncviewer localhost:2222X
``` 
5. Submit job:
```bash
pfeXX:~> qsub -I -X -lselect=1:ncpus=16:mpiprocs=16:model=san,walltime=1:00:00 -q devel
```
6. When finished, before logging off, make sure to kill the VNC server:
```bash
pfeXX:~> vncserver -kill :XX
```
where `XX` is the ID of the server initiated earlier.

### [Batch Jobs](https://www.nas.nasa.gov/hecc/support/kb/sample-pbs-script-for-pleiades_190.html)
See [`sample.pbs`](/nas/sample.pbs).


## [Athena++](https://github.com/PrincetonUniversity/athena/wiki)
### [Configure](https://github.com/PrincetonUniversity/athena/wiki/Configuring)
```bash
pfeXX:~> ./configure.py --prob=streaming_instability -p --eos=isothermal --nghost=3 -hdf5 -h5double -mpi --cxx=icpc -mpi --mpiccmd="icpc -lmpi -lmpi++" --cflag="-xCORE-AVX512"
```

### [Compile](https://github.com/PrincetonUniversity/athena/wiki/Compiling)
#### On Pleiades front end (PFE)
```bash
$ make -j 2
```

## Local Commands
### Visit
- In local directory
  ```bash
  visit -np 4 -sessionfile <fname>
  ```
