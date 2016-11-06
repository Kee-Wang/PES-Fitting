1. About compiler:
        ifort is Intel's commercial compiler.
        mkl is Intel's math kernal library.(This one can be downloaded for free)--`https://software.intel.com/en-us/intel-mkl`
        mkl is located at `cd /opt/intel/mkl/`

2. PGI, pgf90 is another option for compiler. `https://www.pgroup.com/purchase/freepgi.php`
        You need to download a lisence to the `/opt/pgi/` as `lisence.dat` file. A developer free version lisence can be obtained through NVIDIA OpenACC Toolkit programe.
        Ater installation and activation, you can open configured terminal in `/opt/pgi/`, for example, `open pgdbg64.app`. Then you can use commands such as `pgf90` to complie your Fortran codes.
        Or, as Mac user, you can add the following command in the Terminal>>Preference>>Profiles>>(your default terinal)>>Shell>>Startup, check Run command and Run inside shell, and enter: `setenv PATH /opt/pgi/osx86-64/2016/bin:/opt/pgi/osx86-64/2016/mpi/mpich/bin:$PATH;export PATH=/opt/pgi/osx86-64/2016/bin:/opt/pgi/osx86-64/2016/mpi/mpich/bin:$PATH; clear`

3. Use `printenv` to print global variable.
