1. For Fortran compiler:
	-fdefault-real-8 -fdefault-double-8
2. When writting Python code and generate a new file. Always close a file after using it.

3. use zip -r folder.zip folder, and unzip folder to compress and uncompres files

4. Long range interaction might be importance for CO2 hydrate clathrate.

5. Bas's Library fitting routine is weighted while in MSA ther is no weight. The benefit of weighting in Bas's library is that the extremely big data won't affact too much.

6. Molpro could result in problem such as `GLOBAL ERROR fehler on processor   0`

7. About compiler:
	ifort is Intel's commercial compiler.
	mkl is Intel's math kernal library.(This one can be downloaded for free)--`https://software.intel.com/en-us/intel-mkl`
	mkl is located at `cd /opt/intel/mkl/`

8. PGI, pgf90 is another option for compiler. `https://www.pgroup.com/purchase/freepgi.php`
	You need to download a lisence to the `/opt/pgi/` as `lisence.dat` file. A developer free version lisence can be obtained through NVIDIA OpenACC Toolkit programe.
	Ater installation and activation, you can open configured terminal in `/opt/pgi/`, for example, `open pgdbg64.app`. Then you can use commands such as `pgf90` to complie your Fortran codes.
	Or, as Mac user, you can add the following command in the Terminal>>Preference>>Profiles>>(your default terinal)>>Shell>>Startup, check Run command and Run inside shell, and enter: `setenv PATH /opt/pgi/osx86-64/2016/bin:/opt/pgi/osx86-64/2016/mpi/mpich/bin:$PATH;export PATH=/opt/pgi/osx86-64/2016/bin:/opt/pgi/osx86-64/2016/mpi/mpich/bin:$PATH; clear`

9. Use `printenv` to print global variable.

10. Use `tar -vzvf 'test.tar' * ` to compress all current files.
    `tar -czvf 'test.tar.bz2' *` is even better.
	In order to unzip: `tar -xzvf archive.tar.gz`

11. Use `grep -r 'string' ./` to search for string in the all files in current folder.

12. Simply fit the feature if you have only one feature. But for a 4000 * 904 matrix, you have to think about it.
