# sphinx-mac

# Building netcdf

```
cd netcdf
tar -zxvf netcdf-c-4.8.1.tar.gz
cd netcdf-c-4.8.1
#
# Building netcdf
source build_env
#
CC=clang CXX=clang++ ./configure --disable-netcdf-4 --prefix="$PWD"/install_dn --disable-byterange --disable-libxml2
make -j8 && make install
# Optionally
make check # Almost all passed, see below
```

```
PASS: t_dap3a
PASS: test_cvt3
PASS: test_vara
PASS: tst_ncdap3.sh
PASS: testpathcvt.sh
PASS: tst_ber.sh
FAIL: tst_remote3.sh
PASS: tst_formatx.sh
PASS: testurl.sh
PASS: tst_fillmismatch.sh
PASS: tst_zero_len_var.sh
PASS: tst_encode.sh
PASS: test_partvar
PASS: t_misc
PASS: test_nstride_cached
PASS: test_varm3
============================================================================
Testsuite summary for netCDF 4.8.1
============================================================================
# TOTAL: 16
# PASS:  15
# SKIP:  0
# XFAIL: 0
# FAIL:  1
# XPASS: 0
# ERROR: 0
============================================================================
See ncdap_test/test-suite.log
Please report to support-netcdf@unidata.ucar.edu
============================================================================
make[5]: *** [test-suite.log] Error 1
make[4]: *** [check-TESTS] Error 2
make[3]: *** [check-am] Error 2
make[2]: *** [check-recursive] Error 1
make[1]: *** [check] Error 2
make: *** [check-recursive] Error 1
```

# Building sphinx

I could build sphinx-3.1

Modified files: src/configure, sxaccelerate/src/configure, sxaccelerate/src/math/SxComplex.h, src/dft/SxPW.cpp

```
tar -xvf sphinx-3.1.tar.xz
cd sphinx-3.1
#
# Building sphinx
source build_env

# Modify files above
CC=clang CXX=clang++ ./configure --disable-debug --enable-mkl --enable-mklfft --with-mklpath="$INTEL_MKL"/include:"$INTEL_MKL"/lib --with-numlibs=$PWD/../netcdf/netcdf-c-4.8.1/install_dn --disable-openmp --disable-mpi --disable-numlibschecks --prefix="$PWD"/install_dn
make all -j8
make install
```

otool (ldd) output

```
otool -L sxdefectalign
sxdefectalign:
    /Users/fonari/src/sphinx-mac/sphinx-3.1/install_dn/lib/libsxext.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac/sphinx-3.1/install_dn/lib/libsxaddutil.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac/sphinx-3.1/install_dn/lib/libsxexx.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac/sphinx-3.1/install_dn/lib/libsxstruct.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac/sphinx-3.1/install_dn/lib/libsxclassic.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac/sphinx-3.1/install_dn/lib/libsxdft.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac/sphinx-3.1/install_dn/lib/libsxlgpl.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac/sphinx-3.1/install_dn/lib/libsxdirac.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac/sphinx-3.1/install_dn/lib/libsxgeom.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac/sphinx-3.1/install_dn/lib/libsxio2.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac/netcdf/netcdf-c-4.8.1/install_dn/lib/libnetcdf.19.dylib (compatibility version 20.0.0, current version 20.1.0)
    /usr/lib/libz.1.dylib (compatibility version 1.0.0, current version 1.2.12)
    /usr/lib/libcurl.4.dylib (compatibility version 7.0.0, current version 9.0.0)
    /Users/fonari/src/sphinx-mac/sphinx-3.1/install_dn/lib/libsxmath.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac/sphinx-3.1/install_dn/lib/libsxmpi.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac/sphinx-3.1/install_dn/lib/libsxio.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac/sphinx-3.1/install_dn/lib/libsxipc.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac/sphinx-3.1/install_dn/lib/libsxfs.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac/sphinx-3.1/install_dn/lib/libsxutil.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    @rpath/libmkl_intel_lp64.dylib (compatibility version 0.0.0, current version 0.0.0)
    @rpath/libmkl_core.dylib (compatibility version 0.0.0, current version 0.0.0)
    @rpath/libmkl_sequential.dylib (compatibility version 0.0.0, current version 0.0.0)
    /usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 1351.0.0)
    /System/Library/Frameworks/CoreFoundation.framework/Versions/A/CoreFoundation (compatibility version 150.0.0, current version 3208.0.0)
    /usr/lib/libc++.1.dylib (compatibility version 1.0.0, current version 1800.105.0)
```

```
~/src/sphinx-mac/sphinx-3.1/install_dn/bin/sxdefectalign --charge 3 --eps 12.9 --vdef v_elec_relaxed_vga_222_chgm3.cub --vref v_elec_pristine_222.cub --ecut 60 --pos 0.5,0.5,0.5 --relative --average 5.2 --qe
Reading mesh+cell...
mesh=135 x 135 x 135
Reading atoms...
New species Z=31
New species Z=33
Reading potential...
cell defect = [a1={21.2161,0,0},a2={0,21.2161,0},a3={0,0,21.2161}]
Reading mesh+cell...
mesh=135 x 135 x 135
Reading atoms...
New species Z=31
New species Z=33
Reading potential...
cell bulk = [a1={21.2161,0,0},a2={0,21.2161,0},a3={0,0,21.2161}]
Excess Electrons = 3 located at {10.608,10.608,10.608}
Atomic structure specified.
ng=74892
V average: -0.0268551 eV
Averaging (33.0881 points)
Averaging (33.0881 points)
Averaging (33.0881 points)
vAlign=0 eV
+-----------------------------------------------------------------------------
=== Intermediate results (unscreened) ===
Isolated energy       : 3.59048
Periodic energy       : 2.98868
Difference (Hartree)  : -0.601801
Difference (eV)       : -16.3758
+-----------------------------------------------------------------------------
Calculation performed with epsilon = 12.9
+-----------------------------------------------------------------------------
Defect correction (eV): 1.26944 (incl. screening & alignment)
```

Linux run:
```
Running a subprocess at Thu Jan 15 08:37:55 2026
/nfs/working/builds/NB/2026-2/build-029/qe-bin/run_qe sxdefectalign --charge 3 --eps 12.9 --vdef v_elec_relaxed_vga_222_chgm3.cub --vref v_elec_pristine_222.cub --ecut 60 --pos 0.5,0.5,0.5 --relative --average 5.2 --qe
Reading mesh+cell...
mesh=135 x 135 x 135
Reading atoms...
New species Z=31
New species Z=33
Reading potential...
cell defect = [a1={21.2161,0,0},a2={0,21.2161,0},a3={0,0,21.2161}]
Reading mesh+cell...
mesh=135 x 135 x 135
Reading atoms...
New species Z=31
New species Z=33
Reading potential...
cell bulk = [a1={21.2161,0,0},a2={0,21.2161,0},a3={0,0,21.2161}]
Excess Electrons = 3 located at {10.608,10.608,10.608}
Atomic structure specified.
ng=74892
V average: -0.0268551 eV
Averaging (33.0881 points)
Averaging (33.0881 points)
Averaging (33.0881 points)
vAlign=0 eV
+-----------------------------------------------------------------------------
=== Intermediate results (unscreened) ===
Isolated energy       : 3.59048
Periodic energy       : 2.98868
Difference (Hartree)  : -0.601801
Difference (eV)       : -16.3758
+-----------------------------------------------------------------------------
Calculation performed with epsilon = 12.9
+-----------------------------------------------------------------------------
Defect correction (eV): 1.26944 (incl. screening & alignment)

Zero correction energy (eV): 1.26944

All done at Thu Jan 15 08:37:57 2026
```

Results are the same.