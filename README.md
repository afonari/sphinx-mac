# sphinx-mac

# Building netcdf

```
cd netcdf
tar -zxvf netcdf-c-4.8.1.tar.gz
cd netcdf-c-4.8.1
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

```
tar -xvf sphinx-3.1.2.tar.xz
cd sphinx-3.1.2
```
