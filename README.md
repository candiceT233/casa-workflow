# casa-workflow
Code taken from pegasus containers.

## pegasus/casa-nowcast
Container: https://hub.docker.com/r/pegasus/casa-nowcast \
Workflow: https://github.com/pegasus-isi/casa-nowcast-workflow

## pegasus/casa-wind
Container: https://hub.docker.com/r/pegasus/casa-wind \
Workflow: https://github.com/pegasus-isi/casa-wind-workflow


# TODO
Currently the test script `casa-nowcast/test/test_nowcast.sh` seems not working correctly.
### (1) NowcastToWDSS2 gives many netCDF4 warnings
```
starttime: mrtV2_co-nfig00
NetCDF: Unknown file format
NetCDF: Not a valid ID
NetCDF: Not a valid ID
NetCDF: Not a valid ID
NetCDF: Not a valid ID
NetCDF: Not a valid ID
NetCDF: Not a valid ID
NetCDF: Not a valid ID
nsteps: 0
NetCDF: Not a valid ID
NetCDF: Not a valid ID
NetCDF: NC_UNLIMITED size already in use
NetCDF: NC_UNLIMITED in the wrong index
NetCDF: Variable not found
NetCDF: Not a valid ID
NetCDF: Not a valid ID
NetCDF: Not a valid ID
NetCDF: Not a valid ID
NetCDF: Not a valid ID
NetCDF: Not a valid ID
NetCDF: Not a valid ID
NetCDF: Variable not found
nsteps: 1
NetCDF: Not a valid ID
...
```
### (1) merged_netcdf2png step seg-fault
```
./test_nowcast.sh: line 57: 13542 Segmentation fault      (core dumped) $PROGRAM_PATH/netcdf2png/merged_netcdf2png -c $INPUT_PNG -q 235 -z 0,75 -o $OUTPUT_DIR/$OUTPUT_PR_FILE $OUTPUT_DIR/$INPUT_PR_FILE
```
Inside docker container:
```
xmlEscapeEntities : char out of range
Segmentation fault (core dumped)
```

# Dependency In Container 
### libxml2
apt search libxml2
```
libxml2/now 2.9.4+dfsg1-6.1ubuntu1.2 amd64 [installed,local]
  GNOME XML library

libxml2-dev/now 2.9.4+dfsg1-6.1ubuntu1.2 amd64 [installed,local]
  Development files for the GNOME XML library
```
### jansson
ls -l /usr/local/jansson/*/*
```
-rw-r--r-- 1 root root   9288 Nov 14  2019 /usr/local/jansson/include/jansson.h
-rw-r--r-- 1 root root   1183 Nov 14  2019 /usr/local/jansson/include/jansson_config.h
-rw-r--r-- 1 root root 389008 Nov 14  2019 /usr/local/jansson/lib/libjansson.a
-rwxr-xr-x 1 root root    971 Nov 14  2019 /usr/local/jansson/lib/libjansson.la
lrwxrwxrwx 1 root root     19 Nov 14  2019 /usr/local/jansson/lib/libjansson.so -> libjansson.so.4.7.0
lrwxrwxrwx 1 root root     19 Nov 14  2019 /usr/local/jansson/lib/libjansson.so.4 -> libjansson.so.4.7.0
-rwxr-xr-x 1 root root 223936 Nov 14  2019 /usr/local/jansson/lib/libjansson.so.4.7.0
```
### libconfig
apt search libconfig
```
libconfig-dev/now 1.5-0.4 amd64 [installed,local]
  parsing/manipulation of structured config files (development)

libconfig-doc/now 1.5-0.4 all [installed,local]
  parsing/manipulation of structured config files (Documentation)

libconfig9/now 1.5-0.4 amd64 [installed,local]
  parsing/manipulation of structured configuration files
```
### netcdf
apt search netcdf
```
libnetcdf-dev/now 1:4.6.0-2build1 amd64 [installed,local]
  creation, access, and sharing of scientific data

libnetcdf13/now 1:4.6.0-2build1 amd64 [installed,local]
  Interface for scientific data access to large binary data
```
libnetcdf-dev
```
/usr/bin/nc-config
/usr/include/netcdf.h
/usr/include/netcdf_mem.h
/usr/include/netcdf_meta.h
/usr/lib/x86_64-linux-gnu/cmake/netCDF/netCDFConfig.cmake
/usr/lib/x86_64-linux-gnu/cmake/netCDF/netCDFConfigVersion.cmake
/usr/lib/x86_64-linux-gnu/cmake/netCDF/netCDFTargets-none.cmake
/usr/lib/x86_64-linux-gnu/cmake/netCDF/netCDFTargets.cmake
/usr/lib/x86_64-linux-gnu/libnetcdf.settings
/usr/lib/x86_64-linux-gnu/pkgconfig/netcdf.pc
/usr/share/doc/libnetcdf-dev/copyright
/usr/share/man/man1/nc-config.1.gz
/usr/share/man/man3/netcdf.3.gz
/usr/lib/x86_64-linux-gnu/libnetcdf.so
/usr/share/doc/libnetcdf-dev/changelog.Debian.gz
```
libnetcdf13
```
/usr/lib/x86_64-linux-gnu/libnetcdf.so.13
/usr/share/doc/libnetcdf13/changelog.Debian.gz
/usr/share/doc/libnetcdf13/copyright
/usr/share/lintian/overrides/libnetcdf13
```


