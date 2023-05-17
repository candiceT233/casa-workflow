# casa-workflow
Code taken from pegasus containers.

## pegasus/casa-nowcast
Container: https://hub.docker.com/r/pegasus/casa-nowcast \
Workflow: https://github.com/pegasus-isi/casa-nowcast-workflow

## pegasus/casa-wind
Container: https://hub.docker.com/r/pegasus/casa-wind \
Workflow: https://github.com/pegasus-isi/casa-wind-workflow


# TODO: casa-nowcast
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
### (2) merged_netcdf2png step seg-fault
```
./test_nowcast.sh: line 57: 13542 Segmentation fault      (core dumped) $PROGRAM_PATH/netcdf2png/merged_netcdf2png -c $INPUT_PNG -q 235 -z 0,75 -o $OUTPUT_DIR/$OUTPUT_PR_FILE $OUTPUT_DIR/$INPUT_PR_FILE
```
Inside docker container:
```
xmlEscapeEntities : char out of range
Segmentation fault (core dumped)
```