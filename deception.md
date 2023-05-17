# Dependencies in Deception

libxml2
- available in system and spack

jansson
- cannot find `jansson.h` abd `jansson_config.h` in system, installed with spack
- require janson version `spack install jansson@2.9`, other versions agvoe 2.0 give below error
```
./mrtV2: /usr/lib64/libjansson.so.4: no version information available (required by ./mrtV2)
```

libconfig
- cannot find `libconfig.h` in system, installed with spack
```
spack install libconfig ^texinfo@5.0
```
- the system `texinfo@5.1` packages in `/usr` causes below error:
```
WARNING: 'makeinfo' is missing on your system.
         You should only need it if you modified a '.texi' file, or
         any other file indirectly affecting the aspect of the manual.
         You might want to install the Texinfo package:
         <http://www.gnu.org/software/texinfo/>
         The spurious makeinfo call might also be the consequence of
         using a buggy 'make' (AIX, DU, IRIX), in which case you might
         want to install GNU make:
         <http://www.gnu.org/software/make/>
make[2]: *** [libconfig.info] Error 127
make[2]: Leaving directory `/scratch/tang584/spack-stage/spack-stage-libconfig-1.7.2-il2q34eyhsjprxztwnqi24rciqwn42qm/spack-src/doc'
make[1]: *** [all-recursive] Error 1
make[1]: Leaving directory `/scratch/tang584/spack-stage/spack-stage-libconfig-1.7.2-il2q34eyhsjprxztwnqi24rciqwn42qm/spack-src'
make: *** [all] Error 2
```


