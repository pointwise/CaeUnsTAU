# Building the TAU CAE Export Plugin
Copyright 2021 Cadence Design Systems, Inc. All rights reserved worldwide.

To build the **TAU** CAE export plugin you must integrate the source code from 
this repository into your local PluginSDK installation.

This plugin was created with the `mkplugin` options `-c` and `-caeu`.

See [How To Integrate Plugin Code][HowTo] for details.

This plugin statically links against the following third party libraries. Newer
versions may work but have not been tested.
 * [netcdf v4.6.1][netcdf]
 * [hdf5 v1.8.20][hdf5]
 * [zlib v1.2.8][zlib]
   * Usually, zlib is only downloaded for win64 builds. The `libz` library is already installed with most linux and macosx distibutions.

You must add the appropriate include path and link library settings to your TAU plugin project when building with Visual Studio on win64.

This plugin uses a custom `modulelocal.mk` file when building on the linux and macosx platforms.

`modulelocal.mk` expects that `machine` is defined on the `make` command line or in the shell's environment. Where `machine` is one of `macosx` or `linux_x86_64`.

For example,
 * `make machine=value`
 * `setenv machine value`   (on csh)
 * `export machine=value`   (on sh or bash)

`modulelocal.mk` expects the link libraries and headers to be located as shown below. See `PluginSDK/external/README_EXTERNAL.txt` for more details.
```
PluginSDK/
   external/
      common/
         hdf5/
            hdf5.h
            hdf5_hl.h
            H5*.h    <-- many the other hdf5 common headers
         netcdf/
            netcdf.h
         zlib/
            zlib.h
            zconf.h
      linux_x86_64/
         netcdf/
            libnetcdf.lib
         hdf5/
            libhdf5.a
            libhdf5_hl.a
            H5pubconf.h
      macosx/
         netcdf/
            libnetcdf.a
         hdf5/
            libhdf5.a
            libhdf5_hl.a
            H5pubconf.h
      win64/
         netcdf/
            netcdf.lib
         hdf5/
            hdf5.lib
            hdf5_hl.lib
            H5pubconf.h
         zlib/
            zlib.lib
```

[HowTo]: https://github.com/pointwise/How-To-Integrate-Plugin-Code
[netcdf]: https://www.unidata.ucar.edu/downloads/netcdf/
[hdf5]: https://www.hdfgroup.org/
[zlib]: http://www.zlib.net/

## Disclaimer
This file is licensed under the Cadence Public License Version 1.0 (the "License"), a copy of which is found in the LICENSE file, and is distributed "AS IS." 
TO THE MAXIMUM EXTENT PERMITTED BY APPLICABLE LAW, CADENCE DISCLAIMS ALL WARRANTIES AND IN NO EVENT SHALL BE LIABLE TO ANY PARTY FOR ANY DAMAGES ARISING OUT OF OR RELATING TO USE OF THIS FILE. 
Please see the License for the full text of applicable terms.
