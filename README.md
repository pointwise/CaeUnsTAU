## Building the TAU CAE Export Plugin

To build the **TAU** CAE export plugin you must integrate the source code from 
this repository into your local PluginSDK installation.

This plugin was created with the `mkplugin` options `-c` and `-caeu`.

This plugin uses the following custom make files.
 * `modulelocal.mk`

This plugin links against the following third party libraries.
 * [netcdf][netcdf]
 * [hdf5][hdf5]
 * [zlib][zlib]

See [How To Integrate Plugin Code][HowTo] for details.

[HowTo]: https://github.com/pointwise/How-To-Integrate-Plugin-Code
[netcdf]: https://www.unidata.ucar.edu/downloads/netcdf/
[hdf5]: https://www.hdfgroup.org/
[zlib]: http://www.zlib.net/
