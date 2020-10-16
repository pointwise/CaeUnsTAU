/****************************************************************************
 *
 * Copyright (c) 2012-2020 Pointwise, Inc.
 * All rights reserved.
 *
 * This sample Pointwise plugin is not supported by Pointwise, Inc.
 * It is provided freely for demonstration purposes only.
 * SEE THE WARRANTY DISCLAIMER AT THE BOTTOM OF THIS FILE.
 *
 ***************************************************************************/

#ifndef _RTCAEPSUPPORTDATA_H_
#define _RTCAEPSUPPORTDATA_H_

/*------------------------------------*/
/* CaeUnsTAU format item setup data */
/*------------------------------------*/
// mapping between names and ids here should match the mapping used in the 
// defineMarkers function in runtimeWrite.cxx (although this struct should
// exclude unspecified, as it is available by default in pointwise).
//
// Also, due to a lack of clarity in the file-type specs, it is best to ensure
// that the BC ids start at 1 and increment by one. However, their order does
// not matter (as long as it matches the order in runtimeWrite.cxx).
CAEP_BCINFO CaeUnsTAUBCInfo[] = {
    { "symmetry plane", 1 },
    { "axisymmetry axis", 2 },
    { "farfield", 3 },
    { "supersonic inflow", 4 },
    { "supersonic outflow", 5 },
    { "reservoir-pressure inflow", 6 },
    { "exit-pressure outflow", 7 },
    { "dirichlet", 8 },
    { "euler wall", 9 },
    { "sharp edge", 10 },
    { "viscous wall", 11 },
    { "laminar wall", 12 },
    { "turbulent wall", 13 },
    { "engine exhaust", 14 },
    { "engine inflow", 15 },
    { "heat exchanger inflow", 16 },
    { "heat exchanger outflow", 17 },
    { "actuator exhaust", 18 },
    { "actuator inflow", 19 },
    { "chimera", 20 },
    { "actuation", 21 },
    { "mirror plane", 22 },
    { "periodic plane", 23 },
};


// Volume conditions not used by TAU
//CAEP_VCINFO CaeUnsTAUVCInfo[] = {
//    { "viscous-CaeUnsTAU", 200 },
//    { "invisid-CaeUnsTAU", 201 },
//};

const char *CaeUnsTAUFileExt[] = {
    "grid",
    "bmap"
};

#endif /* _RTCAEPSUPPORTDATA_H_ */

/****************************************************************************
 *
 * DISCLAIMER:
 * TO THE MAXIMUM EXTENT PERMITTED BY APPLICABLE LAW, POINTWISE DISCLAIMS
 * ALL WARRANTIES, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED
 * TO, IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE, WITH REGARD TO THIS SCRIPT. TO THE MAXIMUM EXTENT PERMITTED
 * BY APPLICABLE LAW, IN NO EVENT SHALL POINTWISE BE LIABLE TO ANY PARTY
 * FOR ANY SPECIAL, INCIDENTAL, INDIRECT, OR CONSEQUENTIAL DAMAGES
 * WHATSOEVER (INCLUDING, WITHOUT LIMITATION, DAMAGES FOR LOSS OF
 * BUSINESS INFORMATION, OR ANY OTHER PECUNIARY LOSS) ARISING OUT OF THE
 * USE OF OR INABILITY TO USE THIS SCRIPT EVEN IF POINTWISE HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGES AND REGARDLESS OF THE
 * FAULT OR NEGLIGENCE OF POINTWISE.
 *
 ***************************************************************************/
