/****************************************************************************
 *
 * (C) 2021 Cadence Design Systems, Inc. All rights reserved worldwide.
 *
 * This sample source code is not supported by Cadence Design Systems, Inc.
 * It is provided freely for demonstration purposes only.
 * SEE THE WARRANTY DISCLAIMER AT THE BOTTOM OF THIS FILE.
 *
 ***************************************************************************/

#ifndef _RTCAEPINSTANCEDATA_H_
#define _RTCAEPINSTANCEDATA_H_

// struct TAU_DATA is defined in src\plugins\CaeUnsTAU\runtimeWrite.cxx
struct TAU_DATA;

// This macro is included in the CAEP_RTITEM structure declaration:
#define CAEP_RUNTIME_INSTDATADECL   TAU_DATA *tau;

#endif /* _RTCAEPINSTANCEDATA_H_ */

/****************************************************************************
 *
 * This file is licensed under the Cadence Public License Version 1.0 (the
 * "License"), a copy of which is found in the included file named "LICENSE",
 * and is distributed "AS IS." TO THE MAXIMUM EXTENT PERMITTED BY APPLICABLE
 * LAW, CADENCE DISCLAIMS ALL WARRANTIES AND IN NO EVENT SHALL BE LIABLE TO
 * ANY PARTY FOR ANY DAMAGES ARISING OUT OF OR RELATING TO USE OF THIS FILE.
 * Please see the License for the full text of applicable terms.
 *
 ****************************************************************************/
