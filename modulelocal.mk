#**************************************************************************
# 
#  Copyright (c) 2012-2020 Pointwise, Inc.
#  All rights reserved.
# 
#  This sample Pointwise plugin is not supported by Pointwise, Inc.
#  It is provided freely for demonstration purposes only.
#  SEE THE WARRANTY DISCLAIMER AT THE BOTTOM OF THIS FILE.
# 
# modulelocal.mk for src\plugins\CaeUnsTAU plugin
#**************************************************************************

#-----------------------------------------------------------------------
# OPTIONAL, locally defined plugin make targets. If a plugin developer wants
# to extend a plugin's make scheme, they should create a modulelocal.mk file
# in the plugin's base folder. To provide for future SDK upgrades, the standard
# module.mk file should NOT be directly edited.
#
# Uncomment, copy and/or edit the sections below as needed.
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# Adds plugin specific source files to the build.
# These files will be compiled and then linked to the plugin.
# The file paths are relative to the project path macro CaeUnsTAU_LOC.
# For example,
#    myFile.cxx is located in $(CaeUnsTAU_LOC)/myFile.cxx
#    sub/myOtherFile.cxx is located in $(CaeUnsTAU_LOC)/sub/myOtherFile.cxx
#
#CaeUnsTAU_CXXFILES_PRIVATE := \
#    yourfile1.cxx \
#    yourfile2.cxx \
#	$(NULL)

#-----------------------------------------------------------------------
# Adds plugin specific include flags to the build.
#
CaeUnsTAU_INCL_PRIVATE := \
	-I./external/common/netcdf \
	$(NULL)

#-----------------------------------------------------------------------
# Adds plugin specific compile flags to the build.
#
#CaeUnsTAU_CXXFLAGS_PRIVATE := \
#	$(NULL)


#-----------------------------------------------------------------------
# Adds plugin specific -lfile and -Llibpath to build.
# These flags will be added to the link.
#
CaeUnsTAU_LIBS_PRIVATE := \
	-L./external/$(machine)/netcdf \
	-lnetcdf \
	-L./external/$(machine)/hdf5 \
	-lhdf5_hl \
	-lhdf5 \
	-lz \
	$(NULL)

#-----------------------------------------------------------------------
# Adds plugin specific link flags to the build.
#
#CaeUnsTAU_LDFLAGS_PRIVATE := \
#	$(NULL)

#-----------------------------------------------------------------------
# Add any locally defined targets that do NOT produce an output object.
# This would include targets used for cleaning, printing, etc. These targets
# will be automatically added the .PHONY target.
#
#CaeUnsTAU_MAINT_TARGETS_PRIVATE = \
#	$(NULL)

#-----------------------------------------------------------------------
# Sample macro. Prefix with CAE name to prevent conflicts.
#
#CaeUnsTAU_DUMMY = \
#	DUMMY1 \
#	DUMMY2 \
#	$(NULL)

#-----------------------------------------------------------------------
# Sample target. Prefix with CAE name to prevent conflicts.
#
#CaeUnsTAU_sample: CaeUnsTAU_clean CaeUnsTAU
#	@echo "done building CaeUnsTAU"

#****************************************************************************
# 
#  DISCLAIMER:
#  TO THE MAXIMUM EXTENT PERMITTED BY APPLICABLE LAW, POINTWISE DISCLAIMS
#  ALL WARRANTIES, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED
#  TO, IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
#  PURPOSE, WITH REGARD TO THIS SCRIPT. TO THE MAXIMUM EXTENT PERMITTED
#  BY APPLICABLE LAW, IN NO EVENT SHALL POINTWISE BE LIABLE TO ANY PARTY
#  FOR ANY SPECIAL, INCIDENTAL, INDIRECT, OR CONSEQUENTIAL DAMAGES
#  WHATSOEVER (INCLUDING, WITHOUT LIMITATION, DAMAGES FOR LOSS OF
#  BUSINESS INFORMATION, OR ANY OTHER PECUNIARY LOSS) ARISING OUT OF THE
#  USE OF OR INABILITY TO USE THIS SCRIPT EVEN IF POINTWISE HAS BEEN
#  ADVISED OF THE POSSIBILITY OF SUCH DAMAGES AND REGARDLESS OF THE
#  FAULT OR NEGLIGENCE OF POINTWISE.
# 
# ***************************************************************************
