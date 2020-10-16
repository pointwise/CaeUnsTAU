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


 /****************************************************************************
 * CaeUnsTAU implementation of runtimeWrite(), runtimeCreate(), and
 * runtimeDestroy()
 * 
 * more info in SPR 4411
 * 
 * TAU require an external library, netcdf, that has been linked in.
 * The API for netcdf can be found here: 
 * http:// www.unidata.ucar.edu/software/netcdf/docs/netcdf-c/
 * In addition, netcdf links against hdf5 and zlib, so those binaries must be
 * included.
 * 
 * In netcdf, a set of 'dimensions' (referred to by a unique integer 'dimid') is
 * used to create a 'variable' (referred to by a unique integer 'varid'). Both
 * dimensions and variables are given unique string names. In other words, a
 * variable can be thought of as a multi-dimensional array, and a variable gets
 * associated with dimensions that indicate the variable's size.
 * 
 * All functions in this module take as input a pointer to a CAEP_RTITEM 
 * struct (most make use the contained tau struct, defined in
 * rtCaepInstanceData.h).
 * 
 * All functions in this module return true on success and false if an error
 * occurs.
 *
 ***************************************************************************/

#include "apiCAEP.h"
#include "apiCAEPUtils.h"
#include "apiGridModel.h"
#include "apiPWP.h"
#include "runtimeWrite.h"
#include "pwpPlatform.h"

#include <assert.h>

// DO NOT CHANGE THE ORDER OF INCLUSION FOR THE FILES BELOW!
// <string>, <map>, and <netcdf.h> MUST be included in the order they are below
#include <string.h>
#include <string>
#include <map>
#include <netcdf.h>


typedef std::map<std::string, PWP_UINT32>   StringUINT32Map;
typedef std::map<std::string, std::string>  StringStringMap;

struct TAU_DATA {
    // Main dimensions
    unsigned int    numPoints;
    unsigned int    numElements;
    unsigned int    numTets;
    unsigned int    numPyramids;
    unsigned int    numWedges;
    unsigned int    numHexes;
    unsigned int    numBars;
    unsigned int    numSurfaceElements;
    unsigned int    numSurfaceTris;
    unsigned int    numSurfaceQuads;
    int             numSurfBndryMarkers;

    // Marker values
    int             idNc;
    int             idXc;
    int             idYc;
    int             idZc;
    int             idTets;
    int             idPyramids;
    int             idPrisms;
    int             idHexes;
    int             markerId;
    int             surfaceTriId;
    int             surfaceQuadId;
    int             bndryPanelStart;
    int             bndryPanelTriStart;
    int             bndryPanelQuadStart;
    StringUINT32Map mapNameToOutputID;
    StringStringMap mapNameToType;

    // Used to track the indexing of hex and wedge creation during 2D export.
    int             indexHex;
    int             indexWedge;
};

static const char *NetCDFDataModel = "DataModel";
static const char *Thickness = "Thickness";

#define OK(err)   (NC_NOERR == (err))


/****************************************************************************
 * 
 * Initialize global condition information. It fills in the mapNameToOutputID
 * map and the mapNameToType map. For both, the keys are given by a value of
 * condition.name, while the values are given by the corresponding call to
 * condition.tid and condition.type (respectively)
 * 
 ***************************************************************************/
static PWP_BOOL
setupMarkerInfo(CAEP_RTITEM *pRti)
{
    TAU_DATA &tau = *(pRti->tau);
    PWP_BOOL ret = PWP_TRUE;
    StringUINT32Map &mapIds = tau.mapNameToOutputID;
    StringStringMap &mapTypes = tau.mapNameToType;
    PWP_UINT32 nextID = 0;
    PWP_UINT32 numDoms = PwModDomainCount(pRti->model);
    for (PWP_UINT32 i = 0; i < numDoms && ret; ++i) {
        PWGM_HDOMAIN dom = PwModEnumDomains(pRti->model, i);
        PWGM_CONDDATA condition = {0};
        ret = ret && PwDomCondition(dom, &condition);
        if (ret && mapIds.find(condition.name) == mapIds.end()) {
            mapIds[condition.name] = nextID;
            nextID++;
            mapTypes[condition.name] = condition.type;
        }
    }
    // These are set for the 2D export only. The 2D faces will have the boundary
    // condition of symmetry. These are the quads and tris only.
    if ((ret && mapIds.find("Symmetry_1") == mapIds.end()) 
        && CAEPU_RT_DIM_2D(pRti)) {
            mapIds["Symmetry_1"] = nextID;
            nextID++;
            mapTypes["Symmetry_1"] = "Symmetry";
        }
    if ((ret && mapIds.find("Symmetry_2") == mapIds.end()) 
        && CAEPU_RT_DIM_2D(pRti)) {
            mapIds["Symmetry_2"] = nextID;
            mapTypes["Symmetry_2"] = "Symmetry";
        }
    return ret;
}


/****************************************************************************
 * 
 * This function opens the file with a call to nc_create(), and creates the
 * global "type" attribute as requested by the TAU file format.
 * 
 ***************************************************************************/
static PWP_BOOL
startGridFile(CAEP_RTITEM *pRti)
{
    TAU_DATA &tau = *(pRti->tau);
    char file[1024];
    sprintf(file, "%s.grid", pRti->pWriteInfo->fileDest);

    // NetCDF_Classic|NetCDF_64bit|NetCDF4_HDF5|
    //              0|           1|           2|
    PWP_UINT netCDFDataModel = 0;
    PwModGetAttributeUINT(pRti->model, NetCDFDataModel, &netCDFDataModel);
    int dataModel = 0;
    switch (netCDFDataModel) {
    case 1:
        dataModel = NC_64BIT_OFFSET;
        break;
    case 2:
        dataModel = NC_NETCDF4;
        break;
    case 0:
    default:
        dataModel = 0;
        break;
    }

    return OK(nc_create(file, NC_CLOBBER | dataModel, &tau.idNc)) &&
        OK(nc_put_att_text(tau.idNc, NC_GLOBAL, "type", 24,
        "Primary Grid: Tau Format"));
}


/****************************************************************************
 * 
 * Creates the "no_of_points" dimension, and the points_xc, points_yc, and
 * points_zc variables that are based on it. These variables are used to store
 * the xyz coordinates of Pointwise vertices.
 * 
 ***************************************************************************/
static PWP_BOOL
definePointVars(CAEP_RTITEM *pRti)
{
    TAU_DATA &tau = *(pRti->tau);
    int idNc = tau.idNc;
    PWP_UINT32 numPoints = PwModVertexCount(pRti->model);
    int idDim;
    tau.numPoints = numPoints;
    // If in 2D mode, the number of points is doubled.
    // The points are doubled because the 2D faces will be extruded one step
    // in either the positive z direction or the negative z direction. This will
    // implicitly create a second set of points that are the exact same order as
    // the original set, simply offset by the number of original points. The
    // original set has 8 points 0-7 and the new set will have 8 points, 8-15
    if(CAEPU_RT_DIM_2D(pRti))
    {
        numPoints *= 2;
    }
    
    // Create the no_of_points dimension and the XYZ coordinate vars
    return OK(nc_def_dim(idNc, "no_of_points", numPoints, &idDim)) &&
        OK(nc_def_var(idNc, "points_xc", NC_DOUBLE, 1, &idDim, &tau.idXc)) &&
        OK(nc_def_var(idNc, "points_yc", NC_DOUBLE, 1, &idDim, &tau.idYc)) &&
        OK(nc_def_var(idNc, "points_zc", NC_DOUBLE, 1, &idDim, &tau.idZc));
}


/****************************************************************************
 * 
 * Unlike Pointwise, the TAU file-format stores different shapes in different
 * variable arrays. The defineElementVars() function determines the dimensions
 * of these variables by retrieving the counts of different element types for
 * each Pointwise block and computing the sum. The variables are then created.
 * 
 * important translations:
 * tetraeder: tetrahedron
 * pyramid: pyramid
 * prism: wedge
 * hexaeder: hexahedron
 * 
 ***************************************************************************/
static PWP_BOOL
defineElementVars(CAEP_RTITEM *pRti)
{
    TAU_DATA &tau = *(pRti->tau);
    PWP_BOOL ret = PWP_TRUE;
    int idNc = tau.idNc;

    // count number of different element types
    PWP_UINT32 numElements = 0;
    PWP_UINT32 numTets = 0;
    PWP_UINT32 numPyramids = 0;
    PWP_UINT32 numPrisms = 0;
    PWP_UINT32 numHexes = 0;

    PWP_UINT32 numBlocks = PwModBlockCount(pRti->model);
    for (PWP_UINT32 i = 0; i < numBlocks; ++i) {
        PWGM_HBLOCK block = PwModEnumBlocks(pRti->model, i);
        PWGM_ELEMCOUNTS elementCounts;
        numElements += PwBlkElementCount(block, &elementCounts);
        numTets += elementCounts.count[PWGM_ELEMTYPE_TET];
        numPyramids += elementCounts.count[PWGM_ELEMTYPE_PYRAMID];
        numPrisms += elementCounts.count[PWGM_ELEMTYPE_WEDGE];
        numHexes += elementCounts.count[PWGM_ELEMTYPE_HEX];
    }

    // copy counts to pRti for use later.
    tau.numElements = numElements;
    tau.numTets = numTets;
    tau.numPyramids = numPyramids;
    tau.numWedges = numPrisms;
    tau.numHexes = numHexes;

    // note that idElemDim is not used in a variable; this is requested by the
    // file-specs
    int idElemDim = 0;
    ret = OK(nc_def_dim(idNc, "no_of_elements", numElements, &idElemDim));

    // make the two dimensions and variable for each shape type. Shape types
    // with no instances do not get created (as in file-spec)
    if (ret && numTets != 0) {
        // create # of tetraeders dimension
        int idNumTets = 0;
        ret = ret && OK(nc_def_dim(idNc, "no_of_tetraeders", numTets,
            &idNumTets));

        // create points per tetraeder dimension
        int idPtsPerTet = 0;
        ret = ret && OK(nc_def_dim(idNc, "points_per_tetraeder", 4,
            &idPtsPerTet));

        // create the points-of-tetraeders variable based on its two dimensions
        int dimids[2] = { idNumTets, idPtsPerTet };
        ret = ret && OK(nc_def_var(idNc, "points_of_tetraeders", NC_INT, 2,
            dimids, &tau.idTets));
    }
    if (ret && numPyramids != 0) {
        // see comments in 'if (numTets != 0)' section
        int idNumPyr = 0;
        ret = ret && OK(nc_def_dim(idNc, "no_of_pyramids", numPyramids,
            &idNumPyr));

        int idPtsPerPyr = 0;
        ret = ret && OK(nc_def_dim(idNc, "points_per_pyramid", 5,
            &idPtsPerPyr));

        int dimids[2] = { idNumPyr, idPtsPerPyr };
        ret = ret && OK(nc_def_var(idNc, "points_of_pyramids", NC_INT, 2,
            dimids, &tau.idPyramids));
    }
    if (ret && numPrisms != 0) {
        // see comments in 'if (numTets != 0)' section
        int idNumPrism = 0;
        ret = ret && OK(nc_def_dim(idNc, "no_of_prisms", numPrisms,
            &idNumPrism));

        int idPtsPerPrism = 0;
        ret = ret && OK(nc_def_dim(idNc, "points_per_prism", 6,
            &idPtsPerPrism));

        int dimids[2] = { idNumPrism, idPtsPerPrism };
        ret = ret && OK(nc_def_var(idNc, "points_of_prisms", NC_INT, 2, dimids,
            &tau.idPrisms));
    }
    if (ret && numHexes != 0) {
        // see comments in 'if (numTets != 0)' section
        // (hexaeders are 'hexes' in pointwise)
        int idNumHex = 0;
        ret = ret && OK(nc_def_dim(idNc, "no_of_hexaeders", numHexes,
            &idNumHex));

        int idPtsPerHex = 0;
        ret = ret && OK(nc_def_dim(idNc, "points_per_hexaeder", 8,
            &idPtsPerHex));

        int dimids[2] = { idNumHex, idPtsPerHex };
        ret = ret && OK(nc_def_var(idNc, "points_of_hexaeders", NC_INT, 2,
            dimids, &tau.idHexes));
    }

    return ret;
}


/****************************************************************************
 * 
 * The defineOneMarker(CAEP_RTITEM) function
 * 
 * helper for defineMarkers(CAEP_RTITEM) function.
 *
 * Defines a marker attribute based on the condition name and id of a marker.
 * 
 ***************************************************************************/
PWP_BOOL
defineOneMarker(CAEP_RTITEM *pRti, const char *name, const int id)
{
    TAU_DATA &tau = *(pRti->tau);
    PWP_BOOL ret = PWP_TRUE;
    int idNc = tau.idNc;
    char attributeName[10];
    ret = 0 < sprintf(attributeName, "marker_%d", id);
    size_t len = strlen(name) + 1;// add one for null character
    return ret && OK(nc_put_att_text(idNc, NC_GLOBAL, attributeName, len,
        name));
}


/****************************************************************************
 * 
 * The defineMarkers(CAEP_RTITEM) function
 * 
 * Defines the marker variable and marker attributes (attributes are defined
 * using the 'defineOneMarker()' function).
 * 
 ***************************************************************************/
PWP_BOOL
defineMarkers(CAEP_RTITEM *pRti)
{
    TAU_DATA &tau = *(pRti->tau);
    PWP_BOOL ret = PWP_TRUE;
    int idNc = tau.idNc;
    int numMarkers = (int)tau.mapNameToOutputID.size();
    int ncNumMarkers = 0;

    ret = ret && OK(nc_def_dim(idNc, "no_of_markers", numMarkers,
        &ncNumMarkers));

    int dimids[1] = { ncNumMarkers };
    ret = ret && OK(nc_def_var(idNc, "marker", NC_INT, 1, dimids,
        &(tau.markerId)));

    StringUINT32Map &nameMap = tau.mapNameToOutputID;
    StringUINT32Map::iterator iter;
    for (iter = nameMap.begin(); iter != nameMap.end() && ret; ++iter) {
        const char *key = iter->first.c_str();
        PWP_UINT32 value = iter->second;
        ret = ret && defineOneMarker(pRti, key, value);
    }

    return ret;
}


/****************************************************************************
 * 
 * Domain-elements are stored in a similar manner to block-elements in the 
 * TAU file-type. As with block-elements, tris and quads are stored in
 * different variable arrays. This function iterates through the domains and
 * sums their triangle and quadrilateral counts in order to make the
 * "no_of_surfacetriangles" and "no_of_surfacequadrilaterals" dimensions, and
 * their associated variables.
 * 
 ***************************************************************************/
static PWP_BOOL
defineSurfaceElementVars(CAEP_RTITEM *pRti)
{
    TAU_DATA &tau = *(pRti->tau);
    PWP_BOOL ret = PWP_TRUE;
    int idNc = tau.idNc;

    // count number of different surface-element types
    PWP_UINT32 numSurfaceElements = 0;
    PWP_UINT32 numSurfaceTriangles = 0;
    PWP_UINT32 numSurfaceQuadrilaterals = 0;
    PWP_UINT32 numDomains = PwModDomainCount(pRti->model);
    for (PWP_UINT32 i = 0; i < numDomains; ++i) {
        PWGM_HDOMAIN domain = PwModEnumDomains(pRti->model, i);
        PWGM_ELEMCOUNTS elementCounts;
        PwDomElementCount(domain, &elementCounts);
        numSurfaceTriangles += elementCounts.count[PWGM_ELEMTYPE_TRI];
        numSurfaceQuadrilaterals += elementCounts.count[PWGM_ELEMTYPE_QUAD];
    }
    // DO NOT use return value from pwDomElementCount, as it includes bar 
    // elements as surface elements
    numSurfaceElements = numSurfaceTriangles + numSurfaceQuadrilaterals;

    tau.numSurfaceElements = numSurfaceElements;

    // create dimension for total number of surface elements
    int ncNumSurfaceElements = 0;
    ret = ret && OK(nc_def_dim(idNc, "no_of_surfaceelements",
        numSurfaceElements, &ncNumSurfaceElements));

    // create variables (and related dimensions) for each type of surface
    // element
    if (ret && numSurfaceTriangles != 0) {
        // create an id for the no_of_surfacetriangles dimension
        int ncSurfaceTriangles = 0;
        ret = ret && OK(nc_def_dim(idNc, "no_of_surfacetriangles",
            numSurfaceTriangles, &ncSurfaceTriangles));

        // create an id for the points_per_surfacetriangle dimension
        int ncPointsPerTri = 0;
        ret = ret && OK(nc_def_dim(idNc, "points_per_surfacetriangle", 3,
            &ncPointsPerTri));

        // create an id for the points_of_surfacetriangles variable, based on
        // the no_of_surfacetriangles and points_per_surfacetriangle
        // dimensions.
        int dimids[2] = { ncSurfaceTriangles, ncPointsPerTri };
        ret = ret && OK(nc_def_var(idNc, "points_of_surfacetriangles", NC_INT,
            2, dimids, &(tau.surfaceTriId)));
    }
    if (ret && numSurfaceQuadrilaterals != 0) {
        // see 'if (numSurfaceTriangles != 0)' comments
        int ncSurfaceQuads = 0;
        ret = ret && OK(nc_def_dim(idNc, "no_of_surfacequadrilaterals",
            numSurfaceQuadrilaterals, &ncSurfaceQuads));

        int ncPointsPerQuad = 0;
        ret = ret && OK(nc_def_dim(idNc, "points_per_surfacequadrilateral", 4,
            &ncPointsPerQuad));

        int dimids[2] = { ncSurfaceQuads, ncPointsPerQuad };
        ret = ret && OK(nc_def_var(idNc, "points_of_surfacequadrilaterals",
            NC_INT, 2, dimids, &(tau.surfaceQuadId)));
    }

    // Define variable (netCDF multi-dimensional array) that associates
    // domains (boundary panels) with surface elements
    int dimids[1] = { ncNumSurfaceElements };

    // Define variable for boundarymarkers
    ret = ret && OK(nc_def_var(idNc, "boundarymarker_of_surfaces", NC_INT, 1,
        dimids, &(tau.numSurfBndryMarkers)));

    // indexing scheme in the boundarypanel variable and in the 
    // boundarymarker variable saves the domain (boundary-panel) or 
    // boundary condition (boundary-marker) number of all the triangles
    // followed by the domain number of all the quads, so this indexing
    // tells us where these are.
    tau.bndryPanelTriStart = 0;
    tau.bndryPanelQuadStart = numSurfaceTriangles;
    
    return ret;
}


/****************************************************************************
 * 
 * The number of bars is calculated in order to determine the number of
 * surface quads for 2D export only. This function iterates through the domains
 * of each block and calculates the number of bar elements. This value is then
 * returned to the calling function.
 * 
 ***************************************************************************/
static PWP_UINT32
defineDomainElements(CAEP_RTITEM *pRti)
{
    PWP_UINT32 totalBars = 0;
    PWP_UINT32 numDomains = PwModDomainCount(pRti->model);
    for (PWP_UINT32 i = 0; i < numDomains; ++i) {
        PWGM_HDOMAIN domain = PwModEnumDomains(pRti->model, i);
        PWGM_ELEMCOUNTS elementCounts;
        PwDomElementCount(domain, &elementCounts);
        totalBars += elementCounts.count[PWGM_ELEMTYPE_BAR];
    }
    return totalBars;
}


/****************************************************************************
 * 
 * In a 2D Tau export, the tris and quad elements of the blocks will be
 * extruded into prisms and tets respectively. The total number of tris will
 * be double the number of tris in the Pointwise grid. The number of quads will
 * be double as well but with the addition of the extruded bars.
 * 
 ***************************************************************************/
static PWP_BOOL
defineElementVars2D(CAEP_RTITEM *pRti)
{
    TAU_DATA &tau = *(pRti->tau);
    PWP_BOOL ret = PWP_TRUE;
    int idNc = tau.idNc;

    // count number of different surface-element types
    PWP_UINT32 numSurfaceElements = 0;
    PWP_UINT32 numSurfaceTriangles = 0;
    PWP_UINT32 numSurfaceQuadrilaterals = 0;
    PWP_UINT32 numHexes = 0;
    PWP_UINT32 numPrisms = 0;
    PWP_UINT32 numBlocks = PwModBlockCount(pRti->model);

    for (PWP_UINT32 i = 0; i < numBlocks; ++i) {
        PWGM_HBLOCK block = PwModEnumBlocks(pRti->model, i);
        PWGM_ELEMCOUNTS elementCounts;
        PwBlkElementCount(block, &elementCounts);

        // Multiplying each set of elements by two because the 2D grid will
        // be virtually extruded, increasing the number of quads and tris by
        // at least double. Each bar will be extruded to a quad requiring the
        // calculation of all the bar elements.
        numSurfaceTriangles += 2 * elementCounts.count[PWGM_ELEMTYPE_TRI];
        numSurfaceQuadrilaterals += 2 * elementCounts.count[PWGM_ELEMTYPE_QUAD];

        // When extruded a single step, the number of hexes and prisms will each
        // be equal to the number of original quads and tris respectively.
        numHexes += elementCounts.count[PWGM_ELEMTYPE_QUAD]; 
        numPrisms += elementCounts.count[PWGM_ELEMTYPE_TRI];
    }
    // numBars is the total number of bars that when extruded, will become
    // surface quads
    PWP_UINT32 numBars=defineDomainElements(pRti);

    // All bars will become surface quadrilaterals
    numSurfaceQuadrilaterals += numBars;

    // Normal sum would be:
    // numSurfaceTriangles + numSurfaceQuadrilaterals + numBars; But numBars
    // is already added into the numSurfaceQuadrilaterals. This is done because
    // the bars will be turned into quads upon extruding one step
    numSurfaceElements = numSurfaceTriangles + numSurfaceQuadrilaterals;

    tau.numSurfaceElements = numSurfaceElements;
    tau.numSurfaceTris = numSurfaceTriangles;
    tau.numSurfaceQuads = numSurfaceQuadrilaterals;
    tau.numHexes = numHexes;
    // Please note, wedges are the same as prisms in context of TAU.
    tau.numWedges = numPrisms;

    // create dimension for total number of surface elements
    int ncNumSurfaceElements = 0;
    ret = ret && OK(nc_def_dim(idNc, "no_of_surfaceelements",
        numSurfaceElements, &ncNumSurfaceElements));

    // create variables (and related dimensions) for each type of surface
    // element
    if (ret && numSurfaceTriangles != 0) {
        // create an id for the no_of_surfacetriangles dimension
        int ncSurfaceTriangles = 0;
        ret = ret && OK(nc_def_dim(idNc, "no_of_surfacetriangles",
            numSurfaceTriangles, &ncSurfaceTriangles));

        // create an id for the points_per_surfacetriangle dimension
        int ncPointsPerTri = 0;
        ret = ret && OK(nc_def_dim(idNc, "points_per_surfacetriangle", 3,
            &ncPointsPerTri));

        // create an id for the points_of_surfacetriangles variable, based on
        // the no_of_surfacetriangles and points_per_surfacetriangle
        // dimensions.
        int dimids[2] = { ncSurfaceTriangles, ncPointsPerTri };
        ret = ret && OK(nc_def_var(idNc, "points_of_surfacetriangles", NC_INT,
            2, dimids, &(tau.surfaceTriId)));
    }
    if (ret && numSurfaceQuadrilaterals != 0) {
        // see 'if (numSurfaceTriangles != 0)' comments
        int ncSurfaceQuads = 0;
        ret = ret && OK(nc_def_dim(idNc, "no_of_surfacequadrilaterals",
            numSurfaceQuadrilaterals, &ncSurfaceQuads));

        int ncPointsPerQuad = 0;
        ret = ret && OK(nc_def_dim(idNc, "points_per_surfacequadrilateral", 4,
            &ncPointsPerQuad));

        int dimids[2] = { ncSurfaceQuads, ncPointsPerQuad };
        ret = ret && OK(nc_def_var(idNc, "points_of_surfacequadrilaterals",
            NC_INT, 2, dimids, &(tau.surfaceQuadId)));
    }
    if (ret && numHexes != 0) {
        // see comments in 'if (numSurfaceTriangles != 0)' section
        // (hexaeders are 'hexes' in pointwise)
        int idNumHex = 0;
        ret = ret && OK(nc_def_dim(idNc, "no_of_hexaeders", numHexes,
            &idNumHex));

        int idPtsPerHex = 0;
        ret = ret && OK(nc_def_dim(idNc, "points_per_hexaeder", 8,
            &idPtsPerHex));

        int dimids[2] = { idNumHex, idPtsPerHex };
        ret = ret && OK(nc_def_var(idNc, "points_of_hexaeders", NC_INT, 2,
            dimids, &tau.idHexes));
    }
    if (ret && numPrisms != 0) {
        // see comments in 'if (numSurfaceTriangles != 0)' section
        int idNumPrism = 0;
        ret = ret && OK(nc_def_dim(idNc, "no_of_prisms", numPrisms,
            &idNumPrism));

        int idPtsPerPrism = 0;
        ret = ret && OK(nc_def_dim(idNc, "points_per_prism", 6,
            &idPtsPerPrism));

        int dimids[2] = { idNumPrism, idPtsPerPrism };
        ret = ret && OK(nc_def_var(idNc, "points_of_prisms", NC_INT, 2, dimids,
            &tau.idPrisms));
    }

    // Define variable (netCDF multi-dimensional array) that associates
    // domains (bars -> quads) with surface elements
    int dimids[1] = { ncNumSurfaceElements };

    // Define variable for boundarymarkers
    ret = ret && OK(nc_def_var(idNc, "boundarymarker_of_surfaces", NC_INT, 1,
        dimids, &(tau.numSurfBndryMarkers)));

    // indexing scheme in the boundarypanel variable and in the 
    // boundarymarker variable saves the domain (boundary-panel) or 
    // boundary condition (boundary-marker) number of all the triangles
    // followed by the domain number of all the quads, so this indexing
    // tells us where these are.
    tau.bndryPanelTriStart = 0;
    tau.bndryPanelQuadStart = numSurfaceTriangles;
    
    return ret;
}


/****************************************************************************
 * 
 * This function ends the file definition (where netcdf dimensions and netcdf
 * variables are created), so that data can be written out to the file.
 * 
 ***************************************************************************/
static PWP_BOOL
endFileDefinition(CAEP_RTITEM *pRti)
{
    TAU_DATA &tau = *(pRti->tau);
    return OK(nc_enddef(tau.idNc));
}


/****************************************************************************
 * 
 * Uses the points of a pair of quads and creates a hex from the vertices. The
 * index of the hexes are being tracked by the structure tau. The ordering of
 * the points for hexes is the same in Tau as it is in Pointwise.
 *
 ***************************************************************************/
static PWP_BOOL
createHexFromQuads(int offSetQuad[4], PWGM_ELEMDATA data, CAEP_RTITEM *pRti, 
        int idNc ) {
    PWP_BOOL ret = PWP_TRUE;
    TAU_DATA &tau = *(pRti->tau);
    size_t start[2] = { (size_t)tau.indexHex, 0 };
    size_t count[2] = { 1, 8 };
    int hexPoints[8];

    hexPoints[0] = (int)data.index[0];
    hexPoints[1] = (int)data.index[1];
    hexPoints[2] = (int)data.index[2];
    hexPoints[3] = (int)data.index[3];
    // The points originally saved for the creation of the new quad are reversed
    // of the order expected for the TAU format. We are manually reversing the
    // order for the saves so that the hex has the right dimensions.
    hexPoints[4] = (int)offSetQuad[0];
    hexPoints[5] = (int)offSetQuad[1];
    hexPoints[6] = (int)offSetQuad[2];
    hexPoints[7] = (int)offSetQuad[3];

    ret = ret && OK(nc_put_vara_int(idNc, tau.idHexes,
                            start, count, &(hexPoints[0])));
    return ret;
}


/****************************************************************************
 * 
 * Uses the points of a pair of triangles to create a wedge shape. Because of
 * the way triangles are defined in TAU format, the order of the points for the
 * second triangle must be the same direction as the points are defined to
 * maintain the normal facing away from the face.
 * 
 ***************************************************************************/
static PWP_BOOL
createWedgeFromTris(int offSetTri[3], PWGM_ELEMDATA data, CAEP_RTITEM *pRti, 
        int idNc) {
    PWP_BOOL ret = PWP_TRUE;
    TAU_DATA &tau = *(pRti->tau);
    size_t start[2] = { (size_t)tau.indexWedge, 0 };
    size_t count[2] = { 1, 6 };
    int wedgePoints[6],error;

    wedgePoints[0] = (int)data.index[0];
    wedgePoints[1] = (int)data.index[1];
    wedgePoints[2] = (int)data.index[2];
    wedgePoints[3] = (int)offSetTri[0];
    wedgePoints[4] = (int)offSetTri[1];
    wedgePoints[5] = (int)offSetTri[2];

    error=nc_put_vara_int(idNc, tau.idPrisms, start, count, &(wedgePoints[0]));
    ret = ret && OK(error);

    return ret;
}


/****************************************************************************
 * 
 * Uses the points of the original Quad and creates an offset Quad from the
 * vertices. The function then calls createHexFromQuads() to create the
 * associated hex elements. The original quads have their point order reversed
 * but the new quads maintain their original order because the normals must be
 * facing away from the model.
 * 
 ***************************************************************************/
static PWP_BOOL
createQuadElement(PWP_UINT32 idxQuad, PWGM_ELEMDATA data, 
        int idNc, CAEP_RTITEM *pRti) {
    PWP_BOOL ret = PWP_TRUE;
    size_t start[2] = { idxQuad, 0 };
    size_t count[2] = { 1, 4 };
    TAU_DATA &tau = *(pRti->tau);
    int quadPoints[4];
    PWP_UINT32 markerID = tau.mapNameToOutputID["Symmetry_2"];

    // Creation of a new quad will share the same points as the original quad
    // but will be offset by the number of original points.
    quadPoints[0] = (int)data.index[0] + (int)tau.numPoints;
    quadPoints[1] = (int)data.index[1] + (int)tau.numPoints;
    quadPoints[2] = (int)data.index[2] + (int)tau.numPoints;
    quadPoints[3] = (int)data.index[3] + (int)tau.numPoints;

    ret &= OK(nc_put_vara_int(idNc, tau.surfaceQuadId, 
            start, count, &(quadPoints[0])));

    size_t adjustedIdx = idxQuad + tau.bndryPanelQuadStart;

    ret = ret && OK(nc_put_var1_int(idNc,
            tau.numSurfBndryMarkers, &adjustedIdx,
            (int *)&(markerID)));

    ret &= createHexFromQuads(quadPoints, data, pRti, idNc);
    return ret;
}


/****************************************************************************
 * 
 * Uses the points of the original triangle and creates an offset triangle from
 * the vertices. The function then calls createWedgeFromTris() to create the
 * associated wedge elements. The triangle elements behave the same way as the
 * quad elements. Please see the comments in "createQuadElement".
 * 
 ***************************************************************************/
static PWP_BOOL
createTriElement(PWP_UINT32 idxTri, PWGM_ELEMDATA data, 
        int idNc, CAEP_RTITEM *pRti) {
    size_t start[2] = { idxTri, 0 };
    size_t count[2] = { 1, 3 };
    TAU_DATA &tau = *(pRti->tau);
    int triPoints[3];
    PWP_UINT32 markerID = tau.mapNameToOutputID["Symmetry_2"];
    PWP_BOOL ret = PWP_TRUE;

    // Creation of a new tri will share the same points as the original tri
    // but will be offset by the number of original points.
    triPoints[0] = (int)data.index[0] + (int)tau.numPoints;
    triPoints[1] = (int)data.index[1] + (int)tau.numPoints;
    triPoints[2] = (int)data.index[2] + (int)tau.numPoints;

    size_t adjustedIdx = idxTri + tau.bndryPanelTriStart;

    // save marker (Pointwise boundary-condition) information
    ret = ret && OK(nc_put_var1_int(idNc,
        tau.numSurfBndryMarkers, &adjustedIdx,
        (int *)&(markerID)));

    
    ret &= createWedgeFromTris(triPoints,data,pRti,idNc);
    ret &= OK(nc_put_vara_int(idNc, tau.surfaceTriId,
            start, count, &(triPoints[0])));
    return ret;
}


/****************************************************************************
 * 
 * Uses the points of a bar element to create a quad element seperate from the
 * quads being used to create hexes. This will also expand the boundary
 * conditions to the new surface quad.
 * 
 ***************************************************************************/
PWP_BOOL
writeQuadElementFromBar(PWP_UINT32 idxQuad, int idNc, CAEP_RTITEM *pRti) {
    PWP_BOOL ret = PWP_TRUE;
    TAU_DATA &tau = *(pRti->tau);
    int quadPoints[4];
    PWGM_ELEMDATA data = {PWGM_ELEMTYPE_BAR};
    PWGM_CONDDATA boundaryCondition = {""};
    PWP_UINT32 numDomains = PwModDomainCount(pRti->model);

    for (PWP_UINT32 i = 0; i < numDomains && ret; ++i) {

            PWGM_HDOMAIN domain = PwModEnumDomains(pRti->model, i);
            
            ret = ret && PwDomCondition(domain, &boundaryCondition);
            PWP_UINT32 markerID =
                    tau.mapNameToOutputID[boundaryCondition.name];
            PWP_UINT32 numSurfaceElements = PwDomElementCount(domain, 0);

            for (PWP_UINT32 j = 0; j < numSurfaceElements && ret; ++j) {
                PWGM_HELEMENT surfaceElement = PwDomEnumElements(domain, j);

                ret = ret && PwElemDataMod(surfaceElement, &data);

                // assert to make sure int and data.index[0] are compatible for 
                // calls to nc_put_vara_int (on surface-element information)
                assert(sizeof(int) == sizeof(data.index[0]));

                // assert to make sure int and PWP_UINT are compatible for 
                // calls to nc_put_var1_int (on boundarypanel information)
                assert(sizeof(int) == sizeof(PWP_UINT32));

                // assert to make sure int and tid are compatible for 
                // calls to nc_put_var1_int (on marker array)
                assert(sizeof(int) == sizeof(boundaryCondition.tid));

                switch(data.type) {
                case PWGM_ELEMTYPE_BAR:{
                    size_t start[2] = { idxQuad, 0 };
                    size_t count[2] = { 1, 4 };
                    
                    // The ordering of the points for the quads made form bars
                    // will follow the order shown below with the middle two
                    // being the offset points.
                    quadPoints[3] = data.index[0];
                    quadPoints[2] = data.index[0] + tau.numPoints;
                    quadPoints[1] = data.index[1] + tau.numPoints;
                    quadPoints[0] = data.index[1];
                    
                    

                    ret &= OK(nc_put_vara_int(idNc, tau.surfaceQuadId,
                                        start, count, &(quadPoints[0])));

                    size_t adjustedIdx = idxQuad + tau.bndryPanelQuadStart;

                    ret = ret && OK(nc_put_var1_int(idNc,
                            tau.numSurfBndryMarkers, &adjustedIdx,
                            (int *)&(markerID)));
                    idxQuad++;
                                       }
                default:
                    break;
                }
                ret = ret && !pRti->opAborted;
                caeuProgressIncr(pRti);
            }
            ret = ret && !pRti->opAborted;
        }
    return ret;
}


/****************************************************************************
 * 
 * This function simply gets the points for a particular vertex and stores them
 * in an array to be used in calculating the cross product.
 * 
 ***************************************************************************/
static void
getPoints(PWGM_XYZVAL coordinate[3], PWGM_HVERTEX vertex)
{
    PwVertXyzVal(vertex, PWGM_XYZ_X, &(coordinate[0]));
    PwVertXyzVal(vertex, PWGM_XYZ_Y, &(coordinate[1]));
    PwVertXyzVal(vertex, PWGM_XYZ_Z, &(coordinate[2]));
}


/****************************************************************************
 * 
 * This function calculates the vector based on a starting point and an end
 * point and stores it into another array to calculate the cross product and
 * determine the orientation.
 * 
 ***************************************************************************/
static void
createVectors(PWGM_XYZVAL vector[3], PWGM_XYZVAL startCoord[3], 
        PWGM_XYZVAL endCoord[3])
{
    vector[0] = endCoord[0] - startCoord[0];
    vector[1] = endCoord[1] - startCoord[1];
    vector[2] = endCoord[2] - startCoord[2];
}


/****************************************************************************
 * 
 * Calculate the orientation of each domain based on the ordering of the points
 * of the first element of the block. This is necessary to determine the proper
 * order to output the points for the triangles and quads as well as to decide
 * the direction the z component will be incremented. If any domain is not
 * oriented the same way as the first domain, a warning message is displayed to
 * the user and the orientation is assumed to be the direction of the first
 * domain.
 * 
 ***************************************************************************/
static PWP_BOOL
determineOrientationOfDomain(CAEP_RTITEM *pRti)
{
    PWP_UINT32 i;
    int orientation, message;
    PWGM_XYZVAL coordinate0[3], coordinate1[3], coordinate2[3];
    PWGM_XYZVAL vector1[3], vector2[3];
    PWP_UINT32 numBlocks = PwModBlockCount(pRti->model);
    PWGM_HVERTEX vertex;
    PWGM_HELEMENT element;
    PWGM_ELEMDATA data = {PWGM_ELEMTYPE_SIZE};
    PWGM_HBLOCK block;
    PWP_BOOL ret;

    // Using the first block for the direction of extrusion regardless of the
    // orientation of the other blocks.
    block = PwModEnumBlocks(pRti->model, 0);
    element = PwBlkEnumElements(block, 0);
    PwElemDataMod(element, &data);

    vertex = PwModEnumVertices(pRti->model, data.index[0]);
    getPoints(coordinate0, vertex);
    vertex = PwModEnumVertices(pRti->model, data.index[1]);
    getPoints(coordinate1, vertex);

    // if the element is a quad, the vertex to be used is the vertex right next
    // to the 0th point i.e. point 3.
    if (data.vertCnt == 4) {
        vertex = PwModEnumVertices(pRti->model, data.index[3]);
        getPoints(coordinate2, vertex);
    } 
    // if the element is a tri, the vertex to be used is the only other vertex
    // remaining.
    else if (data.vertCnt == 3) {
        vertex = PwModEnumVertices(pRti->model, data.index[2]);
        getPoints(coordinate2, vertex);
    }
    
    createVectors(vector1,coordinate0,coordinate1);
    createVectors(vector2,coordinate0,coordinate2);

    // The block's orientation is determined only by the value of the z comp.
    // This is done using the formula:
    //   <cx, cy, cz> = < ay*bz - az*by, az*bx - ax*bz, ax*by - ay*bx >
    // because the z component is the only required calculation, the third part
    // of the formula is the only part being calculated.
    orientation = (int)vector1[0] * (int)vector2[1] -
        (int)vector1[1] * (int)vector2[0];
    if (orientation > 0) {
        message = 1;
        ret = PWP_TRUE;
    }
    else {
        message = 0;
        ret = PWP_FALSE;
    }
    // This loop performs the same function as above but for each block. This is
    // done to check if the user has domains oriented in a different direction
    // as the first block. This is a problem because the file format requires
    // that the symmetry points be input in the exact same order as the original
    // set and different directions would create 2 pairs of points with only one
    // associated set of original points.
    for ( i=1 ; i < numBlocks ; i++ ) {
        block = PwModEnumBlocks(pRti->model, i);
        element = PwBlkEnumElements(block, 0);
        PwElemDataMod(element, &data);

        vertex = PwModEnumVertices(pRti->model, data.index[0]);
        getPoints(coordinate0, vertex);
        vertex = PwModEnumVertices(pRti->model, data.index[1]);
        getPoints(coordinate1, vertex);

        if (data.vertCnt == 4) {
            vertex = PwModEnumVertices(pRti->model, data.index[3]);
            getPoints(coordinate2, vertex);
        }
        else if (data.vertCnt == 3) {
            vertex = PwModEnumVertices(pRti->model, data.index[2]);
            getPoints(coordinate2, vertex);
        }
        createVectors(vector1,coordinate0,coordinate1);
        createVectors(vector2,coordinate0,coordinate2);
        orientation = (int)vector1[0] * (int)vector2[1] -
            (int)vector1[1] * (int)vector2[0];
        if (orientation > 0) {
            if (message == 0) {
                caeuSendWarningMsg(pRti, "Domains do not have uniform "
                    "orientation. Assuming positive orientation.", 0);
                break;
            }
        }
        else {
            if (message == 1) {
                caeuSendWarningMsg(pRti, "Domains do not have uniform "
                    "orientation. Assuming negative orientation.", 0);
                break;
            }
        }
    }
    return ret;
}


/****************************************************************************
 * 
 * Iterate through the grid's vertices and save their xyz coordinates to their
 * respective TAU variables.
 * 
 ***************************************************************************/
static PWP_BOOL
writePointVars(CAEP_RTITEM *pRti)
{
    TAU_DATA &tau = *(pRti->tau);
    PWP_BOOL ret = PWP_TRUE;
    int idNc = tau.idNc;

    int xcid = tau.idXc;
    int ycid = tau.idYc;
    int zcid = tau.idZc;

    // assert to make sure the value retrieved in PwVertXyzVal and that
    // saved by nc_put_var1_double are compatible.
    assert(sizeof(PWGM_XYZVAL) == sizeof(double));

    PWP_UINT32 numPoints = PwModVertexCount(pRti->model);
    if (caeuProgressBeginStep(pRti, 3 * numPoints)) {
        // originally the xs, ys, and zs were saved in the same loop, however
        // that created a speed bottleneck; putting xs ys and zs in separate
        // loops at least triples performance

        // val is a buffer to save points in so that they can be written out
        // faster.
        PWGM_XYZVAL val[64];
        for (PWP_UINT32 i = 0; i < numPoints && ret; ++i) {
            // retrieve vertex
            PWGM_HVERTEX vertex = PwModEnumVertices(pRti->model, i);

            // retrieve x val from vertex
            ret = ret && PwVertXyzVal(vertex, PWGM_XYZ_X, &(val[i % 64]));

            // if at end of buffer or at last point, clear buffer to memory
            if (i % 64 == 63 || i == numPoints - 1) {
                // save vertex to associated variable
                size_t index[1] = {i / 64 * 64};
                size_t count[1] = {i % 64 + 1};
                ret = ret && OK(nc_put_vara_double(idNc, xcid, index, count,
                    &(val[0])));
            }
            ret = ret && !pRti->opAborted;
            caeuProgressIncr(pRti);
        }
        for (PWP_UINT32 i = 0; i < numPoints && ret; ++i) {
            PWGM_HVERTEX vertex = PwModEnumVertices(pRti->model, i);

            ret = ret && PwVertXyzVal(vertex, PWGM_XYZ_Y, &(val[i % 64]));

            // if at end of buffer or at last point, clear buffer to memory
            if (i % 64 == 63 || i == numPoints - 1) {
                size_t index[1] = {i / 64 * 64};
                size_t count[1] = {i % 64 + 1};
                ret = ret && OK(nc_put_vara_double(idNc, ycid, index, count,
                    &(val[0])));
            }
            ret = ret && !pRti->opAborted;
            caeuProgressIncr(pRti);
        }
        for (PWP_UINT32 i = 0; i < numPoints && ret; ++i) {
            PWGM_HVERTEX vertex = PwModEnumVertices(pRti->model, i);

            ret = ret && PwVertXyzVal(vertex, PWGM_XYZ_Z, &(val[i % 64]));

            // if at end of buffer or at last point, clear buffer to memory
            if (i % 64 == 63 || i == numPoints - 1) {
                size_t index[1] = {i / 64 * 64};
                size_t count[1] = {i % 64 + 1};
                ret = ret && OK(nc_put_vara_double(idNc, zcid, index, count,
                    &(val[0])));
            }
            ret = ret && !pRti->opAborted;
            caeuProgressIncr(pRti);
        }
        caeuProgressEndStep(pRti);
    }
    return ret;
}


/****************************************************************************
 * 
 * Iterate through the vertices in the Pointwise model and save their xyz
 * coordinates to their respective TAU variables. For 2D, this process is
 * repeated, offsetting the z component by the amount defined by the user.
 * 
 ***************************************************************************/
static PWP_BOOL
writePointVars2D(CAEP_RTITEM *pRti)
{
    TAU_DATA &tau = *(pRti->tau);
    PWP_BOOL ret = PWP_TRUE;
    PWP_BOOL direction = determineOrientationOfDomain(pRti);
    PWP_UINT32 i;
    int idNc = tau.idNc;
    PWP_REAL thickness;
    PwModGetAttributeREAL(pRti->model, Thickness, &thickness);
    int xcid = tau.idXc;
    int ycid = tau.idYc;
    int zcid = tau.idZc;

    // assert to make sure the value retrieved in PwVertXyzVal and that
    // saved by nc_put_var1_double are compatible.
    assert(sizeof(PWGM_XYZVAL) == sizeof(double));

    PWP_UINT32 numPoints = PwModVertexCount(pRti->model);
    if (caeuProgressBeginStep(pRti, 3 * numPoints)) {
        // originally the xs, ys, and zs were saved in the same loop, however
        // that created a speed bottleneck; putting xs ys and zs in separate
        // loops at least triples performance

        // val is a buffer to save points in so that they can be written out
        // faster.
        PWGM_XYZVAL val[64];
        PWGM_HVERTEX vertex;
        for (i = 0; i < numPoints && ret; ++i) {
            // retrieve vertex
            vertex = PwModEnumVertices(pRti->model, i);

            // retrieve x val from vertex
            ret = ret && PwVertXyzVal(vertex, PWGM_XYZ_X, &(val[i % 64]));

            // if at end of buffer or at last point, clear buffer to memory
            if (i % 64 == 63 || i == numPoints - 1) {
                // save vertex to associated variable
                size_t index[1] = {i / 64 * 64};
                size_t count[1] = {i % 64 + 1};
                ret = ret && OK(nc_put_vara_double(idNc, xcid, index, count,
                    &(val[0])));
            }
            ret = ret && !pRti->opAborted;
            caeuProgressIncr(pRti);
        }
        for (i = 0; i < numPoints && ret; ++i) {
            vertex = PwModEnumVertices(pRti->model, i);

            ret = ret && PwVertXyzVal(vertex, PWGM_XYZ_Y, &(val[i % 64]));

            // if at end of buffer or at last point, clear buffer to memory
            if (i % 64 == 63 || i == numPoints - 1) {
                size_t index[1] = {i / 64 * 64};
                size_t count[1] = {i % 64 + 1};
                ret = ret && OK(nc_put_vara_double(idNc, ycid, index, count,
                    &(val[0])));
            }
            ret = ret && !pRti->opAborted;
            caeuProgressIncr(pRti);
        }
        for (i = 0; i < numPoints && ret; ++i) {
            vertex = PwModEnumVertices(pRti->model, i);

            ret = ret && PwVertXyzVal(vertex, PWGM_XYZ_Z, &(val[i % 64]));

            // if at end of buffer or at last point, clear buffer to memory
            if (i % 64 == 63 || i == numPoints - 1) {
                size_t index[1] = {i / 64 * 64};
                size_t count[1] = {i % 64 + 1};
                ret = ret && OK(nc_put_vara_double(idNc, zcid, index, count,
                    &(val[0])));
            }
            ret = ret && !pRti->opAborted;
            caeuProgressIncr(pRti);
        }
        // This is for the second set of points that are defined as part of
        // the TAU 2D export. They are duplicates of the original points but
        // will have a z offset.
        for (i = 0; i < numPoints && ret; ++i) {
            // retrieve vertex
            vertex = PwModEnumVertices(pRti->model, i);

            // retrieve x val from vertex
            ret = ret && PwVertXyzVal(vertex, PWGM_XYZ_X, &(val[i % 64]));

            // if at end of buffer or at last point, clear buffer to memory
            if (i % 64 == 63 || i == numPoints - 1) {
                // save vertex to associated variable
                //Adding the numPoints will offset the new point to be the same
                //
                size_t index[1] = {(i / 64 * 64) + numPoints};
                size_t count[1] = {i % 64 + 1};
                ret = ret && OK(nc_put_vara_double(idNc, xcid, index, count,
                    &(val[0])));
            }
            ret = ret && !pRti->opAborted;
            caeuProgressIncr(pRti);
        }
        for (i = 0; i < numPoints && ret; ++i) {
            vertex = PwModEnumVertices(pRti->model, i);

            ret = ret && PwVertXyzVal(vertex, PWGM_XYZ_Y, &(val[i % 64]));

            // if at end of buffer or at last point, clear buffer to memory
            if (i % 64 == 63 || i == numPoints - 1) {
                size_t index[1] = {(i / 64 * 64) + numPoints};
                size_t count[1] = {i % 64 + 1};
                ret = ret && OK(nc_put_vara_double(idNc, ycid, index, count,
                    &(val[0])));
            }
            ret = ret && !pRti->opAborted;
            caeuProgressIncr(pRti);
        }
        for (i = 0; i < numPoints && ret; ++i) {
            vertex = PwModEnumVertices(pRti->model, i);

            ret = ret && PwVertXyzVal(vertex, PWGM_XYZ_Z, &(val[i % 64]));
            
            // The value of 1 is the assumed step size of the extrusion for 2D.
            if (direction) {
                val[i % 64] += thickness;
            }
            // Or, the value will have a negative offset if the user has their
            // domains orientated in the negative z direction.
            else {
                val[i % 64] -= thickness;
            }

            // if at end of buffer or at last point, clear buffer to memory
            if (i % 64 == 63 || i == numPoints - 1) {
                size_t index[1] = {(i / 64 * 64) + numPoints};
                size_t count[1] = {i % 64 + 1};
                ret = ret && OK(nc_put_vara_double(idNc, zcid, index, count,
                    &(val[0])));
            }
            ret = ret && !pRti->opAborted;
            caeuProgressIncr(pRti);
        }
        caeuProgressEndStep(pRti);
    }
    return ret;
}

// These four numbers give the buffer sizes for each shape type. They are used
// only by the writeElementVars() function. Shapes are buffered together,
// because they get saved to different parts of the file, and caching them
// reduces runtime (runtimeWrite() function takes 20% less time with buffering)
#define BUFFER_SIZE_TET 64
#define BUFFER_SIZE_PYR 64
#define BUFFER_SIZE_WEDGE 64
#define BUFFER_SIZE_HEX 64

/****************************************************************************
 * 
 * Iterate through the blocks in the Pointwise model and for each block, iterate
 * through the block-elements and save them to their respective TAU variables,
 * based on their element-type.
 * 
 ***************************************************************************/
static PWP_BOOL
writeElementVars(CAEP_RTITEM *pRti)
{
    TAU_DATA &tau = *(pRti->tau);
    PWP_BOOL ret = PWP_TRUE;
    int idNc = tau.idNc;

    // because the different element types are saved in different arrays, we
    // must track and increment separate indexes for each element type in order
    // to know where the next element needs to be saved.
    
    PWP_UINT32 idxTet = 0;
    PWP_UINT32 idxPyr = 0;
    PWP_UINT32 idxWedge = 0;
    PWP_UINT32 idxHex = 0;

    // caches; See comment above '#define BUFFER_SIZE_TET 64'
    int bufferedTets[BUFFER_SIZE_TET][4];
    int bufferedPyramids[BUFFER_SIZE_PYR][5];
    int bufferedWedges[BUFFER_SIZE_WEDGE][6];
    int bufferedHexes[BUFFER_SIZE_HEX][8];

    PWP_UINT32 numTets = (PWP_UINT32)tau.numTets;
    PWP_UINT32 numPyrs = (PWP_UINT32)tau.numPyramids;
    PWP_UINT32 numWedges = (PWP_UINT32)tau.numWedges;
    PWP_UINT32 numHexes = (PWP_UINT32)tau.numHexes;

    // now we iterate through the blocks, and through each block's block-
    // elements, and output the results to the netcdf file.

    PWP_UINT32 numBlocks = PwModBlockCount(pRti->model);
    if (caeuProgressBeginStep(pRti, tau.numElements)) {

        // Init data to invalid value
        PWGM_ELEMDATA data = {PWGM_ELEMTYPE_SIZE};

        for (PWP_UINT32 i = 0; i < numBlocks && ret; ++i) {
            PWGM_HBLOCK block = PwModEnumBlocks(pRti->model, i);
            PWP_UINT32 numElements = PwBlkElementCount(block, 0);
            for (PWP_UINT32 j = 0; j < numElements && ret; ++j) {
                PWGM_HELEMENT element = PwBlkEnumElements(block, j);

                ret = PwElemDataMod(element, &data);

                // assert to make sure int and data.index[0] are compatible for 
                // calls to nc_put_vara_int
                assert(sizeof(int) == sizeof(data.index[0]));

                // now we enter a switch on the type of the element to
                // determine where and how to save it.

                // note that (thankfully) the indexing scheme for EVERY element
                // type is the same for TAU and CAE, so they don't need to be
                // re-indexed.

                // all cases in the switch are handled in the same manner

                switch(data.type) {
                case PWGM_ELEMTYPE_TET:
                    // copy data to cache
                    for (int k = 0; k < 4; ++k) {
                        bufferedTets[idxTet % BUFFER_SIZE_TET][k] =
                                (int)data.index[k];
                    }
                    idxTet++;
                    // If cache is full, OR we reach the last of this kind of
                    // element, then clear the cache
                    if (0u == idxTet % BUFFER_SIZE_TET || idxTet == numTets) {
                        size_t start[2] = {(idxTet - 1) /BUFFER_SIZE_TET *
                            BUFFER_SIZE_TET, 0};
                        size_t count[2] = {BUFFER_SIZE_TET, 4};
                        // Treat the last buffered group differently (because
                        // it may not use the whole buffer). Note that when
                        // (numTets % BUFFER_SIZE_TET == 0) the if block is not
                        // entered, because that means the whole buffer is
                        // being used (and we don't want to set count[0] to 0).
                        if (idxTet == numTets && numTets %
                                BUFFER_SIZE_TET != 0) {
                            count[0] = numTets % BUFFER_SIZE_TET;
                        }
                        // perform the write and handle errors
                        ret = ret && OK(nc_put_vara_int(idNc, tau.idTets, start,
                            count, &(bufferedTets[0][0])));
                    }
                    break;
                case PWGM_ELEMTYPE_PYRAMID:
                    // see comments in PWGM_ELEMTYPE_TET case
                    for (int k = 0; k < 5; ++k) {
                        bufferedPyramids[idxPyr % BUFFER_SIZE_PYR][k] =
                            (int)data.index[k];
                    }
                    idxPyr++;
                    if (0u == idxPyr % BUFFER_SIZE_PYR || idxPyr == numPyrs) {
                        size_t start[2] = {(idxPyr - 1) / BUFFER_SIZE_PYR *
                            BUFFER_SIZE_PYR, 0};
                        size_t count[2] = {BUFFER_SIZE_PYR, 5};
                        if (idxPyr == numPyrs &&
                                numPyrs % BUFFER_SIZE_PYR != 0) {
                            count[0] = numPyrs % BUFFER_SIZE_PYR;
                        }
                        ret = ret && OK(nc_put_vara_int(idNc, tau.idPyramids,
                            start, count, &(bufferedPyramids[0][0])));
                    }
                    break;
                case PWGM_ELEMTYPE_WEDGE:
                    // see comments in PWGM_ELEMTYPE_TET case
                    for (int k = 0; k < 6; ++k) {
                        bufferedWedges[idxWedge % BUFFER_SIZE_WEDGE][k] =
                                (int)data.index[k];
                    }
                    idxWedge++;
                    if (0u == idxWedge % BUFFER_SIZE_WEDGE || idxWedge ==
                            numWedges) {
                        size_t start[2] = {(idxWedge - 1) / BUFFER_SIZE_WEDGE *
                            BUFFER_SIZE_WEDGE, 0};
                        size_t count[2] = {BUFFER_SIZE_WEDGE, 6};
                        if (idxWedge == numWedges && numWedges %
                                BUFFER_SIZE_WEDGE != 0) {
                            count[0] = numWedges % BUFFER_SIZE_WEDGE;
                        }
                        ret = ret && OK(nc_put_vara_int(idNc, tau.idPrisms,
                            start, count, &(bufferedWedges[0][0])));
                    }
                    break;
                case PWGM_ELEMTYPE_HEX:
                    // see comments in PWGM_ELEMTYPE_TET case
                    for (int k = 0; k < 8; ++k) {
                        bufferedHexes[idxHex % BUFFER_SIZE_HEX][k] =
                                (int)data.index[k];
                    }
                    idxHex++;
                    if (0u == idxHex % BUFFER_SIZE_HEX || idxHex == numHexes) {
                        size_t start[2] = {(idxHex - 1) / BUFFER_SIZE_HEX *
                            BUFFER_SIZE_HEX, 0};
                        size_t count[2] = {BUFFER_SIZE_HEX, 8};
                        if (idxHex == numHexes && numHexes %
                                BUFFER_SIZE_HEX != 0) {
                            count[0] = numHexes % BUFFER_SIZE_HEX;
                        }
                        ret = ret && OK(nc_put_vara_int(idNc, tau.idHexes,
                            start, count, &(bufferedHexes[0][0])));
                    }
                    break;
                default:
                    // this should not occur; somehow an unsupported type was
                    // given.
                    ret = false;
                    break;
                }
                ret = ret && !pRti->opAborted;
                caeuProgressIncr(pRti);
            }
        }
        caeuProgressEndStep(pRti);
    }

    return ret;
}


/****************************************************************************
 * 
 * For each block, iterate through the block-elements and save them to their
 * respective TAU variables, based on their element-type. After saving the
 * element, the program then creates the offset variable and the associated 3D
 * element.
 * 
 ***************************************************************************/
static PWP_BOOL
writeElementVars2D(CAEP_RTITEM *pRti)
{
    TAU_DATA &tau = *(pRti->tau);
    PWP_BOOL ret = PWP_TRUE;
    int idNc = tau.idNc;

    // as with block-elements, we must separately track our indexing into the
    // tri and quad netcdf-variables.

    PWP_UINT32 idxTri = 0;
    PWP_UINT32 idxQuad = 0;
    tau.indexHex = 0;
    tau.indexWedge = 0;

    // now we iterate through the domains, and through each domain's domain-
    // elements, and output the results to the netcdf file.

    PWP_UINT32 numBlocks = PwModBlockCount(pRti->model);
    PWGM_ELEMDATA data = {PWGM_ELEMTYPE_BAR};

    if (caeuProgressBeginStep(pRti, tau.numSurfaceElements-tau.numBars)) {

        // This initializes the variable with a dummy variable which will
        // later be changed. In 3D a domain can not have a bar element.
        
        for (PWP_UINT32 i = 0; i < numBlocks && ret; ++i) {

            PWGM_HBLOCK block = PwModEnumBlocks(pRti->model, i);
            
            PWP_UINT32 markerID =
                    tau.mapNameToOutputID["Symmetry_1"];
            PWP_UINT32 numSurfaceElements = PwBlkElementCount(block, 0);

            for (PWP_UINT32 j = 0; j < numSurfaceElements && ret; ++j) {
                PWGM_HELEMENT surfaceElement = PwBlkEnumElements(block, j);

                ret = ret && PwElemDataMod(surfaceElement, &data);

                // assert to make sure int and data.index[0] are compatible for 
                // calls to nc_put_vara_int (on surface-element information)
                assert(sizeof(int) == sizeof(data.index[0]));

                // assert to make sure int and PWP_UINT are compatible for 
                // calls to nc_put_var1_int (on boundarypanel information)
                assert(sizeof(int) == sizeof(PWP_UINT32));

                switch(data.type) {
                
                case PWGM_ELEMTYPE_TRI: {
                    // save vertex information for current triangle
                    size_t start[2] = { idxTri, 0 };
                    size_t count[2] = { 1, 3 };
                    // unfortunately the point order of surface triangles in
                    // TAU is reversed of the order in Pointwise. The reversal
                    // also must occur for surface quadrilaterals.
                    int triPoints[3] = { (int)data.index[2], (int)data.index[1],
                        (int)data.index[0] };
                    ret = ret && OK(nc_put_vara_int(idNc, tau.surfaceTriId,
                        start, count, &(triPoints[0])));

                    size_t adjustedIdx = idxTri + tau.bndryPanelTriStart;

                    // save marker (Pointwise boundary-condition) information
                    ret = ret && OK(nc_put_var1_int(idNc,
                        tau.numSurfBndryMarkers, &adjustedIdx,
                        (int *)&(markerID)));

                    // increment to next triangle
                    idxTri++;
                    // Create the new triangle element that is extruded from
                    // this element.
                    ret = ret && createTriElement(idxTri,data,idNc,pRti);
                    // increment to next triangle
                    idxTri++;
                    // increment to next wedge
                    tau.indexWedge++;
                    break; }
                case PWGM_ELEMTYPE_QUAD: {
                    // see comments for PWGM_ELEMTYPE_TRI:
                    size_t start[2] = { idxQuad, 0 };
                    size_t count[2] = { 1, 4 };
                    int quadPoints[4] = { (int)data.index[3],
                        (int)data.index[2], (int)data.index[1],
                        (int)data.index[0] };
                    ret = ret && OK(nc_put_vara_int(idNc, tau.surfaceQuadId,
                        start, count, &(quadPoints[0])));

                    size_t adjustedIdx = idxQuad + tau.bndryPanelQuadStart;

                    ret = ret && OK(nc_put_var1_int(idNc,
                        tau.numSurfBndryMarkers, &adjustedIdx,
                        (int *)&(markerID)));
                    
                    // increment to next quad
                    idxQuad++;
                    // Create the new quad element that is offset from this
                    // element.
                    ret = ret && createQuadElement(idxQuad,data,idNc,pRti);
                    // increment to next quad
                    idxQuad++;
                    // increment to next hex
                    tau.indexHex++;
                    break; }
                default:
                    break;
                }
                ret = ret && !pRti->opAborted;
                caeuProgressIncr(pRti);
                caeuProgressIncr(pRti);
            }
            ret = ret && !pRti->opAborted;
        }
        caeuProgressEndStep(pRti);
    }

    // After all of the quads and tris are written, the bars need to be extruded
    // and written to the quad variable. This requires continuing to track the
    // index of the quad.
    ret &= writeQuadElementFromBar(idxQuad, idNc, pRti);

    return ret;
}


/****************************************************************************
 * 
 * Writes out the marker numbers. Due to the lack of clarity in the specs, I
 * have assured that each value is the same as its index.
 * 
 ***************************************************************************/
static PWP_BOOL
writeMarkers(CAEP_RTITEM *pRti)
{
    TAU_DATA &tau = *(pRti->tau);
    PWP_BOOL ret = PWP_TRUE;
    int idNc = tau.idNc;
    for (int ii = 0; ii < (int)tau.mapNameToOutputID.size() && ret; ++ii) {
        size_t index[1] = {(size_t)ii};
        ret = ret && OK(nc_put_var1_int(idNc, tau.markerId, index, &ii));
    }
    return ret;
}


/****************************************************************************
 * 
 * For each domain, iterate through the domain-elements and save them to their
 * respective TAU variables, based on their element-type.
 * 
 ***************************************************************************/
static PWP_BOOL
writeSurfaceElementVars(CAEP_RTITEM *pRti)
{
    TAU_DATA &tau = *(pRti->tau);
    PWP_BOOL ret = PWP_TRUE;
    int idNc = tau.idNc;

    // as with block-elements, we must separately track our indexing into the
    // tri and quad netcdf-variables.

    PWP_UINT32 idxTri = 0;
    PWP_UINT32 idxQuad = 0;

    // now we iterate through the domains, and through each domain's domain-
    // elements, and output the results to the netcdf file.

    PWP_UINT32 numDomains = PwModDomainCount(pRti->model);
    if (caeuProgressBeginStep(pRti, tau.numSurfaceElements)) {

        // This initializes the variable with a dummy variable which will
        // later be changed. In 3D a domain can not have a bar element.
        PWGM_ELEMDATA data = {PWGM_ELEMTYPE_BAR};
        for (PWP_UINT32 i = 0; i < numDomains && ret; ++i) {
            PWGM_HDOMAIN domain = PwModEnumDomains(pRti->model, i);
            PWGM_CONDDATA boundaryCondition = {""};
            ret = ret && PwDomCondition(domain, &boundaryCondition);
            PWP_UINT32 markerID =
                    tau.mapNameToOutputID[boundaryCondition.name];
            PWP_UINT32 numSurfaceElements = PwDomElementCount(domain, 0);
            for (PWP_UINT32 j = 0; j < numSurfaceElements && ret; ++j) {
                PWGM_HELEMENT surfaceElement = PwDomEnumElements(domain, j);

                ret = ret && PwElemDataMod(surfaceElement, &data);

                // assert to make sure int and data.index[0] are compatible for 
                // calls to nc_put_vara_int (on surface-element information)
                assert(sizeof(int) == sizeof(data.index[0]));

                // assert to make sure int and PWP_UINT are compatible for 
                // calls to nc_put_var1_int (on boundarypanel information)
                assert(sizeof(int) == sizeof(PWP_UINT32));

                // assert to make sure int and tid are compatible for 
                // calls to nc_put_var1_int (on marker array)
                assert(sizeof(int) == sizeof(boundaryCondition.tid));

                switch(data.type) {
                case PWGM_ELEMTYPE_TRI: {
                    // save vertex information for current triangle
                    size_t start[2] = { idxTri, 0 };
                    size_t count[2] = { 1, 3 };
                    // unfortunately the point order of surface triangles in
                    // TAU is reversed of the order in Pointwise. The reversal
                    // also must occur for surface quadrilaterals.
                    int triPoints[3] = { (int)data.index[2], (int)data.index[1],
                        (int)data.index[0] };
                    ret = ret && OK(nc_put_vara_int(idNc, tau.surfaceTriId,
                        start, count, &(triPoints[0])));

                    // compute the index into the boundarypanel variable and
                    // the marker variable, based on the offset for 
                    // triangles and the index of the current triangle
                    size_t adjustedIdx = idxTri + tau.bndryPanelTriStart;

                    // save marker (Pointwise boundary-condition) information
                    ret = ret && OK(nc_put_var1_int(idNc,
                        tau.numSurfBndryMarkers, &adjustedIdx,
                        (int *)&(markerID)));

                    // increment to next triangle
                    idxTri++;
                    break; }
                case PWGM_ELEMTYPE_QUAD: {
                    // see comments for PWGM_ELEMTYPE_TRI:
                    size_t start[2] = { idxQuad, 0 };
                    size_t count[2] = { 1, 4 };
                    int quadPoints[4] = { (int)data.index[3],
                        (int)data.index[2], (int)data.index[1],
                        (int)data.index[0] };
                    ret = ret && OK(nc_put_vara_int(idNc, tau.surfaceQuadId,
                        start, count, &(quadPoints[0])));

                    size_t adjustedIdx = idxQuad + tau.bndryPanelQuadStart;

                    ret = ret && OK(nc_put_var1_int(idNc,
                        tau.numSurfBndryMarkers, &adjustedIdx,
                        (int *)&(markerID)));

                    idxQuad++;
                    break; }
                default:// probably a bar
                    break;
                }
                ret = ret && !pRti->opAborted;
                caeuProgressIncr(pRti);
            }
            ret = ret && !pRti->opAborted;
        }
        caeuProgressEndStep(pRti);
    }

    return ret;
}


/****************************************************************************
 * 
 * Close the netcdf file
 * 
 ***************************************************************************/
static PWP_BOOL
closeGridFile(CAEP_RTITEM *pRti)
{
    TAU_DATA &tau = *(pRti->tau);
    int idNc = tau.idNc;
    return OK(nc_close(idNc));
}


/****************************************************************************
 * 
 * Writes out the TAU boundary mapping file.
 * 
 ***************************************************************************/
static PWP_BOOL
writeBmapFile(CAEP_RTITEM *pRti)
{
    TAU_DATA &tau = *(pRti->tau);
    PWP_BOOL ret = PWP_TRUE;

    char strBmap[1024];
    sprintf(strBmap, "%s.%s", pRti->pWriteInfo->fileDest, "bmap");
    FILE *fp = pwpFileOpen(strBmap, pwpWrite | pwpAscii);

    // ensure file pointer is valid
    if (!fp) {
        return PWP_FALSE;
    }

    ret = ret && 0 < fprintf(fp,
            " -----------------------------------------------------\n");
    ret = ret && 0 < fprintf(fp, " BOUNDARY MAPPING\n");
    ret = ret && 0 < fprintf(fp,
            " -----------------------------------------------------\n");

    StringUINT32Map &mNameToId = tau.mapNameToOutputID;
    StringUINT32Map::iterator iterIds;
    StringStringMap &mNameToType = tau.mapNameToType;
    for (iterIds = mNameToId.begin(); iterIds != mNameToId.end() && ret;
            ++iterIds) {
        const char *name = iterIds->first.c_str();
        PWP_UINT32 outputID = iterIds->second;
        const char *type = mNameToType[iterIds->first].c_str();
        ret = ret && 0 < fprintf(fp, "\n");
        ret = ret && 0 < fprintf(fp, "%26s %u\n", "Markers:", outputID);
        ret = ret && 0 < fprintf(fp, "%26s %s\n", "Type:", type);
        ret = ret && 0 < fprintf(fp, "%26s %s\n", "Name:", name);
        // future attributes added here
        ret = ret && 0 < fprintf(fp, "\n");
        ret = ret && 0 < fprintf(fp, " block end\n");
        ret = ret && 0 < fprintf(fp, " ---------------------------\n");
    }

    pwpFileClose(fp);

    return ret;
}


/**************************************/
PWP_BOOL
runtimeWrite(CAEP_RTITEM *pRti, PWGM_HGRIDMODEL model,
             const CAEP_WRITEINFO *pWriteInfo)
{
    PWP_BOOL ret = PWP_TRUE;
    TAU_DATA data = {0};
    pRti->tau = &data; // cppcheck-suppress autoVariables
    // Checks if the user is using 2D or 3D
    PWP_BOOL dim = CAEPU_RT_DIM_2D(pRti);
    PWP_UINT32 cnt;

    if (pRti && model && pWriteInfo) {
        if(dim)
        {
            cnt = 2; /* the # of MAJOR progress steps 2D */
        }
        else
        {
            cnt = 3; /* the # of MAJOR progress steps 3D */
        }
        

        if (caeuProgressInit(pRti, cnt)) {
            
            // These first four functions are the same for 2D and 3D.
            if (!setupMarkerInfo(pRti)) {
                caeuSendErrorMsg(pRti, "Could not setup marker info!", 0);
                ret = false;
            }
            else if (!startGridFile(pRti)) {
                caeuSendErrorMsg(pRti, "Could start grid file info!", 0);
                ret = false;
            }
            else if (!definePointVars(pRti)) {
                caeuSendErrorMsg(pRti, "Could not define point vars!", 0);
                ret = false;
            }
            else if (!defineMarkers(pRti)) {
                caeuSendErrorMsg(pRti, "Could not define markers!", 0);
                ret = false;
            }
            // The 2D and 3D exports now differ.
            if(dim)
            {
                if(!defineElementVars2D(pRti)) {
                    caeuSendErrorMsg(pRti, "Could not define element vars!", 0);
                    ret = false;
                }
                else if (!endFileDefinition(pRti)) {
                    caeuSendErrorMsg(pRti, "Could not create TAU file! If "
                        "writing a large file, ensure that the solver "
                        "attribute 'dataModel' is not set to classic.",
                        0);
                    ret = false;
                }
                // first major progress step
                else if(!writePointVars2D(pRti)) {
                    caeuSendErrorMsg(pRti, "Could not write point vars!", 0);
                    ret = false;
                }
                else if(!writeMarkers(pRti)) {
                    caeuSendErrorMsg(pRti, "Could not write markers!", 0);
                    ret = false;
                }
                // second major progress step
                else if(!writeElementVars2D(pRti)) {
                    caeuSendErrorMsg(pRti, "Could not write element vars!", 0);
                    ret = false;
                }
            }
            else
            {
                if (!defineElementVars(pRti)) {
                    caeuSendErrorMsg(pRti, "Could not define element vars!", 0);
                    ret = false;
                }
                else if (!defineSurfaceElementVars(pRti)) {
                    caeuSendErrorMsg(pRti, "Could not define surface element "
                        "vars!", 0);
                    ret = false;
                }
                else if (!endFileDefinition(pRti)) {
                    caeuSendErrorMsg(pRti, "Could not create TAU file! If "
                        "writing a large file, ensure that the solver "
                        "attribute 'dataModel' is not set to classic.",
                        0);
                    ret = false;
                }
                // first major progress step
                else if (!writePointVars(pRti)) {
                    caeuSendErrorMsg(pRti, "Could not write point vars!", 0);
                    ret = false;
                }
                // second major progress step
                else if (!writeElementVars(pRti)) {
                    caeuSendErrorMsg(pRti, "Could not write element vars!", 0);
                    ret = false;
                }
                else if (!writeMarkers(pRti)) {
                    caeuSendErrorMsg(pRti, "Could not write markers!", 0);
                    ret = false;
                }
                // third major step
                else if (!writeSurfaceElementVars(pRti)) {
                    caeuSendErrorMsg(pRti, "Could not write surface element "
                        "vars!", 0);
                    ret = false;
                }
            }
            // Attempt to close file even if previous step failed.
            if (!closeGridFile(pRti)) {
                ret = false;
            }

            ret = ret && writeBmapFile(pRti);

            caeuProgressEnd(pRti, ret);
            ret = ret && !pRti->opAborted;
        }
    }
    return ret;
}


PWP_BOOL
runtimeCreate(CAEP_RTITEM * /*pRti*/)
{
    PWP_BOOL ret = PWP_TRUE;
    const char *desc = "Controls the underlying NetCDF data model used by TAU "
        "export.";
    const char *NetCDFFormatEnum = "NetCDF_Classic|NetCDF_64bit|NetCDF4_HDF5";

    ret = caeuPublishValueDefinition(NetCDFDataModel, PWP_VALTYPE_ENUM,
            "NetCDF_64bit", "RW", desc, NetCDFFormatEnum) &&

        caeuPublishValueDefinition(Thickness, PWP_VALTYPE_REAL, "1.0",
            "RW", "Offset distance for 2D export", "0 +Inf");

    return ret;
}


PWP_VOID
runtimeDestroy(CAEP_RTITEM * /*pRti*/)
{
}

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
