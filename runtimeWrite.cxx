/****************************************************************************
 *
 * (C) 2021 Cadence Design Systems, Inc. All rights reserved worldwide.
 *
 * This sample source code is not supported by Cadence Design Systems, Inc.
 * It is provided freely for demonstration purposes only.
 * SEE THE WARRANTY DISCLAIMER AT THE BOTTOM OF THIS FILE.
 *
 ***************************************************************************/


 /****************************************************************************
 * TAU requires the following third party libraries.
 *    netcdf: https://www.unidata.ucar.edu/downloads/netcdf/
 *    hdf5  : https://www.hdfgroup.org/
 *    zlib  : http://www.zlib.net/
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
 * More info at WISH http://youtrack:8080/issue/PW-4411
 *
 ***************************************************************************/

#include "apiCAEP.h"
#include "apiCAEPUtils.h"
#include "apiGridModel.h"
#include "apiPWP.h"
#include "runtimeWrite.h"
#include "pwpPlatform.h"


#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <sstream>
#include <string>

#include <hdf5.h>
#include <netcdf.h>
#include <zlib.h>


using StringUINT32Map = std::map<std::string, PWP_UINT32>;
using StringStringMap = std::map<std::string, std::string>;


enum class Orientation { Pos = 1, Neg = -1 };

static const char *NetCDFDataModel = "DataModel";
static const char *Thickness = "Thickness";

#define OK(err)     (NC_NOERR == (err))

// Make sure types are compatible for calls to nc_xxxx() functions
static_assert(sizeof(int) == sizeof(PWP_UINT32), "int size mismatch");
static_assert(sizeof(double) == sizeof(PWGM_XYZVAL), "double size mismatch");

//****************************************************************************
//****************************************************************************
//****************************************************************************

struct TAU_DATA {
    TAU_DATA()
    {
        // nc_inq_libvers() returns "v.m.n other stuff". Trim off " other stuff"
        netcdfVer = nc_inq_libvers();
        netcdfVer  = netcdfVer.substr(0, netcdfVer.find(" "));

        zlibVer = zlibVersion();

        unsigned majnum = 0;
        unsigned minnum = 0;
        unsigned relnum = 0;
        const herr_t err = H5get_libversion(&majnum, &minnum, &relnum);
        if (err >= 0) {
            std::ostringstream os;
            os << majnum << "." << minnum << "." << relnum;
            hdf5Ver = os.str();
        }
    }

    // Main dimensions
    PWP_UINT32      numPoints{ 0 };
    PWP_UINT32      numElements{ 0 };
    PWP_UINT32      numTets{ 0 };
    PWP_UINT32      numPyramids{ 0 };
    PWP_UINT32      numPrisms{ 0 };
    PWP_UINT32      numHexes{ 0 };
    PWP_UINT32      numBars{ 0 };
    PWP_UINT32      numSurfaceElements{ 0 };
    PWP_UINT32      numSurfaceTris{ 0 };
    PWP_UINT32      numSurfaceQuads{ 0 };
    PWP_INT32       numSurfBndryMarkers{ 0 };

    // NetCdf file id values
    int             idNc{ 0 };
    int             idXc{ 0 };
    int             idYc{ 0 };
    int             idZc{ 0 };
    int             idTets{ 0 };
    int             idPyramids{ 0 };
    int             idPrisms{ 0 };
    int             idHexes{ 0 };
    int             idTris{ 0 };
    int             idQuads{ 0 };
    int             idMarker{ 0 };

    // Tracking values for elements created by extrusion during 2D export
    int             extrudedTriStart{ 0 };
    int             extrudedQuadStart{ 0 };

    // Maps user-defined BC name to serialized TAU id
    StringUINT32Map nameToId;

    // Maps user-defined BC name to its TAU BC phsical type name
    StringStringMap nameToType;

    // Used to track the indexing of hex and prism creation during 2D export.
    size_t          indexHex{ 0 };
    size_t          indexPrism{ 0 };

    // These values are for reference only. Not used by plugin.
    std::string     netcdfVer;
    std::string     zlibVer;
    std::string     hdf5Ver;
};


//****************************************************************************
//****************************************************************************
//****************************************************************************

template<size_t NumVerts, size_t BufSz>
class ElemBuffer {
    using BufType = std::unique_ptr<int[][NumVerts]>;

public:
    ElemBuffer(const int idNc, const int idArray) :
        idNc_(idNc),
        idArr_(idArray)
    {
    }

    ~ElemBuffer()
    {
        flush();
    }

    bool add(const PWGM_ELEMDATA &data)
    {
        bool ret = ((BufSz == bufUsed_) ? flush() : true);
        if (ret) {
            for (PWP_UINT32 k = 0; k < NumVerts; ++k) {
                buf_[bufUsed_][k] = (int)data.index[k];
            }
            ++ndx_;
            ++bufUsed_;
        }
        return ret;
    }

    bool flush()
    {
        bool ret = true;
        if (0 != bufUsed_) {
            const size_t start[2]{
                bufFirstNdx_,
                0
            };
            const size_t count[2]{
                bufUsed_,
                NumVerts
            };
            ret = OK(nc_put_vara_int(idNc_, idArr_, start, count, &buf_[0][0]));
            bufUsed_ = 0;
            bufFirstNdx_ = ndx_;
        }
        return ret;
    }

private:
    // netcdf mesh id
    int         idNc_{ 0 };

    // netcdf elemet array id
    int         idArr_{ 0 };

    // the element buffer
    BufType     buf_{ new int[BufSz][NumVerts] };

    // global serialized element index
    PWP_UINT32  ndx_{ 0 };

    // global serialized element index of first item in buf_
    size_t      bufFirstNdx_{ 0 };

    // number of elements waiting in buf_
    PWP_UINT32  bufUsed_{ 0 };
};


static bool
is2D(const CAEP_RTITEM &rti)
{
    return 0 != CAEPU_RT_DIM_2D(&rti);
}


/****************************************************************************
 * 
 * Initialize global condition information. It fills in the nameToId
 * map and the nameToType map. For both, the keys are given by a value of
 * condition.name, while the values are given by the corresponding call to
 * condition.tid and condition.type (respectively)
 * 
 ***************************************************************************/
static bool
setupMarkerInfo(const CAEP_RTITEM &rti)
{
    TAU_DATA &tau = *(rti.tau);
    bool ret = true;
    PWP_UINT32 nextID = 0;
    const PWP_UINT32 NumDoms = PwModDomainCount(rti.model);
    for (PWP_UINT32 i = 0; i < NumDoms && ret; ++i) {
        PWGM_HDOMAIN dom = PwModEnumDomains(rti.model, i);
        PWGM_CONDDATA condition{ 0 };
        ret = ret && PwDomCondition(dom, &condition);
        if (ret && tau.nameToId.find(condition.name) == tau.nameToId.end()) {
            tau.nameToId[condition.name] = nextID++;
            tau.nameToType[condition.name] = condition.type;
        }
    }
    if (ret && is2D(rti)) {
        // For 2D export, the grid is extruded to a 1-cell thick 3D grid.
        // The original 2D base tri/quad elements are placed in the 'Symmetry_1'
        // BC patch. The corresponding extruded top tri/quad elements are placed
        // in the 'Symmetry_2' BC patch. The extruded side quads inherit their
        // BC settings from the corresponding 2D boundary bar elements.
        if (tau.nameToId.find("Symmetry_1") == tau.nameToId.end()) {
            tau.nameToId["Symmetry_1"] = nextID++;
            tau.nameToType["Symmetry_1"] = "Symmetry";
        }
        if (tau.nameToId.find("Symmetry_2") == tau.nameToId.end()) {
            tau.nameToId["Symmetry_2"] = nextID++;
            tau.nameToType["Symmetry_2"] = "Symmetry";
        }
    }
    return ret;
}


static bool
getDataModel(const CAEP_RTITEM &rti, int &dataModel)
{
    // NetCDF_Classic|NetCDF_64bit|NetCDF4_HDF5|
    //              0|           1|           2|
    PWP_UINT val = 0;
    const bool ret = PwModGetAttributeUINT(rti.model, NetCDFDataModel, &val);
    if (ret) {
        switch (val) {
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
    }
    return ret;
}


/****************************************************************************
 * 
 * This function opens the file with a call to nc_create(), and creates the
 * global "type" attribute as requested by the TAU file format.
 * 
 ***************************************************************************/
static bool
startGridFile(const CAEP_RTITEM &rti)
{
    TAU_DATA &tau = *(rti.tau);
    const std::string type("Primary Grid: Tau Format");
    std::string file(rti.pWriteInfo->fileDest);
    file += ".grid";
    int dataModel = 0;
    return getDataModel(rti, dataModel) &&
        OK(nc_create(file.c_str(), NC_CLOBBER | dataModel, &tau.idNc)) &&
        OK(nc_put_att_text(tau.idNc, NC_GLOBAL, "type", type.size(),
            type.c_str()));
}


/****************************************************************************
 * 
 * Creates the "no_of_points" dimension, and the points_xc, points_yc, and
 * points_zc variables that are based on it. These variables are used to store
 * the xyz coordinates of Pointwise vertices.
 * 
 ***************************************************************************/
static bool
definePointVars(const CAEP_RTITEM &rti)
{
    TAU_DATA &tau = *(rti.tau);

    // capture number of points in the grid
    tau.numPoints = PwModVertexCount(rti.model);

    // If 2D, the number of points is doubled because the 2D faces will be
    // extruded one step in either the +Y or -Y direction. This will implicitly
    // create a second set of points that are the exact same order as the
    // original set, simply offset by the number of original points.
    // For instance, given a 2D grid with 8 points indexed 0-7, there will be 8
    // extruded points index as 8-15.
    const PWP_UINT32 NcNumPts = tau.numPoints * (is2D(rti) ? 2 : 1);

    // Create the no_of_points dimension and the XYZ coordinate vars
    int id; // value ignored
    return OK(nc_def_dim(tau.idNc, "no_of_points", NcNumPts, &id)) &&
        OK(nc_def_var(tau.idNc, "points_xc", NC_DOUBLE, 1, &id, &tau.idXc)) &&
        OK(nc_def_var(tau.idNc, "points_yc", NC_DOUBLE, 1, &id, &tau.idYc)) &&
        OK(nc_def_var(tau.idNc, "points_zc", NC_DOUBLE, 1, &id, &tau.idZc));
}


/****************************************************************************
 *
 * Add the variables for an element type
 *
 ***************************************************************************/
static bool
defineElementVar(const int idNc, const char *numElemsLbl, const size_t numElems,
    const char *ptsPerElemLbl, const size_t ptsPerElem, const char *idElemsLbl,
    int &idElems)
{
    bool ret = true;
    // make the two dimensions and variable for each shape type. Shape types
    // with no instances do not get created (as in file-spec)
    if (0 != numElems) {
        int ids[2]{ 0, 0 };
        ret = OK(nc_def_dim(idNc, numElemsLbl, numElems, &ids[0])) &&
            OK(nc_def_dim(idNc, ptsPerElemLbl, ptsPerElem, &ids[1])) &&
            OK(nc_def_var(idNc, idElemsLbl, NC_INT, 2, ids, &idElems));
    }
    return ret;
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
 * prism: prism
 * hexaeder: hexahedron
 * 
 ***************************************************************************/
static bool
defineElementVars(const CAEP_RTITEM &rti)
{
    TAU_DATA &tau = *(rti.tau);

    // count number of different element types
    tau.numElements = 0;
    tau.numTets = 0;
    tau.numPyramids = 0;
    tau.numPrisms = 0;
    tau.numHexes = 0;
    const PWP_UINT32 NumBlocks = PwModBlockCount(rti.model);
    PWGM_ELEMCOUNTS elementCounts;
    for (PWP_UINT32 i = 0; i < NumBlocks; ++i) {
        PWGM_HBLOCK hBlk = PwModEnumBlocks(rti.model, i);
        tau.numElements += PwBlkElementCount(hBlk, &elementCounts);
        tau.numTets += elementCounts.count[PWGM_ELEMTYPE_TET];
        tau.numPyramids += elementCounts.count[PWGM_ELEMTYPE_PYRAMID];
        tau.numPrisms += elementCounts.count[PWGM_ELEMTYPE_WEDGE];
        tau.numHexes += elementCounts.count[PWGM_ELEMTYPE_HEX];
    }

    int id = 0; // don't need this
    return OK(nc_def_dim(tau.idNc, "no_of_elements", tau.numElements, &id)) &&
        defineElementVar(tau.idNc, "no_of_tetraeders", tau.numTets,
            "points_per_tetraeder", 4, "points_of_tetraeders", tau.idTets) &&
        defineElementVar(tau.idNc, "no_of_pyramids", tau.numPyramids,
            "points_per_pyramid", 5, "points_of_pyramids", tau.idPyramids) &&
        defineElementVar(tau.idNc, "no_of_prisms", tau.numPrisms,
            "points_per_prism", 6, "points_of_prisms", tau.idPrisms) &&
        defineElementVar(tau.idNc, "no_of_hexaeders", tau.numHexes,
            "points_per_hexaeder", 8, "points_of_hexaeders", tau.idHexes);
}


/****************************************************************************
 * 
 * The defineMarkers(CAEP_RTITEM) function
 * 
 * Defines the marker variable and marker attributes (attributes are defined
 * using the 'defineOneMarker()' function).
 * 
 ***************************************************************************/
static bool
defineMarkers(const CAEP_RTITEM &rti)
{
    TAU_DATA &tau = *(rti.tau);
    int ncNumMarkers = 0;
    bool ret = OK(nc_def_dim(tau.idNc, "no_of_markers",
            tau.nameToId.size(), &ncNumMarkers)) &&
        OK(nc_def_var(tau.idNc, "marker", NC_INT, 1, &ncNumMarkers,
            &tau.idMarker));
    if (ret) {
        for (auto it = tau.nameToId.cbegin(); tau.nameToId.end() != it; ++it) {
            std::ostringstream attr;
            attr << "marker_" << it->second;
            if (!OK(nc_put_att_text(tau.idNc, NC_GLOBAL, attr.str().c_str(),
                    it->first.size(), it->first.c_str()))) {
                ret = false;
                break;
            }
        }
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
static bool
defineSurfaceElementVars(const CAEP_RTITEM &rti)
{
    TAU_DATA &tau = *(rti.tau);
    const int idNc = tau.idNc;

    // count number of different surface-element types
    PWP_UINT32 numSurfaceTriangles = 0;
    PWP_UINT32 numSurfaceQuadrilaterals = 0;
    PWP_UINT32 numDomains = PwModDomainCount(rti.model);
    for (PWP_UINT32 i = 0; i < numDomains; ++i) {
        PWGM_HDOMAIN domain = PwModEnumDomains(rti.model, i);
        PWGM_ELEMCOUNTS elementCounts;
        PwDomElementCount(domain, &elementCounts);
        numSurfaceTriangles += elementCounts.count[PWGM_ELEMTYPE_TRI];
        numSurfaceQuadrilaterals += elementCounts.count[PWGM_ELEMTYPE_QUAD];
    }
    // DO NOT use return value from pwDomElementCount, as it includes bar 
    // elements as surface elements
    tau.numSurfaceElements = numSurfaceTriangles + numSurfaceQuadrilaterals;

    // create dimension for total number of surface elements
    int ncNumSurfaceElements = 0;
    bool ret = OK(nc_def_dim(idNc, "no_of_surfaceelements",
        tau.numSurfaceElements, &ncNumSurfaceElements));

    // create variables (and related dimensions) for each type of surface
    // element
    if (ret && (0 != numSurfaceTriangles)) {
        // create an id for the no_of_surfacetriangles dimension
        int ncSurfaceTriangles = 0;
        ret = OK(nc_def_dim(idNc, "no_of_surfacetriangles", numSurfaceTriangles,
            &ncSurfaceTriangles));

        // create an id for the points_per_surfacetriangle dimension
        int ncPointsPerTri = 0;
        ret = ret && OK(nc_def_dim(idNc, "points_per_surfacetriangle", 3,
            &ncPointsPerTri));

        // create an id for the points_of_surfacetriangles variable, based on
        // the no_of_surfacetriangles and points_per_surfacetriangle
        // dimensions.
        int dimids[2] = { ncSurfaceTriangles, ncPointsPerTri };
        ret = ret && OK(nc_def_var(idNc, "points_of_surfacetriangles", NC_INT,
            2, dimids, &tau.idTris));
    }
    if (ret && (0 != numSurfaceQuadrilaterals)) {
        // see 'if (numSurfaceTriangles != 0)' comments
        int ncSurfaceQuads = 0;
        ret = OK(nc_def_dim(idNc, "no_of_surfacequadrilaterals",
            numSurfaceQuadrilaterals, &ncSurfaceQuads));

        int ncPointsPerQuad = 0;
        ret = ret && OK(nc_def_dim(idNc, "points_per_surfacequadrilateral", 4,
            &ncPointsPerQuad));

        int dimids[2] = { ncSurfaceQuads, ncPointsPerQuad };
        ret = ret && OK(nc_def_var(idNc, "points_of_surfacequadrilaterals",
            NC_INT, 2, dimids, &(tau.idQuads)));
    }

    // Define variable (netCDF multi-dimensional array) that associates
    // domains (boundary panels) with surface elements
    int dimids[1] = { ncNumSurfaceElements };

    // Define variable for boundarymarkers
    ret = ret && OK(nc_def_var(idNc, "boundarymarker_of_surfaces", NC_INT, 1,
        dimids, &tau.numSurfBndryMarkers));

    // indexing scheme in the boundarypanel variable and in the 
    // boundarymarker variable saves the domain (boundary-panel) or 
    // boundary condition (boundary-marker) number of all the triangles
    // followed by the domain number of all the quads, so this indexing
    // tells us where these are.
    tau.extrudedTriStart = 0;
    tau.extrudedQuadStart = numSurfaceTriangles;
    
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
getNumExtrudedBndryBars(const CAEP_RTITEM &rti)
{
    PWGM_ELEMCOUNTS elementCounts;
    PWP_UINT32 totalBars = 0;
    const PWP_UINT32 numDomains = PwModDomainCount(rti.model);
    for (PWP_UINT32 i = 0; i < numDomains; ++i) {
        PWGM_HDOMAIN domain = PwModEnumDomains(rti.model, i);
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
static bool
defineElementVars2D(const CAEP_RTITEM &rti)
{
    TAU_DATA &tau = *(rti.tau);

    // count number of different surface-element types
    PWP_UINT32 numSurfaceElements = 0;
    PWP_UINT32 numSurfaceTriangles = 0;
    PWP_UINT32 numSurfaceQuadrilaterals = 0;
    PWP_UINT32 numHexes = 0;
    PWP_UINT32 numPrisms = 0;
    PWP_UINT32 numBlocks = PwModBlockCount(rti.model);

    for (PWP_UINT32 i = 0; i < numBlocks; ++i) {
        PWGM_HBLOCK hBlk = PwModEnumBlocks(rti.model, i);
        PWGM_ELEMCOUNTS elementCounts;
        PwBlkElementCount(hBlk, &elementCounts);

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

    // All bars will become surface quadrilaterals
    numSurfaceQuadrilaterals += getNumExtrudedBndryBars(rti);
    numSurfaceElements = numSurfaceTriangles + numSurfaceQuadrilaterals;

    tau.numSurfaceElements = numSurfaceElements;
    tau.numSurfaceTris = numSurfaceTriangles;
    tau.numSurfaceQuads = numSurfaceQuadrilaterals;
    tau.numHexes = numHexes;
    // Please note, prisms are the same as prisms in context of TAU.
    tau.numPrisms = numPrisms;

    const int idNc = tau.idNc;

    // create dimension for total number of surface elements
    int ncNumSurfaceElements = 0;
    bool ret = OK(nc_def_dim(idNc, "no_of_surfaceelements", numSurfaceElements,
        &ncNumSurfaceElements));

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
            2, dimids, &(tau.idTris)));
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
            NC_INT, 2, dimids, &(tau.idQuads)));
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
        dimids, &tau.numSurfBndryMarkers));

    // indexing scheme in the boundarypanel variable and in the 
    // boundarymarker variable saves the domain (boundary-panel) or 
    // boundary condition (boundary-marker) number of all the triangles
    // followed by the domain number of all the quads, so this indexing
    // tells us where these are.
    tau.extrudedTriStart = 0;
    tau.extrudedQuadStart = numSurfaceTriangles;
    
    return ret;
}


/****************************************************************************
 * 
 * This function ends the file definition (where netcdf dimensions and netcdf
 * variables are created), so that data can be written out to the file.
 * 
 ***************************************************************************/
static bool
endFileDefinition(const CAEP_RTITEM &rti)
{
    return OK(nc_enddef(rti.tau->idNc));
}


/****************************************************************************
 * 
 * Uses the points of a pair of quads and creates a hex from the vertices. The
 * index of the hexes are being tracked by the structure tau. The ordering of
 * the points for hexes is the same in Tau as it is in Pointwise.
 *
 ***************************************************************************/
static bool
createHexFromQuad(const CAEP_RTITEM &rti, const PWGM_ELEMDATA &data)
{
    const TAU_DATA &tau = *rti.tau;
    const size_t start[2]{ tau.indexHex, 0 };
    const size_t count[2]{ 1, 8 };
    // The points originally saved for the creation of the new quad are reversed
    // of the order expected for the TAU format. We are manually reversing the
    // order for the saves so that the hex has the right dimensions.
    const int hexPoints[8]{
        (int)data.index[0],
        (int)data.index[1],
        (int)data.index[2],
        (int)data.index[3],
        (int)(data.index[0] + tau.numPoints),
        (int)(data.index[1] + tau.numPoints),
        (int)(data.index[2] + tau.numPoints),
        (int)(data.index[3] + tau.numPoints)
    };

    return OK(nc_put_vara_int(tau.idNc, tau.idHexes, start, count, hexPoints));
}


/****************************************************************************
 * 
 * Uses the points of a pair of triangles to create a prism shape. Because of
 * the way triangles are defined in TAU format, the order of the points for the
 * second triangle must be the same direction as the points are defined to
 * maintain the normal facing away from the face.
 * 
 ***************************************************************************/
static bool
createPrismFromTri(const CAEP_RTITEM &rti, const PWGM_ELEMDATA &data)
{
    const TAU_DATA &tau = *rti.tau;
    const size_t start[2]{ tau.indexPrism, 0 };
    const size_t count[2]{ 1, 6 };
    const int indices[6]{
        (int)data.index[0],
        (int)data.index[1],
        (int)data.index[2],
        (int)(data.index[0] + tau.numPoints),
        (int)(data.index[1] + tau.numPoints),
        (int)(data.index[2] + tau.numPoints)
    };
    return OK(nc_put_vara_int(tau.idNc, tau.idPrisms, start, count, indices));
}


/****************************************************************************
 * 
 * Uses the points of the original Quad and creates an offset Quad from the
 * vertices. The function then calls createHexFromQuad() to create the
 * associated hex element. TAU boundary element normals must point out of the
 * extruded volume. The PW convention is the opposite. As such, the original
 * quad has its point order reversed and the extruded quad can use the original
 * order.
 * 
 ***************************************************************************/
static bool
createExtrudedQuadTopElement(const CAEP_RTITEM &rti, const PWP_UINT32 idxQuad,
    const PWGM_ELEMDATA &data)
{
    // the extruded quad is just an offset of the original indices
    const TAU_DATA &tau = *(rti.tau);
    const size_t start[2]{ idxQuad, 0 };
    const size_t count[2]{ 1, 4 };
    const PWP_UINT32 markerID = tau.nameToId.at("Symmetry_2");
    const size_t adjIdx = idxQuad + tau.extrudedQuadStart;
    const int indices[4]{
        (int)(data.index[0] + tau.numPoints),
        (int)(data.index[1] + tau.numPoints),
        (int)(data.index[2] + tau.numPoints),
        (int)(data.index[3] + tau.numPoints)
    };
    return OK(nc_put_vara_int(tau.idNc, tau.idQuads, start, count, indices)) &&
        OK(nc_put_var1_int(tau.idNc, tau.numSurfBndryMarkers, &adjIdx,
            (int*)&markerID)) &&
        createHexFromQuad(rti, data);
}


/****************************************************************************
 * 
 * Uses the points of the original triangle and creates an offset triangle from
 * the vertices. The function then calls createPrismFromTri() to create the
 * associated prism element. TAU boundary element normals must point out of the
 * extruded volume. The PW convention is the opposite. As such, the original tri
 * has its point order reversed and the extruded tri can use the original order.
 * 
 ***************************************************************************/
static bool
createExtrudedTriTopElement(const CAEP_RTITEM &rti, const PWP_UINT32 idxTri,
    const PWGM_ELEMDATA &data)
{
    // the extruded tri is just an offset of the original indices
    const TAU_DATA &tau = *rti.tau;
    const size_t start[2]{ idxTri, 0 };
    const size_t count[2]{ 1, 3 };
    const PWP_UINT32 markerID = tau.nameToId.at("Symmetry_2");
    const size_t adjIdx = idxTri + tau.extrudedTriStart;
    const int indices[3]{
        (int)(data.index[0] + tau.numPoints),
        (int)(data.index[1] + tau.numPoints),
        (int)(data.index[2] + tau.numPoints)
    };
    return OK(nc_put_vara_int(tau.idNc, tau.idTris, start, count,
            indices)) &&
        OK(nc_put_var1_int(tau.idNc, tau.numSurfBndryMarkers, &adjIdx,
            (int*)&markerID)) &&
        createPrismFromTri(rti, data);
}


/****************************************************************************
 * 
 * Uses the points of a bar element to create a quad element seperate from the
 * quads being used to create hexes. This will also expand the boundary
 * conditions to the new surface quad.
 * 
 ***************************************************************************/
static bool
createExtrudedBndryQuads(CAEP_RTITEM &rti, PWP_UINT32 &idxQuad)
{
    const TAU_DATA &tau = *rti.tau;
    bool ret = true;
    const PWP_UINT32 numDomains = PwModDomainCount(rti.model);
    for (PWP_UINT32 ii = 0; ii < numDomains && ret; ++ii) {
        const PWGM_HDOMAIN hDom = PwModEnumDomains(rti.model, ii);
        PWGM_CONDDATA condData;
        if (!PwDomCondition(hDom, &condData)) {
            ret = false;
            break;
        }
        const PWP_UINT32 markerID = tau.nameToId.at(condData.name);
        const PWP_UINT32 numBars = PwDomElementCount(hDom, nullptr);
        PWGM_ELEMDATA data;
        for (PWP_UINT32 jj = 0; jj < numBars && ret; ++jj) {
            if (!PwElemDataMod(PwDomEnumElements(hDom, jj), &data)) {
                ret = false;
                break;
            }
            if (PWGM_ELEMTYPE_BAR == data.type) {
                // Use the bar's end point indices and their corresponding
                // offset indices to define an extruded side-quad. By TAU
                // convention, the quad normal must point out of the extruded
                // volume.
                const size_t start[2]{ idxQuad, 0 };
                const size_t count[2]{ 1, 4 };
                const int indices[4]{
                    (int)data.index[0],
                    (int)(data.index[0] + tau.numPoints),
                    (int)(data.index[1] + tau.numPoints),
                    (int)data.index[1]
                };
                const size_t adjIdx = idxQuad + tau.extrudedQuadStart;
                if (!OK(nc_put_vara_int(tau.idNc, tau.idQuads, start,
                        count, indices))) {
                    ret = false;
                    break;
                }
                if (!OK(nc_put_var1_int(tau.idNc, tau.numSurfBndryMarkers,
                        &adjIdx, (int*)&markerID))) {
                    ret = false;
                    break;
                }
                ++idxQuad;
            }
            if (!caeuProgressIncr(&rti)) {
                ret = false;
                break;
            }
        }
    }
    return ret;
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
static bool
get2DGridOrientation(CAEP_RTITEM &rti, Orientation &orient)
{
    PWGM_ELEMDATA ed{ PWGM_ELEMTYPE_SIZE };
    PWGM_XYZVAL coords[3][3];

    auto getXYZ = [](const PWGM_HVERTEX &hVert, PWGM_XYZVAL xyz[3])
    {
        return PwVertXyzVal(hVert, PWGM_XYZ_X, &xyz[0]) &&
            PwVertXyzVal(hVert, PWGM_XYZ_Y, &xyz[1]) &&
            PwVertXyzVal(hVert, PWGM_XYZ_Z, &xyz[2]);
    };

    auto setCoords = [&rti, &getXYZ, &ed, &coords](const PWP_UINT32 blkNdx)
    {
        const PWGM_HELEMENT hElem = PwBlkEnumElements(
                PwModEnumBlocks(rti.model, blkNdx), 0);
        return PwElemDataMod(hElem, &ed) && (ed.vertCnt > 2) &&
            getXYZ(ed.vert[0], coords[0]) && getXYZ(ed.vert[1], coords[1]) &&
            getXYZ(ed.vert[ed.vertCnt - 1], coords[2]);
    };

    auto getBlkOrientation = [&setCoords, &coords](const PWP_UINT32 blkNdx,
        Orientation &orient, bool &isYAligned)
    {
        const bool ret = setCoords(blkNdx);
        if (ret) {
            const PWGM_XYZVAL a[3]{
                coords[1][0] - coords[0][0],
                coords[1][1] - coords[0][1],
                coords[1][2] - coords[0][2]
            };
            const PWGM_XYZVAL b[3]{
                coords[2][0] - coords[0][0],
                coords[2][1] - coords[0][1],
                coords[2][2] - coords[0][2]
            };
            // THIS ASSUMES THE GRID IS CONSTRUCTED IN the XZ PLANE:
            // The orientation is determined by the the Y component of
            // (a cross b).
            // The cross product formula:
            //   <cx, cy, cz> = < ay*bz-az*by, az*bx-ax*bz, ax*by-ay*bx >
            enum { X, Y, Z };
            const PWP_REAL cX = (a[Y] * b[Z]) - (a[Z] * b[Y]); // ay*bz-az*by
            const PWP_REAL cY = (a[Z] * b[X]) - (a[X] * b[Z]); // az*bx-ax*bz
            const PWP_REAL cZ = (a[X] * b[Y]) - (a[Y] * b[X]); // ax*by-ay*bx
            // TODO: improve mesh alignment processing
            const PWP_REAL Tol = 1e-6;
            if (std::abs(cX) > Tol) {
                isYAligned = false; // mesh normal is NOT Y-aligned
            }
            if (std::abs(cZ) > Tol) {
                isYAligned = false; // mesh normal is NOT Y-aligned
            }
            orient = (cY > 0.0) ? Orientation::Pos : Orientation::Neg;
        }
        return ret;
    };

    // assume orient of first block is same for all blocks
    bool isYAligned = true; // may be set to false by getBlkOrientation()
    bool ret = getBlkOrientation(0, orient, isYAligned);
    if (ret) {
        const PWP_UINT32 numBlocks = PwModBlockCount(rti.model);
        for (PWP_UINT32 ii = 1; ii < numBlocks; ++ii) {
            Orientation iiOrient;
            if (!getBlkOrientation(ii, iiOrient, isYAligned)) {
                ret = false;
                break;
            }
            if (iiOrient != orient) {
                // Orientation change - issue warning and stop processing
                caeuSendWarningMsg(&rti, "Non-uniform domain orientation.", 0);
                break;
            }
        }
        if (ret && !isYAligned) {
            caeuSendWarningMsg(&rti, "2-D grid normal should be Y-aligned.", 0);
        }
    }
    return ret;
}


// Write out the x, y, or z coordinate values 'passes' times. The coords are
// adjusted by 'offset * pass' (0 to passes-1).
static bool
writePtCoordArray(CAEP_RTITEM &rti, const PWGM_ENUM_XYZ which,
    const PWP_UINT32 passes = 1, const PWP_REAL offset = 0.0)
{
    const TAU_DATA &tau = *(rti.tau);
    bool ret = (0 != passes);
    int cid = 0;
    // to which coord array are we writing?
    switch (which) {
    case PWGM_XYZ_X:
        cid = tau.idXc;
        break;
    case PWGM_XYZ_Y:
        cid = tau.idYc;
        break;
    case PWGM_XYZ_Z:
        cid = tau.idZc;
        break;
    default:
        ret = false;
    }
    if (ret) {
        // For speed, buffer coord values in BufSz chunks
        constexpr size_t BufSz{ 1024 * 1024 };
        std::unique_ptr<PWGM_XYZVAL[]> buf(new PWGM_XYZVAL[BufSz]);
        size_t bufUsed = 0;  // number of coords in buf
        size_t startNdx = 0; // index of first coord in buf
        PWP_UINT32 ndx = 0;  // serialized coord index for all passes

        const PWP_UINT32 NumPts = PwModVertexCount(rti.model);
        for (PWP_UINT32 pass = 0; pass < passes; ++pass) {
            const PWP_REAL PassOffset = offset * pass;
            for (PWP_UINT32 ii = 0; ii < NumPts; ++ii, ++ndx) {
                if (bufUsed < BufSz) {
                    // buf still has room - keep going
                }
                else if (!OK(nc_put_vara_double(tau.idNc, cid, &startNdx,
                        &BufSz, buf.get()))) {
                    // could not flush buf to file
                    ret = false;
                    break;
                }
                else {
                    // buf was flushed to file
                    bufUsed = 0;
                    startNdx = ndx;
                }

                // Append coord to buffer
                const PWGM_HVERTEX vertex = PwModEnumVertices(rti.model, ii);
                if (PwVertXyzVal(vertex, which, &buf[bufUsed])) {
                    buf[bufUsed++] += PassOffset;
                }
                else {
                    ret = false;
                    break;
                }

                if (!caeuProgressIncr(&rti)) {
                    ret = false;
                    break;
                }
            }
        }
        if (ret && (0 != bufUsed)) {
            // flush partial buffer
            ret = OK(nc_put_vara_double(tau.idNc, cid, &startNdx, &bufUsed,
                buf.get()));
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
static bool
writePointVars3D(CAEP_RTITEM &rti)
{
    const PWP_UINT32 NumPoints = PwModVertexCount(rti.model);
    const bool ret = caeuProgressBeginStep(&rti, 3 * NumPoints) &&
        writePtCoordArray(rti, PWGM_XYZ_X) &&
        writePtCoordArray(rti, PWGM_XYZ_Y) &&
        writePtCoordArray(rti, PWGM_XYZ_Z);
    return caeuProgressEndStep(&rti) && ret;
}


/****************************************************************************
 * 
 * Iterate through the vertices in the Pointwise model and save their xyz
 * coordinates to their respective TAU variables. For 2D, this process is
 * repeated, offsetting the Y component by the amount defined by the user.
 * 
 ***************************************************************************/
static bool
writePointVars2D(CAEP_RTITEM &rti)
{
    const PWP_UINT32 NumPoints = PwModVertexCount(rti.model);
    const PWP_UINT32 Passes{ 2 };
    Orientation dir;
    PWP_REAL thickness;
    const bool ret = caeuProgressBeginStep(&rti, 3 * NumPoints) &&
        get2DGridOrientation(rti, dir) &&
        PwModGetAttributeREAL(rti.model, Thickness, &thickness) &&
        writePtCoordArray(rti, PWGM_XYZ_X, Passes) &&
        writePtCoordArray(rti, PWGM_XYZ_Y, Passes, thickness *
            static_cast<PWP_REAL>(dir)) &&
        writePtCoordArray(rti, PWGM_XYZ_Z, Passes);
    return caeuProgressEndStep(&rti) && ret;
}


/****************************************************************************
 * 
 * Iterate through the blocks in the Pointwise model and for each block, iterate
 * through the block-elements and save them to their respective TAU variables,
 * based on their element-type.
 * 
 ***************************************************************************/
static bool
writeElementVars(CAEP_RTITEM &rti)
{
    TAU_DATA &tau = *(rti.tau);
    bool ret = true;
    if (caeuProgressBeginStep(&rti, tau.numElements)) {
        // Cache up to 512K cells per buffer
        ElemBuffer<4, 512 * 1024> tetBuf(tau.idNc, tau.idTets);
        ElemBuffer<5, 512 * 1024> pyrBuf(tau.idNc, tau.idPyramids);
        ElemBuffer<6, 512 * 1024> priBuf(tau.idNc, tau.idPrisms);
        ElemBuffer<8, 512 * 1024> hexBuf(tau.idNc, tau.idHexes);
        PWGM_ELEMDATA data;
        const PWP_UINT32 numBlocks = PwModBlockCount(rti.model);
        for (PWP_UINT32 i = 0; i < numBlocks && ret; ++i) {
            const PWGM_HBLOCK hBlk = PwModEnumBlocks(rti.model, i);
            const PWP_UINT32 numElements = PwBlkElementCount(hBlk, 0);
            for (PWP_UINT32 j = 0; j < numElements && ret; ++j) {
                const PWGM_HELEMENT element = PwBlkEnumElements(hBlk, j);
                if (!PwElemDataMod(element, &data)) {
                    ret = false;
                    break;
                }
                // note that the indexing scheme for EVERY element type is the
                // same for TAU and CAE, so they don't need to be re-indexed.
                switch(data.type) {
                case PWGM_ELEMTYPE_TET:
                    tetBuf.add(data);
                    break;
                case PWGM_ELEMTYPE_PYRAMID:
                    pyrBuf.add(data);
                    break;
                case PWGM_ELEMTYPE_WEDGE:
                    priBuf.add(data);
                    break;
                case PWGM_ELEMTYPE_HEX: {
                    hexBuf.add(data);
                    break; }
                default:
                    // this should not occur; somehow an unsupported type was
                    // given.
                    ret = false;
                    break;
                }
                if (!caeuProgressIncr(&rti)) {
                    ret = false;
                    break;
                }
            }
        }
        if (ret) {
            // No errors so far - capture final flush() status
            ret = tetBuf.flush() && pyrBuf.flush() && priBuf.flush() &&
                hexBuf.flush();
        }
    }
    return caeuProgressEndStep(&rti) && ret;
}


/****************************************************************************
 * 
 * For each block, iterate through the block-elements and save them to their
 * respective TAU variables, based on their element-type. After saving the
 * element, the program then creates the offset variable and the associated 3D
 * element.
 * 
 ***************************************************************************/
static bool
writeElementVars2D(CAEP_RTITEM &rti)
{
    TAU_DATA &tau = *(rti.tau);
    bool ret = true;
    // separately track our indexing into the tri and quad netcdf-variables.
    tau.indexHex = 0;
    tau.indexPrism = 0;
    // iterate through each 2D block's elements
    PWP_UINT32 idxQuad = 0;
    const PWP_UINT32 numBlocks = PwModBlockCount(rti.model);
    if (caeuProgressBeginStep(&rti, tau.numSurfaceElements - tau.numBars)) {
        const int idNc = tau.idNc;
        PWP_UINT32 idxTri = 0;
        PWGM_ELEMDATA data;
        // This initializes the variable with a dummy variable which will
        // later be changed. In 3D a domain can not have a bar element.
        for (PWP_UINT32 i = 0; i < numBlocks && ret; ++i) {
            const PWGM_HBLOCK hBlk = PwModEnumBlocks(rti.model, i);
            const PWP_UINT32 markerID = tau.nameToId["Symmetry_1"];
            const PWP_UINT32 numSurfaceElements = PwBlkElementCount(hBlk, 0);
            for (PWP_UINT32 j = 0; j < numSurfaceElements && ret; ++j) {
                if (!PwElemDataMod(PwBlkEnumElements(hBlk, j), &data)) {
                    ret = false;
                    break;
                }
                switch(data.type) {
                case PWGM_ELEMTYPE_TRI: {
                    // the point order of TAU triangles is the reverse of PW
                    const size_t start[2] = { idxTri, 0 };
                    const size_t count[2] = { 1, 3 };
                    const int indices[3] = {
                        (int)data.index[2],
                        (int)data.index[1],
                        (int)data.index[0]
                    };
                    if (!OK(nc_put_vara_int(idNc, tau.idTris, start, count,
                            indices))) {
                        ret = false;
                        break;
                    }
                    const size_t adjIdx = idxTri + tau.extrudedTriStart;
                    if (!OK(nc_put_var1_int(idNc, tau.numSurfBndryMarkers,
                            &adjIdx, (int*)&markerID))) {
                        ret = false;
                        break;
                    }
                    ++idxTri;
                    // Create the new triangle element that is extruded from
                    // this element.
                    if (!createExtrudedTriTopElement(rti, idxTri, data)) {
                        ret = false;
                        break;
                    }
                    ++idxTri;
                    ++tau.indexPrism;
                    break; }
                case PWGM_ELEMTYPE_QUAD: {
                    // see comments for PWGM_ELEMTYPE_TRI:
                    const size_t start[2] = { idxQuad, 0 };
                    const size_t count[2] = { 1, 4 };
                    const int indices[4] = {
                        (int)data.index[3],
                        (int)data.index[2],
                        (int)data.index[1],
                        (int)data.index[0]
                    };
                    if (!OK(nc_put_vara_int(idNc, tau.idQuads, start, count,
                            indices))) {
                        ret = false;
                        break;
                    }
                    const size_t adjIdx = idxQuad + tau.extrudedQuadStart;
                    if (!OK(nc_put_var1_int(idNc, tau.numSurfBndryMarkers,
                            &adjIdx, (int*)&markerID))) {
                        ret = false;
                        break;
                    }
                    ++idxQuad;
                    if (!createExtrudedQuadTopElement(rti, idxQuad, data)) {
                        ret = false;
                        break;
                    }
                    ++idxQuad;
                    ++tau.indexHex;
                    break; }
                default:
                    break;
                }
                if (!caeuProgressIncr(&rti)) {
                    ret = false;
                    break;
                }
            }
        }
        caeuProgressEndStep(&rti);
    }

    // After all of the quads and tris are written, the bars need to be extruded
    // and written to the quad variable. This requires continuing to track the
    // index of the quad.
    return ret && createExtrudedBndryQuads(rti, idxQuad);
}


/****************************************************************************
 * 
 * Writes out the marker numbers. Due to the lack of clarity in the specs, I
 * have assured that each value is the same as its index.
 * 
 ***************************************************************************/
static bool
writeMarkers(const CAEP_RTITEM &rti)
{
    const TAU_DATA &tau = *(rti.tau);
    bool ret = true;
    for (size_t index = 0; index < tau.nameToId.size(); ++index) {
        const int op = (int)index;
        if (!OK(nc_put_var1_int(tau.idNc, tau.idMarker, &index, &op))) {
            ret = false;
            break;
        }
    }
    return ret;
}


/****************************************************************************
 * 
 * For each domain, iterate through the domain-elements and save them to their
 * respective TAU variables, based on their element-type.
 * 
 ***************************************************************************/
static bool
writeSurfaceElementVars(CAEP_RTITEM &rti)
{
    const TAU_DATA &tau = *(rti.tau);
    // iterate through the domain elements and output to the file
    bool ret = (0 != caeuProgressBeginStep(&rti, tau.numSurfaceElements));
    if (ret) {
        // separately track indexing into the tri and quad netcdf-variables
        PWP_UINT32 idxTri = 0;
        PWP_UINT32 idxQuad = 0;
        const PWP_UINT32 numDoms = PwModDomainCount(rti.model);
        for (PWP_UINT32 ii = 0; ii < numDoms && ret; ++ii) {
            PWGM_CONDDATA condData{ "" };
            const PWGM_HDOMAIN hDom = PwModEnumDomains(rti.model, ii);
            if (!PwDomCondition(hDom, &condData)) {
                ret = false;
                break;
            }
            const PWP_UINT32 markerID = tau.nameToId.at(condData.name);
            const PWP_UINT32 numSurfElems = PwDomElementCount(hDom, 0);
            PWGM_ELEMDATA data;
            for (PWP_UINT32 jj = 0; jj < numSurfElems && ret; ++jj) {
                if (!PwElemDataMod(PwDomEnumElements(hDom, jj), &data)) {
                    ret = false;
                    break;
                }
                switch(data.type) {
                case PWGM_ELEMTYPE_TRI: {
                    // save vertex information for current triangle
                    const size_t start[2]{ idxTri, 0 };
                    const size_t count[2]{ 1, 3 };
                    // unfortunately the point order of surface triangles in
                    // TAU is reversed of the order in Pointwise. The reversal
                    // also must occur for surface quadrilaterals.
                    const int indices[3] = {
                        (int)data.index[2],
                        (int)data.index[1],
                        (int)data.index[0]
                    };
                    if (!OK(nc_put_vara_int(tau.idNc, tau.idTris, start, count,
                            indices))) {
                        ret = false;
                        break;
                    }
                    // compute the index into the boundarypanel variable and
                    // the marker variable, based on the offset for 
                    // triangles and the index of the current triangle
                    const size_t adjIdx = idxTri + tau.extrudedTriStart;
                    // save marker (Pointwise boundary-condition) information
                    if (!OK(nc_put_var1_int(tau.idNc, tau.numSurfBndryMarkers,
                            &adjIdx, (int*)&markerID))) {
                        ret = false;
                        break;
                    }
                    ++idxTri;
                    break; }
                case PWGM_ELEMTYPE_QUAD: {
                    // see comments for PWGM_ELEMTYPE_TRI:
                    const size_t start[2]{ idxQuad, 0 };
                    const size_t count[2]{ 1, 4 };
                    const int indices[4] = {
                        (int)data.index[3],
                        (int)data.index[2],
                        (int)data.index[1],
                        (int)data.index[0]
                    };
                    if (!OK(nc_put_vara_int(tau.idNc, tau.idQuads, start, count,
                            indices))) {
                        ret = false;
                        break;
                    }
                    const size_t adjIdx = idxQuad + tau.extrudedQuadStart;
                    if (!OK(nc_put_var1_int(tau.idNc, tau.numSurfBndryMarkers,
                            &adjIdx, (int*)&markerID))) {
                        ret = false;
                        break;
                    }
                    ++idxQuad;
                    break; }
                default:// probably a bar
                    break;
                }
                if (!caeuProgressIncr(&rti)) {
                    ret = false;
                    break;
                }
            }
        }
    }
    return caeuProgressEndStep(&rti) && ret;
}


/****************************************************************************
 * 
 * Close the netcdf file
 * 
 ***************************************************************************/
static PWP_BOOL
closeGridFile(const CAEP_RTITEM &rti)
{
    TAU_DATA &tau = *(rti.tau);
    const int idNc = tau.idNc;
    return OK(nc_close(idNc));
}


/****************************************************************************
 * 
 * Writes out the TAU boundary mapping file.
 * 
 ***************************************************************************/
static bool
writeBmapFile(const CAEP_RTITEM &rti)
{
    const TAU_DATA &tau = *(rti.tau);
    std::string file(rti.pWriteInfo->fileDest);
    file += ".bmap";
    FILE *fp = pwpFileOpen(file.c_str(), pwpWrite | pwpAscii);
    bool ret = (nullptr != fp);
    if (ret) {
        ret = ret && 0 < fprintf(fp,
                " -----------------------------------------------------\n");
        ret = ret && 0 < fprintf(fp, " BOUNDARY MAPPING\n");
        ret = ret && 0 < fprintf(fp,
                " -----------------------------------------------------\n");

        const StringUINT32Map &mNameToId = tau.nameToId;
        const StringStringMap &mNameToType = tau.nameToType;
        for (auto it = mNameToId.cbegin(); it != mNameToId.end() && ret; ++it) {
            const std::string &name = it->first;
            const PWP_UINT32 outputID = it->second;
            const std::string &type = mNameToType.at(name);
            ret = ret && 0 < fprintf(fp, "\n");
            ret = ret && 0 < fprintf(fp, "%26s %u\n", "Markers:", outputID);
            ret = ret && 0 < fprintf(fp, "%26s %s\n", "Type:", type.c_str());
            ret = ret && 0 < fprintf(fp, "%26s %s\n", "Name:", name.c_str());
            // future attributes added here
            ret = ret && 0 < fprintf(fp, "\n");
            ret = ret && 0 < fprintf(fp, " block end\n");
            ret = ret && 0 < fprintf(fp, " ---------------------------\n");
        }

        pwpFileClose(fp);
    }
    return ret;
}


/**************************************/
PWP_BOOL
runtimeWrite(CAEP_RTITEM *pRti, PWGM_HGRIDMODEL, const CAEP_WRITEINFO *)
{
    // attach runtime instance data to pRti
    TAU_DATA tauData;     // will persist for entire export
    pRti->tau = &tauData; // cppcheck-suppress autoVariables

    PWP_BOOL ret = PWP_TRUE;
    const PWP_UINT32 NumMajorSteps = (is2D(*pRti) ? 2 : 3);
    if (caeuProgressInit(pRti, NumMajorSteps)) {
        // These first four functions are the same for 2D and 3D.
        if (!setupMarkerInfo(*pRti)) {
            caeuSendErrorMsg(pRti, "Could not setup marker info!", 0);
            ret = false;
        }
        else if (!startGridFile(*pRti)) {
            caeuSendErrorMsg(pRti, "Could start grid file info!", 0);
            ret = false;
        }
        else if (!definePointVars(*pRti)) {
            caeuSendErrorMsg(pRti, "Could not define point vars!", 0);
            ret = false;
        }
        else if (!defineMarkers(*pRti)) {
            caeuSendErrorMsg(pRti, "Could not define markers!", 0);
            ret = false;
        }

        // The 2D and 3D exports now differ.
        if (is2D(*pRti)) {
            if (!defineElementVars2D(*pRti)) {
                caeuSendErrorMsg(pRti, "Could not define element vars!", 0);
                ret = false;
            }
            else if (!endFileDefinition(*pRti)) {
                caeuSendErrorMsg(pRti, "Could not create TAU file! If "
                    "writing a large file, ensure that the solver "
                    "attribute 'dataModel' is not set to classic.",
                    0);
                ret = false;
            }
            // first major progress step
            else if (!writePointVars2D(*pRti)) {
                caeuSendErrorMsg(pRti, "Could not write point vars!", 0);
                ret = false;
            }
            else if (!writeMarkers(*pRti)) {
                caeuSendErrorMsg(pRti, "Could not write markers!", 0);
                ret = false;
            }
            // second major progress step
            else if (!writeElementVars2D(*pRti)) {
                caeuSendErrorMsg(pRti, "Could not write element vars!", 0);
                ret = false;
            }
        }
        else {
            if (!defineElementVars(*pRti)) {
                caeuSendErrorMsg(pRti, "Could not define element vars!", 0);
                ret = false;
            }
            else if (!defineSurfaceElementVars(*pRti)) {
                caeuSendErrorMsg(pRti, "Could not define surface element "
                    "vars!", 0);
                ret = false;
            }
            else if (!endFileDefinition(*pRti)) {
                caeuSendErrorMsg(pRti, "Could not create TAU file! If "
                    "writing a large file, ensure that the solver "
                    "attribute 'dataModel' is not set to classic.",
                    0);
                ret = false;
            }
            // first major progress step
            else if (!writePointVars3D(*pRti)) {
                caeuSendErrorMsg(pRti, "Could not write point vars!", 0);
                ret = false;
            }
            // second major progress step
            else if (!writeElementVars(*pRti)) {
                caeuSendErrorMsg(pRti, "Could not write element vars!", 0);
                ret = false;
            }
            else if (!writeMarkers(*pRti)) {
                caeuSendErrorMsg(pRti, "Could not write markers!", 0);
                ret = false;
            }
            // third major step
            else if (!writeSurfaceElementVars(*pRti)) {
                caeuSendErrorMsg(pRti, "Could not write surface element "
                    "vars!", 0);
                ret = false;
            }
        }
        // Attempt to close file even if previous step failed.
        if (!closeGridFile(*pRti)) {
            ret = false;
        }
        else if (ret) {
            ret = writeBmapFile(*pRti);
        }
        caeuProgressEnd(pRti, ret);
    }
    return ret && !pRti->opAborted;
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
 * This file is licensed under the Cadence Public License Version 1.0 (the
 * "License"), a copy of which is found in the included file named "LICENSE",
 * and is distributed "AS IS." TO THE MAXIMUM EXTENT PERMITTED BY APPLICABLE
 * LAW, CADENCE DISCLAIMS ALL WARRANTIES AND IN NO EVENT SHALL BE LIABLE TO
 * ANY PARTY FOR ANY DAMAGES ARISING OUT OF OR RELATING TO USE OF THIS FILE.
 * Please see the License for the full text of applicable terms.
 *
 ****************************************************************************/
