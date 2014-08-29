#error Documentation only

//------------------------------------------------------------------------------
/** @defgroup tools Small tools for mesh generation, file conversion, etc
 *
 *  A random collection of tools which are not directly related to the FE
 *  computations for which InSilico is designed. Most of these tools are
 *  related to pre-processing (mesh generation) and post-processing
 *  (visulisation of data). 
 */

//------------------------------------------------------------------------------
/** \ingroup   tools
 *  \namespace tools
 *  Tools for mesh generation, file conversion and post processing.
 */


//------------------------------------------------------------------------------
/** \ingroup   tools
 *  \namespace tools::analysis
 *  File analysis tools.
 */

/** \ingroup   tools
 *  \namespace tools::analysis::bbox
 *  Computation of the bounding box of a mesh.
 */

/** \ingroup   tools
 *  \namespace tools::analysis::enclosed
 *  Computation of the enclosed volume of a surface.
 */

//------------------------------------------------------------------------------
/** \ingroup   tools
 *  \namespace tools::converter
 *  File conversion tools.
 */

/** \ingroup   tools
 *  \namespace tools::converter::gmsh2smf
 *  Convert msh files from GMSH to smf files.
 */

/** \ingroup   tools
 *  \namespace tools::converter::sgf2vtk
 *  Convert sgf files to vtk files.
 */

/** \ingroup   tools
 *  \namespace tools::converter::sgf2xx
 *  Generic conversion tools from sgf to something.
 */

/** \ingroup   tools
 *  \namespace tools::converter::sgfRefine
 *  Refine a given sgf file.
 */

/** \ingroup   tools
 *  \namespace tools::converter::smf2smf
 *  Read smf and write to smf.
 */

/** \ingroup   tools
 *  \namespace tools::converter::smf2gp
 *  Read smf and write to a gnuplot data file.
 */

/** \ingroup   tools
 *  \namespace tools::converter::smf2vtk
 *  Read smf and write to vtk file.
 */

/** \ingroup   tools
 *  \namespace tools::converter::smf2xx
 *  Generic conversion tools from smf to something.
 */

/** \ingroup   tools
 *  \namespace tools::converter::smfMap
 *  Read smf, apply a coordinate map, write smf.
 */

/** \ingroup   tools
 *  \namespace tools::converter::smfRandom
 *  Read smf, apply a randomised coordinate shift, write smf.
 */

/** \ingroup   tools
 *  \namespace tools::converter::smfReverse
 *  Read smf, reverse the element orientation, write smf.
 */

//------------------------------------------------------------------------------
/** \ingroup   tools
 *  \namespace tools::meshGeneration
 *  Simple mesh generators.
 */

/** \ingroup   tools
 *  \namespace tools::meshGeneration::boundaries2D
 *  Generate surface meshes in 2D.
 */

/** \ingroup   tools
 *  \namespace tools::meshGeneration::boundaries3D
 *  Generate surface meshes in 3D.
 */

/** \ingroup   tools
 *  \namespace tools::meshGeneration::circle
 *  Generate a volume mesh in a circle in 2D.
 */

/** \ingroup   tools
 *  \namespace tools::meshGeneration::unitCube
 *  Make a uniform mesh in \f$ [0,1]^d \f$
 */


