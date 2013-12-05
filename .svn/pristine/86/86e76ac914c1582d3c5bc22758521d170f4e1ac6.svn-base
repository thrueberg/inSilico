# error DOCUMENATION_ONLY
/**  \namespace base::io::sgf
 *    Namespace for reading and writing files in SGF format.
 *
 *    SGF Format Specification
 *    ========================
 *
 *
 *    Description
 *    -----------
 *
 *    The *sgf* format (Simple Grid Format) is used to describe a structured
 *    finite element mesh. It contains of coordinates in a lexicographic order
 *    which directly imply the connectivity of the mesh.
 *    A standard *sgf* file has the following outline
 *    \code{.txt}
 *    # Any number of comment lines which can only
 *    # appear at the beginning of the file and must start with
 *    # the #-character
 *    nElementsX nElementsY nElementsZ
 *    x0  y0  z0
 *    x1  y1  z1
 *    ...
 *    \endcode
 *    The important detail is that the nodes are ordered correctly.
 *    Consider a 2D example of 3x2 uniform grid of the domain (0,3)x(0,2)
 *    \code{.txt}
 *    3  2  0
 *    0. 0. 0.
 *    1. 0. 0.
 *    2. 0. 0.
 *    3. 0. 0.
 *    0. 1. 0.
 *    1. 1. 0.
 *    2. 1. 0.
 *    3. 1. 0.
 *    0. 2. 0.
 *    1. 2. 0.
 *    2. 2. 0.
 *    3. 2. 0.
 *    \endcode
 *
 *    Notes:
 *    -  The first line indicates that we have a 3 x 2 grid in two dimensions,
 *       therefore the last value is 0
 *    -  The grid dimensions imply (3+1) x (2+1) x (0+1) = 4 x 3 x 1 = 12
 *       grid nodes which follow after the head line
 *    -  The nodes are ordered such that the x-coordinate changes first, then the
 *       y-coordinate, then the z-coordinate
 *
 *
 *    ASCII-art element types
 *    -----------------------
 *
 *    A structured grid can have the following element types in the spatial
 *    dimension <DIM>
 *
 *    \code{.txt}
 *
 *     <1>              
 *
 *     0----------1 --> x  
 *
 *
 *    <2>
 *
 *     y
 *     ^
 *     |
 *     2-----------3          
 *     |           |          
 *     |           |          
 *     |           |          
 *     |           |          
 *     |           |          
 *     0-----------1 --> x    
 *
 *
 *     <3> 
 *
 *     z
 *     ^
 *     |
 *     4----------6             
 *     |\         |\            
 *     | \        | \           
 *     |  \       |  \          
 *     |   5----------7         
 *     |   |      |   |         
 *     0---|------2---|--> y    
 *      \  |       \  |         
 *       \ |        \ |         
 *        \|         \|         
 *         1----------3         
 *          \
 *           \
 *           `'
 *             x
 *
 *   \endcode 
 *
 *   Note the lexicographic ordering of the vertex indices as opposed to the
 *   hierarchic ordering in the case of unstructured (i.e. smf format,
 *   see base::io::smf) mesh.
 */
