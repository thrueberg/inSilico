# error DOCUMENATION_ONLY
/**  \namespace base::io::smf
 *    Namespace for reading and writing files in SMF format.
 *
 *    SMF Format Specification
 *    ========================
 *
 *
 *    Description
 *    -----------
 *
 *    The *smf* format (Simple Mesh Format) is used to describe a finite element
 *    mesh. Basically, it contains of coordinates and the connectivity of
 *    elements.
 *    A standard *smf* file has the following outline
 *    \code{.txt}
 *    ! elementShape      shape
 *    #  A useless comment
 *    ! elementNumPoints  numNodes
 *    nNodes  nElements
 *    x0  y0  z0
 *    x1  y1  z1
 *    ...
 *    v11 v12 v13 ... v1N
 *    v21 v22 v23 ... v2N
 *    ...
 *    \endcode
 *    Here the first two lines are mandatory (in any order) as they described the
 *    expected type of element. The other parts of the above file are as follows:
 *    -   `!` marks the begin of a header line, all header lines must come before
 *        the begin of number lines
 *    -   `#` marks the begin of a comment line. Comment lines can appear anywhere
 *        in the header
 *    -   `shape` is the literal associated with base::Shape, see also
 *        base::ShapeName for the existing shape names
 *    -   `numNodes` is the number of nodes per element.
 *    -   `nNodes`   is the total number of nodes in the mesh
 *    -   `nElements` is the total number of elements in the mesh
 *    -   `xi`, `yi`, and `zi` are the x, y, and z-coordinates of node i
 *        (_Note:_  coordinates are always given in *3D*, independent of the
 *         dimension of the problem)
 *    -   `vi1` to  `viN` are the node indices of the element i, `N` obviously
 *        being equal to `numNodes` (_Note:_ the node indices are 0-based).
 *
 *    In addition to the standard specification above, there is an extension which
 *    allows to include external node and/or element files. This is accomplished
 *    by using one or both of the following lines
 *    \code{.txt}
 *    ! externalNodes    nodeFile
 *    ! externalElements elementFile
 *    \endcode
 *    in the header of an smf file. Obviously, `nodeFile` is the name of a file
 *    containing the coordinates and `elementFile` consequently contains the
 *    connectivity. As an exmple consider the following 
 *    \code{.txt}
 *    ! externalNodes    nodes.xyz
 *    nNodes nElements
 *    v11 v12 v13 ... v1N
 *    v21 v22 v23 ... v2N
 *    ...
 *    \endcode
 *    Here, `nodes.xyz` is a text file with node coordinates only and the
 *    connectivity is given inline as in the standard case.
 *    _Note_ that `nNodes` and `nElements`  _always_ appear in the main
 *    smf file are not repeated in an external file.
 *
 *    ASCII-art element types
 *    -----------------------
 *
 *    Find in the following pictures of the currently supported element types.
 *    The notation is `<shape,numNodes>`, where shape is of type base::Shape and
 *    numNodes is the number of nodes of the element.
 *
 *    \code{.txt}
 *
 *     <LINE,2>                <LINE,3>
 *
 *     0----------1 --> x      0----2----1
 *
 *
 *     <TRI,3>                 <TRI,6>
 *
 *     y
 *     ^
 *     |
 *     2                       2
 *     |`\                     |`\
 *     |  `\                   |  `\
 *     |    `\                 5    `4
 *     |      `\               |      `\
 *     |        `\             |        `\
 *     0----------1 --> x      0-----3----1
 *
 *
 *    <QUAD,4>                      <QUAD,9>
 *
 *     y
 *     ^
 *     |
 *     3-----------2                 3-----6-----2
 *     |           |                 |           |
 *     |           |                 |           |
 *     |           |                 7     8     5
 *     |           |                 |           |
 *     |           |                 |           |
 *     0-----------1 --> x           0-----4-----1
 *
 *
 *
 *     <TET,4>                                <TET,10>
 *
 *                       z
 *                     .
 *                   ,/
 *                  /
 *                3                                     3
 *              ,/|`\                                 ,/|`\
 *            ,/  |  `\                             ,/  |  `\
 *          ,/    '.   `\                         ,7    '.   `9
 *        ,/       |     `\                     ,/       8     `\
 *      ,/         |       `\                 ,/         |       `\
 *     0-----------'.--------2 --> y         0--------6--'.--------2
 *      `\.         |      ,/                 `\.         |      ,/
 *         `\.      |    ,/                      `\.      |    ,5
 *            `\.   '. ,/                           `4.   '. ,/
 *               `\. |/                                `\. |/
 *                  `1                                    `1
 *                     `\.
 *                        ` x
 *
 *
 *     <HEX,8>                            <HEX,27>
 *
 *     z
 *     ^
 *     |
 *     4----------7                       4----15----7
 *     |\         |\                      |\         |\
 *     | \        | \                     |12    21  | 14
 *     |  \       |  \                   16  \ 25    19 \
 *     |   5----------6                   |   5----13+---6
 *     |   |      |   |                   |22 |  26  | 24|
 *     0---|------3---|--> y              0---+-11---3   |
 *      \  |       \  |                    \ 17    23 \  18
 *       \ |        \ |                     8 |  20    10|
 *        \|         \|                      \|         \|
 *         1----------2                       1-----9----2
 *          \
 *           \
 *           `'
 *             x
 *
 *   \endcode 
 *
 *   Examples
 *   -------
 *
 *   A simple square \f$ (0,1)^2 \f$ composed of 4 equal squares.
 *    \include quadrat.smf
 *
 */
