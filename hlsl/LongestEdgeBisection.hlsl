// data-structures
struct leb_DiamondParent {
    cbt_Node base, top;
};
leb_DiamondParent leb_DecodeDiamondParent(in const cbt_Node node);

// conforming split / merge
#ifdef CBT_FLAG_WRITE
void leb_Split(in const cbt_Node node);
void leb_Merge(in const cbt_Node node);
#endif

// Subdivision routine O(depth)
float3   leb_DecodeAttributeArray(in const cbt_Node node, in const float3 data);
float3x2 leb_DecodeAttributeArray(in const cbt_Node node, in const float3x2 data);
float3x3 leb_DecodeAttributeArray(in const cbt_Node node, in const float3x3 data);
float3x4 leb_DecodeAttributeArray(in const cbt_Node node, in const float3x4 data);


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

struct leb__SameDepthNeighborIDs {
    uint left, right, edge, node;
};

leb__SameDepthNeighborIDs
leb__CreateSameDepthNeighborIDs(uint left, uint right, uint edge, uint node)
{
    leb__SameDepthNeighborIDs neighborIDs;

    neighborIDs.left = left;
    neighborIDs.right = right;
    neighborIDs.edge = edge;
    neighborIDs.node = node;

    return neighborIDs;
}

leb_DiamondParent
leb__CreateDiamondParent(in const cbt_Node base, in const cbt_Node top)
{
    leb_DiamondParent diamond;

    diamond.base = base;
    diamond.top = top;

    return diamond;
}

/*******************************************************************************
 * GetBitValue -- Returns the value of a bit stored in a 32-bit word
 *
 */
uint leb__GetBitValue(uint bitField, int bitID)
{
    return ((bitField >> bitID) & 1u);
}


/*******************************************************************************
 * SplitNodeIDs -- Updates the IDs of neighbors after one LEB split
 *
 * This code applies the following rules:
 * Split left:
 * LeftID  = 2 * NodeID + 1
 * RightID = 2 * EdgeID + 1
 * EdgeID  = 2 * RightID + 1
 *
 * Split right:
 * LeftID  = 2 * EdgeID
 * RightID = 2 * NodeID
 * EdgeID  = 2 * LeftID
 *
 * The _reserved channel stores NodeID, which is recquired for applying the
 * rules.
 *
 */
leb__SameDepthNeighborIDs 
leb__SplitNeighborIDs(in const leb__SameDepthNeighborIDs  nodeIDs, uint splitBit)
{
    uint b = splitBit;
    uint c = splitBit ^ 1u;
    bool cb = bool(c);
    uint4 idArray = uint4(nodeIDs.left, nodeIDs.right, nodeIDs.edge, nodeIDs.node);

    return leb__CreateSameDepthNeighborIDs(
        (idArray[2 + b] << 1) | uint(cb && bool(idArray[2 + b])),
        (idArray[2 + c] << 1) | uint(cb && bool(idArray[2 + c])),
        (idArray[b    ] << 1) | uint(cb && bool(idArray[b    ])),
        (idArray[3    ] << 1) | b
    );
}

/*******************************************************************************
 * DecodeNodeNeighborIDs -- Decodes the IDs of the cbt_Nodes neighbour to node
 *
 * The IDs are associated to the depth of the input node. As such, they
 * don't necessarily exist in the LEB subdivision.
 *
 */
leb__SameDepthNeighborIDs leb__DecodeSameDepthNeighborIDs(in const cbt_Node node)
{
    int bitID = node.depth > 0 ? node.depth - 1 : 0;
    uint b = leb__GetBitValue(node.id, bitID);
    leb__SameDepthNeighborIDs nodeIDs =
        leb__CreateSameDepthNeighborIDs(0u, 0u, 3u - b, 2u + b);

    for (bitID = node.depth - 2; bitID >= 0; --bitID) {
        nodeIDs = leb__SplitNeighborIDs(nodeIDs, leb__GetBitValue(node.id, bitID));
    }

    return nodeIDs;
}


/*******************************************************************************
 * EdgeNode -- Computes the neighbour of the input node wrt to its longest edge
 *
 */
cbt_Node leb__EdgeNeighbor(in const cbt_Node node)
{
    uint nodeID = leb__DecodeSameDepthNeighborIDs(node).edge;

    return cbt_CreateNode(nodeID, (nodeID == 0u) ? 0 : node.depth);
}


/*******************************************************************************
 * SplitNodeConforming -- Splits a node while producing a conforming LEB
 *
 */
#ifdef CBT_FLAG_WRITE
void leb_SplitNode(in const cbt_Node node)
{
    if (!cbt_IsCeilNode(node)) {
        const uint minNodeID = 1u;
        cbt_Node nodeIterator = node;

        cbt_SplitNode_Fast(nodeIterator);
        nodeIterator = leb__EdgeNeighbor(nodeIterator);

        while (nodeIterator.id > minNodeID) {
            cbt_SplitNode_Fast(nodeIterator);
            nodeIterator = cbt_ParentNode_Fast(nodeIterator);
            
            if (nodeIterator.id > minNodeID) {
                cbt_SplitNode_Fast(nodeIterator);
                nodeIterator = leb__EdgeNeighbor(nodeIterator);
            }
        }
    }
}
#endif // CBT_FLAG_WRITE


/*******************************************************************************
 * DecodeDiamondParent -- Decodes the diamond associated to the Node
 *
 * If the neighbour part does not exist, the parentNode is copied instead.
 *
 */
leb_DiamondParent leb_DecodeDiamondParent(in const cbt_Node node)
{
    cbt_Node parentNode = cbt_ParentNode_Fast(node);
    uint edgeNeighborID = leb__DecodeSameDepthNeighborIDs(parentNode).edge;
    cbt_Node edgeNeighborNode = cbt_CreateNode(
        edgeNeighborID > 0u ? edgeNeighborID : parentNode.id,
        parentNode.depth
    );

    return leb__CreateDiamondParent(parentNode, edgeNeighborNode);
}


/*******************************************************************************
 * HasDiamondParent -- Determines whether a diamond parent is actually stored
 *
 * This procedure checks that the diamond parent is encoded in the CBT.
 * We can perform this test by checking that both the base and top nodes
 * that form the diamond parent are split, i.e., CBT[base] = CBT[top] = 2.
 * This is a crucial operation for implementing the leb_Merge routine.
 *
 */
bool leb__HasDiamondParent(in const leb_DiamondParent diamondParent)
{
    bool canMergeBase = cbt_HeapRead(diamondParent.base) <= 2u;
    bool canMergeTop = cbt_HeapRead(diamondParent.top) <= 2u;

    return canMergeBase && canMergeTop;
}


/*******************************************************************************
 * MergeNode -- Merges a node while producing a conforming LEB
 *
 * This routines makes sure that the children of a diamond (including the
 * input node) all exist in the LEB before calling a merge.
 *
 */
#ifdef CBT_FLAG_WRITE
void
leb_MergeNode(in const cbt_Node node, in const leb_DiamondParent diamondParent)
{
    if ((node.depth > 1) && leb__HasDiamondParent(diamondParent)) {
        cbt_MergeNode_Fast(node);
        // FIXME
        cbt_MergeNode_Fast(cbt_RightChildNode_Fast(diamondParent.top));
    }
}
#endif // CBT_FLAG_WRITE


/*******************************************************************************
 * SplitMatrix3x3 -- Computes a LEB splitting matrix from a split bit
 *
 */
float3x3 leb__SplittingMatrix(uint bitValue)
{
    float b = float(bitValue);
    float c = 1.0f - b;

    return float3x3(
        c   , b   , 0.0f,
        0.5f, 0.0f, 0.5f,
        0.0f, c   , b
    );
}


/*******************************************************************************
 * SquareMatrix3x3 -- Computes the matrix that affects the triangle to the square
 *
 */
float3x3 leb__SquareMatrix(uint bitValue)
{
    float b = float(bitValue);
    float c = 1.0f - b;

    return float3x3(
        c, 0.0f, b,
        b,    c, b,
        b, 0.0f, c
    );
}


/*******************************************************************************
 * WindingMatrix -- Computes the matrix that garantees that triangles have same winding
 *
 */
float3x3 leb__WindingMatrix(uint bitValue)
{
    float b = float(bitValue);
    float c = 1.0f - b;

    return float3x3(
        c, 0.0f, b,
        0, 1.0f, 0,
        b, 0.0f, c
    );
}


/*******************************************************************************
 * DecodeTransformationMatrix -- Computes the splitting matrix associated to a LEB
 * node
 *
 */
float3x3 leb__DecodeTransformationMatrix(in const cbt_Node node)
{
    int bitID = node.depth > 0 ? node.depth - 1 : 0;
    float3x3 xf = leb__SquareMatrix(leb__GetBitValue(node.id, bitID));

    for (bitID = node.depth - 2; bitID >= 0; --bitID) {
        xf = mul(leb__SplittingMatrix(leb__GetBitValue(node.id, bitID)), xf);
    }

    return mul(leb__WindingMatrix((node.depth & 1) ^ 1), xf);
}


/*******************************************************************************
 * DecodeNodeAttributeArray -- Compute the triangle attributes at the input node
 *
 */
float3 leb_DecodeAttributeArray(in const cbt_Node node, in const float3 data)
{
    return mul(leb__DecodeTransformationMatrix(node), data);
}

float3x2 leb_DecodeAttributeArray(in const cbt_Node node, in const float3x2 data)
{
    return mul(leb__DecodeTransformationMatrix(node), data);
}

float3x3 leb_DecodeAttributeArray(in const cbt_Node node, in const float3x3 data)
{
    return mul(leb__DecodeTransformationMatrix(node), data);
}

float3x4 leb_DecodeAttributeArray(in const cbt_Node node, in const float3x4 data)
{
    return mul(leb__DecodeTransformationMatrix(node), data);
}
