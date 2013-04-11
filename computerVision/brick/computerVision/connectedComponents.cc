/**
***************************************************************************
* @file brick/computerVision/connectedComponents.cc
*
* Source file defining symbols declared in connectedComponents.hh.
*
* Copyright (C) 2006,2012-2013 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/


#include <brick/computerVision/connectedComponents.hh>

namespace brick {

  namespace computerVision {

    /// @cond privateCode
    namespace privateCode {

      // This function traverses the blob label graph, starting at a
      // particular node, propagating that node's label so that any blobs
      // which are connected to the starting node inherit the the starting
      // nodes label iff the starting nodes label is less than the current
      // label of the connected blob.
      bool
      propagateLabel(size_t label, size_t node, std::vector<size_t>& labelArray,
                     const std::vector< std::list<size_t> >& neighborsVector)
      {
        // If we already know that the label of this node is not the
        // lowest among all of the nodes it's connected to, then stop
        // now.
        if(labelArray[node] <= label) {
          return false;
        }
        labelArray[node] = label;

        // Make a list of neighbors we need to check.
        std::list<size_t> nodeList;
        std::copy(neighborsVector[node].begin(),
                  neighborsVector[node].end(),
                  std::back_inserter(nodeList));

        // For each neighbor...
        std::list<size_t>::iterator nextNodeIterator = nodeList.begin();
        while(nextNodeIterator != nodeList.end()) {
          // Don't update labels that are already lower than the current label.
          if(labelArray[*nextNodeIterator] <= label) {
            ++nextNodeIterator;
            continue;
          }
          // Propagate the label and add any new neighbors (nodes that
          // are connected to the node to which we're propagating) to our
          // list.
          labelArray[*nextNodeIterator] = label;
          std::copy(neighborsVector[*nextNodeIterator].begin(),
                    neighborsVector[*nextNodeIterator].end(),
                    std::back_inserter(nodeList));
          ++nextNodeIterator;
        }
        return true;
      }

    } // namespace privateCode
    /// @endcond

  } // namespace computerVision
    
} // namespace brick
