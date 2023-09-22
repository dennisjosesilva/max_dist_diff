#pragma once 

#include <morphotree/tree/mtree.hpp>
#include <morphotree/core/box.hpp>

#include <gft.h>

#include <array>

class UltimateAttributeOpening{

    
  public:
    using uint8 = morphotree::uint8;
    using uint32 = morphotree::uint32;
    using MTree = morphotree::MorphologicalTree<uint8>;
    using NodePtr = MTree::NodePtr;
        
    UltimateAttributeOpening(const MTree &ptree, int maxCriterion, std::vector<uint32> attrs_increasing);

    ~UltimateAttributeOpening();

    std::vector<uint8> createRandomColor(std::vector<uint32> img);

    std::vector<uint8> getMaxConstrastImage();

    std::vector<uint32> getAssociatedImage();

    std::vector<uint8> getAssociatedColorImage();
    
    void computeUAO(NodePtr currentNode, int levelNodeNotInNR, bool qPropag);

  private:
    int maxCriterion;
    std::vector<uint32> attrs_increasing;
    MTree maxtree;
    int* maxContrastLUT;
    int* associatedIndexLUT;
};
 