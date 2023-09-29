#include <BasicAttributeComputer.hpp>
#include <algorithm>
#include <cmath>

// =========================================================
// MAIN CLASS
// =========================================================

BasicAttributeComputer::BasicAttributeComputer(
  const Box &domain, const std::vector<uint8> &f)
  : domain_{domain}, f_{f}
{ }

std::vector<BasicAttributes> 
  BasicAttributeComputer::computeAttribute(const MTree &tree)
{
  std::vector<AuxAttributes> auxAttrs(tree.numberOfNodes());
  std::vector<BasicAttributes> attrs(tree.numberOfNodes());

  AuxAttributesComputer auxAttrComputer{domain_, f_, auxAttrs, attrs};
  BasicAttributeIncrementalComputer basicAttrComputer{domain_, f_, auxAttrs, attrs};

  computeAuxAttributes(tree, auxAttrComputer);
  computeBasicAttributes(tree, basicAttrComputer);

  return attrs;
}

void BasicAttributeComputer::computeAuxAttributes(const MTree &tree,
  AuxAttributesComputer &auxAttrComputer)
{
  tree.tranverse([this,&auxAttrComputer](NodePtr node){
    auxAttrComputer.computeAtCNPs(node);

    if (node->parent() != nullptr) 
      auxAttrComputer.mergeToParent(node, node->parent());
  });
}

void BasicAttributeComputer::computeBasicAttributes(const MTree &tree,
  BasicAttributeIncrementalComputer &incrBasicAttrComputer)
{
  tree.tranverse([this, &incrBasicAttrComputer](NodePtr node){
    incrBasicAttrComputer.computeAtCNPs(node);

    if (node->parent() != nullptr)
      incrBasicAttrComputer.mergeToParent(node, node->parent());
  });
}

// =======================================================================
// AUXILIARY ATTRIBUTE COMPUTER 
// =======================================================================
BasicAttributeComputer::AuxAttributesComputer::AuxAttributesComputer(
  const Box &domain, 
  const std::vector<uint8> &f, 
  std::vector<AuxAttributes> &auxAttrs,
  std::vector<BasicAttributes> &attrs)
  : domain_{domain}, f_{f}, auxAttrs_{auxAttrs}, attrs_{attrs}
{}

void BasicAttributeComputer::AuxAttributesComputer::computeAtCNPs(
  NodePtr node)
{
  using morphotree::I32Point;

  AuxAttributes &nodeAuxAttrs = auxAttrs_[node->id()];
  BasicAttributes &nodeAttrs = attrs_[node->id()];

  nodeAttrs.area_          += node->cnps().size();
  nodeAttrs.volume_        += node->cnps().size() * node->level();
  nodeAttrs.level_          = node->level();
  nodeAttrs.varianceLevel_ += node->cnps().size() * (node->level() * node->level());

  // compute auxiliary attributes
  for (uint32 pidx : node->cnps()) {
    I32Point p = domain_.indexToPoint(pidx);

    nodeAuxAttrs.xmin = std::min(nodeAuxAttrs.xmin, p.x());
    nodeAuxAttrs.ymin = std::min(nodeAuxAttrs.ymin, p.y());
    nodeAuxAttrs.xmax = std::max(nodeAuxAttrs.xmax, p.x());
    nodeAuxAttrs.ymax = std::max(nodeAuxAttrs.ymax, p.y());

    nodeAuxAttrs.m10 += p.x();
    nodeAuxAttrs.m01 += p.y();
    nodeAuxAttrs.m11 += (p.x() + p.y());
    nodeAuxAttrs.m20 += (p.x()* p.x());
    nodeAuxAttrs.m02 += (p.y() * p.y());
  }

  // Attributes are finalised.
  nodeAttrs.meanLevel_      = nodeAttrs.volume_ / nodeAttrs.area_;
  nodeAttrs.width_          = nodeAuxAttrs.xmax - nodeAuxAttrs.xmin + 1;
  nodeAttrs.height_         = nodeAuxAttrs.ymax - nodeAuxAttrs.ymin + 1;
  nodeAttrs.rectangularity_ = nodeAttrs.area_ / (nodeAttrs.width_ * nodeAttrs.height_);
  nodeAttrs.ratioWH_        = std::max(nodeAttrs.width_, nodeAttrs.height_) / 
                                std::min(nodeAttrs.width_, nodeAttrs.height_);
}

void BasicAttributeComputer::AuxAttributesComputer::mergeToParent(
  NodePtr node, NodePtr parent)
{
  BasicAttributes &nodeAttrs = attrs_[node->id()];
  BasicAttributes &parentAttrs = attrs_[parent->id()];

  AuxAttributes &nodeAuxAttrs = auxAttrs_[node->id()];
  AuxAttributes &parentAuxAttrs  = auxAttrs_[parent->id()];

  parentAttrs.area_ += nodeAttrs.area_;
  parentAttrs.volume_ += nodeAttrs.volume_;
  parentAttrs.varianceLevel_ += nodeAttrs.varianceLevel_;

  parentAuxAttrs.xmin = std::min(parentAuxAttrs.xmin, nodeAuxAttrs.xmin);
  parentAuxAttrs.ymin = std::min(parentAuxAttrs.ymin, nodeAuxAttrs.ymin);
  parentAuxAttrs.xmax = std::max(parentAuxAttrs.xmax, nodeAuxAttrs.xmax);
  parentAuxAttrs.ymax = std::max(parentAuxAttrs.ymax, nodeAuxAttrs.ymax);

  parentAuxAttrs.m10 += nodeAuxAttrs.m10;
  parentAuxAttrs.m01 += nodeAuxAttrs.m01;
  parentAuxAttrs.m11 += nodeAuxAttrs.m11;
  parentAuxAttrs.m20 += nodeAuxAttrs.m20;
  parentAuxAttrs.m02 += nodeAuxAttrs.m02;
}

// =======================================================================
// INCREMENTAL BASIC ATTRIBUTES COMPUTER 
// =======================================================================
BasicAttributeComputer::BasicAttributeIncrementalComputer
  ::BasicAttributeIncrementalComputer(const Box &domain,
    const std::vector<uint8> &f, const std::vector<AuxAttributes> &auxAttrs,
    std::vector<BasicAttributes> &attrs)
    : domain_{domain}, f_{f}, auxAttrs_{auxAttrs}, attrs_{attrs}
{}

void BasicAttributeComputer::BasicAttributeIncrementalComputer
  ::computeAtCNPs(NodePtr node)
{
  using morphotree::I32Point;

  BasicAttributes      &nodeAttrs    = attrs_[node->id()];
  const AuxAttributes  &nodeAuxAttrs = auxAttrs_[node->id()];

  I32Point centroid{ nodeAuxAttrs.m10 / nodeAttrs.area_, 
    nodeAuxAttrs.m01 / nodeAttrs.area_ };
  
  nodeAttrs.varianceLevel_ += 
    (((nodeAttrs.meanLevel_ - node->level()) * (nodeAttrs.meanLevel_ - node->level()))
    * node->cnps().size());

  for (uint32 pidx : node->cnps()) {
    I32Point p = domain_.indexToPoint(pidx);
    nodeAttrs.moment20_ += ((p.x() - centroid.x()) * (p.x() - centroid.x()));
    nodeAttrs.moment02_ += ((p.y() - centroid.y()) * (p.y() - centroid.y()));
    nodeAttrs.moment11_ += ((p.x() - centroid.x()) * (p.y() - centroid.y()));
  }

  // Attributes are finalised.
  float sumDiffLevel2 = static_cast<float>(nodeAttrs.varianceLevel_);
  float area = static_cast<float>(nodeAttrs.area_);
  nodeAttrs.varianceLevel_ = sumDiffLevel2 / area;

  float moment20 = nodeAttrs.moment20_;
  float moment02 = nodeAttrs.moment02_;
  float moment11 = nodeAttrs.moment11_;

  if (area > 1) {
    float a1 = moment20 + moment02 + 
      std::sqrt(std::pow(moment20 - moment02, 2) + 4 * std::pow(moment11, 2));
    float a2 = moment20 + moment02 - 
      std::sqrt(std::pow(moment20 - moment02, 2) + 4 * std::pow(moment11, 2));
    nodeAttrs.inertia_ = normMoment(area, moment02, 0, 2) 
      + normMoment(area, moment20, 2, 0);  // hu = inertia
    nodeAttrs.orientation_ = 0.5f * std::atan2(2 * moment11, moment20 - moment02); // orientation

    nodeAttrs.lenMajorAxis_ = std::sqrt( (2 * a1) / area ); // length major axis
    nodeAttrs.lenMinorAxis_ = std::sqrt( (2 * a2) / area ); // length minor axis
    nodeAttrs.eccentricity_ = (a2 != 0 ? a1 / a2 : a1 / 0.1); // eccentricity
    nodeAttrs.compactness_ = (1.0f / (2.0f*M_PI)) * (area / (moment20 + moment02)); // compactness
  }
}

float BasicAttributeComputer::BasicAttributeIncrementalComputer
  ::normMoment(float area, float moment, int p, int q) const
{
  return moment / std::pow( area, (p + q + 2.0) / 2.0);
}
