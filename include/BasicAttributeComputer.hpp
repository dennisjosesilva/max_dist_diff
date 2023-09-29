#include <vector>
#include <morphotree/tree/mtree.hpp>
#include <morphotree/core/box.hpp>
#include <morphotree/attributes/attributeComputer.hpp>


class BasicAttributes
{
public:
  float area() const { return area_; }
  float volume() const { return volume_; }
  float level() const { return level_; }
  float meanLevel() const { return meanLevel_; }
  float varianceLevel() const { return varianceLevel_; }
  float width() const { return width_; }
  float height() const { return height_; }
  float rectangularity() const { return rectangularity_; }
  float ratioWH() const { return ratioWH_; }
  float moment02() const { return moment02_; }
  float moment20() const { return moment20_; }
  float moment11() const { return moment11_; }
  float inertia() const { return inertia_; }
  float orientation() const { return orientation_; }
  float lenMajorAxis() const { return lenMajorAxis_; }
  float lenMinorAxis() const { return lenMinorAxis_; }
  float eccentricity() const { return eccentricity_; }
  float compactness() const { return compactness_; }

private:
  float area_;
  float volume_;
  float level_;
  float meanLevel_;
  float varianceLevel_;
  float width_;
  float height_;
  float rectangularity_;
  float ratioWH_;
  float moment02_;
  float moment20_;
  float moment11_;
  float inertia_;
  float orientation_;
  float lenMajorAxis_;
  float lenMinorAxis_;
  float eccentricity_;
  float compactness_;

  friend class BasicAttributeComputer;
};


class BasicAttributeComputer
{
public:
  using uint8 = morphotree::uint8;
  using uint32 = morphotree::uint32;
  using Box = morphotree::Box;
  using MTree = morphotree::MorphologicalTree<uint8>;
  using NodePtr = typename MTree::NodePtr;
  using AttributeComputer = morphotree::AttributeComputer<BasicAttributes, uint8>;

  BasicAttributeComputer(const Box &domain, const std::vector<uint8> &f);

  std::vector<BasicAttributes> computeAttribute(const MTree &tree);

private:
  const Box &domain_;
  const std::vector<uint8> &f_;

private:
  struct AuxAttributes
  {
    int xmax;
    int ymax;
    int xmin;
    int ymin;
    int m10;
    int m01;
    int m20;
    int m02;
    int m11;
  }; 

  class AuxAttributesComputer 
  {
  public:
    AuxAttributesComputer(const Box &domain, const std::vector<uint8> &f,
      std::vector<AuxAttributes> &auxAttrs, std::vector<BasicAttributes> &attrs);

    void computeAtCNPs(NodePtr node);
    void mergeToParent(NodePtr node, NodePtr parent);

  private:
    const std::vector<uint8> &f_;
    const Box &domain_;
    std::vector<AuxAttributes> &auxAttrs_;
    std::vector<BasicAttributes> &attrs_;
  };

  class BasicAttributeIncrementalComputer
  {
  public:
    BasicAttributeIncrementalComputer(const Box &domain, const std::vector<uint8> &f,
      const std::vector<AuxAttributes> &auxAttrs, std::vector<BasicAttributes> &attrs);

    void computeAtCNPs(NodePtr node);
    void mergeToParent(NodePtr node, NodePtr parent);

  private:
    float normMoment(float area, float moment, int p, int q) const;

  private:
    const std::vector<uint8> &f_;
    const Box &domain_;
    const std::vector<AuxAttributes> &auxAttrs_;
    std::vector<BasicAttributes> &attrs_;
  };

private:   
  void computeAuxAttributes(const MTree &tree,
    AuxAttributesComputer &auxAttrComputer);
  
  void computeBasicAttributes(const MTree &tree,
    BasicAttributeIncrementalComputer &incrBasicAttrComputer);
};



std::vector<BasicAttributes> computeBasicAttributes(
  const BasicAttributeComputer::Box &domain,
  const std::vector<BasicAttributeComputer::uint8> &f, 
  const BasicAttributeComputer::MTree &tree);
