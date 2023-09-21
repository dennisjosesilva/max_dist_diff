#include <morphotree/core/box.hpp>
#include <morphotree/core/alias.hpp>
#include <morphotree/tree/ct_builder.hpp>
#include <morphotree/tree/mtree.hpp>
#include <morphotree/adjacency/adjacency4c.hpp>
#include <morphotree/adjacency/adjacency8c.hpp>
#include <morphotree/attributes/extinctionValues/ExtinctionValueLeavesComputer.hpp>
#include <morphotree/attributes/areaComputer.hpp>
#include <morphotree/filtering/extinctionFilter.hpp>

#include <iostream>

#include <MaxDistComputer.hpp>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#define APPDEBUG

int main(int argc, char *argv[])
{
  // import morphotree types
  using morphotree::uint8;
  using morphotree::uint32;
  using morphotree::Box;
  using morphotree::UI32Point;
  using morphotree::MorphologicalTree;
  using morphotree::Adjacency4C;
  using morphotree::Adjacency8C;
  using morphotree::Adjacency;
  using morphotree::buildMaxTree;
  using morphotree::AreaComputer;
  using morphotree::ExtinctionValueLeavesComputer;
  using NodePtr = MorphologicalTree<uint8>::NodePtr;
  using ExtinctionValueComputer = ExtinctionValueLeavesComputer<uint8, uint32>;
  using ExtinctionValueMapType = typename ExtinctionValueComputer::MapType;
  using morphotree::extinctionFilter;
  using morphotree::iextinctionFilter;

  // check number of arguments from the command call
  if (argc < 3) {
    std::cerr << "usage error!\n";
    std::cerr << "usage: max_dist_extiction_values <image> <out_img>\n";
    return -1;
  }

  // read image
  int width, height, nchannels;
  uint8 *data = stbi_load(argv[1], &width, &height, &nchannels, 1);
  //stbi_image_free(data);

  // convert to morphotree image format
  Box domain = Box::fromSize({static_cast<uint32>(width), static_cast<uint32>(height)});
  std::vector<uint8> f(data, data + domain.numberOfPoints());

  // build max-tree
  std::shared_ptr<Adjacency> adj = std::make_shared<Adjacency8C>(domain);
  MorphologicalTree<uint8> maxtree = buildMaxTree(f, adj);

  std::cout << "Number of nodes: " << maxtree.numberOfNodes() << std::endl;

  // Extract leaves
  std::vector<NodePtr> leaves;
  maxtree.tranverse([&leaves](NodePtr node){
    if (node->children().empty()) // node is leaf
      leaves.push_back(node);
  });

  // Exctiction value
  std::vector<uint32> maxDist = computeMaxDistanceAttribute(domain, f, maxtree);
  
  #ifdef APPDEBUG
    //print exctition value 
    ExtinctionValueMapType extVals = ExtinctionValueComputer().compute(maxtree, maxDist);
    for (auto& leafVal : extVals) {
      std::cout << "maxDist[" << leafVal.first << "] = " << leafVal.second << "\n";
    }
  #endif

  // perform extinction filter
  iextinctionFilter(maxtree, maxDist, 15);

  std::cout << "FILTER\n";
  std::cout << "Number of nodes: " << maxtree.numberOfNodes() << std::endl;

  // record filtered image into the disk
  std::vector<uint8> out = maxtree.reconstructImage();
  stbi_write_png(argv[2], domain.width(), domain.height(), 1, out.data(), 0);
  stbi_image_free(data);

  return 0;
}