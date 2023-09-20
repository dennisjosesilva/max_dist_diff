#include "MaxDistComputer.hpp"
#include <unordered_set>

#include <morphotree/adjacency/adjacency.hpp>
#include <morphotree/adjacency/adjacency4c.hpp>

#include <edt_diff.hpp>

#include <iostream>

//#define APPDEBUG

#ifdef APPDEBUG
  #include <sstream>
#endif

// Useful function
std::vector<morphotree::uint32> computeMaxDistanceAttribute(
  const morphotree::Box &domain,
  const std::vector<morphotree::uint8> &f, 
  const MaxDistComputer::MTree &tree)
{
  MaxDistComputer comp{domain, f};
  return comp.computeAttribute(tree);
}


MaxDistComputer::MaxDistComputer(Box domain,
  const std::vector<uint8> &f)
  : domain_{domain}, f_{f}
{}

std::vector<MaxDistComputer::uint32> MaxDistComputer::computeAttribute(
  const MTree &tree) const
{
  using NodePtr = MTree::NodePtr;
  using morphotree::InfAdjacency4C;

  // define useful images 
  std::vector<uint32> maxDist(tree.numberOfNodes(), 0); 
  std::array<std::vector<NodePtr>, 256> levelToNodes = extractLevelMap(tree);
  gft::sImage32 *gftImg = createGFTImage();
  gft::sImage32 *bin = gft::Image32::Create(gftImg);
  gft::sImage32 *root = gft::Image32::Create(gftImg);
  gft::sImage32 *pred = gft::Image32::Create(gftImg);
  gft::sImage32 *cost = gft::Image32::Create(gftImg); // cost = edt^2
  gft::sImage32 *Bedt = gft::Image32::Create(gftImg);
  
  // IFT adjacency 
  gft::sAdjRel* A8 = gft::AdjRel::Neighborhood_8();  

  initPredAndRoot(pred, root);
  
  // define priority queues
  int nb = SQUARE(MIN(gftImg->ncols, gftImg->nrows) / 2.0 + 1);
  gft::sPQueue32 *Q_edt = gft::PQueue32::Create(nb, gftImg->n, cost->data);

  // define variables from incremental contour
  std::vector<std::unordered_set<uint32>> contours(tree.numberOfNodes());
  std::vector<uint8> ncount(domain_.numberOfPoints());
  
  // Contour adjacency
  std::shared_ptr<Adjacency> adj = std::make_shared<InfAdjacency4C>(domain_);
  
  // process the level sets from 255 down to 0
  for (int level=255; level >= 0; level--) {    
    // skip level that does not contain nodes
    const std::vector<NodePtr> nodes = levelToNodes[level];
    if (nodes.empty())
      continue;

    std::cout << "level: " << level << std::endl;

    // removed contour pixels of level "level"
    std::vector<uint32> toRemove;

    // There exist at least one node in "level"
    // So, we have to process them.
    for (NodePtr node : nodes) {
      // process node 
      
      // define Ncontour which will be processed
      std::unordered_set<uint32> &Ncontour = contours[node->id()];
      
      // reuse children contour pixels
      for (NodePtr c : node->children()) {
        for (uint32 pidx : contours[c->id()])
          Ncontour.insert(pidx);
      }

      for (uint32 pidx : node->cnps()) {
         for (uint32 qidx : adj->neighbours(pidx)) {
           if (qidx != Box::UndefinedIndex && f_[pidx] < f_[qidx]) {
              if (ncount[qidx] == 0) {
                // qidx does not have a background neighbour anymore. Remove
                // it from the contour.
                Ncontour.erase(qidx);

                // store removed pixel
                toRemove.push_back(qidx);
            }
          }
        }
      }

      if (!toRemove.empty())
        treeRemoval(toRemove, bin, Q_edt, root, pred, cost, A8);

      // compute new contour points and remove contours 
      // points analysing node CNPs
      for (uint32 pidx : node->cnps()) {
        // incremetally create level-set binary image
        bin->data[pidx] = 1;

        for (uint32 qidx : adj->neighbours(pidx)) {
          if (qidx == Box::UndefinedIndex || f_[pidx] > f_[qidx]) {
            // qidx is background neighbour, thus count it.
            ncount[pidx]++;
          }
          // else if (f_[pidx] < f_[qidx]) {
          //   // pidx was a background pixel of qidx, but it is not anymore
          //   // "remove" pidx from the qidx count.
          //   ncount[qidx]--;

            
          //   if (ncount[qidx] == 0) {
          //     // qidx does not have a background neighbour anymore. Remove
          //     // it from the contour.
          //     Ncontour.erase(qidx);

          //     // store removed pixel
          //     toRemove.push_back(qidx);
          //   }            
          // }          
        } // end loop on pidx neighbours
        
        if (ncount[pidx] > 0) {
            // pidx has at least one background pixel
            // add it to the contour
            Ncontour.insert(pidx);

            // setting up for the diffential ift.
            root->data[pidx] = pidx;
            pred->data[pidx] = NIL;
            cost->data[pidx] = 0;
            gft::PQueue32::FastInsertElem(Q_edt, static_cast<int>(pidx));
          }
          else {  // pidx is pixel of level "level" and it is not a contour pixel
            cost->data[pidx] = INT_MAX;
            insertNeighborsPQueue(pidx, adj, bin, cost, Q_edt);
          }
      } // end loop on node cnps
    } // end of loop on nodes of the levels

    // if there exists contour pixels removed, remove it from
    // ift setting up
    // TODO: Adapt "treeRemoval Function"
    

    // The roots, costs and borders are set. 
    // Compute maximum distance attributes for 
    // the nodes in "level"
    
    #ifdef APPDEBUG
      std::stringstream ss;
      ss << "../out/" << level << ".pgm";
      const char *filename = ss.str().c_str();
      gft::Image32::Write(bin, (char *)filename);
    #endif

    // Run differential EDT
    EDT_DIFF(Q_edt, A8, bin, root, pred, cost, Bedt);

   #ifdef APPDEBUG
      std::stringstream sscost;
      sscost << "../out/cost-" << level << ".pgm";
      const char *filenamecost = sscost.str().c_str();
      gft::Image32::Write(cost, (char *)filenamecost);
    #endif


    // compute the max dist attribute for each node
    for (NodePtr node : nodes) {
      const std::unordered_set<uint32> &NContour = contours[node->id()];
      uint32 maxDistValue = 0;

      // search for maximum distance of N on its contour
      for (uint32 pidx : NContour) {
        if (Bedt->data[pidx] > maxDistValue)
          maxDistValue = Bedt->data[pidx];
      }
      
      maxDist[node->id()] = maxDistValue;
      std::cout << sqrt(maxDist[node->id()]) << " ";
    }
    std::cout << "\n";
  }

  // Clean up memory
  gft::Image32::Destroy(&cost);
  gft::Image32::Destroy(&pred);
  gft::Image32::Destroy(&root);
  gft::PQueue32::Destroy(&Q_edt);
  gft::Image32::Destroy(&gftImg);
  gft::Image32::Destroy(&bin);
  gft::AdjRel::Destroy(&A8);

  return maxDist;
}

std::array<std::vector<MaxDistComputer::NodePtr>, 256>
  MaxDistComputer::extractLevelMap(const MTree &tree) const
{
  std::array<std::vector<NodePtr>, 256> levelToNodes;
  tree.tranverse([&levelToNodes](NodePtr node) {
    levelToNodes[node->level()].push_back(node);
  });

  return levelToNodes;
}

gft::sImage32 *MaxDistComputer::createGFTImage() const
{
  gft::sImage32 *gftImg = gft::Image32::Create(domain_.width(), domain_.height());
  for (uint32 pidx = 0; pidx < domain_.numberOfPoints(); pidx++) {
    gftImg->data[pidx] = static_cast<int>(f_[pidx]);
  }

  return gftImg;
}

void MaxDistComputer::initPredAndRoot(gft::sImage32 *pred, gft::sImage32 *root) const
{
  for (uint32 pidx = 0; pidx < pred->n; pidx++) {
    pred->data[pidx] = NIL;
    root->data[pidx] = pidx;
  }  
}

void MaxDistComputer::insertNeighborsPQueue(
    uint32 pidx, 
    std::shared_ptr<Adjacency> adj,
    gft::sImage32 *bin,
    gft::sImage32 *cost, 
    gft::sPQueue32 *Q) const
{
  for (uint32 qidx : adj->neighbours(pidx)) {
    if (qidx != Box::UndefinedIndex) { // qidx is in the domain
      if (bin->data[qidx] > 0 
        && cost->data[qidx] != INT_MAX 
        && Q->L.elem[qidx].color != GRAY) {
        Q->L.elem[qidx].color = WHITE;
        gft::PQueue32::FastInsertElem(Q, qidx);
      }
    }
  }
}