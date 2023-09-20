#include "MaxDistComputer.hpp"


MaxDistComputer::MaxDistComputer(Box domain,
  const std::vector<uint8> &f)
  : domain_{domain}, f_{f}
{}

std::vector<MaxDistComputer::uint32> MaxDistComputer::computeAttribute(
  const MTree &tree) const
{
  using NodePtr = MTree::NodePtr;

  // define useful images 
  std::vector<uint32> maxDist; 
  std::array<std::vector<NodePtr>, 256> levelToNodes = extractLevelMap(tree);
  gft::sImage32 *gftImg = createGFTImage();
  gft::sImage32 *edt = nullptr;
  gft::sImage32 *bin = gft::Image32::Create(gftImg);
  gft::sImage32 *root = gft::Image32::Create(gftImg);
  gft::sImage32 *pred = gft::Image32::Create(gftImg);
  gft::sImage32 *cost = gft::Image32::Create(gftImg);
  gft::sImage32 *Bedt = gft::Image32::Create(gftImg);
  gft::sPQueue32 *Q = nullptr, *Q_edt = nullptr;

  initPredAndRoot(pred, root);
  
  // define priority queues
  int nb = SQUARE(MIN(gftImg->ncols, gftImg->nrows) / 2.0 + 1);
  gft::sPQueue32 *Q     = gft::PQueue32::Create(nb, gftImg->n, cost->data);
  gft::sPQueue32 *Q_edt = gft::PQueue32::Create(nb, gftImg->n, cost->data);


//////////////////////////////////////////////////////////
  int *boundary = (int *)malloc((gftImg->n+1)*sizeof(int));
  int nboundary = 0;

  int Imin = gft::Image32::GetMinVal(gftImg);
  int Imax = gft::Image32::GetMaxVal(gftImg);
  int T = Imax;
  gft::sAdjRel *A4 = gft::AdjRel::Neighborhood_4();
  gft::sImage32 *min_neighbor= get_min_neighbor(gftImg, A4); //importar essa funcao get_min_neighbor
  //Enquanto: Imax < Imin  ..precisa incluir
  for (NodePtr node: levelToNodes[Imax]){
    for(int p: node->reconstruct()){
      bin->data[p] = 1;

      if(min_neighbor->data[p] < T) {
        boundary[nboundary] = p;
        nboundary++;
        
        root->data[p] = p;
        pred->data[p] = NIL;
        cost->data[p] = 0;
        
        gft::PQueue32::FastInsertElem(Q_edt, p);

      }else {
        cost->data[p] = INT_MAX;
        insert_neighbors_pqueue(p, A4, bin, cost, Q_edt);
      }
    }
  }
  runEDT_DIFF(img, bin, root, pred, cost, Bedt, boundary, nboundary, Q_edt)

////////////////////////////////////////

  for (uint32 pidx = 0; pidx < gftImg->n; pidx++) {
    gft::PQueue32::FastInsertElem(Q, pidx);
  }




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