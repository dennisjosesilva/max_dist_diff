#include "UltimateAttributeOpening.hpp"

    
    UltimateAttributeOpening::UltimateAttributeOpening(morphotree::MorphologicalTree<uint8> maxtree, int maxCriterion, std::vector<morphotree::uint32> attrs_increasing){
      this->maxtree = maxtree;
      this->maxContrastLUT = new int[this->maxtree.numberOfNodes()];
      this->associatedIndexLUT = new int[this->maxtree.numberOfNodes()];
      this->attrs_increasing = attrs_increasing;
      this->maxCriterion = maxCriterion;
      

      for (int id=0; id < maxtree.numberOfNodes(); id++) {
        maxContrastLUT[id] = 0;
        associatedIndexLUT[id] = 0;
      }

      for(NodePtr son: this->maxtree.root()->children()){
		    computeUAO(son, this->maxtree.root()->level(), false);
	    }
    }

    UltimateAttributeOpening::~UltimateAttributeOpening(){
      free(maxContrastLUT);
      free(associatedIndexLUT);
    }

    std::vector<morphotree::uint8> UltimateAttributeOpening::createRandomColor(std::vector<morphotree::uint32> img){
      int max = 0;
      for(int i=0; i < img.size(); i++){
        if(max < img[i])
          max = img[i];
      }

      int r[max+1];
      int g[max+1];
      int b[max+1];
      r[0] = 0;
      g[0] = 0;
      r[0] = 0;
      for(int i= 1; i <= max; i++){
        r[i] = rand() % 256;
        g[i] = rand() % 256;
        b[i] = rand() % 256;
      }
      std::vector<uint8> output(img.size()*3);
      for (int pidx=0; pidx < img.size(); pidx++) {

        uint32 cpidx = pidx * 3;           // (coloured) for 3 channels
        output[cpidx ] = r[ img[pidx] ];
        output[cpidx +1] = g[ img[pidx] ];
        output[cpidx +2] = b[ img[pidx] ];
          
      }
      return output;
    }

    std::vector<morphotree::uint8> UltimateAttributeOpening::getMaxConstrastImage(){
      std::vector<uint8> out(maxtree.numberOfCNPs(), 0);
      
      for (int pidx=0; pidx < maxtree.numberOfCNPs(); pidx++) {
        out[pidx] = maxContrastLUT[maxtree.smallComponent(pidx)->id()];
      }
      return out;
    }

    std::vector<morphotree::uint32> UltimateAttributeOpening::getAssociatedImage(){
      std::vector<uint32> out(maxtree.numberOfCNPs(), 0);
      for (int pidx=0; pidx < maxtree.numberOfCNPs(); pidx++) {
        out[pidx] = associatedIndexLUT[maxtree.smallComponent(pidx)->id()];
      }
      return out;
    }

    std::vector<morphotree::uint8> UltimateAttributeOpening::getAssociatedColorImage(){
      return createRandomColor(getAssociatedImage());
    }

  


    void UltimateAttributeOpening::computeUAO(NodePtr currentNode, int levelNodeNotInNR, bool qPropag){
      NodePtr parentNode = currentNode->parent();
      bool flagPropag = false;
      if (this->attrs_increasing[currentNode->id()] != this->attrs_increasing[currentNode->id()]){
          levelNodeNotInNR = parentNode->level();
      }
      
      if(this->attrs_increasing[currentNode->id()] <= this->maxCriterion){
        int contrast = (int) std::abs( currentNode->level() - levelNodeNotInNR );

        if (this->maxContrastLUT[parentNode->id()] >= contrast) {
            this->maxContrastLUT[currentNode->id()] = this->maxContrastLUT[parentNode->id()];
            this->associatedIndexLUT[currentNode->id()] = this->associatedIndexLUT[parentNode->id()];
        }
        else{
          this->maxContrastLUT[currentNode->id()] = contrast;
          if(qPropag){
            this->associatedIndexLUT[currentNode->id()] = this->associatedIndexLUT[parentNode->id()];
          }
          else{
            this->associatedIndexLUT[currentNode->id()] = this->attrs_increasing[currentNode->id()] + 1;	
          }
          flagPropag = true;
        }
      }

      for(NodePtr son: currentNode->children()){
        this->computeUAO(son, levelNodeNotInNR, flagPropag);
      }
    }

 