#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::cout;
using std::endl;

#include "include.h"

int SPM::M;
int SPM::N;
int SPM::L;
int SPM::dim;
int *SPM::deg;

/**
 * initializes the statics
 * @param L_in the nr of sites
 * @param N_in nr of particles
 */
void SPM::init(int L_in,int N_in){

   L = L_in;
   N = N_in;

   M = 2*L;

   dim = L/2  + 1;

   deg = new int [dim];

   deg[0] = 1;
   deg[L/2] = 1;

   for(int k = 1;k < L/2;++k)
      deg[k] = 2;

}

/**
 * deallocates the static lists
 */
void SPM::clear(){

   delete [] deg;

}

/**
 * constructor, the SPM is completely diagonal in momentum space and two times degenerate in parity, for 0 < k < L/2
 */
SPM::SPM(){

   spm = new double [dim];

}

/**
 * copy constructor
 * @param spm_copy content of this matrix will be copied into the constructed matrix
 */
SPM::SPM(const SPM &spm_copy) {

   spm = new double [dim];

   int inc = 1;

   dcopy_(&dim,spm_copy.spm,&inc,spm,&inc);

}

/**
 * TPM constructor: Creates a SPM initialized on the "bar" of the TPM.
 * @param scale the factor u want the SPM to be scaled with (1/N-1 for normal sp density matrix)
 * @param tpm the TPM out of which the SPM will be initiated.
 */
SPM::SPM(double scale,const TPM &tpm) {

   spm = new double [dim];

   this->bar(scale,tpm);

}

/**
 * destructor
 */
SPM::~SPM(){

   delete [] spm;

}

/**
 * @return nr of particles
 */
int SPM::gN() const{

   return N;

}

/**
 * @return dimension of sp space
 */
int SPM::gM() const{

   return M;

}

/**
 * @return nr of sites
 */
int SPM::gL() const{

   return L;

}

/**
 * @return the pointer to the SPM object
 */
const double *SPM::gspm() const{

   return spm;

}

ostream &operator<<(ostream &output,const SPM &spm_p){

   for(int k = 0;k < spm_p.dim;++k)
      output << k << "\t" << spm_p[k] << endl;

   return output;

}

/**
 * Trace out a set of indices to create the "bar" matrix of a TPM
 * @param scale the factor u want the SPM to be scaled with (1/N-1 for normal sp density matrix)
 * @param tpm the TPM out of which the SPM will be filled
 */
void SPM::bar(double scale,const TPM &tpm){

   for(int k = 0;k <= L/2;++k){

      spm[k] = 0.0;

      //first K = 0
      int k_ = (L - k)%L;

      if(k == k_)
         spm[k] += 2.0*tpm(0,0,0,k,k_,k,k_);
      else
         spm[k] += tpm(0,0,0,k,k_,k,k_) + 3.0 * tpm(1,0,1,k,k_,k,k_);

      //then 0 < K < L/2: 
      for(int K = 1;K < L/2;++K){

         k_ = (K - k + L)%L;

         if(k == k_)
            spm[k] += 2.0*tpm(0,K,0,k,k_,k,k_);
         else
            spm[k] += tpm(0,K,0,k,k_,k,k_) + 3.0 * tpm(1,K,0,k,k_,k,k_);

         //conjugate term
         k_ = ( (L - K) - k + L)%L;

         if(k == k_)
            spm[k] += 2.0*tpm(0,K,0,(L - k)%L,(L - k_)%L,(L - k)%L,(L - k_)%L);
         else
            spm[k] += tpm(0,K,0,(L - k)%L,(L - k_)%L,(L - k)%L,(L - k_)%L) + 3.0 * tpm(1,K,0,(L - k)%L,(L - k_)%L,(L - k)%L,(L - k_)%L);

      }

      //for K == L/2, both positive and negative parity should be added
      k_ = L/2 - k;

      if(k == k_){

         //|0 L/2> is not here, so positive and negative parity should be avaraged.
         if(k != 0 && k_ != 0)
            spm[k] += tpm(0,L/2,0,k,k_,k,k_) + tpm(0,L/2,1,k,k_,k,k_);
         else//only positive parity
            spm[k] += 2.0 * tpm(0,L/2,0,k,k_,k,k_);

      }
      else{

         //|0 L/2> is not here, so positive and negative parity should be avaraged.
         if(k != 0 && k_ != 0){

            //S = 0
            spm[k] += 0.5 * (tpm(0,L/2,0,k,k_,k,k_) + tpm(0,L/2,1,k,k_,k,k_));

            //S = 1
            spm[k] += 1.5 * (tpm(1,L/2,0,k,k_,k,k_) + tpm(1,L/2,1,k,k_,k,k_));

         }
         else//only positive parity
            spm[k] += tpm(0,L/2,0,k,k_,k,k_) + 3.0 * tpm(1,L/2,0,k,k_,k,k_);

      }

      spm[k] *= 0.5 * scale;

   }

}

/**
 * write access to the SPM, change the number on index i
 * @param i row number
 * @return the entry on place i
 */
double &SPM::operator[](int i){

   return spm[i];

}

/**
 * read access to the SPM, view the number on index i
 * @param i row number
 * @return the entry on place i
 */
double SPM::operator[](int i) const {

   return spm[i];

}
