#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::cout;
using std::endl;

#include "include.h"

int SPM::dim;
int *SPM::deg;

/**
 * initializes the statics
 */
void SPM::init(){

   dim = Tools::gL()/2  + 1;

   deg = new int [dim];

   deg[0] = 1;
   deg[Tools::gL()/2] = 1;

   for(int k = 1;k < Tools::gL()/2;++k)
      deg[k] = 2;

}

/**
 * deallocates the static lists
 */
void SPM::clear(){

   delete [] deg;

}

/**
 * constructor, the SPM is completely diagonal in momentum space and two times degenerate in parity, for 0 < k < Tools::gL()/2
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
 * PHM constructor: Creates a SPM initialized on the "bar" of the PHM.
 * @param scale the factor u want the SPM to be scaled with 
 * @param tpm the TPM out of which the SPM will be initiated.
 */
SPM::SPM(double scale,const PHM &phm) {

   spm = new double [dim];

   this->bar(scale,phm);

}

/**
 * destructor
 */
SPM::~SPM(){

   delete [] spm;

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
 * write access to the SPM, change the number on index i
 * @param i row number
 * @return the entry on place i
 */
double &SPM::operator[](int i){

   if(i > Tools::gL()/2)
      return spm[Tools::gL() - i];
   else
      return spm[i];

}

/**
 * read access to the SPM, view the number on index i
 * @param i row number
 * @return the entry on place i
 */
double SPM::operator[](int i) const {

   if(i > Tools::gL()/2)
      return spm[Tools::gL() - i];
   else
      return spm[i];

}

/**
 * Trace out a set of indices to create the "bar" matrix of a TPM
 * @param scale the factor u want the SPM to be scaled with (1/N-1 for normal sp density matrix)
 * @param tpm the TPM out of which the SPM will be filled
 */
void SPM::bar(double scale,const TPM &tpm){

   double ward;
   int K;

   for(int k = 0;k <= Tools::gL()/2;++k){

      spm[k] = 0.0;

      for(int k_ = 0;k_ < Tools::gL();++k_){

         K = (k + k_)%Tools::gL();

         ward = 0.0;

         for(int Z = 0;Z < 2;++Z)
            for(int p = 0;p < 2;++p)
               ward += (2.0*Z + 1)*tpm(Z,K,p,k,k_,k,k_);

         if(k != k_)
            ward /= 4.0*TPM::norm(K,k,k_)*TPM::norm(K,k,k_);
         else
            ward /= 2.0*TPM::norm(K,k,k_)*TPM::norm(K,k,k_);

         spm[k] += ward;

      }

      spm[k] *= 0.5*scale;

   }

}

/**
 * Trace out a set of indices to create the "bar" matrix of a PHM
 * @param scale the factor u want the SPM to be scaled with
 * @param phm the PHM out of which the SPM will be filled
 */
void SPM::bar(double scale,const PHM &phm){

   double ward;
   int K;

   for(int k = 0;k <= Tools::gL()/2;++k){

      spm[k] = 0.0;

      for(int k_ = 0;k_ < Tools::gL();++k_){

         K = (k + k_)%Tools::gL();

         ward = 0.0;

         for(int Z = 0;Z < 2;++Z)
            for(int p = 0;p < 2;++p)
               ward += (2.0*Z + 1)*phm(Z,K,p,k,k_,k,k_);

         ward /= 4.0*PHM::norm(K,k,k_)*PHM::norm(K,k,k_);

         spm[k] += ward;

      }

      spm[k] *= 0.5*scale;

   }

}

/** 
 * This bar function maps a PPHM object directly onto a SPM object, scaling it with a factor scale
 * @param scale the scalefactor
 * @param pphm Input PPHM object
 */
void SPM::bar(double scale,const PPHM &pphm){

   double ward;

   int K_pph;

   for(int k = 0;k <= Tools::gL()/2;++k){

      (*this)[k] = 0.0;

      //first S = 1/2 part
      for(int S_ab = 0;S_ab < 2;++S_ab){//S_ab can be both 0 and 1

         for(int k_a = 0;k_a < Tools::gL();++k_a)
            for(int k_b = 0;k_b < Tools::gL();++k_b){

               ward = 0.0;

               K_pph = (k_a + k_b + k)%Tools::gL();

               for(int pi = 0;pi < 2;++pi)
                  ward += pphm.pph(0,K_pph,pi,S_ab,k_a,k_b,k,S_ab,k_a,k_b,k);

               if(k_a == k_b)
                  ward *= 2.0;

               (*this)[k] += 0.25/(PPHM::norm(K_pph,k_a,k_b,k) * PPHM::norm(K_pph,k_a,k_b,k)) * ward;

            }

      }

      //then S = 3/2 part:
      for(int k_a = 0;k_a < Tools::gL();++k_a)
         for(int k_b = 0;k_b < Tools::gL();++k_b){

            ward = 0.0;

            K_pph = (k_a + k_b + k)%Tools::gL();

            for(int pi = 0;pi < 2;++pi)
               ward += pphm.pph(1,K_pph,pi,1,k_a,k_b,k,1,k_a,k_b,k);

            (*this)[k] += 0.5/(PPHM::norm(K_pph,k_a,k_b,k) * PPHM::norm(K_pph,k_a,k_b,k)) * ward;

         }

      //scaling
      (*this)[k] *= scale;

   }

}

/**
 * @param k the pseudomomentum of the parity state
 * @return the norm corresponding to the state with pseudo-momentum k
 */
double SPM::norm(int k) {

   if(k == 0 || k == Tools::gL()/2)
      return 0.5;
   else
      return 1.0/std::sqrt(2.0);

}
