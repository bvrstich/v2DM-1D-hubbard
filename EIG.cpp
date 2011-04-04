#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::endl;

#include "include.h"

int EIG::M;
int EIG::N;
int EIG::L;
int EIG::dim;

/**
 * initialize the statics
 * @param L_in nr of sites
 * @param N_in nr of particles
 */
void EIG::init(int L_in,int N_in){

   N = N_in;
   L = L_in;

   M = 2*L;
   dim = M*(M - 1);

#ifdef __G_CON
   
   dim += M*M;

#endif

}

/**
 * standard constructor with initialization on the eigenvalues of a SUP object.
 * @param SZ input SUP object that will be destroyed after this function is called. The eigenvectors
 * of the matrix will be stored in the columns of the original SUP matrix.
 */
EIG::EIG(SUP &SZ){

   v_tp = new BlockVector<TPM> * [2];

   for(int i = 0;i < 2;++i)
      v_tp[i] = new BlockVector<TPM>(SZ.tpm(i));

#ifdef __G_CON
   
   v_ph = new BlockVector<PHM>(SZ.phm());

#endif

}

/**
 * Copy constructor\n
 * allocates the memory for the eigenvalues of a SUP object and copies the content of eig_c into it.
 * @param eig_c The input EIG that will be copied into this.
 */
EIG::EIG(const EIG &eig_c){

   v_tp = new BlockVector<TPM> * [2];

   for(int i = 0;i < 2;++i)
      v_tp[i] = new BlockVector<TPM>(eig_c.tpv(i));

#ifdef __G_CON

   v_ph = new BlockVector<PHM>(eig_c.phv());

#endif

}

/**
 * overload equality operator
 * @param eig_c object that will be copied into this.
 */
EIG &EIG::operator=(const EIG &eig_c){

   for(int i = 0;i < 2;++i)
      *v_tp[i] = *eig_c.v_tp[i];

#ifdef __G_CON

   *v_ph = *eig_c.v_ph;

#endif

   return *this;

}

/**
 * Destructor, deallocation of the memory
 */
EIG::~EIG(){

   for(int i = 0;i < 2;++i)
      delete v_tp[i];

   delete [] v_tp;

#ifdef __G_CON
   
   delete v_ph;

#endif

}

/**
 * Diagonalize a SUP matrix when the memory has allready been allocated before
 * @param sup matrix to be diagonalized
 */
void EIG::diagonalize(SUP &sup){

   for(int i = 0;i < 2;++i)
      v_tp[i]->diagonalize(sup.tpm(i));

#ifdef __G_CON
   
   v_ph->diagonalize(sup.phm());

#endif

}

ostream &operator<<(ostream &output,const EIG &eig_p){

   for(int i = 0;i < 2;++i)
      std::cout << eig_p.tpv(i) << std::endl;

#ifdef __G_CON
   
   std::cout << eig_p.phv() << std::endl;

#endif

   return output;

}

/**
 * @return nr of particles
 */
int EIG::gN() const{

   return N;

}

/**
 * @return dimension of sp space
 */
int EIG::gM() const{

   return M;

}

/**
 * @return nr of sites
 */
int EIG::gL() const{

   return L;

}

/** 
 * get the BlockVector<TPM> object containing the eigenvalues of the TPM blocks P and Q
 * @param i == 0, the eigenvalues of the P block will be returned, i == 1, the eigenvalues of the Q block will be returned
 * @return a BlockVector<TPM> object containing the desired eigenvalues
 */
BlockVector<TPM> &EIG::tpv(int i){

   return *v_tp[i];

}

/** 
 * const version\n\n
 * get the BlockVector<TPM> object containing the eigenvalues of the TPM blocks P and Q
 * @param i == 0, the eigenvalues of the P block will be returned, i == 1, the eigenvalues of the Q block will be returned
 * @return a BlockVector<TPM> object containing the desired eigenvalues
 */
const BlockVector<TPM> &EIG::tpv(int i) const{

   return *v_tp[i];

}

#ifdef __G_CON

/** 
 * get the BlockVector<PHM> object containing the eigenvalues of the PHM block G
 * @return a BlockVector<PHM> object containing the desired eigenvalues
 */
BlockVector<PHM> &EIG::phv(){

   return *v_ph;

}

/** 
 * get the BlockVector<PHM> object containing the eigenvalues of the PHM block G, const version
 * @return a BlockVector<PHM> object containing the desired eigenvalues
 */
const BlockVector<PHM> &EIG::phv() const{

   return *v_ph;

}

#endif

/**
 * @return total dimension of the EIG object
 */
int EIG::gdim() const{

   return dim;

}


/**
 * @return the minimal element present in this EIG object.
 * watch out, only works when EIG is filled with the eigenvalues of a diagonalized SUP matrix
 */
double EIG::min() const{

   //lowest eigenvalue of P block
   double ward = v_tp[0]->min();

   //lowest eigenvalue of Q block
   if(ward > v_tp[1]->min())
      ward = v_tp[1]->min();

#ifdef __G_CON

   //lowest eigenvalue of G block
   if(ward > v_ph->min())
      ward = v_ph->min();

#endif

   return ward;

}

/**
 * @return the maximum element present in this EIG object.
 * watch out, only works when EIG is filled with the eigenvalues of a diagonalized SUP matrix
 */
double EIG::max() const{

   //highest eigenvalue of P block
   double ward = v_tp[0]->max();

   //highest eigenvalue of Q block
   if(ward < v_tp[1]->max())
      ward = v_tp[1]->max();

#ifdef __G_CON

   //highest eigenvalue of G block
   if(ward < v_ph->max())
      ward = v_ph->max();

#endif

   return ward;

}

/**
 * @return The deviation of the central path as calculated with the logarithmic barrierfunction, the EIG object is calculated
 * in SUP::center_dev.
 */
double EIG::center_dev() const{

   double sum = v_tp[0]->sum() + v_tp[1]->sum();

   double log_product = v_tp[0]->log_product() + v_tp[1]->log_product();

#ifdef __G_CON

   sum += v_ph->sum();

   log_product += v_ph->log_product();

#endif

   return dim*log(sum/(double)dim) - log_product;

}

/**
 * @return the deviation of the central path measured trough the logarithmic potential barrier (see primal_dual.pdf), when you take a stepsize alpha from
 * the point (S,Z) in the primal dual newton direction (DS,DZ), for which you have calculated the generalized eigenvalues eigen_S and eigen_Z in SUP::line_search.
 * (*this) = eigen_S --> generalized eigenvalues for the DS step
 * @param alpha the stepsize
 * @param eigen_Z --> generalized eigenvalues for the DS step
 * @param c_S = Tr (DS Z)/Tr (SZ): parameter calculated in SUP::line_search
 * @param c_Z = Tr (S DZ)/Tr (SZ): parameter calculated in SUP::line_search
 */
double EIG::centerpot(double alpha,const EIG &eigen_Z,double c_S,double c_Z) const{

   double ward = dim*log(1.0 + alpha*(c_S + c_Z));

   for(int i = 0;i < 2;++i)
      ward -= v_tp[i]->centerpot(alpha) + (eigen_Z.tpv(i)).centerpot(alpha);

#ifdef __G_CON

   ward -= v_ph->centerpot(alpha) + (eigen_Z.phv()).centerpot(alpha);

#endif

   return ward;

}
