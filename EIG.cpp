#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::endl;

#include "include.h"

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

#ifdef __T1_CON
   v_dp = new BlockVector<DPM>(SZ.dpm());
#endif

#ifdef __T2_CON
   v_pph = new BlockVector<PPHM>(SZ.pphm());
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

#ifdef __T1_CON
   v_dp = new BlockVector<DPM>(eig_c.dpv());
#endif

#ifdef __T2_CON
   v_pph = new BlockVector<PPHM>(eig_c.pphv());
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

#ifdef __T1_CON

   *v_dp = *eig_c.v_dp;

#endif

#ifdef __T2_CON

   *v_pph = *eig_c.v_pph;

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

#ifdef __T1_CON
   
   delete v_dp;

#endif

#ifdef __T2_CON
   
   delete v_pph;

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

#ifdef __T1_CON
   
   v_dp->diagonalize(sup.dpm());

#endif

#ifdef __T2_CON
   
   v_pph->diagonalize(sup.pphm());

#endif

}

ostream &operator<<(ostream &output,const EIG &eig_p){

   for(int i = 0;i < 2;++i)
      std::cout << eig_p.tpv(i) << std::endl;

#ifdef __G_CON
   
   std::cout << eig_p.phv() << std::endl;

#endif

#ifdef __T1_CON
   
   std::cout << eig_p.dpv() << std::endl;

#endif

#ifdef __T2_CON
   
   std::cout << eig_p.pphv() << std::endl;

#endif

   return output;

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

#ifdef __T1_CON

/** 
 * get the BlockVector<DPM> object containing the eigenvalues of the DPM block T1 of the SUP matrix
 * @return a BlockVector<DPM> object containing the desired eigenvalues
 */
BlockVector<DPM> &EIG::dpv(){

   return *v_dp;

}

/** 
 * get the BlockVector<DPM> object containing the eigenvalues of the DPM block T1 of the SUP matrix, const version
 * @return a BlockVector<DPM> object containing the desired eigenvalues
 */
const BlockVector<DPM> &EIG::dpv() const{

   return *v_dp;

}

#endif

#ifdef __T2_CON

/** 
 * get the BlockVector<PPHM> object containing the eigenvalues of the PPHM block T2 of the SUP matrix
 * @return a BlockVector<PPHM> object containing the desired eigenvalues
 */
BlockVector<PPHM> &EIG::pphv(){

   return *v_pph;

}

/** 
 * get the BlockVector<PPHM> object containing the eigenvalues of the PPHM block T2 of the SUP matrix, const version
 * @return a BlockVector<PPHM> object containing the desired eigenvalues
 */
const BlockVector<PPHM> &EIG::pphv() const{

   return *v_pph;

}

#endif

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

#ifdef __T1_CON

   //lowest eigenvalue of the T1 block
   if(ward > v_dp->min())
      ward = v_dp->min();

#endif

#ifdef __T2_CON

   //lowest eigenvalue of the T2 block
   if(ward > v_pph->min())
      ward = v_pph->min();

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

#ifdef __T1_CON

   //highest eigenvalue of the T1 block
   if(ward < v_dp->max())
      ward = v_dp->max();

#endif

#ifdef __T2_CON

   //highest eigenvalue of the T2 block
   if(ward < v_pph->max())
      ward = v_pph->max();

#endif

   return ward;

}
