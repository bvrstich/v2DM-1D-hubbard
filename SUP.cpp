#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

/**
 * standard constructor\n
 * Allocates two TPM matrices and optionally a PHM, DPM or PPHM matrix.
 */
SUP::SUP(){

   SZ_tp = new TPM * [2];

   for(int i = 0;i < 2;++i)
      SZ_tp[i] = new TPM();

#ifdef __G_CON
   SZ_ph = new PHM();
#endif

#ifdef __T1_CON
   SZ_dp = new DPM();
#endif

#ifdef __T2_CON
   SZ_pph = new PPHM();
#endif

}

/**
 * standard constructor\n
 * Allocates two TPM matrices and optionally a PHM, DPM or PPHM matrix, then copies the content of
 * input SUP SZ_c into it.
 * @param SZ_c input SUP
 */
SUP::SUP(const SUP &SZ_c){

   SZ_tp = new TPM * [2];

   for(int i = 0;i < 2;++i)
      SZ_tp[i] = new TPM(*SZ_c.SZ_tp[i]);

#ifdef __G_CON
   SZ_ph = new PHM(*SZ_c.SZ_ph);
#endif

#ifdef __T1_CON
   SZ_dp = new DPM(*SZ_c.SZ_dp);
#endif

#ifdef __T2_CON
   SZ_pph = new PPHM(*SZ_c.SZ_pph);
#endif

}

/**
 * Destructor
 */
SUP::~SUP(){

   for(int i = 0;i < 2;++i)
      delete SZ_tp[i];

   delete [] SZ_tp;

#ifdef __G_CON
   
   delete SZ_ph;

#endif

#ifdef __T1_CON

   delete SZ_dp;

#endif

#ifdef __T2_CON
   
   delete SZ_pph;

#endif

}

/**
 * Overload += operator
 * @param SZ_pl The SUP matrix that has to be added to this
 */
SUP &SUP::operator+=(const SUP &SZ_pl){

   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) += (*SZ_pl.SZ_tp[i]);

#ifdef __G_CON
   
   (*SZ_ph) += (*SZ_pl.SZ_ph);

#endif

#ifdef __T1_CON

   (*SZ_dp) += (*SZ_pl.SZ_dp);

#endif

#ifdef __T2_CON

   (*SZ_pph) += (*SZ_pl.SZ_pph);

#endif

   return *this;

}

/**
 * Overload -= operator
 * @param SZ_pl The SUP that will be deducted from this
 */
SUP &SUP::operator-=(const SUP &SZ_pl){

   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) -= (*SZ_pl.SZ_tp[i]);

#ifdef __G_CON
   
   (*SZ_ph) -= (*SZ_pl.SZ_ph);

#endif

#ifdef __T1_CON

   (*SZ_dp) -= (*SZ_pl.SZ_dp);

#endif

#ifdef __T2_CON

   (*SZ_pph) -= (*SZ_pl.SZ_pph);

#endif

   return *this;

}

/**
 * Overload equality operator, copy SZ_c into this
 * @param SZ_c SUP_PQ to be copied into this
 */
SUP &SUP::operator=(const SUP &SZ_c){

   (*SZ_tp[0]) = (*SZ_c.SZ_tp[0]);
   (*SZ_tp[1]) = (*SZ_c.SZ_tp[1]);

#ifdef __G_CON

   (*SZ_ph) = (*SZ_c.SZ_ph);

#endif

#ifdef __T1_CON

   (*SZ_dp) = (*SZ_c.SZ_dp);

#endif

#ifdef __T2_CON

   (*SZ_pph) = (*SZ_c.SZ_pph);

#endif

   return *this;

}

/**
 * overload operator = number, all the blockmatrices in SUP are put equal to the number a.
 * e.g. SZ = 0 makes all the Matrix elements zero.
 * @param a the number
 */
SUP &SUP::operator=(double &a){

   (*SZ_tp[0]) = a;
   (*SZ_tp[1]) = a;

#ifdef __G_CON

   (*SZ_ph) = a;

#endif

#ifdef __T1_CON

   (*SZ_dp) = a;

#endif

#ifdef __T2_CON

   (*SZ_pph) = a;

#endif

   return *this;

}

/**
 * @param i which block you want to have the pointer to.
 * @return pointer to the individual TPM blocks: SZ_tp[i]
 */
TPM &SUP::tpm(int i){

   return *SZ_tp[i];

}

/**
 * const version
 * @param i which block you want to have the pointer to.
 * @return pointer to the individual TPM blocks: SZ_tp[i]
 */
const TPM &SUP::tpm(int i) const{

   return *SZ_tp[i];

}

#ifdef __G_CON

/**
 * @return pointer to the PHM block: SZ_ph
 */
PHM &SUP::phm(){

   return *SZ_ph;

}

/**
 * const version
 * @return pointer to the PHM block: SZ_ph
 */
const PHM &SUP::phm() const{

   return *SZ_ph;

}

#endif

#ifdef __T1_CON

/**
 * @return pointer to the DPM block: SZ_dp
 */
DPM &SUP::dpm(){

   return *SZ_dp;

}

/**
 * const version
 * @return pointer to the DPM block: SZ_dp
 */
const DPM &SUP::dpm() const{

   return *SZ_dp;

}

#endif

#ifdef __T2_CON

/**
 * const version
 * @return pointer to the PPHM block: SZ_pph
 */
const PPHM &SUP::pphm() const{

   return *SZ_pph;

}

/**
 * @return pointer to the PPHM block: SZ_pph
 */
PPHM &SUP::pphm(){

   return *SZ_pph;

}

#endif

/**
 * Initialization of the SUP matrix S, is just u^0: see primal_dual.pdf for more information
 */
void SUP::init_S(){

   (*SZ_tp[0]).unit();

   this->fill();

}

ostream &operator<<(ostream &output,const SUP &SZ_p){

   output << (*SZ_p.SZ_tp[0]) << std::endl;
   output << (*SZ_p.SZ_tp[1]);

#ifdef __G_CON

   output << std::endl;
   output << (*SZ_p.SZ_ph);

#endif

#ifdef __T1_CON

   output << std::endl;
   output << (*SZ_p.SZ_dp);

#endif

#ifdef __T2_CON

   output << std::endl;
   output << (*SZ_p.SZ_pph);

#endif

   return output;

}

/**
 * Fill the SUP matrix with random elements, watch out, not necessarily positive definite
 */
void SUP::fill_Random(){

   SZ_tp[0]->fill_Random();
   SZ_tp[1]->fill_Random();

#ifdef __G_CON

   SZ_ph->fill_Random();

#endif

#ifdef __T1_CON

   SZ_dp->fill_Random();

#endif

#ifdef __T2_CON

   SZ_pph->fill_Random();

#endif

}

/**
 * @param SZ_i input SUP_PQ SZ_i
 * @return inproduct between this and input matrix SZ_i, defined as Tr(this SZ_i)
 */
double SUP::ddot(const SUP &SZ_i) const{

   double ward = 0.0;

   for(int i = 0;i < 2;++i)
      ward += SZ_tp[i]->ddot(*SZ_i.SZ_tp[i]);

#ifdef __G_CON
   
   ward += SZ_ph->ddot(*SZ_i.SZ_ph);

#endif

#ifdef __T1_CON

   ward += SZ_dp->ddot(*SZ_i.SZ_dp);

#endif

#ifdef __T2_CON

   ward += SZ_pph->ddot(*SZ_i.SZ_pph);

#endif

   return ward;

}

/**
 * Invert all the Matrices in SUP and put it in this, watch out, destroys original matrices.
 * Makes use of cholesky decomposition, so only positive definite matrices can be used as input!
 */
void SUP::invert(){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->invert();

#ifdef __G_CON
   
   SZ_ph->invert();

#endif

#ifdef __T1_CON
   
   SZ_dp->invert();

#endif

#ifdef __T2_CON
   
   SZ_pph->invert();

#endif

}

/**
 * Scale all the elements in the blockmatrix with parameter alpha
 * @param alpha the scalefactor
 */
void SUP::dscal(double alpha){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->dscal(alpha);

#ifdef __G_CON
   
   SZ_ph->dscal(alpha);

#endif

#ifdef __T1_CON
   
   SZ_dp->dscal(alpha);

#endif

#ifdef __T2_CON
   
   SZ_pph->dscal(alpha);

#endif

}

/**
 * Take the square root out of the positive semidefinite SUP matrix this and put it in this, watch out, original matrix is destroyed.
 * @param option = 1 positive square root, = -1 negative square root.
 */
void SUP::sqrt(int option){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->sqrt(option);

#ifdef __G_CON

   SZ_ph->sqrt(option);

#endif

#ifdef __T1_CON

   SZ_dp->sqrt(option);

#endif

#ifdef __T2_CON

   SZ_pph->sqrt(option);

#endif

}

/**
 * Multiply symmetric SUP blockmatrix object left en right with symmetric SUP blockmatrix map to 
 * form another symmetric matrix and put it in (*this): this = map*object*map
 * @param map SUP that will be multiplied to the left en to the right of matrix object
 * @param object central SUP
 */
void SUP::L_map(const SUP &map,const SUP &object){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->L_map(map.tpm(i),object.tpm(i));

#ifdef __G_CON

   SZ_ph->L_map(map.phm(),object.phm());

#endif

#ifdef __T1_CON

   SZ_dp->L_map(map.dpm(),object.dpm());

#endif

#ifdef __T2_CON

   SZ_pph->L_map(map.pphm(),object.pphm());

#endif

}

/**
 * add the SUP SZ_p times the constant alpha to this
 * @param alpha the constant to multiply the SZ_p with
 * @param SZ_p the SUP to be multiplied by alpha and added to (*this)
 */
void SUP::daxpy(double alpha,const SUP &SZ_p){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->daxpy(alpha,SZ_p.tpm(i));

#ifdef __G_CON
   
   SZ_ph->daxpy(alpha,SZ_p.phm());

#endif

#ifdef __T1_CON
   
   SZ_dp->daxpy(alpha,SZ_p.dpm());

#endif

#ifdef __T2_CON
   
   SZ_pph->daxpy(alpha,SZ_p.pphm());

#endif

}

/**
 * General matrixproduct between two SUP matrices, act with Matrix::mprod on every block
 * 
 * @param A left hand matrix
 * @param B right hand matrix
 * @return The product AB
 */
SUP &SUP::mprod(const SUP &A,const SUP &B){

   for(int i= 0;i < 2;++i)
      SZ_tp[i]->mprod(A.tpm(i),B.tpm(i));

#ifdef __G_CON

   SZ_ph->mprod(A.phm(),B.phm());

#endif

#ifdef __T1_CON

   SZ_dp->mprod(A.dpm(),B.dpm());

#endif

#ifdef __T2_CON

   SZ_pph->mprod(A.pphm(),B.pphm());

#endif

   return *this;

}

/**
 * Fill the SUP matrix (*this) with a TPM matrix like: this = diag[tpm  Q(tpm)  ( G(tpm) T1(tpm) T2(tpm) ) ]
 * @param tpm input TPM
 */
void SUP::fill(const TPM &tpm){

   *SZ_tp[0] = tpm;
   SZ_tp[1]->Q(1,tpm);

#ifdef __G_CON
   
   SZ_ph->G(tpm);

#endif

#ifdef __T1_CON
   
   SZ_dp->T(tpm);

#endif

#ifdef __T2_CON
   
   SZ_pph->T(tpm);

#endif

}

/**
 * fill the SUP matrix with the TPM matrix stored in the first block:\n\n
 * this = diag[this->tpm(0) Q(this->tpm(0)) ( G(this->tpm(0)) T1(this->tpm(0)) T2(this-tpm(0)) ) ]
 */
void SUP::fill(){

   SZ_tp[1]->Q(1,*SZ_tp[0]);

#ifdef __G_CON

   SZ_ph->G(*SZ_tp[0]);

#endif 

#ifdef __T1_CON

   SZ_dp->T(*SZ_tp[0]);

#endif 

#ifdef __T2_CON

   SZ_pph->T(*SZ_tp[0]);

#endif 

}

/*
 * Output the SUP matrix to a file, associated with the ofstream object output
 * @param output The ofstream object you will dump the SUP in.
 */
void SUP::out(ofstream &output) const{

   output.precision(15);

   //P
   for(int B = 0;B < this->tpm(0).gnr();++B){

      for(int i = 0;i < this->tpm(0).gdim(B);++i)
         for(int j = 0;j < this->tpm(0).gdim(B);++j)
            output << i << "\t" << j << "\t" << this->tpm(0)(B,i,j) << std::endl;

   }

   //Q
   for(int B = 0;B < this->tpm(1).gnr();++B){

      for(int i = 0;i < this->tpm(1).gdim(B);++i)
         for(int j = 0;j < this->tpm(1).gdim(B);++j)
            output << i << "\t" << j << "\t" << this->tpm(1)(B,i,j) << std::endl;

   }

#ifdef __G_CON

   //G
   for(int B = 0;B < this->phm().gnr();++B){

      for(int i = 0;i < this->phm().gdim(B);++i)
         for(int j = 0;j < this->phm().gdim(B);++j)
            output << i << "\t" << j << "\t" << this->phm()(B,i,j) << std::endl;

   }

#endif

#ifdef __T1_CON

   //T1
   for(int B = 0;B < this->dpm().gnr();++B){

      for(int i = 0;i < this->dpm().gdim(B);++i)
         for(int j = 0;j < this->dpm().gdim(B);++j)
            output << i << "\t" << j << "\t" << this->dpm()(B,i,j) << std::endl;

   }

#endif

#ifdef __T2_CON

   //T2
   for(int B = 0;B < this->pphm().gnr();++B){

      for(int i = 0;i < this->pphm().gdim(B);++i)
         for(int j = 0;j < this->pphm().gdim(B);++j)
            output << i << "\t" << j << "\t" << this->pphm()(B,i,j) << std::endl;

   }

#endif

}

/*
 * Input a SUP matrix from a file
 * @param input ifstream object association with input file.
 */
void SUP::in(ifstream &input){

   int I,J;

   //first P
   for(int B = 0;B < (this->tpm(0)).gnr();++B){

      for(int i = 0;i < (this->tpm(0)).gdim(B);++i)
         for(int j = 0;j < (this->tpm(0)).gdim(B);++j)
            input >> I >> J >> (this->tpm(0))(B,i,j);

   }

   //then Q
   for(int B = 0;B < (this->tpm(1)).gnr();++B){

      for(int i = 0;i < (this->tpm(1)).gdim(B);++i)
         for(int j = 0;j < (this->tpm(1)).gdim(B);++j)
            input >> I >> J >> (this->tpm(1))(B,i,j);

   }

#ifdef __G_CON

   //then G
   for(int B = 0;B < (this->phm()).gnr();++B){

      for(int i = 0;i < (this->phm()).gdim(B);++i)
         for(int j = 0;j < (this->phm()).gdim(B);++j)
            input >> I >> J >> (this->phm())(B,i,j);

   }

#endif

#ifdef __T1_CON

   //then T1
   for(int B = 0;B < (this->dpm()).gnr();++B){

      for(int i = 0;i < (this->dpm()).gdim(B);++i)
         for(int j = 0;j < (this->dpm()).gdim(B);++j)
            input >> I >> J >> (this->dpm())(B,i,j);

   }

#endif

#ifdef __T2_CON

   //(finally) T2
   for(int B = 0;B < (this->pphm()).gnr();++B){

      for(int i = 0;i < (this->pphm()).gdim(B);++i)
         for(int j = 0;j < (this->pphm()).gdim(B);++j)
            input >> I >> J >> (this->pphm())(B,i,j);

   }

#endif

}
