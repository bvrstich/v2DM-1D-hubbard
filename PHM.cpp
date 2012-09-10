#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::ofstream;
using std::vector;
using std::cout;
using std::endl;

#include "include.h"

vector< vector<int> > *PHM::ph2s;
int ***PHM::s2ph;

int **PHM::block_char;
int ***PHM::char_block;

/**
 * initialize the static variables and allocate the static lists
 */
void PHM::init(){

   int nr_B = Tools::gL() + 6;

   //allocate some stuff
   ph2s = new vector< vector<int> > [nr_B];

   s2ph = new int ** [nr_B];

   for(int B = 0;B < nr_B;++B){

      s2ph[B] = new int * [Tools::gL()];

      for(int k = 0;k < Tools::gL();++k)
         s2ph[B][k] = new int [Tools::gL()];

   }

   block_char = new int * [nr_B];

   for(int B = 0;B < nr_B;++B)
      block_char[B] = new int [3];

   char_block = new int ** [2];

   for(int S = 0;S < 2;++S){

      char_block[S] = new int * [Tools::gL()/2 + 1];

      for(int K = 0;K <= Tools::gL()/2;++K)
         char_block[S][K] = new int [2];

   }

   //eerst K = 0: positieve en negatieve pariteit:
   int block = 0;

   int ph = 0;

   block_char[block][0] = 0;//S
   block_char[block][1] = 0;//K
   block_char[block][2] = 0;//p

   char_block[0][0][0] = block;

   block_char[block + Tools::gL()/2 + 3][0] = 1;//S
   block_char[block + Tools::gL()/2 + 3][1] = 0;//K
   block_char[block + Tools::gL()/2 + 3][2] = 0;//p

   char_block[1][0][0] = block + Tools::gL()/2 + 3;

   vector<int> v(2);

   //for positive parity: k_a <= k_b
   for(int k_a = 0;k_a < Tools::gL();++k_a)
      for(int k_b = k_a;k_b < Tools::gL();++k_b){

         if( (k_a + k_b)%Tools::gL() == 0 ){

            v[0] = k_a;
            v[1] = k_b;

            ph2s[block].push_back(v);
            ph2s[block + Tools::gL()/2 + 3].push_back(v);

            s2ph[block][k_a][k_b] = ph;
            s2ph[block + Tools::gL()/2 + 3][k_a][k_b] = ph;

            ++ph;

         }

      }

   ++block;

   ph = 0;

   block_char[block][0] = 0;//S
   block_char[block][1] = 0;//K
   block_char[block][2] = 1;//p

   char_block[0][0][1] = block;

   block_char[block + Tools::gL()/2 + 3][0] = 1;//S
   block_char[block + Tools::gL()/2 + 3][1] = 0;//K
   block_char[block + Tools::gL()/2 + 3][2] = 1;//p

   char_block[1][0][1] = block + Tools::gL()/2 + 3;

   //for negative parity: k_a < k_b
   for(int k_a = 1;k_a < Tools::gL();++k_a)
      for(int k_b = k_a + 1;k_b < Tools::gL();++k_b){

         if( (k_a + k_b)%Tools::gL() == 0){

            v[0] = k_a;
            v[1] = k_b;

            ph2s[block].push_back(v);
            ph2s[block + Tools::gL()/2 + 3].push_back(v);

            s2ph[block][k_a][k_b] = ph;
            s2ph[block + Tools::gL()/2 + 3][k_a][k_b] = ph;

            ++ph;

         }

      }

   ++block;

   //now 0 < K < Tools::gL()/2: only keep positive parity blocks
   for(int K = 1;K < Tools::gL()/2;++K){

      ph = 0;

      block_char[block][0] = 0;//S
      block_char[block][1] = K;//K
      block_char[block][2] = 0;//p

      char_block[0][K][0] = block;
      char_block[0][K][1] = block;

      block_char[block + Tools::gL()/2 + 3][0] = 1;//S
      block_char[block + Tools::gL()/2 + 3][1] = K;//K
      block_char[block + Tools::gL()/2 + 3][2] = 0;//p

      char_block[1][K][0] = block + Tools::gL()/2 + 3;
      char_block[1][K][1] = block + Tools::gL()/2 + 3;

      vector<int> v(2);

      for(int k_a = 0;k_a < Tools::gL();++k_a)
         for(int k_b = 0;k_b < Tools::gL();++k_b){

            if( (k_a + k_b)%Tools::gL() == K ){

               v[0] = k_a;
               v[1] = k_b;

               ph2s[block].push_back(v);
               ph2s[block + Tools::gL()/2 + 3].push_back(v);

               s2ph[block][k_a][k_b] = ph;
               s2ph[block + Tools::gL()/2 + 3][k_a][k_b] = ph;

               ++ph;

            }

         }

      ++block;

   }

   //now only K = Tools::gL()/2 is left: first positive parity: 0 <= k_a,k_b <= Tools::gL()/2
   ph = 0;

   block_char[block][0] = 0;//S
   block_char[block][1] = Tools::gL()/2;//K
   block_char[block][2] = 0;//p

   char_block[0][Tools::gL()/2][0] = block;

   block_char[block + Tools::gL()/2 + 3][0] = 1;//S
   block_char[block + Tools::gL()/2 + 3][1] = Tools::gL()/2;//K
   block_char[block + Tools::gL()/2 + 3][2] = 0;//p

   char_block[1][Tools::gL()/2][0] = block + Tools::gL()/2 + 3;

   for(int k_a = 0;k_a <= Tools::gL()/2;++k_a)
      for(int k_b = 0;k_b <= Tools::gL()/2;++k_b){

         if( (k_a + k_b)%Tools::gL() == Tools::gL()/2 ){

            v[0] = k_a;
            v[1] = k_b;

            ph2s[block].push_back(v);
            ph2s[block + Tools::gL()/2 + 3].push_back(v);

            s2ph[block][k_a][k_b] = ph;
            s2ph[block + Tools::gL()/2 + 3][k_a][k_b] = ph;

            ++ph;

         }

      }

   ++block;

   //then negative parity: 0 < k_a,k_b < Tools::gL()/2
   ph = 0;

   block_char[block][0] = 0;//S
   block_char[block][1] = Tools::gL()/2;//K
   block_char[block][2] = 1;//p

   char_block[0][Tools::gL()/2][1] = block;

   block_char[block + Tools::gL()/2 + 3][0] = 1;//S
   block_char[block + Tools::gL()/2 + 3][1] = Tools::gL()/2;//K
   block_char[block + Tools::gL()/2 + 3][2] = 1;//p

   char_block[1][Tools::gL()/2][1] = block + Tools::gL()/2 + 3;

   for(int k_a = 1;k_a < Tools::gL()/2;++k_a)
      for(int k_b = 1;k_b < Tools::gL()/2;++k_b){

         if( (k_a + k_b)%Tools::gL() == Tools::gL()/2 ){

            v[0] = k_a;
            v[1] = k_b;

            ph2s[block].push_back(v);
            ph2s[block + Tools::gL()/2 + 3].push_back(v);

            s2ph[block][k_a][k_b] = ph;
            s2ph[block + Tools::gL()/2 + 3][k_a][k_b] = ph;

            ++ph;

         }

      }

   ++block;

}

/**
 * deallocate the static lists
 */
void PHM::clear(){

   delete [] ph2s;

   for(int B = 0;B < Tools::gL() + 6;++B){

      for(int k = 0;k < Tools::gL();++k)
         delete [] s2ph[B][k];

      delete [] s2ph[B];

   }

   delete [] s2ph;

   for(int B = 0;B < Tools::gL() + 6;++B)
      delete [] block_char[B];

   delete [] block_char;

   for(int S = 0;S < 2;++S){

      for(int K = 0;K <= Tools::gL()/2;++K)
         delete [] char_block[S][K];

      delete [] char_block[S];

   }

   delete [] char_block;

}

/**
 * standard constructor: constructs BlockMatrix object
 */
PHM::PHM() : BlockMatrix(Tools::gL() + 6) {

   //first K = 0: 4 blocks
   //positive parity
   this->setMatrixDim(0,ph2s[0].size(),1);//S = 0
   this->setMatrixDim(Tools::gL()/2 + 3,ph2s[Tools::gL()/2 + 3].size(),3);//S = 1

   //negative parity
   this->setMatrixDim(1,ph2s[1].size(),1);//S = 0
   this->setMatrixDim(Tools::gL()/2 + 4,ph2s[Tools::gL()/2 + 4].size(),3);//S = 1

   //then for 0 < K < Tools::gL()/2
   for(int K = 1;K < Tools::gL()/2;++K){

      this->setMatrixDim(K + 1,ph2s[K + 1].size(),2);//S = 0
      this->setMatrixDim(Tools::gL()/2 + K + 4,ph2s[Tools::gL()/2 + K + 4].size(),6);//S = 1

   }

   //and last for K = Tools::gL()/2: parity positive
   this->setMatrixDim(Tools::gL()/2 + 1,ph2s[Tools::gL()/2 + 1].size(),1);//S = 0
   this->setMatrixDim(Tools::gL() + 4,ph2s[Tools::gL() + 4].size(),3);//S = 1

   //negative parity
   this->setMatrixDim(Tools::gL()/2 + 2,ph2s[Tools::gL()/2 + 2].size(),1);//S = 0
   this->setMatrixDim(Tools::gL() + 5,ph2s[Tools::gL() + 5].size(),3);//S = 1

}

/**
 * copy constructor: constructs BlockMatrix object
 * @param phm_c PHM to be copied into (*this)
 */
PHM::PHM(const PHM &phm_c) : BlockMatrix(phm_c){ }

/**
 * destructor: if counter == 1 the memory for the static lists ph2s en s2ph will be deleted.
 */
PHM::~PHM(){ }

/**
 * access the elements of the matrix in sp mode, 
 * @param B The blockindex of the block you want to access
 * @param k_a first sp momentum index that forms the ph row index i in block B together with k_b
 * @param k_b second sp momentum index that forms the ph row index i in block B together with k_a
 * @param k_c first sp momentum index that forms the ph column index j in block B together with k_d
 * @param k_d second sp momentum index that forms the ph column index j in block B together with k_c
 * @return the number on place PHM(i,j)
 */
double PHM::operator()(int B,int k_a,int k_b,int k_c,int k_d) const{

   return (*this)(block_char[B][0],block_char[B][1],block_char[B][2],k_a,k_b,k_c,k_d);

}

/**
 * access the elements of the matrix in sp mode, 
 * @param S The ph-spin of the block you want to access
 * @param K The ph-momentum of the block you want to access
 * @param p The ph-parity of the block you want to access
 * @param k_a first sp momentum index that forms the ph row index i in block B together with k_b
 * @param k_b second sp momentum index that forms the ph row index i in block B together with k_a
 * @param k_c first sp momentum index that forms the ph column index j in block B together with k_d
 * @param k_d second sp momentum index that forms the ph column index j in block B together with k_c
 * @return the number on place PHM(i,j)
 */
double PHM::operator()(int S,int K,int p,int k_a,int k_b,int k_c,int k_d) const{

   //momentum conservation:
   if( (k_a + k_b)%Tools::gL() != K)
      return 0;

   if( (k_c + k_d)%Tools::gL() != K)
      return 0;

   int K_copy = K;

   int phase_i = get_phase_order(S,K_copy,p,k_a,k_b);

   if(phase_i == 0)
      return 0;

   int phase_j = get_phase_order(S,K,p,k_c,k_d);

   if(phase_j == 0)
      return 0;

   int B = char_block[S][K][p];

   int i = s2ph[B][k_a][k_b];
   int j = s2ph[B][k_c][k_d];

   return phase_i*phase_j*(*this)(B,i,j);

}

/**
 * get the right phase and order of sp indices
 * @param S tp spin
 * @param K tp momentum
 * @param p tp parity
 * @param k_a first sp index
 * @param k_b second sp index
 * @return the phase
 */
int PHM::get_phase_order(int S,int &K,int p,int &k_a,int &k_b){

   int phase = 1;

   if(K == 0){

      //if k_a == 0 or Tools::gL()/2, only positive parity is present.
      if(k_a == 0 || k_a == Tools::gL()/2){

         if(p == 1)
            return 0;

      }
      else if(k_a > k_b){

         int hulp = k_a;
         k_a = k_b;
         k_b = hulp;

         if(p == 1)
            phase *= -1;

      }

   }
   else if(K == Tools::gL()/2){

      //again, if k_a == 0 or Tools::gL()/2, only positive parity is present.
      if(k_a == 0 || k_a == Tools::gL()/2){

         if(p == 1)
            return 0;

      }
      else if(k_a > Tools::gL()/2){//switch

         k_a = Tools::gL() - k_a;
         k_b = Tools::gL() - k_b;

         if(p == 1)
            phase *= -1;

      }

   }
   else if(K > Tools::gL()/2){

      K = Tools::gL() - K;
      k_a = (Tools::gL() - k_a)%Tools::gL();
      k_b = (Tools::gL() - k_b)%Tools::gL();

      if(p == 1)
         phase *= -1;

   }

   return phase;

}

ostream &operator<<(ostream &output,const PHM &phm_p){

   int S,K,p;

   for(int B = 0;B < phm_p.gnr();++B){

      S = phm_p.block_char[B][0];
      K = phm_p.block_char[B][1];
      p = phm_p.block_char[B][2];

      output << "S =\t" << S << "\tK =\t" << K << "\tparity =\t" << p << "\tdimension =\t" << phm_p.gdim(B) << "\tdegeneracy =\t" << phm_p.gdeg(B) << std::endl;
      output << std::endl;

      for(int i = 0;i < phm_p.gdim(B);++i)
         for(int j = 0;j < phm_p.gdim(B);++j){

            output << i << "\t" << j << "\t|\t" << phm_p.ph2s[B][i][0] << "\t" << phm_p.ph2s[B][i][1]

               << "\t" << phm_p.ph2s[B][j][0] << "\t" << phm_p.ph2s[B][j][1] << "\t" << phm_p(B,i,j) << endl;

         }

      output << endl;

   }

   return output;

}

/**
 * The G map, maps a TPM object on a PHM object.
 * @param tpm input TPM
 */
void PHM::G(const TPM &tpm){

   //construct the SPM corresponding to the TPM
   SPM spm(1.0/(Tools::gN() - 1.0),tpm);

   int k_a,k_b,k_c,k_d;
   int k_a_,k_b_,k_d_;

   int S,K,p;
   int K_;

   int psign;

   double ward;

   for(int B = 0;B < gnr();++B){

      S = block_char[B][0];
      K = block_char[B][1];
      p = block_char[B][2];

      if(K == 0 || K == Tools::gL()/2){

         //sign of parity
         psign = 1 - 2*p;

         for(int i = 0;i < gdim(B);++i){

            k_a = ph2s[B][i][0];
            k_b = ph2s[B][i][1];

            //transform k_a and k_b to tpm sp-momentum:
            k_a_ = (Tools::gL() - k_a)%Tools::gL();
            k_b_ = (Tools::gL() - k_b)%Tools::gL();

            for(int j = i;j < gdim(B);++j){

               k_c = ph2s[B][j][0];
               k_d = ph2s[B][j][1];

               //transform k_d to tpm sp-momentum:
               k_d_ = (Tools::gL() - k_d)%Tools::gL();

               //first term: k_a - k_d

               //the K_ is the tp momentum in the tpm
               K_ = (k_a + k_d_)%Tools::gL();

               ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  for(int par = 0;par < 2;++par)
                     ward -= Tools::g6j(0,0,Z,S) * (2.0*Z + 1.0) * tpm(Z,K_,par,k_a,k_d_,k_c,k_b_);

               //now for the norms
               ward /= 2.0 * TPM::norm(K_,k_a,k_d_) * TPM::norm(K_,k_c,k_b_);

               if(k_a == k_d_)
                  ward *= std::sqrt(2.0);

               if(k_b_ == k_c)
                  ward *= std::sqrt(2.0);

               (*this)(B,i,j) = ward;

               //second term: -k_a -k_d
               ward = 0.0;

               K_ = (k_a_ + k_d_)%Tools::gL();

               for(int Z = 0;Z < 2;++Z)
                  for(int par = 0;par < 2;++par)
                     ward -= Tools::g6j(0,0,Z,S) * (2.0*Z + 1.0) * tpm(Z,K_,par,k_a_,k_d_,k_c,k_b);

               //now for the norms
               ward /= 2.0 * TPM::norm(K_,k_a_,k_d_) * TPM::norm(K_,k_c,k_b);

               if(k_a_ == k_d_)
                  ward *= std::sqrt(2.0);

               if(k_b == k_c)
                  ward *= std::sqrt(2.0);

               (*this)(B,i,j) += psign * ward;

               //now the G norm
               (*this)(B,i,j) *= PHM::norm(K,k_a,k_b) * PHM::norm(K,k_c,k_d);

            }

            (*this)(B,i,i) += spm[k_a];

         }

      }
      else{

         for(int i = 0;i < gdim(B);++i){

            k_a = ph2s[B][i][0];

            //transform k_b to tpm sp-momentum:
            k_b_ = (-ph2s[B][i][1] + Tools::gL())%Tools::gL();

            for(int j = i;j < gdim(B);++j){

               k_c = ph2s[B][j][0];

               //transform k_d to tpm sp-momentum:
               k_d_ = (-ph2s[B][j][1] + Tools::gL())%Tools::gL();

               //the K_ is the tp momentum in the tpm
               K_ = (k_a + k_d_)%Tools::gL();

               (*this)(B,i,j) = 0;

               for(int Z = 0;Z < 2;++Z)
                  for(int par = 0;par < 2;++par)
                     (*this)(B,i,j) -= Tools::g6j(0,0,Z,S) * (2.0*Z + 1.0) * tpm(Z,K_,par,k_a,k_d_,k_c,k_b_);

               //now for the norms
               (*this)(B,i,j) /= 4.0 * TPM::norm(K_,k_a,k_d_) * TPM::norm(K_,k_c,k_b_);

               if(k_a == k_d_)
                  (*this)(B,i,j) *= std::sqrt(2.0);

               if(k_b_ == k_c)
                  (*this)(B,i,j) *= std::sqrt(2.0);

            }

            (*this)(B,i,i) += spm[k_a];

         }

      }

   }

   this->symmetrize();

}

double PHM::norm(int K,int k_a,int k_b){

   if(K == 0){

      if(k_a == 0 || k_a == Tools::gL()/2)
         return 0.5;
      else
         return 1.0/std::sqrt(2.0);

   }
   else if(K < Tools::gL()/2)
      return 1.0/std::sqrt(2.0);
   else if(K == Tools::gL()/2){

      if(k_a == 0 || k_a == Tools::gL()/2)
         return 0.5;
      else
         return 1.0/std::sqrt(2.0);

   }
   else
      return 1.0/std::sqrt(2.0);

}

/**
 * The bar function that maps a PPHM object onto a PHM object by tracing away the first pair of incdices of the PPHM
 * @param pphm Input PPHM object
 */
void PHM::bar(const PPHM &pphm){

   int k_a,k_b,k_c,k_d;
   int k_a_,k_b_;

   double ward,hard;

   int psign;

   int Z,K,p;

   int K_pph;

   for(int B = 0;B < gnr();++B){//loop over the blocks PHM

      Z = block_char[B][0];
      K = block_char[B][1];
      p = block_char[B][2];

      psign = 1 - 2*p;

      for(int i = 0;i < gdim(B);++i){

         k_a = ph2s[B][i][0];
         k_b = ph2s[B][i][1];

         k_a_ = (Tools::gL() - k_a)%Tools::gL();
         k_b_ = (Tools::gL() - k_b)%Tools::gL();

         for(int j = i;j < gdim(B);++j){

            k_c = ph2s[B][j][0];
            k_d = ph2s[B][j][1];

            //init
            (*this)(B,i,j) = 0.0;

            if(K == 0 || K == Tools::gL()/2){

               //first the contribution of the S = 1/2 block of the PPHM matrix:
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_de = 0;S_de < 2;++S_de){

                     ward = 0.0;

                     for(int k_l = 0;k_l < Tools::gL();++k_l){

                        K_pph = (k_l + k_a_ + k_b_)%Tools::gL();

                        hard = 0.0;

                        for(int pi = 0;pi < 2;++pi)
                           hard += pphm(0,K_pph,pi,S_ab,k_l,k_a_,k_b_,S_de,k_l,k_c,k_d);

                        if(k_l == k_a_)
                           hard *= std::sqrt(2.0);

                        if(k_l == k_c)
                           hard *= std::sqrt(2.0);

                        ward += 0.5/ ( PPHM::norm(K_pph,k_l,k_a_,k_b_) * PPHM::norm(K_pph,k_l,k_c,k_d) ) * hard;

                     }

                     (*this)(B,i,j) += 2.0 * psign * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) ) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * ward;

                  }

               //then the S = 3/2 contribution: this only occurs when Z = 1
               if(Z == 1){

                  ward = 0.0;

                  for(int k_l = 0;k_l < Tools::gL();++k_l){

                     K_pph = (k_l + k_a_ + k_b_)%Tools::gL();

                     hard = 0.0;

                     for(int pi = 0;pi < 2;++pi)
                        hard += pphm(1,K_pph,pi,1,k_l,k_a_,k_b_,1,k_l,k_c,k_d);

                     ward += 0.5/ ( PPHM::norm(K_pph,k_l,k_a_,k_b_) * PPHM::norm(K_pph,k_l,k_c,k_d) ) * hard;

                  }

                  (*this)(B,i,j) += psign * ward * 4.0/3.0;

               }

            }

            //first the contribution of the S = 1/2 block of the PPHM matrix:
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_de = 0;S_de < 2;++S_de){

                  ward = 0.0;

                  for(int k_l = 0;k_l < Tools::gL();++k_l){

                     K_pph = (k_l + k_a + k_b)%Tools::gL();

                     hard = 0.0;

                     for(int pi = 0;pi < 2;++pi)
                        hard += pphm(0,K_pph,pi,S_ab,k_l,k_a,k_b,S_de,k_l,k_c,k_d);

                     if(k_l == k_a)
                        hard *= std::sqrt(2.0);

                     if(k_l == k_c)
                        hard *= std::sqrt(2.0);

                     ward += 0.5/ ( PPHM::norm(K_pph,k_l,k_a,k_b) * PPHM::norm(K_pph,k_l,k_c,k_d) ) * hard;

                  }

                  (*this)(B,i,j) += 2.0 * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) ) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * ward;

               }

            //then the S = 3/2 contribution: this only occurs when Z = 1
            if(Z == 1){

               ward = 0.0;

               for(int k_l = 0;k_l < Tools::gL();++k_l){

                  K_pph = (k_l + k_a + k_b)%Tools::gL();

                  hard = 0.0;

                  for(int pi = 0;pi < 2;++pi)
                     hard += pphm(1,K_pph,pi,1,k_l,k_a,k_b,1,k_l,k_c,k_d);

                  ward += 0.5/ ( PPHM::norm(K_pph,k_l,k_a,k_b) * PPHM::norm(K_pph,k_l,k_c,k_d) ) * hard;

               }

               (*this)(B,i,j) += ward * 4.0/3.0;

            }

            (*this)(B,i,j) *= PHM::norm(K,k_a,k_b) * PHM::norm(K,k_c,k_d);

         }
      }

   }

   this->symmetrize();

}
