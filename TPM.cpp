#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::cout;
using std::endl;
using std::ios;

#include "include.h"

//def of some statics
vector< vector<int> > *TPM::t2s;
int ***TPM::s2t;

int **TPM::block_char;
int ***TPM::char_block;

/**
 * static function that initializes the static variables and allocates and fill the static lists
 */
void TPM::init(){

   int nr_B = Tools::gL() + 4;//nr of blocks

   //allocate
   t2s = new vector< vector<int> > [nr_B];

   s2t = new int ** [nr_B];

   for(int B = 0;B < nr_B;++B){

      s2t[B] = new int * [Tools::gL()];

      for(int k_a = 0;k_a < Tools::gL();++k_a)
         s2t[B][k_a] = new int [Tools::gL()];

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

   int block = 0;

   int tp;

   vector<int> v(2);

   //loop over the pseudo-momentum

   //first the K = 0 blocks:

   //first S = 0 (only parity = +1)
   block_char[block][0] = 0;
   block_char[block][1] = 0;
   block_char[block][2] = 0;//parity baby! (0 means +1, 1 means -1)

   char_block[0][0][0] = block;

   tp = 0;

   for(int k_a = 0;k_a < Tools::gL();++k_a)
      for(int k_b = k_a;k_b < Tools::gL();++k_b){

         if( (k_a + k_b)%Tools::gL() == 0 ){

            v[0] = k_a;
            v[1] = k_b;

            t2s[block].push_back(v);

            s2t[block][k_a][k_b]  = tp;
            s2t[block][k_b][k_a]  = tp;

            ++tp;

         }

      }

   //then S = 1: block shifted with Tools::gL()/2 + 2 (only parity = -1)
   block_char[block + Tools::gL()/2 + 2][0] = 1;
   block_char[block + Tools::gL()/2 + 2][1] = 0;
   block_char[block + Tools::gL()/2 + 2][2] = 1;//parity baby! (0 means +1, 1 means -1)

   char_block[1][0][1] = block + Tools::gL()/2 + 2;

   tp = 0;

   for(int k_a = 0;k_a < Tools::gL();++k_a)
      for(int k_b = k_a + 1;k_b < Tools::gL();++k_b){

         if( (k_a + k_b)%Tools::gL() == 0 ){

            v[0] = k_a;
            v[1] = k_b;

            t2s[block + Tools::gL()/2 + 2].push_back(v);

            s2t[block + Tools::gL()/2 + 2][k_a][k_b]  = tp;
            s2t[block + Tools::gL()/2 + 2][k_b][k_a]  = tp;

            ++tp;

         }

      }

   ++block;


   //then 0 < K < Tools::gL()/2
   for(int K = 1;K < Tools::gL()/2;++K){

      //first S = 0
      block_char[block][0] = 0;
      block_char[block][1] = K;
      block_char[block][2] = 0;//parity baby! (0 means +1, 1 means -1)

      //both parities refer to the same block
      char_block[0][K][0] = block;
      char_block[0][K][1] = block;

      tp = 0;

      for(int k_a = 0;k_a < Tools::gL();++k_a)
         for(int k_b = k_a;k_b < Tools::gL();++k_b){

            if( (k_a + k_b)%Tools::gL() == K ){

               v[0] = k_a;
               v[1] = k_b;

               t2s[block].push_back(v);

               s2t[block][k_a][k_b]  = tp;
               s2t[block][k_b][k_a]  = tp;

               ++tp;

            }

         }

      //then S = 1: block shifted with Tools::gL()/2 + 2
      block_char[block + Tools::gL()/2 + 2][0] = 1;
      block_char[block + Tools::gL()/2 + 2][1] = K;
      block_char[block + Tools::gL()/2 + 2][2] = 0;//parity baby! (0 means +1, 1 means -1)

      //both parities refer to the same block
      char_block[1][K][0] = block + Tools::gL()/2 + 2;
      char_block[1][K][1] = block + Tools::gL()/2 + 2;

      tp = 0;

      for(int k_a = 0;k_a < Tools::gL();++k_a)
         for(int k_b = k_a + 1;k_b < Tools::gL();++k_b){

            if( (k_a + k_b)%Tools::gL() == K ){

               v[0] = k_a;
               v[1] = k_b;

               t2s[block + Tools::gL()/2 + 2].push_back(v);

               s2t[block + Tools::gL()/2 + 2][k_a][k_b]  = tp;
               s2t[block + Tools::gL()/2 + 2][k_b][k_a]  = tp;

               ++tp;

            }

         }

      ++block;

   }

   //now only the L/2 pseudo momentum block is left.
   //The sp-momenta here only go from 0 -> L/2

   //first postive parity:
   //first S = 0;
   block_char[block][0] = 0;
   block_char[block][1] = Tools::gL()/2;
   block_char[block][2] = 0;//parity baby! (0 means +1, 1 means -1)

   char_block[0][Tools::gL()/2][0] = block;

   tp = 0;

   //the momenta only go to L/2
   for(int k_a = 0;k_a <= Tools::gL()/2;++k_a)
      for(int k_b = k_a;k_b <= Tools::gL()/2;++k_b){

         if( (k_a + k_b)%Tools::gL() == Tools::gL()/2 ){

            v[0] = k_a;
            v[1] = k_b;

            t2s[block].push_back(v);

            s2t[block][k_a][k_b]  = tp;
            s2t[block][k_b][k_a]  = tp;

            ++tp;

         }

      }

   //then S = 1: block shifted with Tools::gL()/2 + 2
   block_char[block + Tools::gL()/2 + 2][0] = 1;
   block_char[block + Tools::gL()/2 + 2][1] = Tools::gL()/2;
   block_char[block + Tools::gL()/2 + 2][2] = 0;//parity baby! (0 means +1, 1 means -1)

   char_block[1][Tools::gL()/2][0] = block + Tools::gL()/2 + 2;

   tp = 0;

   for(int k_a = 0;k_a <= Tools::gL()/2;++k_a)
      for(int k_b = k_a + 1;k_b <= Tools::gL()/2;++k_b){

         if( (k_a + k_b)%Tools::gL() == Tools::gL()/2 ){

            v[0] = k_a;
            v[1] = k_b;

            t2s[block + Tools::gL()/2 + 2].push_back(v);

            s2t[block + Tools::gL()/2 + 2][k_a][k_b]  = tp;
            s2t[block + Tools::gL()/2 + 2][k_b][k_a]  = tp;

            ++tp;

         }

      }

   ++block;

   //then negative parity: this block hasn't got the |0 L/2> basisvector

   //first S = 0;
   block_char[block][0] = 0;
   block_char[block][1] = Tools::gL()/2;
   block_char[block][2] = 1;//parity baby! (0 means +1, 1 means -1)

   char_block[0][Tools::gL()/2][1] = block;

   tp = 0;

   //the only difference is that k_a starts from 1, and the momenta only go to Tools::gL()/2
   for(int k_a = 1;k_a < Tools::gL()/2;++k_a)
      for(int k_b = k_a;k_b < Tools::gL()/2;++k_b){

         if( (k_a + k_b)%Tools::gL() == Tools::gL()/2 ){

            v[0] = k_a;
            v[1] = k_b;

            t2s[block].push_back(v);

            s2t[block][k_a][k_b]  = tp;
            s2t[block][k_b][k_a]  = tp;

            ++tp;

         }

      }

   //then S = 1: block shifted with Tools::gL()/2 + 2
   block_char[block + Tools::gL()/2 + 2][0] = 1;
   block_char[block + Tools::gL()/2 + 2][1] = Tools::gL()/2;
   block_char[block + Tools::gL()/2 + 2][2] = 1;//parity baby! (0 means +1, 1 means -1)

   char_block[1][Tools::gL()/2][1] = block + Tools::gL()/2 + 2;

   tp = 0;

   for(int k_a = 1;k_a < Tools::gL()/2;++k_a)
      for(int k_b = k_a + 1;k_b < Tools::gL()/2;++k_b){

         if( (k_a + k_b)%Tools::gL() == Tools::gL()/2 ){

            v[0] = k_a;
            v[1] = k_b;

            t2s[block + Tools::gL()/2 + 2].push_back(v);

            s2t[block + Tools::gL()/2 + 2][k_a][k_b]  = tp;
            s2t[block + Tools::gL()/2 + 2][k_b][k_a]  = tp;

            ++tp;

         }

      }

}

/**
 * deallocates the statics lists
 */
void TPM::clear(){

   delete [] t2s;

   for(int B = 0;B < Tools::gL() + 4;++B){

      for(int k_a = 0;k_a < Tools::gL();++k_a)
         delete [] s2t[B][k_a];

      delete [] s2t[B];

   }

   delete [] s2t;

   for(int B = 0;B < Tools::gL() + 4;++B)
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
 * standard constructor for a spinsymmetrical, translationally invariant tp matrix, with parity (k <--> -k) taken into account: 
 * constructs BlockMatrix object with Tools::gL() + 4 blocks, Tools::gL()/2 + 2 for S = 0 and S = 1,
 */
TPM::TPM() : BlockMatrix(Tools::gL() + 4) {

   //K = 0
   this->setMatrixDim(0,t2s[0].size(),1);
   this->setMatrixDim(Tools::gL()/2 + 2,t2s[Tools::gL()/2 + 2].size(),3);

   //0 < K < Tools::gL()/2
   for(int B = 1;B < Tools::gL()/2;++B){

      this->setMatrixDim(B,t2s[B].size(),2);
      this->setMatrixDim(B + Tools::gL()/2 + 2,t2s[B + Tools::gL()/2 + 2].size(),6);

   }

   //K = Tools::gL()/2: 4 more blocks (positive and negative parity in S = 0/1)
   this->setMatrixDim(Tools::gL()/2,t2s[Tools::gL()/2].size(),1);
   this->setMatrixDim(Tools::gL()/2 + 1,t2s[Tools::gL()/2 + 1].size(),1);

   this->setMatrixDim(Tools::gL() + 2,t2s[Tools::gL() + 2].size(),3);
   this->setMatrixDim(Tools::gL() + 3,t2s[Tools::gL() + 3].size(),3);

}

/**
 * copy constructor for a spinsymmetrical, translationally invariant tp matrix, with parity (k <--> -k) taken into account: 
 * constructs BlockMatrix object with Tools::gL() + 4 blocks, Tools::gL()/2 + 2 for S = 0 and S = 1,
 * @param tpm_c The TPM object to be copied into (*this)
 */
TPM::TPM(const TPM &tpm_c) : BlockMatrix(tpm_c){ }

/**
 * destructor
 */
TPM::~TPM(){ }

/**
 * access the elements of the the blocks in sp mode, the symmetry or antisymmetry of the blocks is automatically accounted for:\n\n
 * Antisymmetrical for S = 1, symmetrical in the sp orbitals for S = 0\n\n
 * @param B The blockindex
 * @param k_a first sp momentum index that forms the tp row index i of block B, together with k_b
 * @param k_b second sp momentum index that forms the tp row index i of block B, together with k_a
 * @param k_c first sp momentum index that forms the tp column index j of block B, together with k_d
 * @param k_d second sp momentum index that forms the tp column index j of block B, together with k_c
 * @return the number on place TPM(B,i,j) with the right phase.
 */

double TPM::operator()(int B,int k_a,int k_b,int k_c,int k_d) const{

   return (*this)(block_char[B][0],block_char[B][1],block_char[B][2],k_a,k_b,k_c,k_d);

}

/**
 * access the elements of the the blocks in sp mode, the symmetry or antisymmetry of the blocks is automatically accounted for:\n\n
 * Antisymmetrical for S = 1, symmetrical in the sp orbitals for S = 0\n\n
 * @param S The two-particle spinquantumnumber
 * @param K The two-particle momentum
 * @param p The two-particle parity
 * @param k_a first sp momentum index that forms the tp row index i of block B(S,K,p), together with k_b
 * @param k_b second sp momentum index that forms the tp row index i of block B(S,K,p), together with k_a
 * @param k_c first sp momentum index that forms the tp column index j of block B(S,K,p), together with k_d
 * @param k_d second sp momentum index that forms the tp column index j of block B(S,K,p), together with k_c
 * @return the number on place TPM(B,i,j) with the right phase.
 */
double TPM::operator()(int S,int K,int p,int k_a,int k_b,int k_c,int k_d) const{

   //momentum checks out
   if( (k_a + k_b)%Tools::gL() != K)
      return 0;

   if( (k_c + k_d)%Tools::gL() != K)
      return 0;

   int copy_K = K;

   int phase_i = get_phase_order(S,K,p,k_a,k_b);

   if(phase_i == 0)
      return 0;

   int phase_j = get_phase_order(S,copy_K,p,k_c,k_d);

   if(phase_j == 0)
      return 0;

   int B = char_block[S][K][p];

   int i = s2t[B][k_a][k_b];
   int j = s2t[B][k_c][k_d];

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
int TPM::get_phase_order(int S,int &K,int p,int &k_a,int &k_b){

   int phase = 1;

   //for the K = 0 , S = 0 blocks only one parity type is present
   if(K == 0){

      if(S == 0){

         if(p == 1)
            return 0;

      }
      else{

         if(p == 0)
            return 0;

      }

   }
   else if(K == Tools::gL()/2){//for K = L/2 and k_a k_b = 0 L/2 , only positive parity is present

      if(k_a == 0 || k_b == 0){

         if(p == 1)
            return 0;

      }
      else if(k_a > Tools::gL()/2){

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

   if(S == 1){

      if(k_a > k_b)
         phase *= -1;
      else if(k_a == k_b)
         return 0;

   }

   return phase;

}

ostream &operator<<(ostream &output,const TPM &tpm_p){

   int S,K,p;

   for(int B = 0;B < tpm_p.gnr();++B){

      S = tpm_p.block_char[B][0];
      K = tpm_p.block_char[B][1];
      p = tpm_p.block_char[B][2];

      output << "S =\t" << S << "\tK =\t" << K << "\tp = \t" << p << "\tdimension =\t" << tpm_p.gdim(B) << "\tdegeneracy =\t" << tpm_p.gdeg(B) << std::endl;
      output << std::endl;

      for(int i = 0;i < tpm_p.gdim(B);++i)
         for(int j = 0;j < tpm_p.gdim(B);++j){

            output << i << "\t" << j << "\t|\t" << tpm_p.t2s[B][i][0] << "\t" << tpm_p.t2s[B][i][1]

               << "\t" << tpm_p.t2s[B][j][0] << "\t" << tpm_p.t2s[B][j][1] << "\t" << tpm_p(B,i,j) << endl;

         }

      output << std::endl;

   }

   return output;

}

/**
 * Output to a file with no parity included.
 * @param filename name and location of the file you want to print to.
 */
void TPM::out_sp(const char *filename) const{

   ofstream output(filename);
   output.precision(15);

   for(int B = 0;B < gnr();++B){

      for(int i = 0;i < gdim(B);++i)
         for(int j = 0;j < gdim(B);++j)
            output << block_char[B][0] << "\t" << block_char[B][1] << "\t" << block_char[B][2] << "\t" << t2s[B][i][0] << "\t" << t2s[B][i][1] 

               << "\t" << t2s[B][j][0] << "\t" << t2s[B][j][1] << "\t" << (*this)(B,i,j) << endl;

   }

}

/**
 * construct the spinsymmetrical hubbard hamiltonian in momentum space with on site repulsion U
 * @param U onsite repulsion term
 */
void TPM::hubbard(double U){

   int k_a,k_b,k_c,k_d;//sp momentum 

   double ward = 1.0/(Tools::gN() - 1.0);

   int S,K,p;

   for(int B = 0;B < gnr();++B){

      S = block_char[B][0];
      K = block_char[B][1];
      p = block_char[B][2];

      for(int i = 0;i < gdim(B);++i){

         k_a = t2s[B][i][0];
         k_b = t2s[B][i][1];

         for(int j = i;j < gdim(B);++j){

            k_c = t2s[B][j][0];
            k_d = t2s[B][j][1];

            //init
            (*this)(B,i,j) = 0;

            //hopping (kinetic energy):
            if(i == j)
               (*this)(B,i,i) = -2.0 * ward * ( cos( 2.0 * k_a * 3.141592653589793238462 / (double) Tools::gL())
               
                     + cos( 2.0 * k_b * 3.141592653589793238462 / (double) Tools::gL()) );

            //on-site repulsion
            if(S == 0){

               double ward = 2.0*U / (double) Tools::gL();

               if(k_a == k_b)
                  ward /= std::sqrt(2.0);

               if(k_c == k_d)
                  ward /= std::sqrt(2.0);

               if(K != Tools::gL()/2)
                  (*this)(B,i,j) += ward;
               else//K = Tools::gL()/2
                  if(p == 0){//positive parity

                     if(k_a == 0)
                        ward /= std::sqrt(2.0);
                     if(k_c == 0)
                        ward /= std::sqrt(2.0);

                     (*this)(B,i,j) += 2.0*ward;

                  }

            }

         }
      }

   }

   this->symmetrize();

}

/**
 * The spincoupled Q map
 * @param option = 1, regular Q map , = -1 inverse Q map
 * @param tpm_d the TPM of which the Q map is taken and saved in this.
 */
void TPM::Q(int option,const TPM &tpm_d){

   double a = 1;
   double b = 1.0/(Tools::gN()*(Tools::gN() - 1.0));
   double c = 1.0/(Tools::gN() - 1.0);

   this->Q(option,a,b,c,tpm_d);

}

/**
 * The spincoupled Q-like map: see primal-dual.pdf for more info (form: Q^S(A,B,C)(TPM) )
 * @param option = 1, regular Q-like map , = -1 inverse Q-like map
 * @param A factor in front of the two particle piece of the map
 * @param B factor in front of the no particle piece of the map
 * @param C factor in front of the single particle piece of the map
 * @param tpm_d the TPM of which the Q-like map is taken and saved in this.
 */
void TPM::Q(int option,double A,double B,double C,const TPM &tpm_d){

   //for inverse
   if(option == -1){

      B = (B*A + 2.0 * B * C * Tools::gL() - 2.0*C*C)/( A * (C*(2.0 * Tools::gL() - 2.0) -  A) 
            * ( A + 2.0*B*Tools::gL()*(2*Tools::gL() - 1.0) - 2.0*C*(2.0*Tools::gL() - 1.0) ) );

      C = C/(A*(C*(2.0*Tools::gL() - 2.0) - A));
      A = 1.0/A;

   }

   SPM spm(C,tpm_d);

   //de trace*2 omdat mijn definitie van trace in berekeningen over alle (alpha,beta) loopt
   double ward = B*tpm_d.trace()*2.0;

   int k_a,k_b;

   for(int B = 0;B < gnr();++B){

      for(int i = 0;i < gdim(B);++i){

         k_a = t2s[B][i][0];
         k_b = t2s[B][i][1];

         //tp part is only nondiagonal part
         for(int j = i;j < gdim(B);++j)
            (*this)(B,i,j) = A * tpm_d(B,i,j);

         (*this)(B,i,i) += ward - spm[k_a] - spm[k_b];

      }

   }

   this->symmetrize();

}

/**
 * initialize this onto the unitmatrix with trace nr of pairs 
 */
void TPM::unit(){

   double ward = Tools::gN()*(Tools::gN() - 1.0)/(2.0*Tools::gL()*(2.0*Tools::gL() - 1.0));

   for(int B = 0;B < gnr();++B){

      for(int i = 0;i < gdim(B);++i){

         (*this)(B,i,i) = ward;

         for(int j = i + 1;j < gdim(B);++j)
            (*this)(B,i,j) = (*this)(B,j,i) = 0.0;

      }
   }

}

/**
 * orthogonal projection onto the space of traceless matrices
 */
void TPM::proj_Tr(){

   double ward = (2.0 * this->trace())/(2.0*Tools::gL()*(2.0*Tools::gL() - 1));

   this->min_unit(ward);

}

/**
 * Deduct the unitmatrix times a constant (scale) from this.\n\n
 * this -= scale* 1
 * @param scale the constant
 */

void TPM::min_unit(double scale){

   for(int B = 0;B < gnr();++B)
      for(int i = 0;i < gdim(B);++i)
         (*this)(B,i,i) -= scale;

}

/**
 * Collaps a SUP matrix S onto a TPM matrix like this:\n\n
 * sum_i Tr (S u^i)f^i = this
 * @param option = 0, project onto full symmetric matrix space, = 1 project onto traceless symmetric matrix space
 * @param S input SUP
 */
void TPM::collaps(int option,const SUP &S){

   *this = S.tpm(0);

   TPM hulp;

   hulp.Q(1,S.tpm(1));

   *this += hulp;

#ifdef __G_CON

   hulp.G(S.phm());

   *this += hulp;

#endif

#ifdef __T1_CON

   hulp.T(S.dpm());

   *this += hulp;

#endif

#ifdef __T2_CON

   hulp.T(S.pphm());

   *this += hulp;

#endif

   if(option == 1)
      this->proj_Tr();

}

/**
 * @return The expectation value of the total spin for the TPM.
 */
double TPM::spin() const{

   double ward = 0.0;

   int S;

   for(int B = 0;B < gnr();++B){

      S = block_char[B][0];

      if(S == 0){

         for(int i = 0;i < gdim(B);++i)
            ward += -1.5 * (Tools::gN() - 2.0)/(Tools::gN() - 1.0) * gdeg(B) * (*this)(B,i,i);

      }
      else{

         for(int i = 0;i < this->gdim(B);++i)
            ward += gdeg(B) * ( -1.5 * (Tools::gN() - 2.0)/(Tools::gN() - 1.0) + 2.0 ) * (*this)(B,i,i);

      }

   }

   return ward;

}

/**
 * Fill a TPM object from a file.
 * @param input The ifstream object, corresponding to the file containing the TPM
 */
void TPM::in(ifstream &input){

   double block,dim,deg;
   int I,J;

   for(int B = 0;B < gnr();++B){

      input >> block >> dim >> deg;

      for(int i = 0;i < gdim(B);++i)
         for(int j = 0;j < gdim(B);++j)
            input >> I >> J >> (*this)(B,i,j);

   }

}

/**
 * @param K the tp momentum
 * @param k_a the first sp momentum index
 * @param k_b the second sp momentum index
 * @return the norm of the wavefunction.
 */
double TPM::norm(int K,int k_a,int k_b){

   if(K == 0)
      return 0.5;
   else if(K < Tools::gL()/2)
      return 1.0/std::sqrt(2.0);
   else if(K == Tools::gL()/2){

      if(k_a == 0 || k_b == 0)
         return 0.5;
      else
         return 1.0/std::sqrt(2.0);

   }
   else
      return 1.0/std::sqrt(2.0);

}

/**
 * The G down map, maps a PHM object onto a TPM object using the G map.
 * @param phm input PHM
 */
void TPM::G(const PHM &phm){

   double ward;

   SPM spm(1.0/(Tools::gN() - 1.0),phm);

   int k_a,k_b,k_c,k_d;
   int k_a_,k_b_,k_c_,k_d_;

   int S,K,K_ph,p;

   int sign,psign;

   for(int B = 0;B < gnr();++B){

      S = block_char[B][0];
      K = block_char[B][1];
      p = block_char[B][2];

      sign = 1 - 2*S;
      psign = 1 - 2*p;

      for(int i = 0;i < gdim(B);++i){

         k_a = t2s[B][i][0];
         k_b = t2s[B][i][1];

         k_a_ = (Tools::gL() - k_a)%Tools::gL();
         k_b_ = (Tools::gL() - k_b)%Tools::gL();

         //tp part is only nondiagonal part
         for(int j = i;j < gdim(B);++j){

            k_c = t2s[B][j][0];
            k_d = t2s[B][j][1];

            k_c_ = (Tools::gL() - k_c)%Tools::gL();
            k_d_ = (Tools::gL() - k_d)%Tools::gL();

            (*this)(B,i,j) = 0.0;

            //four ph exchange terms:
            if(K == 0 || K == Tools::gL()/2){

               //1)
               K_ph = (k_a_ + k_d_)%Tools::gL();

               ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  for(int pi = 0;pi < 2;++pi)
                     ward -= (2.0*Z + 1.0) * Tools::g6j(0,0,S,Z) * phm(Z,K_ph,pi,k_a_,k_d_,k_c,k_b);

               ward /= 2.0 * PHM::norm(K_ph,k_a_,k_d_) * PHM::norm(K_ph,k_c,k_b);

               (*this)(B,i,j) += psign*ward;

               //2)
               K_ph = (k_b_ + k_c_)%Tools::gL();

               ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  for(int pi = 0;pi < 2;++pi)
                     ward -= (2.0*Z + 1.0) * Tools::g6j(0,0,S,Z) * phm(Z,K_ph,pi,k_b_,k_c_,k_d,k_a);

               ward /= 2.0 * PHM::norm(K_ph,k_b_,k_c_) * PHM::norm(K_ph,k_d,k_a);

               (*this)(B,i,j) += psign*ward;

               //3)
               K_ph = (k_b_ + k_d_)%Tools::gL();

               ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  for(int pi = 0;pi < 2;++pi)
                     ward -= (2.0*Z + 1.0) * Tools::g6j(0,0,S,Z) * phm(Z,K_ph,pi,k_b_,k_d_,k_c,k_a);

               ward /= 2.0 * PHM::norm(K_ph,k_b_,k_d_) * PHM::norm(K_ph,k_c,k_a);

               (*this)(B,i,j) += psign*sign*ward;

               //4)
               K_ph = (k_a_ + k_c_)%Tools::gL();

               ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  for(int pi = 0;pi < 2;++pi)
                     ward -= (2.0*Z + 1.0) * Tools::g6j(0,0,S,Z) * phm(Z,K_ph,pi,k_a_,k_c_,k_d,k_b);

               ward /= 2.0 * PHM::norm(K_ph,k_a_,k_c_) * PHM::norm(K_ph,k_d,k_b);

               (*this)(B,i,j) += psign*sign*ward;

            }

            //four regular ph terms:
            //1)
            K_ph = (k_a + k_d_)%Tools::gL();

            ward = 0.0;

            for(int Z = 0;Z < 2;++Z)
               for(int pi = 0;pi < 2;++pi)
                  ward -= (2.0*Z + 1.0) * Tools::g6j(0,0,S,Z) * phm(Z,K_ph,pi,k_a,k_d_,k_c,k_b_);

            ward /= 2.0 * PHM::norm(K_ph,k_a,k_d_) * PHM::norm(K_ph,k_c,k_b_);

            (*this)(B,i,j) += ward;

            //2)
            K_ph = (k_b + k_c_)%Tools::gL();

            ward = 0.0;

            for(int Z = 0;Z < 2;++Z)
               for(int pi = 0;pi < 2;++pi)
                  ward -= (2.0*Z + 1.0) * Tools::g6j(0,0,S,Z) * phm(Z,K_ph,pi,k_b,k_c_,k_d,k_a_);

            ward /= 2.0 * PHM::norm(K_ph,k_b,k_c_) * PHM::norm(K_ph,k_d,k_a_);

            (*this)(B,i,j) += ward;

            //3)
            K_ph = (k_b + k_d_)%Tools::gL();

            ward = 0.0;

            for(int Z = 0;Z < 2;++Z)
               for(int pi = 0;pi < 2;++pi)
                  ward -= (2.0*Z + 1.0) * Tools::g6j(0,0,S,Z) * phm(Z,K_ph,pi,k_b,k_d_,k_c,k_a_);

            ward /= 2.0 * PHM::norm(K_ph,k_b,k_d_) * PHM::norm(K_ph,k_c,k_a_);

            (*this)(B,i,j) += sign*ward;

            //4)
            K_ph = (k_a + k_c_)%Tools::gL();

            ward = 0.0;

            for(int Z = 0;Z < 2;++Z)
               for(int pi = 0;pi < 2;++pi)
                  ward -= (2.0*Z + 1.0) * Tools::g6j(0,0,S,Z) * phm(Z,K_ph,pi,k_a,k_c_,k_d,k_b_);

            ward /= 2.0 * PHM::norm(K_ph,k_a,k_c_) * PHM::norm(K_ph,k_d,k_b_);

            (*this)(B,i,j) += sign*ward;

            //tp-norm:
            (*this)(B,i,j) *= TPM::norm(K,k_a,k_b) * TPM::norm(K,k_c,k_d);

            if(k_a == k_b)
               (*this)(B,i,j) /= std::sqrt(2.0);

            if(k_c == k_d)
               (*this)(B,i,j) /= std::sqrt(2.0);

         }

         //sp term is diagonal
         (*this)(B,i,i) += spm[k_a] + spm[k_b];

      }

   }

   this->symmetrize();

}

/**
 * Construct a spincoupled, translationally invariant and parity symmetric TPM matrix out of a spincoupled, translationally invariant and parity symmetric DPM matrix.
 * For the definition and derivation see symmetry.pdf
 * @param dpm input DPM
 */
void TPM::bar(const DPM &dpm){

   int k_a,k_b,k_c,k_d;
   int k_c_,k_d_;

   double ward,hard;

   int K,K_dp,p,psign;

   //first the S = 0 part, easiest:
   for(int B = 0;B < Tools::gL()/2 + 2;++B){

      K = block_char[B][1];
      p = block_char[B][2];

      psign = 1 - 2*p;

      for(int i = 0;i < gdim(B);++i){

         k_a = t2s[B][i][0];
         k_b = t2s[B][i][1];

         for(int j = i;j < gdim(B);++j){

            k_c = t2s[B][j][0];
            k_d = t2s[B][j][1];

            k_c_ = (Tools::gL() - k_c)%Tools::gL(); 
            k_d_ = (Tools::gL() - k_d)%Tools::gL();

            (*this)(B,i,j) = 0.0;

            //only total S = 1/2 can remain because cannot couple to S = 3/2 with intermediate S = 0
            if(K == 0 || K == Tools::gL()/2){

               for(int k = 0;k < Tools::gL();++k){

                  K_dp = (k_a + k_b + k)%Tools::gL();

                  hard = 0.0;

                  for(int pi = 0;pi < 2;++pi)
                     hard += dpm(0,K_dp,pi,0,k_a,k_b,k,0,k_c_,k_d_,k) / ( DPM::norm(0,K_dp,pi,0,k_a,k_b,k) * DPM::norm(0,K_dp,pi,0,k_c_,k_d_,k) );

                  (*this)(B,i,j) += psign * hard;

               }

            }

            for(int k = 0;k < Tools::gL();++k){

               K_dp = (k_a + k_b + k)%Tools::gL();

               hard = 0.0;

               for(int pi = 0;pi < 2;++pi)
                  hard += dpm(0,K_dp,pi,0,k_a,k_b,k,0,k_c,k_d,k)/ ( DPM::norm(0,K_dp,pi,0,k_a,k_b,k) * DPM::norm(0,K_dp,pi,0,k_c,k_d,k) );

               (*this)(B,i,j) += hard;

            }

            (*this)(B,i,j) *= TPM::norm(K,k_a,k_b) * TPM::norm(K,k_c,k_d);

         }
      }

   }

   //then the S = 1 part:
   for(int B = Tools::gL()/2 + 2;B < gnr();++B){

      K = block_char[B][1];
      p = block_char[B][2];

      psign = 1 - 2*p;

      for(int i = 0;i < gdim(B);++i){

         k_a = t2s[B][i][0];
         k_b = t2s[B][i][1];

         for(int j = i;j < gdim(B);++j){

            k_c = t2s[B][j][0];
            k_d = t2s[B][j][1];

            k_c_ = (Tools::gL() - k_c)%Tools::gL();
            k_d_ = (Tools::gL() - k_d)%Tools::gL();

            (*this)(B,i,j) = 0.0;

            if(K == 0 || K == Tools::gL()/2){

               for(int Z = 0;Z < 2;++Z){//loop over the dpm blocks: S = 1/2 and 3/2 = Z + 1/2

                  ward = 0.0;

                  for(int k = 0;k < Tools::gL();++k){

                     K_dp = (k + k_a + k_b)%Tools::gL();

                     hard = 0.0;

                     for(int pi = 0;pi < 2;++pi)
                        hard += dpm(Z,K_dp,pi,1,k_a,k_b,k,1,k_c_,k_d_,k) / (DPM::norm(Z,K_dp,pi,1,k_a,k_b,k) * DPM::norm(Z,K_dp,pi,1,k_c_,k_d_,k) );

                     ward += 0.5 * psign * hard;

                  }

                  (*this)(B,i,j) += ward * (2 * (Z + 0.5) + 1.0)/3.0;

               }

            }

            for(int Z = 0;Z < 2;++Z){//loop over the dpm blocks: S = 1/2 and 3/2 = Z + 1/2

               ward = 0.0;

               for(int k = 0;k < Tools::gL();++k){

                  K_dp = (k + k_a + k_b)%Tools::gL();

                  hard = 0.0;

                  for(int pi = 0;pi < 2;++pi)
                     hard += dpm(Z,K_dp,pi,1,k_a,k_b,k,1,k_c,k_d,k) / (DPM::norm(Z,K_dp,pi,1,k_a,k_b,k) * DPM::norm(Z,K_dp,pi,1,k_c,k_d,k) );

                  ward += 0.5 * hard;

               }

               (*this)(B,i,j) += ward * (2 * (Z + 0.5) + 1.0)/3.0;

            }

            (*this)(B,i,j) *= TPM::norm(K,k_a,k_b) * TPM::norm(K,k_c,k_d);

         }
      }

   }

   this->symmetrize();

}

/** 
 * The T1-down map that maps a DPM on TPM. This is just a Q-like map using the TPM::bar (dpm) as input.
 * @param dpm the input DPM matrix
 */
void TPM::T(const DPM &dpm){

   TPM tpm;
   tpm.bar(dpm);

   double a = 1;
   double b = 1.0/(3.0*Tools::gN()*(Tools::gN() - 1.0));
   double c = 0.5/(Tools::gN() - 1.0);

   this->Q(1,a,b,c,tpm);

}

/**
 * The bar function that maps a PPHM object onto a TPM object by tracing away the last pair of incdices of the PPHM
 * @param pphm Input PPHM object
 */
void TPM::bar(const PPHM &pphm){

   int k_a,k_b,k_c,k_d;

   int k_a_,k_b_;

   int Z,K,p;
   int psign;

   int K_pph;

   double ward,hard;

   for(int B = 0;B < gnr();++B){//loop over the tp blocks

      Z = block_char[B][0];//spin of the TPM - block
      K = block_char[B][1];//momentum of the TPM - block
      p = block_char[B][2];//parity of the TPM - block

      psign = 1 - 2*p;

      for(int i = 0;i < gdim(B);++i){

         k_a = t2s[B][i][0];
         k_b = t2s[B][i][1];

         k_a_ = (Tools::gL() - k_a)%Tools::gL();
         k_b_ = (Tools::gL() - k_b)%Tools::gL();

         for(int j = i;j < gdim(B);++j){

            k_c = t2s[B][j][0];
            k_d = t2s[B][j][1];

            (*this)(B,i,j) = 0.0;

            for(int S = 0;S < 2;++S){//loop over three particle spin: 1/2 and 3/2

               ward = 0.0;

               if(K == 0 || K == Tools::gL()/2){

                  for(int k_l = 0;k_l < Tools::gL();++k_l){

                     hard = 0.0;

                     K_pph = (k_a_ + k_b_ + k_l)%Tools::gL();

                     for(int pi = 0;pi < 2;++pi)
                        hard += pphm(S,K_pph,pi,Z,k_a_,k_b_,k_l,Z,k_c,k_d,k_l);

                     ward += psign * 0.5/( PPHM::norm(K_pph,k_a_,k_b_,k_l) * PPHM::norm(K_pph,k_c,k_d,k_l) )* hard;

                  }

               }

               for(int k_l = 0;k_l < Tools::gL();++k_l){

                  hard = 0.0;

                  K_pph = (k_a + k_b + k_l)%Tools::gL();

                  for(int pi = 0;pi < 2;++pi)
                     hard += pphm(S,K_pph,pi,Z,k_a,k_b,k_l,Z,k_c,k_d,k_l);

                  ward += 0.5/( PPHM::norm(K_pph,k_a,k_b,k_l) * PPHM::norm(K_pph,k_c,k_d,k_l) )* hard;

               }

               (*this)(B,i,j) += (2*(S + 0.5) + 1.0) * ward;

            }

            (*this)(B,i,j) *= TPM::norm(K,k_a,k_b) * TPM::norm(K,k_c,k_d)/(2.0*Z + 1.0);

         }
      }

   }

   this->symmetrize();

}

/**
 * The spincoupled T2-down map that maps a PPHM on a TPM object.
 * @param pphm Input PPHM object
 */
void TPM::T(const PPHM &pphm){

   //first make the bar tpm
   TPM tpm;
   tpm.bar(pphm);

   //then make the bar phm
   PHM phm;
   phm.bar(pphm);

   //also make the bar spm with the correct scale factor
   SPM spm;
   spm.bar(0.5/(Tools::gN() - 1.0),pphm);

   int k_a,k_b,k_c,k_d;
   int k_a_,k_b_,k_c_,k_d_;

   int K_ph;

   double ward;

   int sign;
   int psign;

   int S,K,p;

   for(int B = 0;B < gnr();++B){//loop over the blocks

      S = block_char[B][0];
      K = block_char[B][1];
      p = block_char[B][2];

      sign = 1 - 2*S;
      psign = 1 - 2*p;

      for(int i = 0;i < gdim(B);++i){

         k_a = t2s[B][i][0];
         k_b = t2s[B][i][1];

         //and for access to the phm elements:
         k_a_ = (Tools::gL() - k_a)%Tools::gL();
         k_b_ = (Tools::gL() - k_b)%Tools::gL();

         for(int j = i;j < gdim(B);++j){

            k_c = t2s[B][j][0];
            k_d = t2s[B][j][1];

            //and for access to the phm elements:
            k_c_ = (Tools::gL() - k_c)%Tools::gL();
            k_d_ = (Tools::gL() - k_d)%Tools::gL();

            (*this)(B,i,j) = 0.0;

            //four ph exchange terms:
            if(K == 0 || K == Tools::gL()/2){

               //1)
               K_ph = (k_a_ + k_d_)%Tools::gL();

               ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  for(int pi = 0;pi < 2;++pi)
                     ward -= (2.0*Z + 1.0) * Tools::g6j(0,0,S,Z) * phm(Z,K_ph,pi,k_a_,k_d_,k_c,k_b);

               ward /= 2.0 * PHM::norm(K_ph,k_a_,k_d_) * PHM::norm(K_ph,k_c,k_b);

               (*this)(B,i,j) += psign*ward;

               //2)
               K_ph = (k_b_ + k_c_)%Tools::gL();

               ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  for(int pi = 0;pi < 2;++pi)
                     ward -= (2.0*Z + 1.0) * Tools::g6j(0,0,S,Z) * phm(Z,K_ph,pi,k_b_,k_c_,k_d,k_a);

               ward /= 2.0 * PHM::norm(K_ph,k_b_,k_c_) * PHM::norm(K_ph,k_d,k_a);

               (*this)(B,i,j) += psign*ward;

               //3)
               K_ph = (k_b_ + k_d_)%Tools::gL();

               ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  for(int pi = 0;pi < 2;++pi)
                     ward -= (2.0*Z + 1.0) * Tools::g6j(0,0,S,Z) * phm(Z,K_ph,pi,k_b_,k_d_,k_c,k_a);

               ward /= 2.0 * PHM::norm(K_ph,k_b_,k_d_) * PHM::norm(K_ph,k_c,k_a);

               (*this)(B,i,j) += psign*sign*ward;

               //4)
               K_ph = (k_a_ + k_c_)%Tools::gL();

               ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  for(int pi = 0;pi < 2;++pi)
                     ward -= (2.0*Z + 1.0) * Tools::g6j(0,0,S,Z) * phm(Z,K_ph,pi,k_a_,k_c_,k_d,k_b);

               ward /= 2.0 * PHM::norm(K_ph,k_a_,k_c_) * PHM::norm(K_ph,k_d,k_b);

               (*this)(B,i,j) += psign*sign*ward;

            }

            //four regular ph terms:
            //1)
            K_ph = (k_a + k_d_)%Tools::gL();

            ward = 0.0;

            for(int Z = 0;Z < 2;++Z)
               for(int pi = 0;pi < 2;++pi)
                  ward -= (2.0*Z + 1.0) * Tools::g6j(0,0,S,Z) * phm(Z,K_ph,pi,k_a,k_d_,k_c,k_b_);

            ward /= 2.0 * PHM::norm(K_ph,k_a,k_d_) * PHM::norm(K_ph,k_c,k_b_);

            (*this)(B,i,j) += ward;

            //2)
            K_ph = (k_b + k_c_)%Tools::gL();

            ward = 0.0;

            for(int Z = 0;Z < 2;++Z)
               for(int pi = 0;pi < 2;++pi)
                  ward -= (2.0*Z + 1.0) * Tools::g6j(0,0,S,Z) * phm(Z,K_ph,pi,k_b,k_c_,k_d,k_a_);

            ward /= 2.0 * PHM::norm(K_ph,k_b,k_c_) * PHM::norm(K_ph,k_d,k_a_);

            (*this)(B,i,j) += ward;

            //3)
            K_ph = (k_b + k_d_)%Tools::gL();

            ward = 0.0;

            for(int Z = 0;Z < 2;++Z)
               for(int pi = 0;pi < 2;++pi)
                  ward -= (2.0*Z + 1.0) * Tools::g6j(0,0,S,Z) * phm(Z,K_ph,pi,k_b,k_d_,k_c,k_a_);

            ward /= 2.0 * PHM::norm(K_ph,k_b,k_d_) * PHM::norm(K_ph,k_c,k_a_);

            (*this)(B,i,j) += sign*ward;

            //4)
            K_ph = (k_a + k_c_)%Tools::gL();

            ward = 0.0;

            for(int Z = 0;Z < 2;++Z)
               for(int pi = 0;pi < 2;++pi)
                  ward -= (2.0*Z + 1.0) * Tools::g6j(0,0,S,Z) * phm(Z,K_ph,pi,k_a,k_c_,k_d,k_b_);

            ward /= 2.0 * PHM::norm(K_ph,k_a,k_c_) * PHM::norm(K_ph,k_d,k_b_);

            (*this)(B,i,j) += sign*ward;

            //tp-norm:
            (*this)(B,i,j) *= TPM::norm(K,k_a,k_b) * TPM::norm(K,k_c,k_d);

            if(k_a == k_b)
               (*this)(B,i,j) /= std::sqrt(2.0);

            if(k_c == k_d)
               (*this)(B,i,j) /= std::sqrt(2.0);

            //finally the tp part
            (*this)(B,i,j) += tpm(B,i,j);
           
         }

         //sp part is diagonal for translationaly invariance
         (*this)(B,i,i) += spm[k_a] + spm[k_b];

      }

   }

   this->symmetrize();

}

/**
 * Construct the right hand side of the Tools::gN()ewton equation for the determination of the search direction, 
 * the gradient of the potential:
 * @param t scaling factor of the potential
 * @param ham Hamiltonian of the current problem
 * @param P SUP matrix containing the inverse of the constraint matrices (carrier space matrices).
 */
void TPM::constr_grad(double t,const TPM &ham,const SUP &P){

   //eerst P conditie 
   *this = P.tpm(0);

   //de Q conditie toevoegen
   TPM hulp;

   hulp.Q(1,P.tpm(1));

   *this += hulp;

#ifdef __G_CON

   hulp.G(P.phm());

   *this += hulp;

#endif

#ifdef __T1_CON

   hulp.T(P.dpm());

   *this += hulp;

#endif

#ifdef __T2_CON

   hulp.T(P.pphm());

   *this +=hulp;

#endif

   this->dscal(t);

   *this -= ham;

   this->proj_Tr();

}

/**
 * solve the Tools::gN()ewton equations for the determination of the search direction,
 * @param t scaling factor of the potential
 * @param P SUP matrix containing the inverse of the constraint matrices (carrier space matrices).
 * @param b right hand side (the gradient constructed int TPM::constr_grad)
 * @return nr of iterations needed to converge to the desired accuracy
 */
int TPM::solve(double t,const SUP &P,TPM &b){

   int iter = 0;

   //delta = 0
   *this = 0;

   //residu:
   TPM r(b);

   //norm van het residu
   double rr = r.ddot(r);

   //enkele variabelen
   double rr_old,ward;

   TPM Hb;

   while(rr > 1.0e-10){ 

      Hb.H(t,b,P);

      ward = rr/b.ddot(Hb);

      //delta += ward*b
      this->daxpy(ward,b);

      //r -= ward*Hb
      r.daxpy(-ward,Hb);

      //nieuwe variabelen berekenen en oude overdragen
      rr_old = rr;
      rr = r.ddot(r);

      //nieuwe b nog:
      b.dscal(rr/rr_old);

      b += r;

      ++iter;

   }

   return iter;

}

/**
 * perform a line search what step size in along the Tools::gN()ewton direction is ideal.
 * @param t potential scaling factor
 * @param P SUP matrix containing the inverse of the constraints (carrier space matrices)
 * @param ham Hamiltonian of the problem
 * @return the steplength
 */
double TPM::line_search(double t,SUP &P,const TPM &ham){

   double tolerance = 1.0e-5*t;

   if(tolerance < 1.0e-12)
      tolerance = 1.0e-12;

   //neem de wortel uit P
   P.sqrt(1);

   //maak eerst een SUP van delta
   SUP S_delta;

   S_delta.fill(*this);

   //hulpje om dingskes in te steken:
   SUP hulp;

   hulp.L_map(P,S_delta);

   EIG eigen(hulp);

   double a = 0;

   double b = -1.0/eigen.min();

   double c(0);

   double ham_delta = ham.ddot(*this);

   while(b - a > tolerance){

      c = (b + a)/2.0;

      if( (ham_delta - t*eigen.lsfunc(c)) < 0.0)
         a = c;
      else
         b = c;

   }

   return c;

}

/**
 * The hessian-map of the Tools::gN()ewton system:
 * @param t potential scaling factor
 * @param b the TPM on which the hamiltonian will work, the image will be put in (*this)
 * @param P the SUP matrix containing the constraints, (can be seen as the metric).
 */
void TPM::H(double t,const TPM &b,const SUP &P){

   //eerst de P conditie:

   this->L_map(P.tpm(0),b);

   TPM hulp;

   //maak Q(b)
   TPM Q_b;
   Q_b.Q(1,b);

   //stop Q(rdm)^{-1}Q(b)Q(rdm)^{-1} in hulp
   hulp.L_map(P.tpm(1),Q_b);

   //maak Q(hulp) en stop in Q_b
   Q_b.Q(1,hulp);

   //en tel op bij this
   *this += Q_b;

#ifdef __G_CON

   //hulpje voor het PHM stuk
   PHM hulp_ph;
   PHM G_b;

   //stop G(b) in G_b
   G_b.G(b);

   //bereken G(rdm)^{-1}G(b)G(rdm)^{-1} en stop in hulp_ph
   hulp_ph.L_map(P.phm(),G_b);

   //tenslotte nog de antisymmetrische G hierop:
   hulp.G(hulp_ph);

   //en optellen bij this
   *this += hulp;

#endif
   
#ifdef __T1_CON

   //hulpjes voor het DPM stuk
   DPM hulp_dp;
   DPM T1_b;

   //stop T1(b) in T1_b
   T1_b.T(b);

   hulp_dp.L_map(P.dpm(),T1_b);

   hulp.T(hulp_dp);

   *this += hulp;

#endif

#ifdef __T2_CON

   PPHM hulp_pph;
   PPHM T2_b;

   T2_b.T(b);

   hulp_pph.L_map(P.pphm(),T2_b);

   hulp.T(hulp_pph);

   *this+=hulp;

#endif

   //nog schalen met t:
   this->dscal(t);

   //en projecteren op spoorloze ruimte
   this->proj_Tr();

}

/**
 * perform a line search what step size in along the Newton direction is ideal, this one is used for extrapolation.
 * @param t potential scaling factor
 * @param rdm TPM containing the current approximation of the rdm.
 * @param ham Hamiltonian of the problem
 * @return the steplength
 */
double TPM::line_search(double t,const TPM &rdm,const TPM &ham){

   SUP P;

   P.fill(rdm);

   P.invert();

   return this->line_search(t,P,ham);

}
