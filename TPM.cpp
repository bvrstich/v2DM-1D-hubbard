#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::endl;
using std::ios;

#include "include.h"

//def of some statics
vector< vector<int> > *TPM::t2s;
int ***TPM::s2t;

int **TPM::block_char;
int ***TPM::char_block;

double **TPM::_6j;

int TPM::M;
int TPM::N;
int TPM::L;

double TPM::Sa = 1;
double TPM::Sc = 0;

/**
 * static function that initializes the static variables and allocates and fill the static lists
 * @param L_in nr of sites
 * @param N_in nr of particles
 */
void TPM::init(int L_in,int N_in){

   L = L_in;
   N = N_in;

   M = L*2;

   int nr_B = L + 4;//nr of blocks

   //allocate
   t2s = new vector< vector<int> > [nr_B];

   s2t = new int ** [nr_B];

   for(int B = 0;B < nr_B;++B){

      s2t[B] = new int * [L];

      for(int k_a = 0;k_a < L;++k_a)
         s2t[B][k_a] = new int [L];

   }

   block_char = new int * [nr_B];

   for(int B = 0;B < nr_B;++B)
      block_char[B] = new int [3];

   char_block = new int ** [2];

   for(int S = 0;S < 2;++S){

      char_block[S] = new int * [L/2 + 1];

      for(int K = 0;K <= L/2;++K)
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

   for(int k_a = 0;k_a < L;++k_a)
      for(int k_b = k_a;k_b < L;++k_b){

         if( (k_a + k_b)%L == 0 ){

            v[0] = k_a;
            v[1] = k_b;

            t2s[block].push_back(v);

            s2t[block][k_a][k_b]  = tp;
            s2t[block][k_b][k_a]  = tp;

            ++tp;

         }

      }

   //then S = 1: block shifted with L/2 + 2 (only parity = -1)
   block_char[block + L/2 + 2][0] = 1;
   block_char[block + L/2 + 2][1] = 0;
   block_char[block + L/2 + 2][2] = 1;//parity baby! (0 means +1, 1 means -1)

   char_block[1][0][1] = block + L/2 + 2;

   tp = 0;

   for(int k_a = 0;k_a < L;++k_a)
      for(int k_b = k_a + 1;k_b < L;++k_b){

         if( (k_a + k_b)%L == 0 ){

            v[0] = k_a;
            v[1] = k_b;

            t2s[block + L/2 + 2].push_back(v);

            s2t[block + L/2 + 2][k_a][k_b]  = tp;
            s2t[block + L/2 + 2][k_b][k_a]  = tp;

            ++tp;

         }

      }

   ++block;


   //then 0 < K < L/2
   for(int K = 1;K < L/2;++K){

      //first S = 0
      block_char[block][0] = 0;
      block_char[block][1] = K;
      block_char[block][2] = 0;//parity baby! (0 means +1, 1 means -1)

      //both parities refer to the same block
      char_block[0][K][0] = block;
      char_block[0][K][1] = block;

      tp = 0;

      for(int k_a = 0;k_a < L;++k_a)
         for(int k_b = k_a;k_b < L;++k_b){

            if( (k_a + k_b)%L == K ){

               v[0] = k_a;
               v[1] = k_b;

               t2s[block].push_back(v);

               s2t[block][k_a][k_b]  = tp;
               s2t[block][k_b][k_a]  = tp;

               ++tp;

            }

         }

      //then S = 1: block shifted with L/2 + 2
      block_char[block + L/2 + 2][0] = 1;
      block_char[block + L/2 + 2][1] = K;
      block_char[block + L/2 + 2][2] = 0;//parity baby! (0 means +1, 1 means -1)

      //both parities refer to the same block
      char_block[1][K][0] = block + L/2 + 2;
      char_block[1][K][1] = block + L/2 + 2;

      tp = 0;

      for(int k_a = 0;k_a < L;++k_a)
         for(int k_b = k_a + 1;k_b < L;++k_b){

            if( (k_a + k_b)%L == K ){

               v[0] = k_a;
               v[1] = k_b;

               t2s[block + L/2 + 2].push_back(v);

               s2t[block + L/2 + 2][k_a][k_b]  = tp;
               s2t[block + L/2 + 2][k_b][k_a]  = tp;

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
   block_char[block][1] = L/2;
   block_char[block][2] = 0;//parity baby! (0 means +1, 1 means -1)

   char_block[0][L/2][0] = block;

   tp = 0;

   //the momenta only go to L/2!
   for(int k_a = 0;k_a <= L/2;++k_a)
      for(int k_b = k_a;k_b <= L/2;++k_b){

         if( (k_a + k_b)%L == L/2 ){

            v[0] = k_a;
            v[1] = k_b;

            t2s[block].push_back(v);

            s2t[block][k_a][k_b]  = tp;
            s2t[block][k_b][k_a]  = tp;

            ++tp;

         }

      }

   //then S = 1: block shifted with L/2 + 2
   block_char[block + L/2 + 2][0] = 1;
   block_char[block + L/2 + 2][1] = L/2;
   block_char[block + L/2 + 2][2] = 0;//parity baby! (0 means +1, 1 means -1)

   char_block[1][L/2][0] = block + L/2 + 2;

   tp = 0;

   for(int k_a = 0;k_a <= L/2;++k_a)
      for(int k_b = k_a + 1;k_b <= L/2;++k_b){

         if( (k_a + k_b)%L == L/2 ){

            v[0] = k_a;
            v[1] = k_b;

            t2s[block + L/2 + 2].push_back(v);

            s2t[block + L/2 + 2][k_a][k_b]  = tp;
            s2t[block + L/2 + 2][k_b][k_a]  = tp;

            ++tp;

         }

      }

   ++block;

   //then negative parity: this block hasn't got the |0 L/2> basisvector

   //first S = 0;
   block_char[block][0] = 0;
   block_char[block][1] = L/2;
   block_char[block][2] = 1;//parity baby! (0 means +1, 1 means -1)

   char_block[0][L/2][1] = block;

   tp = 0;

   //the only difference is that k_a starts from 1, and the momenta only go to L/2
   for(int k_a = 1;k_a < L/2;++k_a)
      for(int k_b = k_a;k_b < L/2;++k_b){

         if( (k_a + k_b)%L == L/2 ){

            v[0] = k_a;
            v[1] = k_b;

            t2s[block].push_back(v);

            s2t[block][k_a][k_b]  = tp;
            s2t[block][k_b][k_a]  = tp;

            ++tp;

         }

      }

   //then S = 1: block shifted with L/2 + 2
   block_char[block + L/2 + 2][0] = 1;
   block_char[block + L/2 + 2][1] = L/2;
   block_char[block + L/2 + 2][2] = 1;//parity baby! (0 means +1, 1 means -1)

   char_block[1][L/2][1] = block + L/2 + 2;

   tp = 0;

   for(int k_a = 1;k_a < L/2;++k_a)
      for(int k_b = k_a + 1;k_b < L/2;++k_b){

         if( (k_a + k_b)%L == L/2 ){

            v[0] = k_a;
            v[1] = k_b;

            t2s[block + L/2 + 2].push_back(v);

            s2t[block + L/2 + 2][k_a][k_b]  = tp;
            s2t[block + L/2 + 2][k_b][k_a]  = tp;

            ++tp;

         }

      }

   //allocate 6j
   _6j = new double * [2];

   for(int S = 0;S < 2;++S)
      _6j[S] = new double [2]; 

   //initialize
   _6j[0][0] = -0.5;
   _6j[0][1] = 0.5;
   _6j[1][0] = 0.5;
   _6j[1][1] = 1.0/6.0;

   init_overlap();

}

/**
 * initialize the overlapmatrix parameters
 */
void TPM::init_overlap(){

   Sa += 1.0;
   Sc += (2.0*N - M)/((N - 1.0)*(N - 1.0));

}

/**
 * deallocates the statics lists
 */
void TPM::clear(){

   delete [] t2s;

   for(int B = 0;B < L + 4;++B){

      for(int k_a = 0;k_a < L;++k_a)
         delete [] s2t[B][k_a];

      delete [] s2t[B];

   }

   delete [] s2t;

   for(int B = 0;B < L + 4;++B)
      delete [] block_char[B];

   delete [] block_char;

   for(int S = 0;S < 2;++S){

      for(int K = 0;K <= L/2;++K)
         delete [] char_block[S][K];

      delete [] char_block[S];

   }

   delete [] char_block;

   for(int S = 0;S < 2;++S)
      delete [] _6j[S];

   delete [] _6j;

}

/**
 * standard constructor for a spinsymmetrical, translationally invariant tp matrix, with parity (k <--> -k) taken into account: 
 * constructs BlockMatrix object with L + 4 blocks, L/2 + 2 for S = 0 and S = 1,
 */
TPM::TPM() : BlockMatrix(L + 4) {

   //K = 0
   this->setMatrixDim(0,t2s[0].size(),1);
   this->setMatrixDim(L/2 + 2,t2s[L/2 + 2].size(),3);

   //0 < K < L/2
   for(int B = 1;B < L/2;++B){

      this->setMatrixDim(B,t2s[B].size(),2);
      this->setMatrixDim(B + L/2 + 2,t2s[B + L/2 + 2].size(),6);

   }

   //K = L/2: 4 more blocks (positive and negative parity in S = 0/1)
   this->setMatrixDim(L/2,t2s[L/2].size(),1);
   this->setMatrixDim(L/2 + 1,t2s[L/2 + 1].size(),1);

   this->setMatrixDim(L + 2,t2s[L + 2].size(),3);
   this->setMatrixDim(L + 3,t2s[L + 3].size(),3);

}

/**
 * copy constructor for a spinsymmetrical, translationally invariant tp matrix, with parity (k <--> -k) taken into account: 
 * constructs BlockMatrix object with L + 4 blocks, L/2 + 2 for S = 0 and S = 1,
 * @param tpm_c The TPM object to be copied into (*this)
 */
TPM::TPM(const TPM &tpm_c) : BlockMatrix(tpm_c){

}

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
   if( (k_a + k_b)%(M/2) != K)
      return 0;

   if( (k_c + k_d)%(M/2) != K)
      return 0;

   //for the K = 0 , S = 0 blocks only one parity type is present
   if(K == 0 && S == 0 && p == 1)
      return 0;

   if(K == 0 && S == 1 && p == 0)
      return 0;

   int B = char_block[S][K][p];

   if(S == 0){

      int i = s2t[B][k_a][k_b];
      int j = s2t[B][k_c][k_d];

      return (*this)(B,i,j);

   }
   else{

      if( (k_a == k_b) || (k_c == k_d) )
         return 0;
      else{

         int i = s2t[B][k_a][k_b];
         int j = s2t[B][k_c][k_d];

         int phase = 1;

         if(k_a > k_b)
            phase *= -1;
         if(k_c > k_d)
            phase *= -1;

         return phase*(*this)(B,i,j);

      }

   }

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
 * @return number of particles
 */
int TPM::gN() const{

   return N;

}

/**
 * @return number of sp orbitals
 */
int TPM::gM() const{

   return M;

}

/**
 * @return number of sites
 */
int TPM::gL() const{

   return L;

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

   double ward = 1.0/(N - 1.0);

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
               (*this)(B,i,i) = -2.0 * ward * ( cos( 4.0 * k_a * 3.141592653589793238462 / (double) M)  + cos( 4.0 * k_b * 3.141592653589793238462 / (double) M) );

            //on-site repulsion
            if(S == 0){

               double ward = 4.0*U / (double) M;

               if(k_a == k_b)
                  ward /= std::sqrt(2.0);

               if(k_c == k_d)
                  ward /= std::sqrt(2.0);

               if(K != L/2)
                  (*this)(B,i,j) += ward;
               else//K = L/2
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
   double b = 1.0/(N*(N - 1.0));
   double c = 1.0/(N - 1.0);

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

      B = (B*A + B*C*M - 2.0*C*C)/( A * (C*(M - 2.0) -  A) * ( A + B*M*(M - 1.0) - 2.0*C*(M - 1.0) ) );
      C = C/(A*(C*(M - 2.0) - A));
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
 * initialize this onto the unitmatrix with trace N*(N - 1)/2
 */
void TPM::unit(){

   double ward = N*(N - 1.0)/(M*(M - 1.0));

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

   double ward = (2.0 * this->trace())/(M*(M - 1));

   this->min_unit(ward);

}

/**
 * Primal hessian map:\n\n
 * Hb = D_1 b D_1 + D_2 Q(b) D_2 + D_3 G(b) D_3 + D_4 T1(b) D_4 + D_5 T2(b) D5 \n\n
 * with D_1, D_2, D_3, D_4 and D_5 the P, Q, G, T1 and T2 blocks of the SUP D. 
 * @param b TPM domain matrix, hessian will act on it and the image will be put in this
 * @param D SUP matrix that defines the structure of the hessian map. (see primal-dual.pdf for more info)
 */

void TPM::H(const TPM &b,const SUP &D){

   this->L_map(D.tpm(0),b);

   //maak Q(b)
   TPM Qb;
   Qb.Q(1,b);

   TPM hulp;

   hulp.L_map(D.tpm(1),Qb);

   Qb.Q(1,hulp);

   *this += Qb;

   this->proj_Tr();

}

/**
 * Implementation of a linear conjugate gradient algoritm for the solution of the primal Newton equations\n\n
 * H(*this) =  b\n\n 
 * in which H represents the hessian map.
 * @param b righthandside of the equation
 * @param D SUP matrix that defines the structure of the hessian
 * @return return number of iterations needed to converge to the desired accuracy
 */

int TPM::solve(TPM &b,const SUP &D){

   *this = 0;

   //de r initialiseren op b
   TPM r(b);

   double rr = r.ddot(r);
   double rr_old,ward;

   //nog de Hb aanmaken ook, maar niet initialiseren:
   TPM Hb;

   int cg_iter = 0;

   while(rr > 1.0e-15){

      ++cg_iter;

      Hb.H(b,D);

      ward = rr/b.ddot(Hb);

      //delta += ward*b
      this->daxpy(ward,b);

      //r -= ward*Hb
      r.daxpy(-ward,Hb);

      //nieuwe r_norm maken
      rr_old = rr;
      rr = r.ddot(r);

      //eerst herschalen van b
      b.dscal(rr/rr_old);

      //dan r er bijtellen
      b += r;

   }

   return cg_iter;

}

/**
 * ( Overlapmatrix of the U-basis ) - map, maps a TPM onto a different TPM, this map is actually a Q-like map
 * for which the paramaters a,b and c are calculated in primal_dual.pdf. Since it is a Q-like map the inverse
 * can be taken as well.
 * @param option = 1 direct overlapmatrix-map is used , = -1 inverse overlapmatrix map is used
 * @param tpm_d the input TPM
 */

void TPM::S(int option,const TPM &tpm_d){

   this->Q(option,Sa,0.0,Sc,tpm_d);

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
            ward += -1.5 * (N - 2.0)/(N - 1.0) * (*this)(B,i,i);

      }
      else{

         for(int i = 0;i < this->gdim(B);++i)
            ward += 3.0 * ( -1.5 * (N - 2.0)/(N - 1.0) + 2.0 ) * (*this)(B,i,i);

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
