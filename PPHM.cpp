#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::ofstream;
using std::cout;
using std::endl;

#include "include.h"

vector< vector<int> > *PPHM::pph2s;
int *****PPHM::s2pph;

double **PPHM::_6j;

int **PPHM::block_char;
int ***PPHM::char_block;

int PPHM::M;
int PPHM::L;
int PPHM::N;

/**
 * initialize the static lists and static variables
 * @param L_in nr of sites
 * @param N_in nr of particles
 */
void PPHM::init(int L_in,int N_in){

   L = L_in;
   N = N_in;

   M = 2*L;

   //allocate first
   pph2s = new vector< vector<int> > [L + 6];

   s2pph = new int **** [L + 6];

   for(int B = 0;B < L + 6;++B){

      s2pph[B] = new int *** [2];

      for(int S_ab = 0;S_ab < 2;++S_ab){

         s2pph[B][S_ab] = new int ** [L];

         for(int k_a = 0;k_a < L;++k_a){

            s2pph[B][S_ab][k_a] = new int * [L];

            for(int k_b = 0;k_b < L;++k_b)
               s2pph[B][S_ab][k_a][k_b] = new int [L];

         }
      }
   }

   block_char = new int * [L + 6];

   for(int B = 0;B < L + 6;++B)
      block_char[B] = new int [3];

   char_block = new int ** [2];

   for(int S = 0;S < 2;++S){

      char_block[S] = new int * [L/2 + 1];

      for(int K = 0;K <= L/2;++K)
         char_block[S][K] = new int [2];

   }

   //initialize
   int block = 0;

   int pph = 0;

   vector<int> v(4);

   //first K = 0, positive parity

   //S = 1/2
   block_char[block][0] = 0;//S
   block_char[block][1] = 0;//K
   block_char[block][2] = 0;//p

   char_block[0][0][0] = block;

   //S_ab = 0

   //first three terms that only have positive parity

   //1)
   v[0] = 0;//S_ab
   v[1] = 0;
   v[2] = 0;
   v[3] = 0;

   pph2s[block].push_back(v);

   s2pph[block][0][0][0][0] = pph;

   ++pph;

   //2)
   v[0] = 0;//S_ab
   v[1] = 0;
   v[2] = L/2;
   v[3] = L/2;

   pph2s[block].push_back(v);

   s2pph[block][0][0][L/2][L/2] = pph;

   ++pph;

   //3)
   v[0] = 0;//S_ab
   v[1] = L/2;
   v[2] = L/2;
   v[3] = 0;

   pph2s[block].push_back(v);

   s2pph[block][0][L/2][L/2][0] = pph;

   ++pph;

   //then the states which have a k = 0
   for(int k_a = 1;k_a < L/2;++k_a){

      int k_b = L - k_a;

      v[0] = 0;//S_ab
      v[1] = k_a;
      v[2] = k_b;
      v[3] = 0;

      pph2s[block].push_back(v);

      s2pph[block][0][k_a][k_b][0] = pph;

      ++pph;

      v[0] = 0;//S_ab
      v[1] = 0;
      v[2] = k_a;
      v[3] = k_b;

      pph2s[block].push_back(v);

      s2pph[block][0][0][k_a][k_b] = pph;

      ++pph;

   }

   //then the rest
   for(int k_a = 1;k_a < L;++k_a)
      for(int k_b = k_a;k_b < L;++k_b)
         for(int k_c = 1;k_c < L;++k_c){

            if(k_a + k_b + k_c == L){// no % !

               v[0] = 0;
               v[1] = k_a;
               v[2] = k_b;
               v[3] = k_c;

               pph2s[block].push_back(v);

               s2pph[block][0][k_a][k_b][k_c] = pph;

               ++pph;

            }

         }

   //S_ab = 1

   //first one that has only positive parity:
   v[0] = 1;
   v[1] = 0;
   v[2] = L/2;
   v[3] = L/2;

   pph2s[block].push_back(v);

   s2pph[block][1][0][L/2][L/2] = pph;

   ++pph;

   //then those terms with a k = 0
   for(int k_a = 1;k_a < L/2;++k_a){

      int k_b = L - k_a;

      v[0] = 1;
      v[1] = 0;
      v[2] = k_a;
      v[3] = k_b;

      pph2s[block].push_back(v);

      s2pph[block][1][0][k_a][k_b] = pph;

      ++pph;

   }

   //then the rest:
   for(int k_a = 1;k_a < L;++k_a)
      for(int k_b = k_a + 1;k_b < L;++k_b)
         for(int k_c = 1;k_c < L;++k_c){

            if(k_a + k_b + k_c == L){// no % !

               v[0] = 1;
               v[1] = k_a;
               v[2] = k_b;
               v[3] = k_c;

               pph2s[block].push_back(v);

               s2pph[block][1][k_a][k_b][k_c] = pph;

               ++pph;

            }

         }

   pph = 0;

   //S = 3/2
   block_char[block + L/2 + 3][0] = 1;//S
   block_char[block + L/2 + 3][1] = 0;//K
   block_char[block + L/2 + 3][2] = 0;//p

   char_block[1][0][0] = block + L/2 + 3;

   //only S_ab = 1 for S = 3/2 of course

   //first one that has only positive parity:
   v[0] = 1;
   v[1] = 0;
   v[2] = L/2;
   v[3] = L/2;

   pph2s[block + L/2 + 3].push_back(v);

   s2pph[block + L/2 + 3][1][0][L/2][L/2] = pph;

   ++pph;

   //then those terms with a k = 0
   for(int k_a = 1;k_a < L/2;++k_a){

      int k_b = L - k_a;

      v[0] = 1;
      v[1] = 0;
      v[2] = k_a;
      v[3] = k_b;

      pph2s[block + L/2 + 3].push_back(v);

      s2pph[block + L/2 + 3][1][0][k_a][k_b] = pph;

      ++pph;

   }

   //then the rest:
   for(int k_a = 1;k_a < L;++k_a)
      for(int k_b = k_a + 1;k_b < L;++k_b)
         for(int k_c = 1;k_c < L;++k_c){

            if(k_a + k_b + k_c == L){// no % !

               v[0] = 1;
               v[1] = k_a;
               v[2] = k_b;
               v[3] = k_c;

               pph2s[block + L/2 + 3].push_back(v);

               s2pph[block + L/2 + 3][1][k_a][k_b][k_c] = pph;

               ++pph;

            }

         }

   ++block;

   pph = 0;

   //K = 0 : negative parity

   //S = 1/2
   block_char[block][0] = 0;//S
   block_char[block][1] = 0;//K
   block_char[block][2] = 1;//p

   char_block[0][0][1] = block;

   //S_ab = 0: terms that have a zero
   for(int k_a = 1;k_a < L/2;++k_a){

      int k_b = L - k_a;

      v[0] = 0;
      v[1] = 0;
      v[2] = k_a;
      v[3] = k_b;

      pph2s[block].push_back(v);

      s2pph[block][0][0][k_a][k_b] = pph;

      ++pph;

   }

   //then the rest:
   for(int k_a = 1;k_a < L;++k_a)
      for(int k_b = k_a;k_b < L;++k_b)
         for(int k_c = 1;k_c < L;++k_c){

            if(k_a + k_b + k_c == L){// no % !

               v[0] = 0;
               v[1] = k_a;
               v[2] = k_b;
               v[3] = k_c;

               pph2s[block].push_back(v);

               s2pph[block][0][k_a][k_b][k_c] = pph;

               ++pph;

            }

         }

   //S_ab = 1: terms that have a zero
   for(int k_a = 1;k_a < L/2;++k_a){

      int k_b = L - k_a;

      v[0] = 1;
      v[1] = k_a;
      v[2] = k_b;
      v[3] = 0;

      pph2s[block].push_back(v);

      s2pph[block][1][k_a][k_b][0] = pph;

      ++pph;

      v[0] = 1;
      v[1] = 0;
      v[2] = k_a;
      v[3] = k_b;

      pph2s[block].push_back(v);

      s2pph[block][1][0][k_a][k_b] = pph;

      ++pph;

   }

   //then the rest:
   for(int k_a = 1;k_a < L;++k_a)
      for(int k_b = k_a + 1;k_b < L;++k_b)
         for(int k_c = 1;k_c < L;++k_c){

            if(k_a + k_b + k_c == L){// no % !

               v[0] = 1;
               v[1] = k_a;
               v[2] = k_b;
               v[3] = k_c;

               pph2s[block].push_back(v);

               s2pph[block][1][k_a][k_b][k_c] = pph;

               ++pph;

            }

         }

   pph = 0;

   //S = 3/2 K = 0 negative parity
   block_char[block + L/2 + 3][0] = 1;//S
   block_char[block + L/2 + 3][1] = 0;//K
   block_char[block + L/2 + 3][2] = 1;//p

   char_block[1][0][1] = block + L/2 + 3;

   //terms with a zero
   for(int k_a = 1;k_a < L/2;++k_a){

      int k_b = L - k_a;

      v[0] = 1;
      v[1] = k_a;
      v[2] = k_b;
      v[3] = 0;

      pph2s[block + L/2 + 3].push_back(v);

      s2pph[block + L/2 + 3][1][k_a][k_b][0] = pph;

      ++pph;

      v[0] = 1;
      v[1] = 0;
      v[2] = k_a;
      v[3] = k_b;

      pph2s[block + L/2 + 3].push_back(v);

      s2pph[block + L/2 + 3][1][0][k_a][k_b] = pph;

      ++pph;

   }

   //then the rest:
   for(int k_a = 1;k_a < L;++k_a)
      for(int k_b = k_a + 1;k_b < L;++k_b)
         for(int k_c = 1;k_c < L;++k_c){

            if(k_a + k_b + k_c == L){// no % !

               v[0] = 1;
               v[1] = k_a;
               v[2] = k_b;
               v[3] = k_c;

               pph2s[block + L/2 + 3].push_back(v);

               s2pph[block + L/2 + 3][1][k_a][k_b][k_c] = pph;

               ++pph;

            }

         }

   ++block;

   //now for 0 < K < L/2, it should be easier
   for(int K = 1;K < L/2;++K){

      pph = 0;

      //first S = 1/2
      block_char[block][0] = 0;//S
      block_char[block][1] = K;//K
      block_char[block][2] = 0;//p

      //both positive and negative parity refer to the same block
      char_block[0][K][0] = block;
      char_block[0][K][1] = block;

      for(int S_ab = 0;S_ab < 2;++S_ab){

         for(int k_a = 0;k_a < L;++k_a)
            for(int k_b = k_a + S_ab;k_b < L;++k_b)
               for(int k_c = 0;k_c < L;++k_c){

                  if( (k_a + k_b + k_c)%L == K ){

                     v[0] = S_ab;
                     v[1] = k_a;
                     v[2] = k_b;
                     v[3] = k_c;

                     pph2s[block].push_back(v);

                     s2pph[block][S_ab][k_a][k_b][k_c] = pph;

                     ++pph;

                  }

               }

      }

      pph = 0;

      //then S = 3/2
      block_char[block + L/2 + 3][0] = 1;//S
      block_char[block + L/2 + 3][1] = K;//K
      block_char[block + L/2 + 3][2] = 0;//p

      //both positive and negative parity refer to the same block
      char_block[1][K][0] = block + L/2 + 3;
      char_block[1][K][1] = block + L/2 + 3;

      //only S_ab = 1 possible here
      for(int k_a = 0;k_a < L;++k_a)
         for(int k_b = k_a + 1;k_b < L;++k_b)
            for(int k_c = 0;k_c < L;++k_c){

               if( (k_a + k_b + k_c)%L == K ){

                  v[0] = 1;
                  v[1] = k_a;
                  v[2] = k_b;
                  v[3] = k_c;

                  pph2s[block + L/2 + 3].push_back(v);

                  s2pph[block + L/2 + 3][1][k_a][k_b][k_c] = pph;

                  ++pph;

               }

            }
      
      ++block;

   }

   pph = 0;

   //finally the K = L/2 blocks, first positive parity

   //S = 1/2
   block_char[block][0] = 0;//S
   block_char[block][1] = L/2;//K
   block_char[block][2] = 0;//p

   char_block[0][L/2][0] = block;

   //S_ab = 0: first three terms that only have positive parity

   //1)
   v[0] = 0;//S_ab
   v[1] = 0;
   v[2] = 0;
   v[3] = L/2;

   pph2s[block].push_back(v);

   s2pph[block][0][0][0][L/2] = pph;

   ++pph;

   //2)
   v[0] = 0;//S_ab
   v[1] = 0;
   v[2] = L/2;
   v[3] = 0;

   pph2s[block].push_back(v);

   s2pph[block][0][0][L/2][0] = pph;

   ++pph;

   //3)
   v[0] = 0;//S_ab
   v[1] = L/2;
   v[2] = L/2;
   v[3] = L/2;

   pph2s[block].push_back(v);

   s2pph[block][0][L/2][L/2][L/2] = pph;

   ++pph;

   //then the states containing k = L/2
   for(int k_a = 1;k_a < L/2;++k_a){

      int k_b = L - k_a;

      v[0] = 0;
      v[1] = k_a;
      v[2] = k_b;
      v[3] = L/2;

      pph2s[block].push_back(v);

      s2pph[block][0][k_a][k_b][L/2] = pph;

      ++pph;

      v[0] = 0;
      v[1] = k_a;
      v[2] = L/2;
      v[3] = k_b;

      pph2s[block].push_back(v);

      s2pph[block][0][k_a][L/2][k_b] = pph;

      ++pph;

   }

   //the states not containing k = L/2 can be split up in two groups:

   //first k_a,k_b,k_c < L/2
   for(int k_a = 0;k_a < L/2;++k_a)
      for(int k_b = k_a;k_b < L/2;++k_b)
         for(int k_c = 0;k_c < L/2;++k_c){

            if(k_a + k_b + k_c == L/2){

               v[0] = 0;
               v[1] = k_a;
               v[2] = k_b;
               v[3] = k_c;

               pph2s[block].push_back(v);

               s2pph[block][0][k_a][k_b][k_c] = pph;

               ++pph;

            }

         }

   //second k_a <= k_b and k_c > L/2 with sum = L + L/2
   for(int k_a = 1;k_a < L;++k_a){

      if(k_a == L/2)
         k_a++;

      for(int k_b = k_a;k_b < L;++k_b){

         if(k_b == L/2)
            k_b++;

         for(int k_c = L/2 + 1;k_c < L;++k_c){

            if(k_a + k_b + k_c == L + L/2){

               v[0] = 0;
               v[1] = k_a;
               v[2] = k_b;
               v[3] = k_c;

               pph2s[block].push_back(v);

               s2pph[block][0][k_a][k_b][k_c] = pph;

               ++pph;

            }

         }
      }
   }

   //S_ab = 1 here only one term that has only positive parity

   v[0] = 1;//S_ab
   v[1] = 0;
   v[2] = L/2;
   v[3] = 0;

   pph2s[block].push_back(v);

   s2pph[block][1][0][L/2][0] = pph;

   ++pph;

   //then the states containing k = L/2: only |k_a L/2 k_b> because for S_ab = 1 |k_a k_b L/2> has negative parity automatically
   for(int k_a = 1;k_a < L/2;++k_a){

      int k_b = L - k_a;

      v[0] = 1;
      v[1] = k_a;
      v[2] = L/2;
      v[3] = k_b;

      pph2s[block].push_back(v);

      s2pph[block][1][k_a][L/2][k_b] = pph;

      ++pph;

   }

   //the states not containing k = L/2 can be split up in two groups:

   //first k_a,k_b,k_c < L/2
   for(int k_a = 0;k_a < L/2;++k_a)
      for(int k_b = k_a + 1;k_b < L/2;++k_b)
         for(int k_c = 0;k_c < L/2;++k_c){

            if(k_a + k_b + k_c == L/2){

               v[0] = 1;
               v[1] = k_a;
               v[2] = k_b;
               v[3] = k_c;

               pph2s[block].push_back(v);

               s2pph[block][1][k_a][k_b][k_c] = pph;

               ++pph;

            }

         }

   //second k_a < k_b and k_c > L/2 with sum = L + L/2
   for(int k_a = 1;k_a < L;++k_a){

      if(k_a == L/2)
         k_a++;

      for(int k_b = k_a + 1;k_b < L;++k_b){

         if(k_b == L/2)
            k_b++;

         for(int k_c = L/2 + 1;k_c < L;++k_c){

            if(k_a + k_b + k_c == L + L/2){

               v[0] = 1;
               v[1] = k_a;
               v[2] = k_b;
               v[3] = k_c;

               pph2s[block].push_back(v);

               s2pph[block][1][k_a][k_b][k_c] = pph;

               ++pph;

            }

         }
      }
   }

   pph = 0;

   //now S = 3/2 with positive parity and K = L/2
   block_char[block + L/2 + 3][0] = 1;//S
   block_char[block + L/2 + 3][1] = L/2;//K
   block_char[block + L/2 + 3][2] = 0;//p

   char_block[1][L/2][0] = block + L/2 + 3;

   v[0] = 1;//S_ab
   v[1] = 0;
   v[2] = L/2;
   v[3] = 0;

   pph2s[block + L/2 + 3].push_back(v);

   s2pph[block + L/2 + 3][1][0][L/2][0] = pph;

   ++pph;

   //then the states containing k = L/2: only |k_a L/2 k_b> because for S_ab = 1 |k_a k_b L/2> has negative parity automatically
   for(int k_a = 1;k_a < L/2;++k_a){

      int k_b = L - k_a;

      v[0] = 1;
      v[1] = k_a;
      v[2] = L/2;
      v[3] = k_b;

      pph2s[block + L/2 + 3].push_back(v);

      s2pph[block + L/2 + 3][1][k_a][L/2][k_b] = pph;

      ++pph;

   }

   //the states not containing k = L/2 can be split up in two groups:

   //first k_a,k_b,k_c < L/2
   for(int k_a = 0;k_a < L/2;++k_a)
      for(int k_b = k_a + 1;k_b < L/2;++k_b)
         for(int k_c = 0;k_c < L/2;++k_c){

            if(k_a + k_b + k_c == L/2){

               v[0] = 1;
               v[1] = k_a;
               v[2] = k_b;
               v[3] = k_c;

               pph2s[block + L/2 + 3].push_back(v);

               s2pph[block + L/2 + 3][1][k_a][k_b][k_c] = pph;

               ++pph;

            }

         }

   //second k_a < k_b and k_c > L/2 with sum = L + L/2
   for(int k_a = 1;k_a < L;++k_a){

      if(k_a == L/2)
         k_a++;

      for(int k_b = k_a + 1;k_b < L;++k_b){

         if(k_b == L/2)
            k_b++;

         for(int k_c = L/2 + 1;k_c < L;++k_c){

            if(k_a + k_b + k_c == L + L/2){

               v[0] = 1;
               v[1] = k_a;
               v[2] = k_b;
               v[3] = k_c;

               pph2s[block + L/2 + 3].push_back(v);

               s2pph[block + L/2 + 3][1][k_a][k_b][k_c] = pph;

               ++pph;

            }

         }
      }
   }

   ++block;

   pph = 0;

   //now negative parity, K = L/2
   
   //S = 1/2
   block_char[block][0] = 0;//S
   block_char[block][1] = L/2;//K
   block_char[block][2] = 1;//p

   char_block[0][L/2][1] = block;

   //S_ab = 0:

   //for terms with k = L/2: only |k_a L/2 k_b> remains, |k_a k_b L/2> only has positive parity with S_ab = 0
   for(int k_a = 1;k_a < L/2;++k_a){

      int k_b = L - k_a;

      v[0] = 0;
      v[1] = k_a;
      v[2] = L/2;
      v[3] = k_b;

      pph2s[block].push_back(v);

      s2pph[block][0][k_a][L/2][k_b] = pph;

      ++pph;

   }

   //the states not containing k = L/2 can be split up in two groups:

   //first k_a,k_b,k_c < L/2
   for(int k_a = 0;k_a < L/2;++k_a)
      for(int k_b = k_a;k_b < L/2;++k_b)
         for(int k_c = 0;k_c < L/2;++k_c){

            if(k_a + k_b + k_c == L/2){

               v[0] = 0;
               v[1] = k_a;
               v[2] = k_b;
               v[3] = k_c;

               pph2s[block].push_back(v);

               s2pph[block][0][k_a][k_b][k_c] = pph;

               ++pph;

            }

         }

   //second k_a < k_b and k_c > L/2 with sum = L + L/2
   for(int k_a = 1;k_a < L;++k_a){

      if(k_a == L/2)
         k_a++;

      for(int k_b = k_a;k_b < L;++k_b){

         if(k_b == L/2)
            k_b++;

         for(int k_c = L/2 + 1;k_c < L;++k_c){

            if(k_a + k_b + k_c == L + L/2){

               v[0] = 0;
               v[1] = k_a;
               v[2] = k_b;
               v[3] = k_c;

               pph2s[block].push_back(v);

               s2pph[block][0][k_a][k_b][k_c] = pph;

               ++pph;

            }

         }
      }
   }

   //S_ab = 1:

   //for terms with k = L/2: two terms now
   for(int k_a = 1;k_a < L/2;++k_a){

      int k_b = L - k_a;

      v[0] = 1;
      v[1] = k_a;
      v[2] = k_b;
      v[3] = L/2;

      pph2s[block].push_back(v);

      s2pph[block][1][k_a][k_b][L/2] = pph;

      ++pph;

      v[0] = 1;
      v[1] = k_a;
      v[2] = L/2;
      v[3] = k_b;

      pph2s[block].push_back(v);

      s2pph[block][1][k_a][L/2][k_b] = pph;

      ++pph;

   }

   //the states not containing k = L/2 can be split up in two groups:

   //first k_a,k_b,k_c < L/2
   for(int k_a = 0;k_a < L/2;++k_a)
      for(int k_b = k_a + 1;k_b < L/2;++k_b)
         for(int k_c = 0;k_c < L/2;++k_c){

            if(k_a + k_b + k_c == L/2){

               v[0] = 1;
               v[1] = k_a;
               v[2] = k_b;
               v[3] = k_c;

               pph2s[block].push_back(v);

               s2pph[block][1][k_a][k_b][k_c] = pph;

               ++pph;

            }

         }

   //second k_a < k_b and k_c > L/2 with sum = L + L/2
   for(int k_a = 1;k_a < L;++k_a){

      if(k_a == L/2)
         k_a++;

      for(int k_b = k_a + 1;k_b < L;++k_b){

         if(k_b == L/2)
            k_b++;

         for(int k_c = L/2 + 1;k_c < L;++k_c){

            if(k_a + k_b + k_c == L + L/2){

               v[0] = 1;
               v[1] = k_a;
               v[2] = k_b;
               v[3] = k_c;

               pph2s[block].push_back(v);

               s2pph[block][1][k_a][k_b][k_c] = pph;

               ++pph;

            }

         }
      }
   }

   pph = 0;

   //S = 3/2, negative parity, K = L/2
   block_char[block + L/2 + 3][0] = 1;//S
   block_char[block + L/2 + 3][1] = L/2;//K
   block_char[block + L/2 + 3][2] = 1;//p

   char_block[1][L/2][1] = block + L/2 + 3;

   //only S_ab = 1:

   //for terms with k = L/2: the two terms can occur
   for(int k_a = 1;k_a < L/2;++k_a){

      int k_b = L - k_a;

      v[0] = 1;
      v[1] = k_a;
      v[2] = k_b;
      v[3] = L/2;

      pph2s[block + L/2 + 3].push_back(v);

      s2pph[block + L/2 + 3][1][k_a][k_b][L/2] = pph;

      ++pph;

      v[0] = 1;
      v[1] = k_a;
      v[2] = L/2;
      v[3] = k_b;

      pph2s[block + L/2 + 3].push_back(v);

      s2pph[block + L/2 + 3][1][k_a][L/2][k_b] = pph;

      ++pph;

   }

   //the states not containing k = L/2 can be split up in two groups:

   //first k_a,k_b,k_c < L/2
   for(int k_a = 0;k_a < L/2;++k_a)
      for(int k_b = k_a + 1;k_b < L/2;++k_b)
         for(int k_c = 0;k_c < L/2;++k_c){

            if(k_a + k_b + k_c == L/2){

               v[0] = 1;
               v[1] = k_a;
               v[2] = k_b;
               v[3] = k_c;

               pph2s[block + L/2 + 3].push_back(v);

               s2pph[block + L/2 + 3][1][k_a][k_b][k_c] = pph;

               ++pph;

            }

         }

   //second k_a < k_b and k_c > L/2 with sum = L + L/2
   for(int k_a = 1;k_a < L;++k_a){

      if(k_a == L/2)
         k_a++;

      for(int k_b = k_a + 1;k_b < L;++k_b){

         if(k_b == L/2)
            k_b++;

         for(int k_c = L/2 + 1;k_c < L;++k_c){

            if(k_a + k_b + k_c == L + L/2){

               v[0] = 1;
               v[1] = k_a;
               v[2] = k_b;
               v[3] = k_c;

               pph2s[block + L/2 + 3].push_back(v);

               s2pph[block + L/2 + 3][1][k_a][k_b][k_c] = pph;

               ++pph;

            }

         }
      }
   }

}

/**
 * deallocate the static lists
 */
void PPHM::clear(){

   delete [] pph2s;

   for(int B = 0;B < L + 6;++B){

      for(int S_ab = 0;S_ab < 2;++S_ab){

         for(int k_a = 0;k_a < L;++k_a){

            for(int k_b = 0;k_b < L;++k_b)
               delete [] s2pph[B][S_ab][k_a][k_b];

            delete [] s2pph[B][S_ab][k_a];

         }

         delete [] s2pph[B][S_ab];

      }

      delete [] s2pph[B];

   }

   delete [] s2pph;

   for(int B = 0;B < L + 6;++B)
      delete [] block_char[B];

   delete [] block_char;

   for(int S = 0;S < 2;++S){

      for(int K = 0;K <= L/2;++K)
         delete [] char_block[S][K];

      delete [] char_block[S];

   }

   delete [] char_block;

}

/**
 * standard constructor
 */
PPHM::PPHM() : BlockMatrix(L + 6) {

   //first K = 0: 4 blocks

   //positive parity
   this->setMatrixDim(0,pph2s[0].size(),2);//S = 1/2
   this->setMatrixDim(L/2 + 3,pph2s[L/2 + 3].size(),4);//S = 3/2

   //negative parity
   this->setMatrixDim(1,pph2s[1].size(),2);//S = 1/2
   this->setMatrixDim(L/2 + 4,pph2s[L/2 + 4].size(),4);//S = 3/2

   //then for 0 < K < L/2
   for(int K = 1;K < L/2;++K){

      this->setMatrixDim(K + 1,pph2s[K + 1].size(),4);//S = 1/2
      this->setMatrixDim(L/2 + K + 4,pph2s[L/2 + K + 4].size(),8);//S = 3/2

   }

   //and last for K = L/2: parity positive
   this->setMatrixDim(L/2 + 1,pph2s[L/2 + 1].size(),2);//S = 1/2
   this->setMatrixDim(L + 4,pph2s[L + 4].size(),4);//S = 3/2

   //negative parity
   this->setMatrixDim(L/2 + 2,pph2s[L/2 + 2].size(),2);//S = 1/2
   this->setMatrixDim(L + 5,pph2s[L + 5].size(),4);//S = 3/2

}

/**
 * copy constructor: constructs BlockMatrix object, and copies the content of the pphm_c blocks into it,
 * @param pphm_c PPHM object to be copied into (*this)
 */
PPHM::PPHM(const PPHM &pphm_c) : BlockMatrix(pphm_c) { }

/**
 * Destructor
 */
PPHM::~PPHM(){ }

/** 
 * @return nr of particles
 */
int PPHM::gN() const{

   return N;

}

/**
 * @return nr of sp orbitals
 */
int PPHM::gM() const{

   return M;

}

/**
 * @return nr of sites
 */
int PPHM::gL() const{

   return L;

}

/**
 * @param block block index
 * @return the spin of the block.
 */
int PPHM::gS(int block) const{

   return block_char[block][0];

}

/**
 * @param block block index
 * @return the momentum of the block.
 */
int PPHM::gK(int block) const{

   return block_char[block][1];

}

/**
 * @param block block index
 * @return the pph parity of the block.
 */
int PPHM::gp(int block) const{

   return block_char[block][2];

}

/**
 * access the elements of the matrix in sp mode, special symmetry and antisymmetry relations are automatically accounted for:\n\n
 * @param B The  blockindex
 * @param S_ab The intermediate spinquantumnumber of k_a and k_b.
 * @param k_a first sp index that forms the pph row index i together with k_b, k_c and S_ab in block B
 * @param k_b second sp index that forms the pph row index i together with k_a, k_c and S_ab in block B
 * @param k_c third sp index that forms the pph row index i together with k_a, k_b and S_ab in block B
 * @param S_de The intermediate spinquantumnumber of k_d and k_e.
 * @param k_d first sp index that forms the pph column index j together with k_e, k_z and S_de in block B
 * @param k_e second sp index that forms the pph column index j together with k_d, k_z and S_de in block B
 * @param k_z third sp index that forms the pph column index j together with k_d, k_e and S_de in block B
 * @return the number on place PPHM(B,i,j) with the right phase.
 */
/*
   double PPHM::operator()(int B,int S_ab,int k_a,int k_b,int k_c,int S_de,int k_d,int k_e,int k_z) const {

   int S = block_char[B][0];
   int K = block_char[B][1];

   return (*this)(S,K,S_ab,k_a,k_b,k_c,S_de,k_d,k_e,k_z);

   }
 */
/**
 * access the elements of the matrix in sp mode, special symmetry and antisymmetry relations are automatically accounted for:\n\n
 * @param S The pphm-spin index, when == 0 then access the block S = 1/2, for spinindex == 1 we access the S = 3/2.
 * @param K The pphm-momentum index
 * @param S_ab The intermediate spinquantumnumber of k_a and k_b.
 * @param k_a first sp index that forms the pph row index i together with k_b, k_c and S_ab in block B
 * @param k_b second sp index that forms the pph row index i together with k_a, k_c and S_ab in block B
 * @param k_c third sp index that forms the pph row index i together with k_a, k_b and S_ab in block B
 * @param S_de The intermediate spinquantumnumber of k_d and k_e.
 * @param k_d first sp index that forms the pph column index j together with k_e, k_z and S_de in block B
 * @param k_e second sp index that forms the pph column index j together with k_d, k_z and S_de in block B
 * @param k_z third sp index that forms the pph column index j together with k_d, k_e and S_de in block B
 * @return the number on place PPHM(B,i,j) with the right phase.
 */
/*
   double PPHM::operator()(int S,int K,int S_ab,int k_a,int k_b,int k_c,int S_de,int k_d,int k_e,int k_z) const {

//check the momentum
if( (k_a + k_b + k_c)%(M/2) != K)
return 0;

if( (k_d + k_e + k_z)%(M/2) != K)
return 0;

//associated blockindex
int B = (M/2)*S + K;

int i,j;

int phase_i = get_inco(S,K,S_ab,k_a,k_b,k_c,i);

if(phase_i == 0)
return 0;

int phase_j = get_inco(S,K,S_de,k_d,k_e,k_z,j);

if(phase_j == 0)
return 0;

return phase_i*phase_j* (*this)(B,i,j);

}
 */
/** 
 * Member function that gets the pph-index and phase corresponding to the sp indices S, K, S_ab, k_a, k_b, k_c.
 * @param S spin-index of the state, 0 -> S = 1/2, 1 -> S = 3/2
 * @param K momentum-index of the state, 0 <= K < M/2
 * @param S_ab intermediate spincoupling of k_a and k_b. = 0 or 1
 * @param k_a first sp orbital
 * @param k_b second sp orbital
 * @param k_c third sp orbital
 * @param i the corresponding pph index will be stored in this int after calling the function
 * @return the phase needed to get to a normal ordering of indices that corresponds to a pph index i
 */
/*
   int PPHM::get_inco(int S,int K,int S_ab,int k_a,int k_b,int k_c,int &i) const{

//block index:
int B = (M/2)*S + K;

if(S == 0){//S = 1/2

if(S_ab == 0){//symmetric in spatial sp's

if(k_a <= k_b)
i = s2pph[B][0][k_a][k_b][k_c];
else
i = s2pph[B][0][k_b][k_a][k_c];

return 1;

}
else{//antisymmetric in spatial sp's

if(k_a == k_b)
return 0;

if(k_a < k_b){

i = s2pph[B][1][k_a][k_b][k_c];

return 1;

}
else{

i = s2pph[B][1][k_b][k_a][k_c];

return -1;

}

}

}
else{//S = 3/2

if(S_ab == 0)//no possibile for S = 3/2
return 0;

if(k_a == k_b)//no possibile for S = 3/2
return 0;

if(k_a < k_b){

i = s2pph[B][0][k_a][k_b][k_c];

return 1;

}
else{

i = s2pph[B][0][k_b][k_a][k_c];

return -1;

}

}

}
 */
/**
 * The spincoupled, translationally invariant T2 map, maps a TPM onto a PPHM object. See notes for more info
 * be aware that the k_c and k_z in the T2 notation become -k_c and -k_z in TPM space (remember the G-map)
 * @param tpm input TPM matrix
 */
/*
   void PPHM::T(const TPM &tpm){

   SPM spm(1.0/(N - 1.0),tpm);

   int k_a,k_b,k_c,k_d,k_e,k_z;
   int S_ab,S_de;

   double norm_ab,norm_de;
   int sign_ab,sign_de;

//first the S = 1/2 blocks, these should be the most difficult ones.
for(int B = 0;B < M/2;++B){

for(int i = 0;i < gdim(B);++i){

S_ab = pph2s[B][i][0];

k_a = pph2s[B][i][1];
k_b = pph2s[B][i][2];

//change to tp-notation:
k_c = (-pph2s[B][i][3] + M/2)%(M/2);

sign_ab = 1 - 2*S_ab;

norm_ab = 1.0;

if(k_a == k_b)
norm_ab /= std::sqrt(2.0);

for(int j = i;j < gdim(B);++j){

S_de = pph2s[B][j][0];

k_d = pph2s[B][j][1];
k_e = pph2s[B][j][2];

//change to tp-notation:
k_z = (-pph2s[B][j][3] + M/2)%(M/2);

sign_de = 1 - 2*S_de;

norm_de = 1.0;

if(k_d == k_e)
norm_de /= std::sqrt(2.0);

//start the map: init
(*this)(B,i,j) = 0.0;

//sp term becomes diagonal here:
if(i == j)
(*this)(B,i,j) += spm[k_c];

//tp(1)
if(k_c == k_z)
if(S_ab == S_de)
(*this)(B,i,j) += tpm(S_ab,(k_a + k_b)%(M/2),k_a,k_b,k_d,k_e);

//tp(2)
if(k_a == k_d){

double ward = 0.0;

for(int J = 0;J < 2;++J)
for(int Z = 0;Z < 2;++Z)
ward += (2*J + 1.0) * (2*Z + 1.0) * _6j[J][S_ab] * _6j[J][S_de] * _6j[J][Z] * tpm(Z,(k_c + k_e)%(M/2),k_c,k_e,k_z,k_b);

ward *= norm_ab * norm_de * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) );

if(k_c == k_e)
   ward *= std::sqrt(2.0);

if(k_z == k_b)
   ward *= std::sqrt(2.0);

   (*this)(B,i,j) -= ward;

   }

//tp(3)
if(k_b == k_d){

   double ward = 0.0;

   for(int J = 0;J < 2;++J)
      for(int Z = 0;Z < 2;++Z)
         ward += (2*J + 1.0) * (2*Z + 1.0) * _6j[J][S_ab] * _6j[J][S_de] * _6j[J][Z] * tpm(Z,(k_c + k_e)%(M/2),k_c,k_e,k_z,k_a);

   ward *= norm_ab * norm_de * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) );

   if(k_c == k_e)
      ward *= std::sqrt(2.0);

   if(k_z == k_a)
      ward *= std::sqrt(2.0);

   (*this)(B,i,j) -= sign_ab * ward;

}

//tp(4)
if(k_a == k_e){

   double ward = 0.0;

   for(int J = 0;J < 2;++J)
      for(int Z = 0;Z < 2;++Z)
         ward += (2*J + 1.0) * (2*Z + 1.0) * _6j[J][S_ab] * _6j[J][S_de] * _6j[J][Z] * tpm(Z,(k_c + k_d)%(M/2),k_c,k_d,k_z,k_b);

   ward *= norm_ab * norm_de * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) );

   if(k_c == k_d)
      ward *= std::sqrt(2.0);

   if(k_z == k_b)
      ward *= std::sqrt(2.0);

   (*this)(B,i,j) -= sign_de * ward;

}

//tp(5)
if(k_b == k_e){

   double ward = 0.0;

   for(int J = 0;J < 2;++J)
      for(int Z = 0;Z < 2;++Z)
         ward += (2*J + 1.0) * (2*Z + 1.0) * _6j[J][S_ab] * _6j[J][S_de] * _6j[J][Z] * tpm(Z,(k_c + k_d)%(M/2),k_c,k_d,k_z,k_a);

   ward *= norm_ab * norm_de * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) );

   if(k_c == k_d)
      ward *= std::sqrt(2.0);

   if(k_z == k_a)
      ward *= std::sqrt(2.0);

   (*this)(B,i,j) -= sign_ab * sign_de * ward;

}

}

}

}

//the easier S = 3/2 part:
for(int B = M/2;B < M;++B){

   for(int i = 0;i < gdim(B);++i){

      k_a = pph2s[B][i][1];
      k_b = pph2s[B][i][2];

      //change to correct sp-momentum
      k_c = (-pph2s[B][i][3] + M/2)%(M/2);

      for(int j = i;j < gdim(B);++j){

         k_d = pph2s[B][j][1];
         k_e = pph2s[B][j][2];

         //change to correct sp-momentum
         k_z = (-pph2s[B][j][3] + M/2)%(M/2);

         //init
         (*this)(B,i,j) = 0.0;

         //sp part is diagonal
         if(i == j)
            (*this)(B,i,j) += spm[k_c];

         //tp(1)
         if(k_c == k_z)
            (*this)(B,i,j) += tpm(1,(k_a + k_b)%(M/2),k_a,k_b,k_d,k_e);

         //tp(2)
         if(k_a == k_d){

            double ward = 0.0;

            for(int Z = 0;Z < 2;++Z)
               ward += (2*Z + 1.0) * _6j[1][Z] * tpm(Z,(k_c + k_e)%(M/2),k_c,k_e,k_z,k_b);

            if(k_c == k_e)
               ward *= std::sqrt(2.0);

            if(k_z == k_b)
               ward *= std::sqrt(2.0);

            (*this)(B,i,j) -= ward;

         }

         //tp(3)
         if(k_b == k_d){

            double ward = 0.0;

            for(int Z = 0;Z < 2;++Z)
               ward += (2*Z + 1.0) * _6j[1][Z] * tpm(Z,(k_c + k_e)%(M/2),k_c,k_e,k_z,k_a);

            if(k_c == k_e)
               ward *= std::sqrt(2.0);

            if(k_z == k_a)
               ward *= std::sqrt(2.0);

            (*this)(B,i,j) += ward;

         }

         //tp(5)
         if(k_b == k_e){

            double ward = 0.0;

            for(int Z = 0;Z < 2;++Z)
               ward += (2*Z + 1.0) * _6j[1][Z] * tpm(Z,(k_c + k_d)%(M/2),k_c,k_d,k_z,k_a);

            if(k_c == k_d)
               ward *= std::sqrt(2.0);

            if(k_z == k_a)
               ward *= std::sqrt(2.0);

            (*this)(B,i,j) -= ward;

         }

      }

   }

}

this->symmetrize();

}
*/
ostream &operator<<(ostream &output,const PPHM &pphm_p){

   for(int B = 0;B < pphm_p.gnr();++B){

      output << "(" << pphm_p.gS(B) << "," << pphm_p.gK(B) << "," << pphm_p.gp(B) << ")\t" << pphm_p.gdim(B) << "\t" << pphm_p.gdeg(B) << std::endl;
      output << std::endl;

      for(int i = 0;i < pphm_p.gdim(B);++i)
         for(int j = 0;j < pphm_p.gdim(B);++j){

            output << i << "\t" << j << "\t|\t" << 

               pphm_p.pph2s[B][i][0] << "\t" << pphm_p.pph2s[B][i][1] << "\t" << pphm_p.pph2s[B][i][2] << "\t" << pphm_p.pph2s[B][i][3] << 

               "\t" << pphm_p.pph2s[B][j][0] << "\t" << pphm_p.pph2s[B][j][1] << "\t" << pphm_p.pph2s[B][j][2] << "\t" << pphm_p.pph2s[B][j][3] 

               << "\t" << pphm_p(B,i,j) << endl;

         }

      output << endl;

   }

   return output;

}

/**
 * static test function that prints the basis constructed
 */
void PPHM::print_basis(){

    for(int B = 0;B < L + 6;++B){

      cout << endl;
      cout << B << "\t" << gdim(B) << "\t" << gdeg(B) << "\t(" << block_char[B][0] << "," << block_char[B][1] << "," << 1 - 2*block_char[B][2] << ")\t" << endl;
      cout << endl;

      for(unsigned int i = 0;i < pph2s[B].size();++i)
         cout << i << "\t|\t" << pph2s[B][i][0] << "\t" << pph2s[B][i][1] << "\t" << pph2s[B][i][2] << "\t" << pph2s[B][i][3] << endl;

   }

}

/**
 * test function that print the total dimension of the object
 */
int PPHM::total_dim(){

   int tmp = 0;

   for(int B = 0;B < gnr();++B)
      tmp = gdim(B)*gdeg(B);

   return tmp;

}
