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

vector< vector<int> > *DPM::dp2s;
int *****DPM::s2dp;

double **DPM::_6j;

int **DPM::block_char;
int ***DPM::char_block;

int DPM::M;
int DPM::N;
int DPM::L;

/**
 * allocate and initialize the static lists and variables
 * @param L_in nr of sites
 * @param N_in nr of particles
 */
void DPM::init(int L_in,int N_in){

   L = L_in;
   N = N_in;

   M = 2*L;

   //allocate the lists
   dp2s = new vector< vector<int> > [L + 6];

   s2dp = new int **** [L + 6];

   for(int B = 0;B < L + 6;++B){

      s2dp[B] = new int *** [2];

      for(int S_ab = 0;S_ab < 2;++S_ab){

         s2dp[B][S_ab] = new int ** [L];

         for(int k_a = 0;k_a < L;++k_a){

            s2dp[B][S_ab][k_a] = new int * [L];

            for(int k_b = 0;k_b < L;++k_b)
               s2dp[B][S_ab][k_a][k_b] = new int [L];

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

   //initialize the lists

   //start with K = 0
   int block = 0;

   int dp = 0;

   vector<int> v(4);

   //positive parity

   //S = 1/2: first part S_ab = 0, k_a == k_b != k_c
   block_char[block][0] = 0;//S = 1/2
   block_char[block][1] = 0;//K = 0
   block_char[block][2] = 0;//p = positive,

   char_block[0][0][0] = block;

   for(int k_a = 1;k_a <= L/2;++k_a){

      for(int k_c = 0;k_c < k_a;++k_c){

         if( (2*k_a + k_c)%L == 0 ){

            v[0] = 0;//S_ab
            v[1] = k_a;
            v[2] = k_a;
            v[3] = k_c;

            dp2s[block].push_back(v);

            s2dp[block][0][k_a][k_a][k_c] = dp;

            ++dp;

         }

      }

      for(int k_c = k_a + 1;k_c < L;++k_c){

         if( (2*k_a + k_c)%L == 0 ){

            v[0] = 0;//S_ab
            v[1] = k_a;
            v[2] = k_a;
            v[3] = k_c;

            dp2s[block].push_back(v);

            s2dp[block][0][k_a][k_a][k_c] = dp;

            ++dp;

         }

      }

   }

   //now S_ab = 0, k_a < k_b < k_c, but start with k_a == 0
   for(int k_b = 1;k_b < L/2;++k_b){

      int k_c = L - k_b;

      v[0] = 0;//S_ab
      v[1] = 0;
      v[2] = k_b;
      v[3] = k_c;

      dp2s[block].push_back(v);

      s2dp[block][0][0][k_b][k_c] = dp;

      ++dp;

   }

   //then, for both S_ab = 0/1, k_a < k_b < k_c
   for(int S_ab = 0;S_ab < 2;++S_ab){

      for(int k_a = 1;k_a < L;++k_a)
         for(int k_b = k_a + 1;k_b < L;++k_b)
            for(int k_c = k_b + 1;k_c < L;++k_c){

               if(k_a + k_b + k_c == L){//not modulo! special case

                  v[0] = S_ab;
                  v[1] = k_a;
                  v[2] = k_b;
                  v[3] = k_c;

                  dp2s[block].push_back(v);

                  s2dp[block][S_ab][k_a][k_b][k_c] = dp;

                  ++dp;

               }

            }

   }

   //now S = 3/2: only witouth k_a = 0 in positive parity
   dp = 0;

   block_char[block + L/2 + 3][0] = 1;//S = 3/2
   block_char[block + L/2 + 3][1] = 0;//K = 0
   block_char[block + L/2 + 3][2] = 0;//p = positive,

   char_block[1][0][0] = block + L/2 + 3;

   //only S_ab == 1!
   for(int k_a = 1;k_a < L;++k_a)
      for(int k_b = k_a + 1;k_b < L;++k_b)
         for(int k_c = k_b + 1;k_c < L;++k_c){

            if(k_a + k_b + k_c == L){//not modulo! special case

               v[0] = 1;
               v[1] = k_a;
               v[2] = k_b;
               v[3] = k_c;

               dp2s[block + L/2 + 3].push_back(v);

               s2dp[block + L/2 + 3][1][k_a][k_b][k_c] = dp;

               ++dp;

            }

         }

   ++block;

   //Allright! Now, K=0 S=1/2: negative parity
   dp = 0;

   block_char[block][0] = 0;//S = 1/2
   block_char[block][1] = 0;//K = 0
   block_char[block][2] = 1;//p = negative,

   char_block[0][0][1] = block;

   //start again with S_ab = 0 and k_a == k_b != k_c
   for(int k_a = 1;k_a < L/2;++k_a){

      for(int k_c = 0;k_c < k_a;++k_c){

         if( (2*k_a + k_c)%L == 0 ){

            v[0] = 0;
            v[1] = k_a;
            v[2] = k_a;
            v[3] = k_c;

            dp2s[block].push_back(v);

            s2dp[block][0][k_a][k_a][k_c] = dp;

            ++dp;

         }

      }

      for(int k_c = k_a + 1;k_c < L;++k_c){

         if( (2*k_a + k_c)%L == 0 ){

            v[0] = 0;
            v[1] = k_a;
            v[2] = k_a;
            v[3] = k_c;

            dp2s[block].push_back(v);

            s2dp[block][0][k_a][k_a][k_c] = dp;

            ++dp;

         }

      }

   }

   //then S_ab == 0 k_a < k_b < k_c, without 0
   for(int k_a = 1;k_a < L;++k_a)
      for(int k_b = k_a + 1;k_b < L;++k_b)
         for(int k_c = k_b + 1;k_c < L;++k_c){

            if(k_a + k_b + k_c == L){//not modulo! special case

               v[0] = 0;
               v[1] = k_a;
               v[2] = k_b;
               v[3] = k_c;

               dp2s[block].push_back(v);

               s2dp[block][0][k_a][k_b][k_c] = dp;

               ++dp;

            }

         }

   //then the negative parity with 0, which are S_ab = 1 states!
   for(int k_b = 1;k_b < L/2;++k_b){

      int k_c = L - k_b;

      v[0] = 1;
      v[1] = 0;
      v[2] = k_b;
      v[3] = k_c;

      dp2s[block].push_back(v);

      s2dp[block][1][0][k_b][k_c] = dp;

      ++dp;

   }

   //then the rest with S_ab == 1: k_a < k_b < k_c without 0
   for(int k_a = 1;k_a < L;++k_a)
      for(int k_b = k_a + 1;k_b < L;++k_b)
         for(int k_c = k_b + 1;k_c < L;++k_c){

            if(k_a + k_b + k_c == L){//not modulo! special case

               v[0] = 1;
               v[1] = k_a;
               v[2] = k_b;
               v[3] = k_c;

               dp2s[block].push_back(v);

               s2dp[block][1][k_a][k_b][k_c] = dp;

               ++dp;

            }

         }

   //Now S = 3/2 , K = 0 , with negative parity
   dp = 0;

   block_char[block + L/2 + 3][0] = 1;//S = 3/2
   block_char[block + L/2 + 3][1] = 0;//K = 0
   block_char[block + L/2 + 3][2] = 1;//p = negative,

   char_block[1][0][1] = block + L/2 + 3;

   //first with a k=0: 
   for(int k_b = 1;k_b < L/2;++k_b){

      int k_c = L - k_b;

      v[0] = 1;
      v[1] = 0;
      v[2] = k_b;
      v[3] = k_c;

      dp2s[block + L/2 + 3].push_back(v);

      s2dp[block + L/2 + 3][1][0][k_b][k_c] = dp;

      ++dp;

   }

   //then the rest k_a < k_b < k_c
   for(int k_a = 1;k_a < L;++k_a)
      for(int k_b = k_a + 1;k_b < L;++k_b)
         for(int k_c = k_b + 1;k_c < L;++k_c){

            if(k_a + k_b + k_c == L){//not modulo! special case

               v[0] = 1;
               v[1] = k_a;
               v[2] = k_b;
               v[3] = k_c;

               dp2s[block + L/2 + 3].push_back(v);

               s2dp[block + L/2 + 3][1][k_a][k_b][k_c] = dp;

               ++dp;

            }

         }

   ++block;

   //Now 0 < K < L/2: degeneracy between parity + and -

   for(int K = 1;K < L/2;++K){

      //first S = 1/2:
      dp = 0;

      block_char[block][0] = 0;//S = 1/2
      block_char[block][1] = K;
      block_char[block][2] = 0;

      char_block[0][K][0] = block;
      char_block[0][K][1] = block;

      //first S_ab = 0 and k_a = k_b != k_c
      for(int k_a = 0;k_a < L;++k_a){

         for(int k_c = 0;k_c < k_a;++k_c){

            if( (2*k_a + k_c)%L == K ){

               v[0] = 0;
               v[1] = k_a;
               v[2] = k_a;
               v[3] = k_c;

               dp2s[block].push_back(v);

               s2dp[block][0][k_a][k_a][k_c] = dp;

               ++dp;

            }

         }

         for(int k_c = k_a + 1;k_c < L;++k_c){

            if( (2*k_a + k_c)%L == K ){

               v[0] = 0;
               v[1] = k_a;
               v[2] = k_a;
               v[3] = k_c;

               dp2s[block].push_back(v);

               s2dp[block][0][k_a][k_a][k_c] = dp;

               ++dp;

            }

         }

      }

      //then, for S_ab = 0/1, k_a < k_b < k_c
      for(int S_ab = 0;S_ab < 2;++S_ab){

         for(int k_a = 0;k_a < L;++k_a)
            for(int k_b = k_a + 1;k_b < L;++k_b)
               for(int k_c = k_b + 1;k_c < L;++k_c){

                  if( (k_a + k_b + k_c)%L == K ){

                     v[0] = S_ab;
                     v[1] = k_a;
                     v[2] = k_b;
                     v[3] = k_c;

                     dp2s[block].push_back(v);

                     s2dp[block][S_ab][k_a][k_b][k_c] = dp;

                     ++dp;

                  }

               }

      }

      //now S = 3/2
      dp = 0;

      block_char[block + L/2 + 3][0] = 1;//S = 3/2
      block_char[block + L/2 + 3][1] = K;
      block_char[block + L/2 + 3][2] = 0;

      char_block[1][K][0] = block + L/2 + 3;
      char_block[1][K][1] = block + L/2 + 3;

      for(int k_a = 0;k_a < L;++k_a)
         for(int k_b = k_a + 1;k_b < L;++k_b)
            for(int k_c = k_b + 1;k_c < L;++k_c){

               if( (k_a + k_b + k_c)%L == K ){

                  v[0] = 1;
                  v[1] = k_a;
                  v[2] = k_b;
                  v[3] = k_c;

                  dp2s[block + L/2 + 3].push_back(v);

                  s2dp[block + L/2 + 3][1][k_a][k_b][k_c] = dp;

                  ++dp;

               }

            }

      ++block;

   }

   //now the last ones: K = L/2 positive parity
   dp = 0;

   block_char[block][0] = 0;//S = 1/2
   block_char[block][1] = L/2;
   block_char[block][2] = 0;

   char_block[0][L/2][0] = block;

   //again first S_ab = 0, k_a == k_b != k_c

   //for positive parity one extra state: |0 0 L/2>
   v[0] = 0;//S_ab
   v[1] = 0;//k_a
   v[2] = 0;//k_b
   v[3] = L/2;//k_c

   dp2s[block].push_back(v);

   s2dp[block][0][0][0][L/2] = dp;

   ++dp;

   for(int k_a = 1;k_a < L/2;++k_a){

      for(int k_c = 0;k_c < k_a;++k_c){

         if( (2*k_a + k_c)%L == L/2 ){

            v[0] = 0;//S_ab
            v[1] = k_a;//k_a
            v[2] = k_a;//k_b
            v[3] = k_c;//k_c

            dp2s[block].push_back(v);

            s2dp[block][0][k_a][k_a][k_c] = dp;

            ++dp;

         }

      }

      for(int k_c = k_a + 1;k_c < L;++k_c){

         if( (2*k_a + k_c)%L == L/2 ){

            v[0] = 0;//S_ab
            v[1] = k_a;//k_a
            v[2] = k_a;//k_b
            v[3] = k_c;//k_c

            dp2s[block].push_back(v);

            s2dp[block][0][k_a][k_a][k_c] = dp;

            ++dp;

         }

      }

   }

   //then, for S_ab = 0,1: k_a < k_b < k_c: this is quite involved.
   for(int S_ab = 0;S_ab < 2;++S_ab){

      //first with k_a = 0, k_b < k_c < L/2
      for(int k_b = 1;k_b < L/2;++k_b)
         for(int k_c = k_b + 1;k_c < L/2;++k_c){

            if( (k_b + k_c) == L/2 ){

               v[0] = S_ab;//S_ab
               v[1] = 0;//k_a
               v[2] = k_b;//k_b
               v[3] = k_c;//k_c

               dp2s[block].push_back(v);

               s2dp[block][S_ab][0][k_b][k_c] = dp;

               ++dp;

            }

         }

      //then k_a < k_b < k_c
      for(int k_a = 1;k_a < L/2;++k_a){

         for(int k_b = k_a + 1;k_b < L/2;++k_b)
            for(int k_c = k_b + 1;k_c < L/2;++k_c){

               if( (k_a + k_b + k_c) == L/2){

                  v[0] = S_ab;//S_ab
                  v[1] = k_a;//k_a
                  v[2] = k_b;//k_b
                  v[3] = k_c;//k_c

                  dp2s[block].push_back(v);

                  s2dp[block][S_ab][k_a][k_b][k_c] = dp;

                  ++dp;

               }

            }

         //only S_ab == 0 has k's equal to L/2 for positive parity
         for(int k_b = L/2 + S_ab;k_b < L;++k_b)
            for(int k_c = k_b + 1;k_c < L;++k_c){

               if( (k_a + k_b + k_c)%L == L/2 ){

                  v[0] = S_ab;//S_ab
                  v[1] = k_a;//k_a
                  v[2] = k_b;//k_b
                  v[3] = k_c;//k_c

                  dp2s[block].push_back(v);

                  s2dp[block][S_ab][k_a][k_b][k_c] = dp;

                  ++dp;

               }

            }

      }

   }
   
   //then K = L/2, S = 3/2 positive parity:
   dp = 0;

   block_char[block + L/2 + 3][0] = 1;//S = 3/2
   block_char[block + L/2 + 3][1] = L/2;
   block_char[block + L/2 + 3][2] = 0;

   char_block[1][L/2][0] = block + L/2 + 3;

   //first with k_a = 0, k_b < k_c < L/2
   for(int k_b = 1;k_b < L/2;++k_b)
      for(int k_c = k_b + 1;k_c < L/2;++k_c){

         if( (k_b + k_c) == L/2 ){

            v[0] = 1;//S_ab
            v[1] = 0;//k_a
            v[2] = k_b;//k_b
            v[3] = k_c;//k_c

            dp2s[block + L/2 + 3].push_back(v);

            s2dp[block + L/2 + 3][1][0][k_b][k_c] = dp;

            ++dp;

         }

      }

   //then k_a < k_b < k_c
   for(int k_a = 1;k_a < L/2;++k_a){

      for(int k_b = k_a + 1;k_b < L/2;++k_b)
         for(int k_c = k_b + 1;k_c < L/2;++k_c){

            if( (k_a + k_b + k_c) == L/2){

               v[0] = 1;//S_ab
               v[1] = k_a;//k_a
               v[2] = k_b;//k_b
               v[3] = k_c;//k_c

               dp2s[block + L/2 + 3].push_back(v);

               s2dp[block + L/2 + 3][1][k_a][k_b][k_c] = dp;

               ++dp;

            }

         }

      //only negative parity has k's equal to L/2
      for(int k_b = L/2 + 1;k_b < L;++k_b)
         for(int k_c = k_b + 1;k_c < L;++k_c){

            if( (k_a + k_b + k_c)%L == L/2 ){

               v[0] = 1;//S_ab
               v[1] = k_a;//k_a
               v[2] = k_b;//k_b
               v[3] = k_c;//k_c

               dp2s[block + L/2 + 3].push_back(v);

               s2dp[block + L/2 + 3][1][k_a][k_b][k_c] = dp;

               ++dp;

            }

         }

   }

   block++;

   //K = L/2, negative parity, start with S = 1/2
   dp = 0;

   block_char[block][0] = 0;//S = 1/2
   block_char[block][1] = L/2;
   block_char[block][2] = 1;//negative parity

   char_block[0][L/2][1] = block;

   //again first S_ab = 0, k_a == k_b != k_c
   for(int k_a = 1;k_a < L/2;++k_a){

      for(int k_c = 0;k_c < k_a;++k_c){

         if( (2*k_a + k_c)%L == L/2 ){

            v[0] = 0;//S_ab
            v[1] = k_a;//k_a
            v[2] = k_a;//k_b
            v[3] = k_c;//k_c

            dp2s[block].push_back(v);

            s2dp[block][0][k_a][k_a][k_c] = dp;

            ++dp;

         }

      }

      for(int k_c = k_a + 1;k_c < L;++k_c){

         if( (2*k_a + k_c)%L == L/2 ){

            v[0] = 0;//S_ab
            v[1] = k_a;//k_a
            v[2] = k_a;//k_b
            v[3] = k_c;//k_c

            dp2s[block].push_back(v);

            s2dp[block][0][k_a][k_a][k_c] = dp;

            ++dp;

         }

      }

   }

   //then k_a < k_b < k_c:
   for(int S_ab = 0;S_ab < 2;++S_ab){

      //first with k_a = 0, k_b < k_c < L/2
      for(int k_b = 1;k_b < L/2;++k_b)
         for(int k_c = k_b + 1;k_c < L/2;++k_c){

            if( (k_b + k_c) == L/2 ){

               v[0] = S_ab;//S_ab
               v[1] = 0;//k_a
               v[2] = k_b;//k_b
               v[3] = k_c;//k_c

               dp2s[block].push_back(v);

               s2dp[block][S_ab][0][k_b][k_c] = dp;

               ++dp;

            }

         }

      //then k_a < k_b < k_c
      for(int k_a = 1;k_a < L/2;++k_a){

         for(int k_b = k_a + 1;k_b < L/2;++k_b)
            for(int k_c = k_b + 1;k_c < L/2;++k_c){

               if( (k_a + k_b + k_c) == L/2){

                  v[0] = S_ab;//S_ab
                  v[1] = k_a;//k_a
                  v[2] = k_b;//k_b
                  v[3] = k_c;//k_c

                  dp2s[block].push_back(v);

                  s2dp[block][S_ab][k_a][k_b][k_c] = dp;

                  ++dp;

               }

            }

         //only S_ab == 1 has k's equal to L/2 for negative parity
         for(int k_b = L/2 + (1 - S_ab);k_b < L;++k_b)
            for(int k_c = k_b + 1;k_c < L;++k_c){

               if( (k_a + k_b + k_c)%L == L/2 ){

                  v[0] = S_ab;//S_ab
                  v[1] = k_a;//k_a
                  v[2] = k_b;//k_b
                  v[3] = k_c;//k_c

                  dp2s[block].push_back(v);

                  s2dp[block][S_ab][k_a][k_b][k_c] = dp;

                  ++dp;

               }

            }

      }

   }

   //at last we have the last block, K = L/2, S = 3/2, negative parity
   dp = 0;

   block_char[block + L/2 + 3][0] = 1;//S = 3/2
   block_char[block + L/2 + 3][1] = L/2;
   block_char[block + L/2 + 3][2] = 1;

   char_block[1][L/2][1] = block + L/2 + 3;

   //first with k_a = 0, k_b < k_c < L/2
   for(int k_b = 1;k_b < L/2;++k_b)
      for(int k_c = k_b + 1;k_c < L/2;++k_c){

         if( (k_b + k_c) == L/2 ){

            v[0] = 1;//S_ab
            v[1] = 0;//k_a
            v[2] = k_b;//k_b
            v[3] = k_c;//k_c

            dp2s[block + L/2 + 3].push_back(v);

            s2dp[block + L/2 + 3][1][0][k_b][k_c] = dp;

            ++dp;

         }

      }

   //then k_a < k_b < k_c
   for(int k_a = 1;k_a < L/2;++k_a){

      for(int k_b = k_a + 1;k_b < L/2;++k_b)
         for(int k_c = k_b + 1;k_c < L/2;++k_c){

            if( (k_a + k_b + k_c) == L/2){

               v[0] = 1;//S_ab
               v[1] = k_a;//k_a
               v[2] = k_b;//k_b
               v[3] = k_c;//k_c

               dp2s[block + L/2 + 3].push_back(v);

               s2dp[block + L/2 + 3][1][k_a][k_b][k_c] = dp;

               ++dp;

            }

         }

      //only negative parity has k's equal to L/2
      for(int k_b = L/2;k_b < L;++k_b)
         for(int k_c = k_b + 1;k_c < L;++k_c){

            if( (k_a + k_b + k_c)%L == L/2 ){

               v[0] = 1;//S_ab
               v[1] = k_a;//k_a
               v[2] = k_b;//k_b
               v[3] = k_c;//k_c

               dp2s[block + L/2 + 3].push_back(v);

               s2dp[block + L/2 + 3][1][k_a][k_b][k_c] = dp;

               ++dp;

            }

         }

   }

}

/**
 * deallocates the static lists
 */
void DPM::clear(){

   delete [] dp2s;

   for(int B = 0;B < L + 6;++B){

      for(int S_ab = 0;S_ab < 2;++S_ab){

         for(int k_a = 0;k_a < L;++k_a){

            for(int k_b = 0;k_b < L;++k_b)
               delete [] s2dp[B][S_ab][k_a][k_b];

            delete [] s2dp[B][S_ab][k_a];

         }

         delete [] s2dp[B][S_ab];

      }

      delete [] s2dp[B];

   }

   delete [] s2dp;

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
 * standard constructor: constructs BlockMatrix object with 
 * @param M nr of sp orbitals
 * @param N nr of particles
 */
DPM::DPM() : BlockMatrix(M) {
   /*
      this->N = N;
      this->M = M;

      if(counter == 0)
      set_dimlist(M);

   //set the dimension and the degeneracies of the blocks
   for(int B = 0;B < M/2;++B){

   this->setMatrixDim(B,dimlist[B],2);
   this->setMatrixDim(B + M/2,dimlist[B + M/2],4);

   }

   if(counter == 0)//make the lists
   construct_lists();

   ++counter;
    */
}

/**
 * copy constructor: constructs BlockMatrix object with 2 * (M/2) blocks, on (M/2) for S=1/2 and (M/2) for S=3/2, and copies the content of the dpm_c blocks into it,
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and dp basis.
 * @param dpm_c DPM to be copied into (*this)
 */
DPM::DPM(const DPM &dpm_c) : BlockMatrix(dpm_c) {
   /*
      this->N = dpm_c.gN();
      this->M = dpm_c.gM();

      if(counter == 0)
      construct_lists();

      ++counter;
    */
}

/**
 * destructor: if counter == 1 the memory for the static lists dp2s en s2dp will be deleted.
 */
DPM::~DPM(){
   /*
      if(counter == 1){

   //first delete S = 1/2 blocks
   for(int B = 0;B < M/2;++B){

   for(int S_ab = 0;S_ab < 2;++S_ab){

   for(int k_a = 0;k_a < M/2;++k_a){

   for(int k_b = 0;k_b < M/2;++k_b)
   delete [] s2dp[B][S_ab][k_a][k_b];

   delete [] s2dp[B][S_ab][k_a];

   }

   delete [] s2dp[B][S_ab];

   }

   delete [] s2dp[B];

   }


   //then the S = 3/2 part
   for(int B = M/2;B < M;++B){

   for(int k_a = 0;k_a < M/2;++k_a){

   for(int k_b = 0;k_b < M/2;++k_b)
   delete [] s2dp[B][0][k_a][k_b];

   delete [] s2dp[B][0][k_a];

   }

   delete [] s2dp[B][0];

   delete [] s2dp[B];

   }

   delete [] s2dp;

   //now delete dp2s 
   for(int B = 0;B < gnr();++B){

   for(int i = 0;i < gdim(B);++i)
   delete [] dp2s[B][i];

   delete [] dp2s[B];

   }

   delete [] dp2s;

   for(int B = 0;B < gnr();++B)
   delete [] block_char[B];

   delete [] block_char;

   for(int S = 0;S < 2;++S)
   delete [] _6j[S];

   delete [] _6j;

   }

   --counter;
   */
}

/** 
 * Function that allocates and initializes the lists needed in the program, called when the first DPM object is constructed,
 */
/*
   void DPM::construct_lists(){

//first allocation
dp2s = new int ** [gnr()];//two total spinblocks

for(int B = 0;B < gnr();++B){

dp2s[B] = new int * [gdim(B)];//dimension of the blocks

for(int i = 0;i < gdim(B);++i)
dp2s[B][i] = new int [4];//amount of information stored for a blockindex: (S_ab,a,b,c)

}

s2dp = new int **** [gnr()];//nr of blocks

for(int B = 0;B < M/2;++B){//first loop over the S = 1/2 blocks

s2dp[B] = new int *** [2];//for the S = 1/2, we have that S_ab can be 0 or 1

for(int S_ab = 0;S_ab < 2;++S_ab){//loop and allocate

s2dp[B][S_ab] = new int ** [M/2];

for(int k_a = 0;k_a < M/2;++k_a){

s2dp[B][S_ab][k_a] = new int * [M/2];

for(int k_b = 0;k_b < M/2;++k_b)
s2dp[B][S_ab][k_a][k_b] = new int [M/2];

}

}

}

for(int B = M/2;B < M;++B){//then loop over the S = 3/2 blocks

s2dp[B] = new int *** [1];//for the S = 3/2, we have that S_ab can be only 1

s2dp[B][0] = new int ** [M/2];

for(int k_a = 0;k_a < M/2;++k_a){//loop and allocate

s2dp[B][0][k_a] = new int * [M/2];

for(int k_b = 0;k_b < M/2;++k_b)
s2dp[B][0][k_a][k_b] = new int [M/2];

}

}

//allocate the block_char list
block_char = new int * [M];

for(int B = 0;B < M;++B)
block_char[B] = new int [2];

//initialize:
for(int i = 0;i < M/2;++i){

block_char[i][0] = 0;//S (by 0 I mean 1/2)
block_char[i][1] = i;//K

block_char[i + M/2][0] = 1;//S (by 1 I mean 3/2)
block_char[i + M/2][1] = i;//K

}

//initialize the lists
int teller;

//first S == 1/2
for(int B = 0;B < M/2;++B){

   teller = 0;

   for(int k_a = 0;k_a < M/2;++k_a){//S_ab == 0 and k_a == k_b != c

      for(int k_c = 0;k_c < k_a;++k_c){

         if( (2*k_a + k_c)%(M/2) == block_char[B][1]){

            s2dp[B][0][k_a][k_a][k_c] = teller;

            dp2s[B][teller][0] = 0;//S_ab

            dp2s[B][teller][1] = k_a;
            dp2s[B][teller][2] = k_a;
            dp2s[B][teller][3] = k_c;

            ++teller;

         }

      }

      for(int k_c = k_a + 1;k_c < M/2;++k_c){

         if( (2*k_a + k_c)%(M/2) == block_char[B][1]){

            s2dp[B][0][k_a][k_a][k_c] = teller;

            dp2s[B][teller][0] = 0;//S_ab

            dp2s[B][teller][1] = k_a;
            dp2s[B][teller][2] = k_a;
            dp2s[B][teller][3] = k_c;

            ++teller;

         }

      }

   }

   //S=1/2 S_ab=0 a != b != c
   for(int k_a = 0;k_a < M/2;++k_a)
      for(int k_b = k_a + 1;k_b < M/2;++k_b)
         for(int k_c = k_b + 1;k_c < M/2;++k_c){

            if( (k_a + k_b + k_c)%(M/2) == block_char[B][1] ){

               s2dp[B][0][k_a][k_b][k_c] = teller;

               dp2s[B][teller][0] = 0;//S_ab

               dp2s[B][teller][1] = k_a;
               dp2s[B][teller][2] = k_b;
               dp2s[B][teller][3] = k_c;

               ++teller;

            }

         }

   //S == 0, S_ab == 1, a != b != c
   for(int k_a = 0;k_a < M/2;++k_a)
      for(int k_b = k_a + 1;k_b < M/2;++k_b)
         for(int k_c = k_b + 1;k_c < M/2;++k_c){

            if( (k_a + k_b + k_c)%(M/2) == block_char[B][1] ){

               s2dp[B][1][k_a][k_b][k_c] = teller;

               dp2s[B][teller][0] = 1;//S_ab

               dp2s[B][teller][1] = k_a;
               dp2s[B][teller][2] = k_b;
               dp2s[B][teller][3] = k_c;

               ++teller;

            }

         }

}

for(int B = M/2;B < M;++B){//now for the S=3/2 blocks:

   teller = 0;

   //S == 3/2, S_ab == 1, a != b != c
   for(int k_a = 0;k_a < M/2;++k_a)
      for(int k_b = k_a + 1;k_b < M/2;++k_b)
         for(int k_c = k_b + 1;k_c < M/2;++k_c){

            if( (k_a + k_b + k_c)%(M/2) == block_char[B][1] ){

               s2dp[B][0][k_a][k_b][k_c] = teller;//watch out! 0 S_ab index means S_ab == 1 here!

               dp2s[B][teller][0] = 1;//S_ab

               dp2s[B][teller][1] = k_a;
               dp2s[B][teller][2] = k_b;
               dp2s[B][teller][3] = k_c;

               ++teller;

            }

         }

}

//allocate
_6j = new double * [2];

for(int S = 0;S < 2;++S)
_6j[S] = new double [2];

//initialize
_6j[0][0] = -0.5;
_6j[0][1] = 0.5;
_6j[1][0] = 0.5;
_6j[1][1] = 1.0/6.0;

}
*/
/**
 * static function that calculates the dimensions of the different blocks.
 * @param M the dimension of sp space
 */
/*
   void DPM::set_dimlist(int M){

   dimlist =  new int [M];

   int teller;

   for(int B = 0;B < M/2;++B){//loop over the blocks

   teller = 0;

   for(int k_a = 0;k_a < M/2;++k_a){//S_ab == 0 and k_a == k_b != c

   for(int k_c = 0;k_c < k_a;++k_c){

   if( (2*k_a + k_c)%(M/2) == B)
   ++teller;

   }

   for(int k_c = k_a + 1;k_c < M/2;++k_c){

   if( (2*k_a + k_c)%(M/2) == B)
   ++teller;

   }

   }

//S=1/2 S_ab=0 a != b != c
for(int k_a = 0;k_a < M/2;++k_a)
for(int k_b = k_a + 1;k_b < M/2;++k_b)
for(int k_c = k_b + 1;k_c < M/2;++k_c)
if( (k_a + k_b + k_c)%(M/2) == B )
++teller;

//S == 1/2, S_ab == 1, a != b != c
for(int k_a = 0;k_a < M/2;++k_a)
for(int k_b = k_a + 1;k_b < M/2;++k_b)
for(int k_c = k_b + 1;k_c < M/2;++k_c)
if( (k_a + k_b + k_c)%(M/2) == B)
++teller;

dimlist[B] = teller;

}

for(int B = M/2;B < M;++B){

teller = 0;

//S == 3/2, S_ab == 1, a != b != c
for(int k_a = 0;k_a < M/2;++k_a)
for(int k_b = k_a + 1;k_b < M/2;++k_b)
for(int k_c = k_b + 1;k_c < M/2;++k_c)
if( (k_a + k_b + k_c)%(M/2) == B - M/2 )
++teller;

dimlist[B] = teller;

}

}
 */
/**
 * @return number of particles
 */
int DPM::gN() const {

   return N;

}

/**
 * @return number of single particle oribals
 */
int DPM::gM() const{

   return M;

}

/**
 * @return number of sites
 */
int DPM::gL() const{

   return L;

}

/**
 * access the elements of the matrix in sp mode, special symmetry and antisymmetry relations are automatically accounted for:\n\n
 * DPM(B,S_ab,k_a,k_b,k_c,S_de,k_d,k_e,k_z) = sum_S_ac (some terms dependent on spin) DPM(B,S_ac,k_a,k_c,k_b,S_de,k_d,k_e,k_z) etc...
 * @param B The block index
 * @param S_ab The intermediate spinquantumnumber of k_a and k_b.
 * @param k_a first sp index that forms the dp row index i of block B together with k_b, k_c and S_ab
 * @param k_b second sp index that forms the dp row index i of block B together with k_a, k_c and S_ab
 * @param k_c third sp index that forms the dp row index i of block B together with k_a, k_b and S_ab
 * @param S_de The intermediate spinquantumnumber of k_d and k_e.
 * @param k_d first sp index that forms the dp column index j of block B together with k_e, k_z and S_de
 * @param k_e second sp index that forms the dp column index j of block B together with k_d, k_z and S_de
 * @param k_z third sp index that forms the dp column index j of block B together with k_d, k_e and S_de
 * @return the number on place DPM(B,i,j) with the right phase and forefactor.
 */
/*
   double DPM::operator()(int B,int S_ab,int k_a,int k_b,int k_c,int S_de,int k_d,int k_e,int k_z) const {

   int S = block_char[B][0];
   int K = block_char[B][1];

   return (*this)(S,K,S_ab,k_a,k_b,k_c,S_de,k_d,k_e,k_z);

   }
 */
/**
 * access the elements of the matrix in sp mode, special symmetry and antisymmetry relations are automatically accounted for:\n\n
 * DPM(S,K,S_ab,k_a,k_b,k_c,S_de,k_d,k_e,k_z) = sum_S_ac (some terms dependent on spin) DPM(S,K,S_ac,k_a,k_c,k_b,S_de,k_d,k_e,k_z) etc...
 * @param S dp-spin block quantumnumber: S = 0 means S == 1/2 and S = 1 means S == 3/2 (for simplicity)
 * @param K dp-momentum block quantumnumber: K = 0 -> M-1
 * @param S_ab The intermediate spinquantumnumber of k_a and k_b.
 * @param k_a first sp index that forms the dp row index i together with k_b, k_c and S_ab, with dp spin and momentum SK
 * @param k_b second sp index that forms the dp row index i together with k_a, k_c and S_ab, with dp spin and momentum SK
 * @param k_c third sp index that forms the dp row index i together with k_a, k_b and S_ab, with dp spin and momentum SK
 * @param S_de The intermediate spinquantumnumber of k_d and k_e.
 * @param k_d first sp index that forms the dp column index j together with k_e, k_z and S_de, with dp spin and momentum SK
 * @param k_e second sp index that forms the dp column index j together with k_d, k_z and S_de, with dp spin and momentum SK
 * @param k_z third sp index that forms the dp column index j together with k_d, k_e and S_de, with dp spin and momentum SK
 * @return the number on place DPM(B,i,j) with the right phase and forefactor.
 */
/*
   double DPM::operator()(int S,int K,int S_ab,int k_a,int k_b,int k_c,int S_de,int k_d,int k_e,int k_z) const {

//check if the momentum is correct:
if( (k_a + k_b + k_c)%(M/2) != K)
return 0.0;

if( (k_d + k_e + k_z)%(M/2) != K)
return 0.0;

//blockindex:
int B = S*(M/2) + K;

int *i = new int [2];
double *coef_i = new double [2];

int dim_i = get_inco(B,S_ab,k_a,k_b,k_c,i,coef_i);

if(dim_i == 0){

delete [] i;
delete [] coef_i;

return 0.0;

}

int *j = new int [2];
double *coef_j = new double [2];

int dim_j = get_inco(B,S_de,k_d,k_e,k_z,j,coef_j);

if(dim_j == 0){

delete [] i;
delete [] j;

delete [] coef_i;
delete [] coef_j;

return 0.0;

}

double ward = 0.0;

for(int I = 0;I < dim_i;++I)
for(int J = 0;J < dim_j;++J)
ward += coef_i[I] * coef_j[J] * (*this)(B,i[I],j[J]);

delete [] i;
delete [] j;

delete [] coef_i;
delete [] coef_j;

return ward;

}
 */
/** 
 * Static member function that gets the dp-indices and their coefficients of the (s and t )-p indices S_ab,a,b,c.
 * @param B block index of the state
 * @param S_ab intermediate spincoupling of a and b.
 * @param k_a first sp-momentum orbital
 * @param k_b second sp-momentum orbital
 * @param k_c third sp-momentum orbital
 * @param i pointer of dim 1 or 2 containing the indices occuring in the expansion of this particular dp state in the normal basis (a==b,c a < b < c).
 * @param coef pointer of dim 1 or 2 containing the coefficients occuring in the expansion.
 * @return the number of terms in the expansion (1 or 2), also the dim of pointers i and coef. When zero is returned this is not a valid element.
 */
/*
   int DPM::get_inco(int B,int S_ab,int k_a,int k_b,int k_c,int *i,double *coef) const{

//they cannot all be equal
if(k_a == k_b && k_b == k_c)
return 0;

if(B < M/2){//spin 1/2 block:

//if normal basis:
if(k_a == k_b){

if(S_ab == 1)//spin has to be zero for k_a == k_b
return 0;

i[0] = s2dp[B][0][k_a][k_b][k_c];
coef[0] = 1;

return 1;

}
else if (k_a < k_b && k_b < k_c){

i[0] = s2dp[B][S_ab][k_a][k_b][k_c];
coef[0] = 1;

return 1;

}
else{//anomal basis:

int min,max,phase;

//first order a and b for code saving reasons
if(k_a < k_b){

min = k_a;
max = k_b;

phase = 1;

}
else{

min = k_b;
max = k_a;

phase = 1 - 2*S_ab;

if(k_c > max){//we still have one simple dim = 1 term left: b < a < c

i[0] = s2dp[B][S_ab][k_b][k_a][k_c];
coef[0] = phase;

return 1;

}

}

//now we have four possibilities left:
//don't forget to multiply every result by phase to get the right k_a and k_b for min and max!
// 1) k_c < min < max
// 2) k_c == min < max
// 3) min < k_c < max
// 4) min < max == c
if(k_c < min){//k_c < min < max

//the S_ca == 0 part:
i[0] = s2dp[B][0][k_c][min][max];
coef[0] = phase * std::sqrt(2.0*S_ab + 1) * (1 - 2*S_ab) * _6j[0][S_ab];

//the S_ca == 1 part:
i[1] = s2dp[B][1][k_c][min][max];
coef[1] = phase * std::sqrt(2.0*S_ab + 1) * (1 - 2*S_ab) * std::sqrt(3.0) * _6j[1][S_ab];

return 2;

}
else if(k_c == min){//k_c == min < max: this will also be a 1 dim list, because S_ac can only be 0 if k_a == k_c.

   i[0] = s2dp[B][0][k_c][min][max];
   coef[0] = std::sqrt(2.0) * phase * std::sqrt(2.0*S_ab + 1) * (1 - 2*S_ab) * _6j[0][S_ab];

   return 1;

}
else if(k_c < max){//min < k_c < max

   //S_ac == 0 part:
   i[0] = s2dp[B][0][min][k_c][max];
   coef[0] = phase * std::sqrt(2.0*S_ab + 1.0) * (1 - 2*S_ab) * _6j[0][S_ab];

   //S_ac == 1 part:
   i[1] = s2dp[B][1][min][k_c][max];
   coef[1] = - phase * std::sqrt(2.0*S_ab + 1.0) * (1 - 2*S_ab) * std::sqrt(3.0) * _6j[1][S_ab];

   return 2;

}
else{// min < k_c == max: also a 1 dim list, S_bc can only be 0 if k_b == k_c

   i[0] = s2dp[B][0][max][k_c][min];
   coef[0] = phase * std::sqrt(2.0) * std::sqrt(2.0*S_ab + 1.0) *_6j[0][S_ab];

   return 1;

}

}

}
else{//spin 3/2 block, totally antisymmetrical in the spatial sp orbs.

   //only S_ab == 1 can couple to 3/2's.
   if(S_ab == 0)
      return 0;

   //if any of the sp orbs are equal, antisymmetry leads to zero:
   if(k_a == k_b || k_b == k_c || k_c == k_a)
      return 0;

   if(k_a < k_b){

      if(k_b < k_c){//k_a < k_b < k_c

         i[0] = s2dp[B][0][k_a][k_b][k_c];
         coef[0] = 1;

      }
      else if(k_c < k_a){//k_c < k_a < k_b

         i[0] = s2dp[B][0][k_c][k_a][k_b];
         coef[0] = 1;

      }
      else{//k_a < k_c < k_b

         i[0] = s2dp[B][0][k_a][k_c][k_b];
         coef[0] = -1;

      }

   }
   else{//k_b < k_a

      if(k_a < k_c){//k_b < k_a < k_c

         i[0] = s2dp[B][0][k_b][k_a][k_c];
         coef[0] = -1;

      }
      else if(k_c < k_b){//k_c < k_b < k_a

         i[0] = s2dp[B][0][k_c][k_b][k_a];
         coef[0] = -1;

      }
      else{//k_b < k_c < k_a

         i[0] = s2dp[B][0][k_b][k_c][k_a];
         coef[0] = 1;

      }

   }

   return 1;

}

}
*/
/**
 * The spincoupled, translationally invariant T1-like (generalized T1) map: maps a TPM object (tpm) on a DPM object (*this)
 * @param A term before the tp part of the map
 * @param B term before the np part of the map
 * @param C term before the sp part of the map
 * @param tpm input TPM
 */
/*
   void DPM::T(double A,double B,double C,const TPM &tpm) {

//make sp matrix out of tpm
SPM spm(C,tpm);

double ward = 2.0*B*tpm.trace();

int k_a,k_b,k_c,k_d,k_e,k_z;
int S_ab,S_de;

int sign_ab,sign_de;

double norm_ab,norm_de;

double hard;

//start with the S = 1/2 blocks, these are the most difficult:
for(int B = 0;B < M/2;++B){

for(int i = 0;i < gdim(B);++i){

S_ab = dp2s[B][i][0];

k_a = dp2s[B][i][1];
k_b = dp2s[B][i][2];
k_c = dp2s[B][i][3];

sign_ab = 1 - 2*S_ab;

norm_ab = 1.0;

if(k_a == k_b)
norm_ab /= std::sqrt(2.0);

for(int j = i;j < gdim(B);++j){

S_de = dp2s[B][j][0];

k_d = dp2s[B][j][1];
k_e = dp2s[B][j][2];
k_z = dp2s[B][j][3];

sign_de = 1 - 2*S_de;

norm_de = 1.0;

if(k_d == k_e)
norm_de /= std::sqrt(2.0);

hard = std::sqrt( (2*S_ab + 1.0) * (2*S_de + 1.0) ) * _6j[S_ab][S_de];

//init
(*this)(B,i,j) = 0.0;

//the np + sp part
if(i == j)
(*this)(B,i,j) = ward - spm[k_a] - spm[k_b] - spm[k_c];

//other parts are a bit more difficult.

//tp(1)
if(k_c == k_z)
if(S_ab == S_de)
(*this)(B,i,j) += A * tpm(S_ab,(k_a + k_b)%(M/2),k_a,k_b,k_d,k_e);

//tp(2)
if(k_b == k_z){

if(k_a == k_c)
(*this)(B,i,j) += std::sqrt(2.0) * A * norm_ab * sign_ab * sign_de * hard * tpm(S_de,(k_a + k_c)%(M/2),k_a,k_c,k_d,k_e);
else
(*this)(B,i,j) += A * norm_ab * sign_ab * sign_de * hard * tpm(S_de,(k_a + k_c)%(M/2),k_a,k_c,k_d,k_e);

}

//tp(3)
if(k_a == k_z){

   if(k_b == k_c)
      (*this)(B,i,j) += std::sqrt(2.0) * A * norm_ab * sign_de * hard * tpm(S_de,(k_b + k_c)%(M/2),k_b,k_c,k_d,k_e);
   else
      (*this)(B,i,j) += A * norm_ab * sign_de * hard * tpm(S_de,(k_b + k_c)%(M/2),k_b,k_c,k_d,k_e);

}

//tp(4)
if(k_c == k_e){

   if(k_d == k_z)
      (*this)(B,i,j) += std::sqrt(2.0) * A * norm_de * sign_ab * sign_de * hard * tpm(S_ab,(k_a + k_b)%(M/2),k_a,k_b,k_d,k_z);
   else
      (*this)(B,i,j) += A * norm_de * sign_ab * sign_de * hard * tpm(S_ab,(k_a + k_b)%(M/2),k_a,k_b,k_d,k_z);

}

//tp(5)
if(k_b == k_e){

   double hulp = 0.0;

   //sum over intermediate spin
   for(int Z = 0;Z < 2;++Z)
      hulp += (2*Z + 1.0) * _6j[Z][S_ab] * _6j[Z][S_de] * tpm(Z,(k_a + k_c)%(M/2),k_a,k_c,k_d,k_z);

   //correct for norms of the tpm
   if(k_a == k_c)
      hulp *= std::sqrt(2.0);

   if(k_d == k_z)
      hulp *= std::sqrt(2.0);

   (*this)(B,i,j) += A * norm_ab * norm_de * sign_ab * sign_de * std::sqrt( (2*S_ab + 1.0) * (2*S_de + 1.0) ) * hulp;

}

//tp(6)
if(k_a == k_e){

   double hulp = 0.0;

   //sum over intermediate spin
   for(int Z = 0;Z < 2;++Z)
      hulp += (2*Z + 1.0) * _6j[Z][S_ab] * _6j[Z][S_de] * tpm(Z,(k_b + k_c)%(M/2),k_b,k_c,k_d,k_z);

   if(k_b == k_c)
      hulp *= std::sqrt(2.0);

   if(k_d == k_z)
      hulp *= std::sqrt(2.0);

   (*this)(B,i,j) += A * sign_de * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * hulp;

}

//tp(7)
if(k_c == k_d){

   if(k_e == k_z)
      (*this)(B,i,j) += std::sqrt(2.0) * A * norm_de * sign_ab * hard * tpm(S_ab,(k_a + k_b)%(M/2),k_a,k_b,k_e,k_z);
   else
      (*this)(B,i,j) += A * norm_de * sign_ab * hard * tpm(S_ab,(k_a + k_b)%(M/2),k_a,k_b,k_e,k_z);

}

//tp(8)
if(k_b == k_d){

   double hulp = 0.0;

   //sum over intermediate spin
   for(int Z = 0;Z < 2;++Z)
      hulp += (2*Z + 1.0) * _6j[Z][S_ab] * _6j[Z][S_de] * tpm(Z,(k_a + k_c)%(M/2),k_a,k_c,k_e,k_z);

   if(k_a == k_c)
      hulp *= std::sqrt(2.0);

   if(k_e == k_z)
      hulp *= std::sqrt(2.0);

   (*this)(B,i,j) += A * sign_ab * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * hulp;

}

//tp(9)
if(k_a == k_d){

   double hulp = 0.0;

   //sum over intermediate spin
   for(int Z = 0;Z < 2;++Z)
      hulp += (2*Z + 1.0) * _6j[Z][S_ab] * _6j[Z][S_de] * tpm(Z,(k_b + k_c)%(M/2),k_b,k_c,k_e,k_z);

   if(k_b == k_c)
      hulp *= std::sqrt(2.0);

   if(k_e == k_z)
      hulp *= std::sqrt(2.0);

   (*this)(B,i,j) += A * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * hulp;

}

}
}

}

//then the S = 3/2 blocks, this should be easy, totally antisymmetrical 
for(int B = M/2;B < M;++B){

   for(int i = 0;i < gdim(B);++i){

      k_a = dp2s[B][i][1];
      k_b = dp2s[B][i][2];
      k_c = dp2s[B][i][3];

      for(int j = i;j < gdim(B);++j){

         k_d = dp2s[B][j][1];
         k_e = dp2s[B][j][2];
         k_z = dp2s[B][j][3];

         (*this)(B,i,j) = 0.0;

         //np + sp part:
         if(i == j)
            (*this)(B,i,j) = ward - spm[k_a] - spm[k_b] - spm[k_c];

         //tp(1)
         if(k_c == k_z)
            (*this)(B,i,j) += A * tpm(1,(k_a + k_b)%(M/2),k_a,k_b,k_d,k_e);

         //tp(2)
         if(k_b == k_z)
            (*this)(B,i,j) -= A * tpm(1,(k_a + k_c)%(M/2),k_a,k_c,k_d,k_e);

         //tp(4)
         if(k_c == k_e)
            (*this)(B,i,j) -= A * tpm(1,(k_a + k_b)%(M/2),k_a,k_b,k_d,k_z);

         //tp(5)
         if(k_b == k_e)
            (*this)(B,i,j) += A * tpm(1,(k_a + k_c)%(M/2),k_a,k_c,k_d,k_z);

         //tp(7)
         if(k_c == k_d)
            (*this)(B,i,j) += A * tpm(1,(k_a + k_b)%(M/2),k_a,k_b,k_e,k_z);

         //tp(8)
         if(k_b == k_d)
            (*this)(B,i,j) -= A * tpm(1,(k_a + k_c)%(M/2),k_a,k_c,k_e,k_z);

         //tp(9)
         if(k_a == k_d)
            (*this)(B,i,j) += A * tpm(1,(k_b + k_c)%(M/2),k_b,k_c,k_e,k_z);

      }
   }

}

this->symmetrize();

}
*/
/**
 * The T1-map: maps a TPM object (tpm) on a DPM object (*this). 
 * @param tpm input TPM
 */
/*
   void DPM::T(const TPM &tpm){

   double a = 1.0;
   double b = 1.0/(N*(N - 1.0));
   double c = 1.0/(N - 1.0);

   this->T(a,b,c,tpm);

   }
 */
/** 
 * The hat function maps a TPM object tpm to a DPM object (*this) so that bar(this) = tpm,
 * The inverse of the TPM::bar function. It is a T1-like map.
 * @param tpm input TPM
 */
/*
   void DPM::hat(const TPM &tpm){

   double a = 1.0/(M - 4.0);
   double b = 1.0/((M - 4.0)*(M - 3.0)*(M - 2.0));
   double c = 1.0/((M - 4.0)*(M - 3.0));

   this->T(a,b,c,tpm);

   }

   ostream &operator<<(ostream &output,const DPM &dpm_p){

   for(int B = 0;B < dpm_p.gnr();++B){

   output << B << "\t(" << dpm_p.gS(B) << "," << dpm_p.gK(B) << ")\t" << dpm_p.gdim(B) << "\t" << dpm_p.gdeg(B) << std::endl;
   output << std::endl;

   for(int i = 0;i < dpm_p.gdim(B);++i)
   for(int j = 0;j < dpm_p.gdim(B);++j){

   output << i << "\t" << j << "\t|\t" << 

   dpm_p.dp2s[B][i][0] << "\t" << dpm_p.dp2s[B][i][1] << "\t" << dpm_p.dp2s[B][i][2] << "\t" << dpm_p.dp2s[B][i][3] << 

   "\t" << dpm_p.dp2s[B][j][0] << "\t" << dpm_p.dp2s[B][j][1] << "\t" << dpm_p.dp2s[B][j][2] << "\t" << dpm_p.dp2s[B][j][3] << "\t" << dpm_p(B,i,j) << endl;

   }

   output << endl;

   }

   return output;

   }
 */
/**
 * @return the dp spin-index corresponding to the blockindex (0 for 1/2 and 1 for 3/2)
 */
/*
   int DPM::gS(int block) const{

   return block_char[block][0];

   }
 */
/**
 * @return the dp momentum-index corresponding to the blockindex
 */
/*
   int DPM::gK(int block) const{

   return block_char[block][1];

   }
 */
/**
 * Output to file, to be read by the spin_pd program.
 * @param filename output file
 */
/*
   void DPM::out_sp(const char *filename) const{

   ofstream output(filename);
   output.precision(15);

   for(int B = 0;B < gnr();++B){

   for(int i = 0;i < gdim(B);++i)
   for(int j = i;j < gdim(B);++j)
   output << block_char[B][0] << "\t" << dp2s[B][i][0] << "\t" << dp2s[B][i][1] << "\t" << dp2s[B][i][2] << "\t" << dp2s[B][i][3] << "\t"

   << dp2s[B][j][0] << "\t" << dp2s[B][j][1] << "\t" << dp2s[B][j][2] << "\t" << dp2s[B][j][3] << "\t" << (*this)(B,i,j) << endl;

   }

   }
 */
