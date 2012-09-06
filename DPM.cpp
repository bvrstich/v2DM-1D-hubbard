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

int **DPM::block_char;
int ***DPM::char_block;

/**
 * allocate and initialize the static lists and variables
 */
void DPM::init(){

   //allocate the lists
   dp2s = new vector< vector<int> > [Tools::gL() + 6];

   s2dp = new int **** [Tools::gL() + 6];

   for(int B = 0;B < Tools::gL() + 6;++B){

      s2dp[B] = new int *** [2];

      for(int S_ab = 0;S_ab < 2;++S_ab){

         s2dp[B][S_ab] = new int ** [Tools::gL()];

         for(int k_a = 0;k_a < Tools::gL();++k_a){

            s2dp[B][S_ab][k_a] = new int * [Tools::gL()];

            for(int k_b = 0;k_b < Tools::gL();++k_b)
               s2dp[B][S_ab][k_a][k_b] = new int [Tools::gL()];

         }
      }
   }

   block_char = new int * [Tools::gL() + 6];

   for(int B = 0;B < Tools::gL() + 6;++B)
      block_char[B] = new int [3];

   char_block = new int ** [2];

   for(int S = 0;S < 2;++S){

      char_block[S] = new int * [Tools::gL()/2 + 1];

      for(int K = 0;K <= Tools::gL()/2;++K)
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

   for(int k_a = 1;k_a <= Tools::gL()/2;++k_a){

      for(int k_c = 0;k_c < k_a;++k_c){

         if( (2*k_a + k_c)%Tools::gL() == 0 ){

            v[0] = 0;//S_ab
            v[1] = k_a;
            v[2] = k_a;
            v[3] = k_c;

            dp2s[block].push_back(v);

            s2dp[block][0][k_a][k_a][k_c] = dp;

            ++dp;

         }

      }

      for(int k_c = k_a + 1;k_c < Tools::gL();++k_c){

         if( (2*k_a + k_c)%Tools::gL() == 0 ){

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

   //now S_ab = 0, k_a != k_b != k_c (transitive) , but with k_c == 0: only positive parity contribution!
   for(int k_a = 1;k_a < Tools::gL()/2;++k_a){

      int k_b = Tools::gL() - k_a;

      v[0] = 0;//S_ab
      v[1] = k_a;
      v[2] = k_b;
      v[3] = 0;

      dp2s[block].push_back(v);

      s2dp[block][0][k_a][k_b][0] = dp;

      ++dp;

   }

   //then, for both S_ab = 0/1, k_a < k_b < k_c: no zero momentum
   for(int S_ab = 0;S_ab < 2;++S_ab){

      for(int k_a = 1;k_a < Tools::gL();++k_a)
         for(int k_b = k_a + 1;k_b < Tools::gL();++k_b)
            for(int k_c = k_b + 1;k_c < Tools::gL();++k_c){

               if(k_a + k_b + k_c == Tools::gL()){//not modulo! special case

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

   block_char[block + Tools::gL()/2 + 3][0] = 1;//S = 3/2
   block_char[block + Tools::gL()/2 + 3][1] = 0;//K = 0
   block_char[block + Tools::gL()/2 + 3][2] = 0;//p = positive,

   char_block[1][0][0] = block + Tools::gL()/2 + 3;

   //only S_ab == 1!
   for(int k_a = 1;k_a < Tools::gL();++k_a)
      for(int k_b = k_a + 1;k_b < Tools::gL();++k_b)
         for(int k_c = k_b + 1;k_c < Tools::gL();++k_c){

            if(k_a + k_b + k_c == Tools::gL()){//not modulo! special case

               v[0] = 1;
               v[1] = k_a;
               v[2] = k_b;
               v[3] = k_c;

               dp2s[block + Tools::gL()/2 + 3].push_back(v);

               s2dp[block + Tools::gL()/2 + 3][1][k_a][k_b][k_c] = dp;

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
   for(int k_a = 1;k_a < Tools::gL()/2;++k_a){

      for(int k_c = 0;k_c < k_a;++k_c){

         if( (2*k_a + k_c)%Tools::gL() == 0 ){

            v[0] = 0;
            v[1] = k_a;
            v[2] = k_a;
            v[3] = k_c;

            dp2s[block].push_back(v);

            s2dp[block][0][k_a][k_a][k_c] = dp;

            ++dp;

         }

      }

      for(int k_c = k_a + 1;k_c < Tools::gL();++k_c){

         if( (2*k_a + k_c)%Tools::gL() == 0 ){

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
   for(int k_a = 1;k_a < Tools::gL();++k_a)
      for(int k_b = k_a + 1;k_b < Tools::gL();++k_b)
         for(int k_c = k_b + 1;k_c < Tools::gL();++k_c){

            if(k_a + k_b + k_c == Tools::gL()){//not modulo! special case

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
   for(int k_a = 1;k_a < Tools::gL()/2;++k_a){

      int k_b = Tools::gL() - k_a;

      v[0] = 1;
      v[1] = k_a;
      v[2] = k_b;
      v[3] = 0;

      dp2s[block].push_back(v);

      s2dp[block][1][k_a][k_b][0] = dp;

      ++dp;

   }

   //then the rest with S_ab == 1: k_a < k_b < k_c without 0
   for(int k_a = 1;k_a < Tools::gL();++k_a)
      for(int k_b = k_a + 1;k_b < Tools::gL();++k_b)
         for(int k_c = k_b + 1;k_c < Tools::gL();++k_c){

            if(k_a + k_b + k_c == Tools::gL()){//not modulo! special case

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

   block_char[block + Tools::gL()/2 + 3][0] = 1;//S = 3/2
   block_char[block + Tools::gL()/2 + 3][1] = 0;//K = 0
   block_char[block + Tools::gL()/2 + 3][2] = 1;//p = negative,

   char_block[1][0][1] = block + Tools::gL()/2 + 3;

   //first with a k=0: 
   for(int k_b = 1;k_b < Tools::gL()/2;++k_b){

      int k_c = Tools::gL() - k_b;

      v[0] = 1;
      v[1] = 0;
      v[2] = k_b;
      v[3] = k_c;

      dp2s[block + Tools::gL()/2 + 3].push_back(v);

      s2dp[block + Tools::gL()/2 + 3][1][0][k_b][k_c] = dp;

      ++dp;

   }

   //then the rest k_a < k_b < k_c
   for(int k_a = 1;k_a < Tools::gL();++k_a)
      for(int k_b = k_a + 1;k_b < Tools::gL();++k_b)
         for(int k_c = k_b + 1;k_c < Tools::gL();++k_c){

            if(k_a + k_b + k_c == Tools::gL()){//not modulo! special case

               v[0] = 1;
               v[1] = k_a;
               v[2] = k_b;
               v[3] = k_c;

               dp2s[block + Tools::gL()/2 + 3].push_back(v);

               s2dp[block + Tools::gL()/2 + 3][1][k_a][k_b][k_c] = dp;

               ++dp;

            }

         }

   ++block;

   //Now 0 < K < Tools::gL()/2: degeneracy between parity + and -

   for(int K = 1;K < Tools::gL()/2;++K){

      //first S = 1/2:
      dp = 0;

      block_char[block][0] = 0;//S = 1/2
      block_char[block][1] = K;
      block_char[block][2] = 0;

      char_block[0][K][0] = block;
      char_block[0][K][1] = block;

      //first S_ab = 0 and k_a = k_b != k_c
      for(int k_a = 0;k_a < Tools::gL();++k_a){

         for(int k_c = 0;k_c < k_a;++k_c){

            if( (2*k_a + k_c)%Tools::gL() == K ){

               v[0] = 0;
               v[1] = k_a;
               v[2] = k_a;
               v[3] = k_c;

               dp2s[block].push_back(v);

               s2dp[block][0][k_a][k_a][k_c] = dp;

               ++dp;

            }

         }

         for(int k_c = k_a + 1;k_c < Tools::gL();++k_c){

            if( (2*k_a + k_c)%Tools::gL() == K ){

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

         for(int k_a = 0;k_a < Tools::gL();++k_a)
            for(int k_b = k_a + 1;k_b < Tools::gL();++k_b)
               for(int k_c = k_b + 1;k_c < Tools::gL();++k_c){

                  if( (k_a + k_b + k_c)%Tools::gL() == K ){

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

      block_char[block + Tools::gL()/2 + 3][0] = 1;//S = 3/2
      block_char[block + Tools::gL()/2 + 3][1] = K;
      block_char[block + Tools::gL()/2 + 3][2] = 0;

      char_block[1][K][0] = block + Tools::gL()/2 + 3;
      char_block[1][K][1] = block + Tools::gL()/2 + 3;

      for(int k_a = 0;k_a < Tools::gL();++k_a)
         for(int k_b = k_a + 1;k_b < Tools::gL();++k_b)
            for(int k_c = k_b + 1;k_c < Tools::gL();++k_c){

               if( (k_a + k_b + k_c)%Tools::gL() == K ){

                  v[0] = 1;
                  v[1] = k_a;
                  v[2] = k_b;
                  v[3] = k_c;

                  dp2s[block + Tools::gL()/2 + 3].push_back(v);

                  s2dp[block + Tools::gL()/2 + 3][1][k_a][k_b][k_c] = dp;

                  ++dp;

               }

            }

      ++block;

   }

   //now the last ones: K = Tools::gL()/2 positive parity
   dp = 0;

   block_char[block][0] = 0;//S = 1/2
   block_char[block][1] = Tools::gL()/2;
   block_char[block][2] = 0;

   char_block[0][Tools::gL()/2][0] = block;

   //again first S_ab = 0, k_a == k_b != k_c

   //for positive parity one extra state: |0 0 Tools::gL()/2>
   v[0] = 0;//S_ab
   v[1] = 0;//k_a
   v[2] = 0;//k_b
   v[3] = Tools::gL()/2;//k_c

   dp2s[block].push_back(v);

   s2dp[block][0][0][0][Tools::gL()/2] = dp;

   ++dp;

   for(int k_a = 1;k_a < Tools::gL()/2;++k_a){

      for(int k_c = 0;k_c < k_a;++k_c){

         if( (2*k_a + k_c)%Tools::gL() == Tools::gL()/2 ){

            v[0] = 0;//S_ab
            v[1] = k_a;//k_a
            v[2] = k_a;//k_b
            v[3] = k_c;//k_c

            dp2s[block].push_back(v);

            s2dp[block][0][k_a][k_a][k_c] = dp;

            ++dp;

         }

      }

      for(int k_c = k_a + 1;k_c < Tools::gL();++k_c){

         if( (2*k_a + k_c)%Tools::gL() == Tools::gL()/2 ){

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

   //then, for S_ab == 0 and positive parity, terms for which k_c = Tools::gL()/2
   for(int k_a = 1;k_a < Tools::gL()/2;++k_a){

      int k_b = Tools::gL() - k_a;

      v[0] = 0;
      v[1] = k_a;
      v[2] = k_b;
      v[3] = Tools::gL()/2;

      dp2s[block].push_back(v);

      s2dp[block][0][k_a][k_b][Tools::gL()/2] = dp;

      ++dp;

   }

   //for both S_ab = 0 and 1: k_a < k_b < k_c: this is quite involved.
   for(int S_ab = 0;S_ab < 2;++S_ab){

      //first with k_a = 0, k_b < k_c < Tools::gL()/2
      for(int k_b = 1;k_b < Tools::gL()/2;++k_b)
         for(int k_c = k_b + 1;k_c < Tools::gL()/2;++k_c){

            if( (k_b + k_c) == Tools::gL()/2 ){

               v[0] = S_ab;//S_ab
               v[1] = 0;//k_a
               v[2] = k_b;//k_b
               v[3] = k_c;//k_c

               dp2s[block].push_back(v);

               s2dp[block][S_ab][0][k_b][k_c] = dp;

               ++dp;

            }

         }

      //then k_a < k_b < k_c with no zero's
      for(int k_a = 1;k_a < Tools::gL()/2;++k_a){

         for(int k_b = k_a + 1;k_b < Tools::gL()/2;++k_b)
            for(int k_c = k_b + 1;k_c < Tools::gL()/2;++k_c){

               if( (k_a + k_b + k_c) == Tools::gL()/2){

                  v[0] = S_ab;//S_ab
                  v[1] = k_a;//k_a
                  v[2] = k_b;//k_b
                  v[3] = k_c;//k_c

                  dp2s[block].push_back(v);

                  s2dp[block][S_ab][k_a][k_b][k_c] = dp;

                  ++dp;

               }

            }

         for(int k_b = Tools::gL()/2 + 1;k_b < Tools::gL();++k_b)
            for(int k_c = k_b + 1;k_c < Tools::gL();++k_c){

               if( (k_a + k_b + k_c)%Tools::gL() == Tools::gL()/2 ){

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

   //then K = Tools::gL()/2, S = 3/2 positive parity:
   dp = 0;

   block_char[block + Tools::gL()/2 + 3][0] = 1;//S = 3/2
   block_char[block + Tools::gL()/2 + 3][1] = Tools::gL()/2;
   block_char[block + Tools::gL()/2 + 3][2] = 0;

   char_block[1][Tools::gL()/2][0] = block + Tools::gL()/2 + 3;

   //first with k_a = 0, k_b < k_c < Tools::gL()/2
   for(int k_b = 1;k_b < Tools::gL()/2;++k_b)
      for(int k_c = k_b + 1;k_c < Tools::gL()/2;++k_c){

         if( (k_b + k_c) == Tools::gL()/2 ){

            v[0] = 1;//S_ab
            v[1] = 0;//k_a
            v[2] = k_b;//k_b
            v[3] = k_c;//k_c

            dp2s[block + Tools::gL()/2 + 3].push_back(v);

            s2dp[block + Tools::gL()/2 + 3][1][0][k_b][k_c] = dp;

            ++dp;

         }

      }

   //then k_a < k_b < k_c
   for(int k_a = 1;k_a < Tools::gL()/2;++k_a){

      for(int k_b = k_a + 1;k_b < Tools::gL()/2;++k_b)
         for(int k_c = k_b + 1;k_c < Tools::gL()/2;++k_c){

            if( (k_a + k_b + k_c) == Tools::gL()/2){

               v[0] = 1;//S_ab
               v[1] = k_a;//k_a
               v[2] = k_b;//k_b
               v[3] = k_c;//k_c

               dp2s[block + Tools::gL()/2 + 3].push_back(v);

               s2dp[block + Tools::gL()/2 + 3][1][k_a][k_b][k_c] = dp;

               ++dp;

            }

         }

      //only negative parity has k's equal to Tools::gL()/2
      for(int k_b = Tools::gL()/2 + 1;k_b < Tools::gL();++k_b)
         for(int k_c = k_b + 1;k_c < Tools::gL();++k_c){

            if( (k_a + k_b + k_c)%Tools::gL() == Tools::gL()/2 ){

               v[0] = 1;//S_ab
               v[1] = k_a;//k_a
               v[2] = k_b;//k_b
               v[3] = k_c;//k_c

               dp2s[block + Tools::gL()/2 + 3].push_back(v);

               s2dp[block + Tools::gL()/2 + 3][1][k_a][k_b][k_c] = dp;

               ++dp;

            }

         }

   }

   block++;

   //K = Tools::gL()/2, negative parity, start with S = 1/2
   dp = 0;

   block_char[block][0] = 0;//S = 1/2
   block_char[block][1] = Tools::gL()/2;
   block_char[block][2] = 1;//negative parity

   char_block[0][Tools::gL()/2][1] = block;

   //again first S_ab = 0, k_a == k_b != k_c
   for(int k_a = 1;k_a < Tools::gL()/2;++k_a){

      for(int k_c = 0;k_c < k_a;++k_c){

         if( (2*k_a + k_c)%Tools::gL() == Tools::gL()/2 ){

            v[0] = 0;//S_ab
            v[1] = k_a;//k_a
            v[2] = k_a;//k_b
            v[3] = k_c;//k_c

            dp2s[block].push_back(v);

            s2dp[block][0][k_a][k_a][k_c] = dp;

            ++dp;

         }

      }

      for(int k_c = k_a + 1;k_c < Tools::gL();++k_c){

         if( (2*k_a + k_c)%Tools::gL() == Tools::gL()/2 ){

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

      //first with k_a = 0, k_b < k_c < Tools::gL()/2
      for(int k_b = 1;k_b < Tools::gL()/2;++k_b)
         for(int k_c = k_b + 1;k_c < Tools::gL()/2;++k_c){

            if( (k_b + k_c) == Tools::gL()/2 ){

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
      for(int k_a = 1;k_a < Tools::gL()/2;++k_a){

         for(int k_b = k_a + 1;k_b < Tools::gL()/2;++k_b)
            for(int k_c = k_b + 1;k_c < Tools::gL()/2;++k_c){

               if( (k_a + k_b + k_c) == Tools::gL()/2){

                  v[0] = S_ab;//S_ab
                  v[1] = k_a;//k_a
                  v[2] = k_b;//k_b
                  v[3] = k_c;//k_c

                  dp2s[block].push_back(v);

                  s2dp[block][S_ab][k_a][k_b][k_c] = dp;

                  ++dp;

               }

            }

         //only S_ab == 1 has k's equal to Tools::gL()/2 for negative parity
         if(S_ab == 1){

            int k_b = Tools::gL() - k_a;

            v[0] = 1;
            v[1] = k_a;
            v[2] = k_b;
            v[3] = Tools::gL()/2;

            dp2s[block].push_back(v);

            s2dp[block][1][k_a][k_b][Tools::gL()/2] = dp;

            ++dp;

         }

         for(int k_b = Tools::gL()/2 + 1;k_b < Tools::gL();++k_b)
            for(int k_c = k_b + 1;k_c < Tools::gL();++k_c){

               if( (k_a + k_b + k_c)%Tools::gL() == Tools::gL()/2 ){

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

   //at last we have the last block, K = Tools::gL()/2, S = 3/2, negative parity
   dp = 0;

   block_char[block + Tools::gL()/2 + 3][0] = 1;//S = 3/2
   block_char[block + Tools::gL()/2 + 3][1] = Tools::gL()/2;
   block_char[block + Tools::gL()/2 + 3][2] = 1;

   char_block[1][Tools::gL()/2][1] = block + Tools::gL()/2 + 3;

   //first with k_a = 0, k_b < k_c < Tools::gL()/2
   for(int k_b = 1;k_b < Tools::gL()/2;++k_b)
      for(int k_c = k_b + 1;k_c < Tools::gL()/2;++k_c){

         if( (k_b + k_c) == Tools::gL()/2 ){

            v[0] = 1;//S_ab
            v[1] = 0;//k_a
            v[2] = k_b;//k_b
            v[3] = k_c;//k_c

            dp2s[block + Tools::gL()/2 + 3].push_back(v);

            s2dp[block + Tools::gL()/2 + 3][1][0][k_b][k_c] = dp;

            ++dp;

         }

      }

   //then k_a < k_b < k_c
   for(int k_a = 1;k_a < Tools::gL()/2;++k_a){

      for(int k_b = k_a + 1;k_b < Tools::gL()/2;++k_b)
         for(int k_c = k_b + 1;k_c < Tools::gL()/2;++k_c){

            if( (k_a + k_b + k_c) == Tools::gL()/2){

               v[0] = 1;//S_ab
               v[1] = k_a;//k_a
               v[2] = k_b;//k_b
               v[3] = k_c;//k_c

               dp2s[block + Tools::gL()/2 + 3].push_back(v);

               s2dp[block + Tools::gL()/2 + 3][1][k_a][k_b][k_c] = dp;

               ++dp;

            }

         }

      //only negative parity has k's equal to Tools::gL()/2
      for(int k_b = Tools::gL()/2;k_b < Tools::gL();++k_b)
         for(int k_c = k_b + 1;k_c < Tools::gL();++k_c){

            if( (k_a + k_b + k_c)%Tools::gL() == Tools::gL()/2 ){

               v[0] = 1;//S_ab
               v[1] = k_a;//k_a
               v[2] = k_b;//k_b
               v[3] = k_c;//k_c

               dp2s[block + Tools::gL()/2 + 3].push_back(v);

               s2dp[block + Tools::gL()/2 + 3][1][k_a][k_b][k_c] = dp;

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

   for(int B = 0;B < Tools::gL() + 6;++B){

      for(int S_ab = 0;S_ab < 2;++S_ab){

         for(int k_a = 0;k_a < Tools::gL();++k_a){

            for(int k_b = 0;k_b < Tools::gL();++k_b)
               delete [] s2dp[B][S_ab][k_a][k_b];

            delete [] s2dp[B][S_ab][k_a];

         }

         delete [] s2dp[B][S_ab];

      }

      delete [] s2dp[B];

   }

   delete [] s2dp;

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
 * standard constructor: constructs BlockMatrix object with 
 */
DPM::DPM() : BlockMatrix(Tools::gL() + 6) {

   //first K = 0: 4 blocks

   //positive parity
   this->setMatrixDim(0,dp2s[0].size(),2);//S = 1/2
   this->setMatrixDim(Tools::gL()/2 + 3,dp2s[Tools::gL()/2 + 3].size(),4);//S = 3/2

   //negative parity
   this->setMatrixDim(1,dp2s[1].size(),2);//S = 1/2
   this->setMatrixDim(Tools::gL()/2 + 4,dp2s[Tools::gL()/2 + 4].size(),4);//S = 3/2

   //then for 0 < K < Tools::gL()/2
   for(int K = 1;K < Tools::gL()/2;++K){

      this->setMatrixDim(K + 1,dp2s[K + 1].size(),4);//S = 1/2
      this->setMatrixDim(Tools::gL()/2 + K + 4,dp2s[Tools::gL()/2 + K + 4].size(),8);//S = 3/2

   }

   //and last for K = Tools::gL()/2: parity positive
   this->setMatrixDim(Tools::gL()/2 + 1,dp2s[Tools::gL()/2 + 1].size(),2);//S = 1/2
   this->setMatrixDim(Tools::gL() + 4,dp2s[Tools::gL() + 4].size(),4);//S = 3/2

   //negative parity
   this->setMatrixDim(Tools::gL()/2 + 2,dp2s[Tools::gL()/2 + 2].size(),2);//S = 1/2
   this->setMatrixDim(Tools::gL() + 5,dp2s[Tools::gL() + 5].size(),4);//S = 3/2

}

/**
 * copy constructor: constructs BlockMatrix object
 * @param dpm_c DPM to be copied into (*this)
 */
DPM::DPM(const DPM &dpm_c) : BlockMatrix(dpm_c) { }

/**
 * destructor 
 */
DPM::~DPM(){ }

/**
 * @param block the blockindex
 * @return the dp spin-index corresponding to the blockindex (0 for 1/2 and 1 for 3/2)
 */
int DPM::gS(int block) const{

   return block_char[block][0];

}

/**
 * @param block the blockindex
 * @return the dp momentum-index corresponding to the blockindex
 */
int DPM::gK(int block) const{

   return block_char[block][1];

}

/**
 * @param block the blockindex
 * @return the dp parity corresponding to the blockindex
 */
int DPM::gp(int block) const{

   return block_char[block][2];

}


ostream &operator<<(ostream &output,const DPM &dpm_p){

   for(int B = 0;B < dpm_p.gnr();++B){

      output << B << "\t(" << dpm_p.gS(B) << "," << dpm_p.gK(B) << "," << 1 - 2*dpm_p.gp(B) << ")\t" << dpm_p.gdim(B) << "\t" << dpm_p.gdeg(B) << std::endl;
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

/**
 * static test function that prints the basis constructed
 */
void DPM::print_basis(){

   for(int B = 0;B < Tools::gL() + 6;++B){

      cout << endl;
      cout << B << "\t(" << block_char[B][0] << "," << block_char[B][1] << "," << 1 - 2*block_char[B][2] << ")\t" << endl;
      cout << endl;

      for(unsigned int i = 0;i < dp2s[B].size();++i)
         cout << i << "\t|\t" << dp2s[B][i][0] << "\t" << dp2s[B][i][1] << "\t" << dp2s[B][i][2] << "\t" << dp2s[B][i][3] << endl;

   }

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
double DPM::operator()(int B,int S_ab,int k_a,int k_b,int k_c,int S_de,int k_d,int k_e,int k_z) const {

   return (*this)(block_char[B][0],block_char[B][1],block_char[B][2],S_ab,k_a,k_b,k_c,S_de,k_d,k_e,k_z);

}

/**
 * access the elements of the matrix in sp mode, special symmetry and antisymmetry relations are automatically accounted for:\n\n
 * DPM(S,K,S_ab,k_a,k_b,k_c,S_de,k_d,k_e,k_z) = sum_S_ac (some terms dependent on spin) DPM(S,K,S_ac,k_a,k_c,k_b,S_de,k_d,k_e,k_z) etc...
 * @param S dp-spin block quantumnumber: S = 0 means S == 1/2 and S = 1 means S == 3/2 (for simplicity)
 * @param K dp-momentum block quantumnumber
 * @param p dp-parity block quantumnumber
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
double DPM::operator()(int S,int K,int p,int S_ab,int k_a,int k_b,int k_c,int S_de,int k_d,int k_e,int k_z) const {

   //check if the momentum is correct:
   if( (k_a + k_b + k_c)%Tools::gL() != K)
      return 0.0;

   if( (k_d + k_e + k_z)%Tools::gL() != K)
      return 0.0;

   //they cannot all be equal
   if(k_a == k_b && k_b == k_c)
      return 0.0;

   if(k_d == k_e && k_e == k_z)
      return 0.0;

   int K_copy = K;

   int pphase_i = get_phase_order(S,K_copy,p,S_ab,k_a,k_b,k_c);

   if(pphase_i == 0)
      return 0.0;

   int *i = new int [2];
   double *coef_i = new double [2];

   int dim_i = get_inco(S,K_copy,p,S_ab,k_a,k_b,k_c,i,coef_i);

   if(dim_i == 0){

      delete [] i;
      delete [] coef_i;

      return 0.0;

   }

   int pphase_j = get_phase_order(S,K,p,S_de,k_d,k_e,k_z);

   if(pphase_j == 0){

      delete [] i;
      delete [] coef_i;

      return 0.0;

   }

   int *j = new int [2];
   double *coef_j = new double [2];

   int dim_j = get_inco(S,K,p,S_de,k_d,k_e,k_z,j,coef_j);

   if(dim_j == 0){

      delete [] i;
      delete [] j;

      delete [] coef_i;
      delete [] coef_j;

      return 0.0;

   }

   int B = char_block[S][K][p];

   double ward = 0.0;

   for(int I = 0;I < dim_i;++I)
      for(int J = 0;J < dim_j;++J)
         ward += pphase_i * pphase_j * coef_i[I] * coef_j[J] * (*this)(B,i[I],j[J]);

   delete [] i;
   delete [] j;

   delete [] coef_i;
   delete [] coef_j;

   return ward;

}

/**
 * function that transforms the sp momenta, en dp momentum, if they or not correct (i.e > Tools::gL()/2)
 * @param S dp spin
 * @param K dp momentum
 * @param p dp parity
 * @param k_a first momentum
 * @param k_b first momentum
 * @param k_c first momentum
 * @return the parity phase and the parity transformed momenta
 */
int DPM::get_phase_order(int S,int &K,int p,int &S_ab,int &k_a,int &k_b,int &k_c){

   if(S == 0){//S = 1/2

      if(K == 0){

         //check if there are two equal
         int flag = 0;
         int index;

         if(k_a == k_b){

            flag = 1;
            index = k_a;

         }
         else if(k_b == k_c){

            flag = 1;
            index = k_b;

         }
         else if (k_c == k_a){

            flag = 1;
            index = k_c;

         }

         if(flag == 1){//there are two equal, with index as the double momentum

            if(index > Tools::gL()/2){//transform!

               k_a = (Tools::gL()-k_a)%Tools::gL();
               k_b = (Tools::gL()-k_b)%Tools::gL();
               k_c = (Tools::gL()-k_c)%Tools::gL();

               return (1 - 2*p);

            }
            else if(index == Tools::gL()/2){//only positive parity!

               if(p == 0)
                  return 1;
               else
                  return 0;

            }
            else
               return 1;

         }
         else{//all different

            if(k_a == 0){

               if(p == 0){//positive parity

                  if(S_ab == 0){

                     if(k_c < k_b){//switcheroo!

                        k_a = k_c;
                        k_c = 0;

                        return -1;

                     }
                     else{

                        k_a = k_b;
                        k_b = k_c;
                        k_c = 0;

                        return -1;

                     }

                  }
                  else{//S_ab == 1

                     if(k_c < k_b){

                        k_a = k_c;
                        k_c = 0;

                        S_ab = 0;

                        return 1;

                     }
                     else{

                        k_a = k_b;
                        k_b = k_c;
                        k_c = 0;

                        S_ab = 0;

                        return 1;

                     }

                  }

               }
               else{//negative parity

                  if(S_ab == 0){

                     if(k_c < k_b){

                        k_a = k_c;
                        k_c = 0;

                        S_ab = 1;

                        return 1;

                     }
                     else{

                        k_a = k_b;
                        k_b = k_c;
                        k_c = 0;

                        S_ab = 1;

                        return -1;

                     }

                  }
                  else{//S_ab == 1

                     if(k_c < k_b){//switcheroo

                        k_a = k_c;
                        k_c = 0;

                        return 1;

                     }
                     else{

                        k_a = k_b;
                        k_b = k_c;
                        k_c = 0;

                        return -1;

                     }

                  }

               }

            }
            else if(k_b == 0){

               if(p == 0){//positive parity

                  if(S_ab == 0){

                     if(k_c < k_a){

                        k_b = k_a;
                        k_a = k_c;
                        k_c = 0;

                        return -1;

                     }
                     else{

                        k_b = k_c;
                        k_c = 0;

                        return -1;

                     }

                  }
                  else{//S_ab == 1

                     if(k_c < k_a){

                        k_b = k_a;
                        k_a = k_c;
                        k_c = 0;

                        S_ab = 0;

                        return 1;

                     }
                     else{

                        k_b = k_c;
                        k_c = 0;

                        S_ab = 0;

                        return 1;

                     }

                  }

               }
               else{//negative parity

                  if(S_ab == 0){

                     if(k_c < k_a){

                        k_b = k_a;
                        k_a = k_c;
                        k_c = 0;

                        S_ab = 1;

                        return 1;

                     }
                     else{

                        k_b = k_c;
                        k_c = 0;

                        S_ab = 1;

                        return -1;

                     }

                  }
                  else{//S_ab == 1

                     if(k_c < k_a){

                        k_b = k_a;
                        k_a = k_c;
                        k_c = 0;

                        return -1;

                     }
                     else{

                        k_b = k_c;
                        k_c = 0;

                        return 1;

                     }

                  }

               }

            }
            else if(k_c == 0){

               if(p == 0){//positive parity

                  if(S_ab == 0){

                     if(k_b < k_a){

                        k_c = k_b;
                        k_b = k_a;
                        k_a = k_c;
                        k_c = 0;

                        return 1;

                     }
                     else
                        return 1;

                  }
                  else//S_ab = 1 cannot have positive parity
                     return 0;

               }
               else{//negative parity

                  if(S_ab == 0)//cannot have negative parity
                     return 0;
                  else{

                     if(k_b < k_a){

                        k_c = k_b;
                        k_b = k_a;
                        k_a = k_c;
                        k_c = 0;

                        return -1;

                     }
                     else
                        return 1;

                  }

               }

            }
            else{//none of them are zero

               if(k_a + k_b + k_c == 2*Tools::gL()){//transform!

                  k_a = Tools::gL() - k_a;
                  k_b = Tools::gL() - k_b;
                  k_c = Tools::gL() - k_c;

                  return (1 - 2*p);

               }
               else
                  return 1;

            }

         }

      }
      else if(K == Tools::gL()/2){

         //check if there are two equal
         int flag = 0;
         int index;

         if(k_a == k_b){

            flag = 1;
            index = k_a;

         }
         else if(k_b == k_c){

            flag = 1;
            index = k_b;

         }
         else if (k_c == k_a){

            flag = 1;
            index = k_c;

         }

         if(flag == 1){//there are two equal, with index as the double momentum

            if(index > Tools::gL()/2){

               k_a = (Tools::gL() - k_a)%Tools::gL();
               k_b = (Tools::gL() - k_b)%Tools::gL();
               k_c = (Tools::gL() - k_c)%Tools::gL();

               return 1 - 2*p;

            }
            else if(index == 0){//only positive parity

               if(p == 0)
                  return 1;
               else
                  return 0;

            }
            else
               return 1;

         }
         else{//all different

            if(k_a == Tools::gL()/2){

               if(p == 0){//positive parity

                  if(S_ab == 0){

                     if(k_c < k_b){

                        k_a = k_c;
                        k_c = Tools::gL()/2;

                        return -1;

                     }
                     else{

                        k_a = k_b;
                        k_b = k_c;
                        k_c = Tools::gL()/2;

                        return -1;

                     }

                  }
                  else{//S_ab == 1

                     if(k_c < k_b){

                        k_a = k_c;
                        k_c = Tools::gL()/2;

                        S_ab = 0;

                        return 1;

                     }
                     else{

                        k_a = k_b;
                        k_b = k_c;
                        k_c = Tools::gL()/2;

                        S_ab = 0;

                        return 1;

                     }

                  }

               }
               else{//negative parity

                  if(S_ab == 0){

                     if(k_c < k_b){

                        k_a = k_c;
                        k_c = Tools::gL()/2;

                        S_ab = 1;

                        return 1;

                     }
                     else{

                        k_a = k_b;
                        k_b = k_c;
                        k_c = Tools::gL()/2;

                        S_ab = 1;

                        return -1;

                     }

                  }
                  else{//S_ab == 1

                     if(k_c < k_b){

                        k_a = k_c;
                        k_c = Tools::gL()/2;

                        return 1;

                     }
                     else{

                        k_a = k_b;
                        k_b = k_c;
                        k_c = Tools::gL()/2;

                        return -1;

                     }

                  }

               }

            }
            else if(k_b == Tools::gL()/2){

               if(p == 0){//positive parity

                  if(S_ab == 0){

                     if(k_c < k_a){

                        k_b = k_a;
                        k_a = k_c;
                        k_c = Tools::gL()/2;

                        return -1;

                     }
                     else{

                        k_b = k_c;
                        k_c = Tools::gL()/2;

                        return -1;

                     }

                  }
                  else{//S_ab == 1

                     if(k_c < k_a){

                        k_b = k_a;
                        k_a = k_c;
                        k_c = Tools::gL()/2;

                        S_ab = 0;

                        return 1;

                     }
                     else{

                        k_b = k_c;
                        k_c = Tools::gL()/2;

                        S_ab = 0;

                        return 1;

                     }

                  }

               }
               else{//negative parity

                  if(S_ab == 0){

                     if(k_c < k_a){

                        k_b = k_a;
                        k_a = k_c;
                        k_c = Tools::gL()/2;

                        S_ab = 1;

                        return 1;

                     }
                     else{

                        k_b = k_c;
                        k_c = Tools::gL()/2;

                        S_ab = 1;

                        return -1;

                     }

                  }
                  else{//S_ab == 1

                     if(k_c < k_a){

                        k_b = k_a;
                        k_a = k_c;
                        k_c = Tools::gL()/2;

                        return -1;

                     }
                     else{

                        k_b = k_c;
                        k_c = Tools::gL()/2;

                        return 1;

                     }

                  }

               }

            }
            else if(k_c == Tools::gL()/2){

               if(p == 0){//positive parity

                  if(S_ab == 0){

                     if(k_b < k_a){

                        k_c = k_a;
                        k_a = k_b;
                        k_b = k_c;
                        k_c = Tools::gL()/2;

                        return 1;

                     }
                     else
                        return 1;

                  }
                  else//S_ab == 1 does not exist
                     return 0;


               }
               else{//negative parity

                  if(S_ab == 0)//does not exist for negative parity
                     return 0;
                  else{//S_ab == 1

                     if(k_b < k_a){

                        k_c = k_a;
                        k_a = k_b;
                        k_b = k_c;
                        k_c = Tools::gL()/2;

                        return -1;

                     }
                     else
                        return 1;

                  }

               }

            }
            else{//none of the indices equals Tools::gL()/2

               if(k_a == 0 || k_b == 0 || k_c == 0){

                  if(k_a + k_b + k_c == Tools::gL()/2)
                     return 1;
                  else{//parity conversed states

                     k_a = (Tools::gL() - k_a)%Tools::gL();
                     k_b = (Tools::gL() - k_b)%Tools::gL();
                     k_c = (Tools::gL() - k_c)%Tools::gL();

                     return 1 - 2*p;

                  }

               }
               else{//none of the indices equals 0: 

                  if(k_a + k_b + k_c == Tools::gL()/2)
                     return 1;
                  else if(k_a + k_b + k_c == Tools::gL() + Tools::gL()/2){

                     //if there are two momenta < Tools::gL()/2, its a parity conversed state
                     if( (k_a < Tools::gL()/2 && k_b < Tools::gL()/2) || (k_a < Tools::gL()/2 && k_c < Tools::gL()/2) || (k_b < Tools::gL()/2 && k_c < Tools::gL()/2) ){

                        k_a = Tools::gL() - k_a;
                        k_b = Tools::gL() - k_b;
                        k_c = Tools::gL() - k_c;

                        return 1 - 2*p;

                     }
                     else
                        return 1;

                  }
                  else{//k_a + k_b + k_c = 2*Tools::gL() + Tools::gL()/2: conversed state

                     k_a = Tools::gL() - k_a;
                     k_b = Tools::gL() - k_b;
                     k_c = Tools::gL() - k_c;

                     return 1 - 2*p;

                  }

               }

            }

         }

      }
      else if (K > Tools::gL()/2){

         K = Tools::gL() - K;

         k_a = (Tools::gL() - k_a)%Tools::gL();
         k_b = (Tools::gL() - k_b)%Tools::gL();
         k_c = (Tools::gL() - k_c)%Tools::gL();

         return (1 - 2*p);

      }
      else
         return 1;

   }
   else{//S = 3/2: totally antisymmetric

      if(K == 0){

         if(k_a == 0 || k_b == 0 || k_c == 0){

            if(p == 0)//no positive parity possible due to antisymmetry
               return 0;
            else
               return 1;

         }
         else{//no zero's

            if(k_a + k_b + k_c == 2*Tools::gL()){//transform!

               k_a = Tools::gL() - k_a;
               k_b = Tools::gL() - k_b;
               k_c = Tools::gL() - k_c;

               return 1 - 2*p;

            }
            else
               return 1;

         }

      }
      else if(K == Tools::gL()/2){

         if(k_a == Tools::gL()/2 || k_b == Tools::gL()/2 || k_c == Tools::gL()/2){

            if(p == 0)//no positive parity possible due to antisymmetry
               return 0;
            else
               return 1;

         }
         else{//no Tools::gL()/2's

            if(k_a == 0 || k_b == 0 || k_c == 0){

               if(k_a + k_b + k_c == Tools::gL()/2)
                  return 1;
               else{//parity conversed states

                  k_a = (Tools::gL() - k_a)%Tools::gL();
                  k_b = (Tools::gL() - k_b)%Tools::gL();
                  k_c = (Tools::gL() - k_c)%Tools::gL();

                  return 1 - 2*p;

               }

            }
            else{//none of the indices equals 0: 

               if(k_a + k_b + k_c == Tools::gL()/2)
                  return 1;
               else if(k_a + k_b + k_c == Tools::gL() + Tools::gL()/2){

                  //if there are two momenta < Tools::gL()/2, its a parity conversed state
                  if( (k_a < Tools::gL()/2 && k_b < Tools::gL()/2) || (k_a < Tools::gL()/2 && k_c < Tools::gL()/2) || (k_b < Tools::gL()/2 && k_c < Tools::gL()/2) ){

                     k_a = Tools::gL() - k_a;
                     k_b = Tools::gL() - k_b;
                     k_c = Tools::gL() - k_c;

                     return 1 - 2*p;

                  }
                  else
                     return 1;

               }
               else{//k_a + k_b + k_c = 2*Tools::gL() + Tools::gL()/2: conversed state

                  k_a = Tools::gL() - k_a;
                  k_b = Tools::gL() - k_b;
                  k_c = Tools::gL() - k_c;

                  return 1 - 2*p;

               }

            }

         }

      }
      else if(K > Tools::gL()/2){

         K = Tools::gL() - K;

         k_a = (Tools::gL() - k_a)%Tools::gL();
         k_b = (Tools::gL() - k_b)%Tools::gL();
         k_c = (Tools::gL() - k_c)%Tools::gL();

         return 1 - 2*p;

      }
      else
         return 1;

   }

}


/** 
 * Static member function that gets the dp-indices and their coefficients of the (s and t )-p indices S_ab,a,b,c.
 * @param S dp-spin
 * @param K dp-momentum
 * @param p dp-parity
 * @param S_ab intermediate spincoupling of a and b.
 * @param k_a first sp-momentum orbital
 * @param k_b second sp-momentum orbital
 * @param k_c third sp-momentum orbital
 * @param i pointer of dim 1 or 2 containing the indices occuring in the expansion of this particular dp state in the normal basis (a==b,c a < b < c).
 * @param coef pointer of dim 1 or 2 containing the coefficients occuring in the expansion.
 * @return the number of terms in the expansion (1 or 2), also the dim of pointers i and coef. When zero is returned this is not a valid element.
 */
int DPM::get_inco(int S,int K,int p,int S_ab,int k_a,int k_b,int k_c,int *i,double *coef){

   int B = char_block[S][K][p];

   if(S == 0){//spin 1/2 block:

      //if normal basis:
      if(K == 0 && k_c == 0){

         i[0] = s2dp[char_block[S][K][p]][S_ab][k_a][k_b][k_c];
         coef[0] = 1;

         return 1;

      }

      if(K == Tools::gL()/2 && k_c == Tools::gL()/2){

         i[0] = s2dp[char_block[S][K][p]][S_ab][k_a][k_b][k_c];
         coef[0] = 1;

         return 1;

      }

      if(k_a == k_b){

         if(S_ab == 1)//spin has to be zero for k_a == k_b
            return 0;

         i[0] = s2dp[char_block[S][K][p]][0][k_a][k_a][k_c];
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
            coef[0] = phase * std::sqrt(2.0*S_ab + 1) * (1 - 2*S_ab) * Tools::g6j(0,0,0,S_ab);

            //the S_ca == 1 part:
            i[1] = s2dp[B][1][k_c][min][max];
            coef[1] = phase * std::sqrt(2.0*S_ab + 1) * (1 - 2*S_ab) * std::sqrt(3.0) * Tools::g6j(0,0,1,S_ab);

            return 2;

         }
         else if(k_c == min){//k_c == min < max: this will also be a 1 dim list, because S_ac can only be 0 if k_a == k_c.

            i[0] = s2dp[B][0][k_c][min][max];
            coef[0] = std::sqrt(2.0) * phase * std::sqrt(2.0*S_ab + 1) * (1 - 2*S_ab) * Tools::g6j(0,0,0,S_ab);

            return 1;

         }
         else if(k_c < max){//min < k_c < max

            //S_ac == 0 part:
            i[0] = s2dp[B][0][min][k_c][max];
            coef[0] = phase * std::sqrt(2.0*S_ab + 1.0) * (1 - 2*S_ab) * Tools::g6j(0,0,0,S_ab);

            //S_ac == 1 part:
            i[1] = s2dp[B][1][min][k_c][max];
            coef[1] = - phase * std::sqrt(2.0*S_ab + 1.0) * (1 - 2*S_ab) * std::sqrt(3.0) * Tools::g6j(0,0,1,S_ab);

            return 2;

         }
         else{// min < k_c == max: also a 1 dim list, S_bc can only be 0 if k_b == k_c

            i[0] = s2dp[B][0][max][k_c][min];
            coef[0] = phase * std::sqrt(2.0) * std::sqrt(2.0*S_ab + 1.0) *Tools::g6j(0,0,0,S_ab);

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

            i[0] = s2dp[B][1][k_a][k_b][k_c];
            coef[0] = 1;

         }
         else if(k_c < k_a){//k_c < k_a < k_b

            i[0] = s2dp[B][1][k_c][k_a][k_b];
            coef[0] = 1;

         }
         else{//k_a < k_c < k_b

            i[0] = s2dp[B][1][k_a][k_c][k_b];
            coef[0] = -1;

         }

      }
      else{//k_b < k_a

         if(k_a < k_c){//k_b < k_a < k_c

            i[0] = s2dp[B][1][k_b][k_a][k_c];
            coef[0] = -1;

         }
         else if(k_c < k_b){//k_c < k_b < k_a

            i[0] = s2dp[B][1][k_c][k_b][k_a];
            coef[0] = -1;

         }
         else{//k_b < k_c < k_a

            i[0] = s2dp[B][1][k_b][k_c][k_a];
            coef[0] = 1;

         }

      }

      return 1;

   }

}

/**
 * The spincoupled, translationally invariant and parity symmetric T1-like (generalized T1) map: maps a TPM object (tpm) on a DPM object (*this)
 * @param A term before the tp part of the map
 * @param B term before the np part of the map
 * @param C term before the sp part of the map
 * @param tpm input TPM
 */
void DPM::T(double A,double B,double C,const TPM &tpm) {

   //make sp matrix out of tpm
   SPM spm(C,tpm);

   double ward = 2.0*B*tpm.trace();

   int k_a,k_b,k_c,k_d,k_e,k_z;
   int S_ab,S_de;

   int k_d_,k_e_,k_z_;

   int K,p;

   int K_tp;

   int sign_ab,sign_de;

   int psign;

   double norm_ab,norm_de;

   double hard,tard,kard;

   //start with the S = 1/2 blocks, these are the most difficult:
   for(int B = 0;B < Tools::gL()/2 + 3;++B){

      K = block_char[B][1];
      p = block_char[B][2];

      psign = 1 - 2*p;

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

            k_d_ = (Tools::gL() - k_d)%Tools::gL();
            k_e_ = (Tools::gL() - k_e)%Tools::gL();
            k_z_ = (Tools::gL() - k_z)%Tools::gL();

            sign_de = 1 - 2*S_de;

            norm_de = 1.0;

            if(k_d == k_e)
               norm_de /= std::sqrt(2.0);

            hard = std::sqrt( (2*S_ab + 1.0) * (2*S_de + 1.0) ) * Tools::g6j(0,0,S_ab,S_de);

            //init
            (*this)(B,i,j) = 0.0;

            //the np + sp part
            if(i == j)
               (*this)(B,i,j) = ward - spm[k_a] - spm[k_b] - spm[k_c];

            kard = 0.0;

            if(K == 0 || K == Tools::gL()/2){//difference between parity + and - blocks: the exchange terms

               //tp(1)
               if(S_ab == S_de){

                  if(k_c == k_z_){

                     K_tp = (k_a + k_b)%Tools::gL();

                     tard = 0.0;

                     for(int pi = 0;pi < 2;++pi)
                        tard += tpm(S_ab,K_tp,pi,k_a,k_b,k_d_,k_e_);

                     kard += A * psign * tard / ( 2.0 * TPM::norm(K_tp,k_a,k_b) * TPM::norm(K_tp,k_d_,k_e_) );

                  }

               }

               //tp(2)
               if(k_b == k_z_){

                  K_tp = (k_a + k_c)%Tools::gL();

                  tard = 0.0;

                  for(int pi = 0;pi < 2;++pi)
                     tard += tpm(S_de,K_tp,pi,k_a,k_c,k_d_,k_e_);

                  if(k_a == k_c)
                     tard *= std::sqrt(2.0);

                  kard += A * psign * norm_ab * sign_ab * sign_de * hard * tard / (2.0 * TPM::norm(K_tp,k_a,k_c) * TPM::norm(K_tp,k_d_,k_e_) );

               }

               //tp(3)
               if(k_a == k_z_){

                  K_tp = (k_b + k_c)%Tools::gL();

                  tard = 0.0;

                  for(int pi = 0;pi < 2;++pi)
                     tard += tpm(S_de,K_tp,pi,k_b,k_c,k_d_,k_e_);

                  if(k_b == k_c)
                     tard *= std::sqrt(2.0);

                  kard += A * psign * norm_ab * sign_de * hard * tard / (2.0 * TPM::norm(K_tp,k_b,k_c) * TPM::norm(K_tp,k_d_,k_e_) ) ;

               }

               //tp(4)
               if(k_c == k_e_){

                  K_tp = (k_a + k_b)%Tools::gL();

                  tard = 0.0;

                  for(int pi = 0;pi < 2;++pi)
                     tard += tpm(S_ab,K_tp,pi,k_a,k_b,k_d_,k_z_);

                  if(k_d == k_z)
                     tard *= std::sqrt(2.0);

                  kard += A * psign * norm_de * sign_ab * sign_de * hard * tard / (2.0 * TPM::norm(K_tp,k_a,k_b) * TPM::norm(K_tp,k_d_,k_z_) ) ;

               }

               //tp(5)
               if(k_b == k_e_){

                  K_tp = (k_a + k_c)%Tools::gL();

                  tard = 0.0;

                  //sum over intermediate spin
                  for(int pi = 0;pi < 2;++pi)
                     for(int Z = 0;Z < 2;++Z)
                        tard += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,K_tp,pi,k_a,k_c,k_d_,k_z_);

                  //correct for norms of the tpm
                  if(k_a == k_c)
                     tard *= std::sqrt(2.0);

                  if(k_d == k_z)
                     tard *= std::sqrt(2.0);

                  kard += A * psign * norm_ab * norm_de * sign_ab * sign_de * std::sqrt( (2*S_ab + 1.0) * (2*S_de + 1.0) ) * tard 

                     / (2.0 * TPM::norm(K_tp,k_a,k_c) * TPM::norm(K_tp,k_d_,k_z_) ) ;

               }

               //tp(6)
               if(k_a == k_e_){

                  tard = 0.0;

                  K_tp = (k_b + k_c)%Tools::gL();

                  //sum over intermediate spin
                  for(int pi = 0;pi < 2;++pi)
                     for(int Z = 0;Z < 2;++Z)
                        tard += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,K_tp,pi,k_b,k_c,k_d_,k_z_);

                  if(k_b == k_c)
                     tard *= std::sqrt(2.0);

                  if(k_d == k_z)
                     tard *= std::sqrt(2.0);

                  kard += A * psign * sign_de * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * tard

                     / (2.0 * TPM::norm(K_tp,k_b,k_c) * TPM::norm(K_tp,k_d_,k_z_) ) ;

               }

               //tp(7)
               if(k_c == k_d_){

                  tard = 0.0;

                  K_tp = (k_a + k_b)%Tools::gL();

                  for(int pi = 0;pi < 2;++pi)
                     tard += tpm(S_ab,K_tp,pi,k_a,k_b,k_e_,k_z_);

                  if(k_e == k_z)
                     tard *= std::sqrt(2.0);

                  kard += A * psign * norm_de * sign_ab * hard * tard / (2.0 * TPM::norm(K_tp,k_a,k_b) * TPM::norm(K_tp,k_e,k_z) );

               }

               //tp(8)
               if(k_b == k_d_){

                  tard = 0.0;

                  K_tp = (k_a + k_c)%Tools::gL();

                  //sum over intermediate spin
                  for(int pi = 0;pi < 2;++pi)
                     for(int Z = 0;Z < 2;++Z)
                        tard += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,K_tp,pi,k_a,k_c,k_e_,k_z_);

                  if(k_a == k_c)
                     tard *= std::sqrt(2.0);

                  if(k_e == k_z)
                     tard *= std::sqrt(2.0);

                  kard += A * psign * sign_ab * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * tard

                     / (2.0 * TPM::norm(K_tp,k_a,k_c) * TPM::norm(K_tp,k_e_,k_z_) );

               }

               //tp(9)
               if(k_a == k_d_){

                  tard = 0.0;

                  K_tp = (k_b + k_c)%Tools::gL();

                  //sum over intermediate spin
                  for(int pi = 0;pi < 2;++pi)
                     for(int Z = 0;Z < 2;++Z)
                        tard += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,K_tp,pi,k_b,k_c,k_e_,k_z_);

                  if(k_b == k_c)
                     tard *= std::sqrt(2.0);

                  if(k_e == k_z)
                     tard *= std::sqrt(2.0);

                  kard += A * psign * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * tard

                     / (2.0 * TPM::norm(K_tp,k_b,k_c) * TPM::norm(K_tp,k_e_,k_z_));

               }

            }

            //tp(1)
            if(k_c == k_z)
               if(S_ab == S_de){

                  K_tp = (k_a + k_b)%Tools::gL();

                  tard = 0.0;

                  for(int pi = 0;pi < 2;++pi)
                     tard += tpm(S_ab,K_tp,pi,k_a,k_b,k_d,k_e);

                  kard += A * tard / ( 2.0 * TPM::norm(K_tp,k_a,k_b) * TPM::norm(K_tp,k_d,k_e) ) ;

               }

            //tp(2)
            if(k_b == k_z){

               K_tp = (k_a + k_c)%Tools::gL();

               tard = 0.0;

               for(int pi = 0;pi < 2;++pi)
                  tard += tpm(S_de,K_tp,pi,k_a,k_c,k_d,k_e);

               if(k_a == k_c)
                  tard *= std::sqrt(2.0);

               kard += A * norm_ab * sign_ab * sign_de * hard * tard / (2.0 * TPM::norm(K_tp,k_a,k_c) * TPM::norm(K_tp,k_d,k_e) );

            }

            //tp(3)
            if(k_a == k_z){

               K_tp = (k_b + k_c)%Tools::gL();

               tard = 0.0;

               for(int pi = 0;pi < 2;++pi)
                  tard += tpm(S_de,K_tp,pi,k_b,k_c,k_d,k_e);

               if(k_b == k_c)
                  tard *= std::sqrt(2.0);

               kard += A * norm_ab * sign_de * hard * tard / (2.0 * TPM::norm(K_tp,k_b,k_c) * TPM::norm(K_tp,k_d,k_e) ) ;

            }

            //tp(4)
            if(k_c == k_e){

               K_tp = (k_a + k_b)%Tools::gL();

               tard = 0.0;

               for(int pi = 0;pi < 2;++pi)
                  tard += tpm(S_ab,K_tp,pi,k_a,k_b,k_d,k_z);

               if(k_d == k_z)
                  tard *= std::sqrt(2.0);

               kard += A * norm_de * sign_ab * sign_de * hard * tard / (2.0 * TPM::norm(K_tp,k_a,k_b) * TPM::norm(K_tp,k_d,k_z) ) ;

            }

            //tp(5)
            if(k_b == k_e){

               K_tp = (k_a + k_c)%Tools::gL();

               tard = 0.0;

               //sum over intermediate spin
               for(int pi = 0;pi < 2;++pi)
                  for(int Z = 0;Z < 2;++Z)
                     tard += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,K_tp,pi,k_a,k_c,k_d,k_z);

               //correct for norms of the tpm
               if(k_a == k_c)
                  tard *= std::sqrt(2.0);

               if(k_d == k_z)
                  tard *= std::sqrt(2.0);

               kard += A * norm_ab * norm_de * sign_ab * sign_de * std::sqrt( (2*S_ab + 1.0) * (2*S_de + 1.0) ) * tard 

                  / (2.0 * TPM::norm(K_tp,k_a,k_c) * TPM::norm(K_tp,k_d,k_z) ) ;

            }

            //tp(6)
            if(k_a == k_e){

               tard = 0.0;

               K_tp = (k_b + k_c)%Tools::gL();

               //sum over intermediate spin
               for(int pi = 0;pi < 2;++pi)
                  for(int Z = 0;Z < 2;++Z)
                     tard += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,K_tp,pi,k_b,k_c,k_d,k_z);

               if(k_b == k_c)
                  tard *= std::sqrt(2.0);

               if(k_d == k_z)
                  tard *= std::sqrt(2.0);

               kard += A * sign_de * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * tard

                  / (2.0 * TPM::norm(K_tp,k_b,k_c) * TPM::norm(K_tp,k_d,k_z) ) ;

            }

            //tp(7)
            if(k_c == k_d){

               tard = 0.0;

               K_tp = (k_a + k_b)%Tools::gL();

               for(int pi = 0;pi < 2;++pi)
                  tard += tpm(S_ab,K_tp,pi,k_a,k_b,k_e,k_z);

               if(k_e == k_z)
                  tard *= std::sqrt(2.0);

               kard += A * norm_de * sign_ab * hard * tard / (2.0 * TPM::norm(K_tp,k_a,k_b) * TPM::norm(K_tp,k_e,k_z) );

            }

            //tp(8)
            if(k_b == k_d){

               tard = 0.0;

               K_tp = (k_a + k_c)%Tools::gL();

               //sum over intermediate spin
               for(int pi = 0;pi < 2;++pi)
                  for(int Z = 0;Z < 2;++Z)
                     tard += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,K_tp,pi,k_a,k_c,k_e,k_z);

               if(k_a == k_c)
                  tard *= std::sqrt(2.0);

               if(k_e == k_z)
                  tard *= std::sqrt(2.0);

               kard += A * sign_ab * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * tard

                  / (2.0 * TPM::norm(K_tp,k_a,k_c) * TPM::norm(K_tp,k_e,k_z) );

            }

            //tp(9)
            if(k_a == k_d){

               tard = 0.0;

               K_tp = (k_b + k_c)%Tools::gL();

               //sum over intermediate spin
               for(int pi = 0;pi < 2;++pi)
                  for(int Z = 0;Z < 2;++Z)
                     tard += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,K_tp,pi,k_b,k_c,k_e,k_z);

               if(k_b == k_c)
                  tard *= std::sqrt(2.0);

               if(k_e == k_z)
                  tard *= std::sqrt(2.0);

               kard += A * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * tard

                  / (2.0 * TPM::norm(K_tp,k_b,k_c) * TPM::norm(K_tp,k_e,k_z));

            }

            //finally the DPM norm
            kard *= DPM::norm(0,K,p,S_ab,k_a,k_b,k_c) * DPM::norm(0,K,p,S_de,k_d,k_e,k_z);

            (*this)(B,i,j) += kard;

         }
      }

   }

   //then the S = 3/2 blocks, this should be easy, totally antisymmetrical 
   for(int B = Tools::gL()/2 + 3;B < gnr();++B){

      K = block_char[B][1];
      p = block_char[B][2];

      psign = 1 - 2*p;

      for(int i = 0;i < gdim(B);++i){

         k_a = dp2s[B][i][1];
         k_b = dp2s[B][i][2];
         k_c = dp2s[B][i][3];

         for(int j = i;j < gdim(B);++j){

            k_d = dp2s[B][j][1];
            k_e = dp2s[B][j][2];
            k_z = dp2s[B][j][3];

            k_d_ = (Tools::gL() - k_d)%Tools::gL();
            k_e_ = (Tools::gL() - k_e)%Tools::gL();
            k_z_ = (Tools::gL() - k_z)%Tools::gL();

            (*this)(B,i,j) = 0.0;

            //np + sp part:
            if(i == j)
               (*this)(B,i,j) = ward - spm[k_a] - spm[k_b] - spm[k_c];

            kard = 0.0;

            if(K == 0 || K == Tools::gL()/2){//difference between parity + and - blocks

               //tp(1)
               if(k_c == k_z_){

                  K_tp = (k_a + k_b)%Tools::gL();

                  tard = 0.0;

                  for(int pi = 0;pi < 2;++pi)
                     tard += tpm(1,K_tp,pi,k_a,k_b,k_d_,k_e_);

                  kard += A * psign * tard / ( 2.0 * TPM::norm(K_tp,k_a,k_b) * TPM::norm(K_tp,k_d_,k_e_) );

               }

               //tp(2)
               if(k_b == k_z_){

                  K_tp = (k_a + k_c)%Tools::gL();

                  tard = 0.0;

                  for(int pi = 0;pi < 2;++pi)
                     tard += tpm(1,K_tp,pi,k_a,k_c,k_d_,k_e_);

                  kard -= A * psign * tard / (2.0 * TPM::norm(K_tp,k_a,k_c) * TPM::norm(K_tp,k_d_,k_e_) );

               }

               //tp(3):
               if(k_a == k_z_){

                  K_tp = (k_b + k_c)%Tools::gL();

                  tard = 0.0;

                  for(int pi = 0;pi < 2;++pi)
                     tard += tpm(1,K_tp,pi,k_b,k_c,k_d_,k_e_);

                  kard += A * psign * tard / (2.0 * TPM::norm(K_tp,k_b,k_c) * TPM::norm(K_tp,k_d_,k_e_) );

               }

               //tp(4)
               if(k_c == k_e_){

                  K_tp = (k_a + k_b)%Tools::gL();

                  tard = 0.0;

                  for(int pi = 0;pi < 2;++pi)
                     tard += tpm(1,K_tp,pi,k_a,k_b,k_d_,k_z_);

                  kard -= A * psign * tard / (2.0 * TPM::norm(K_tp,k_a,k_b) * TPM::norm(K_tp,k_d_,k_z_) ) ;

               }

               //tp(5)
               if(k_b == k_e_){

                  K_tp = (k_a + k_c)%Tools::gL();

                  tard = 0.0;

                  for(int pi = 0;pi < 2;++pi)
                     tard += tpm(1,K_tp,pi,k_a,k_c,k_d_,k_z_);

                  kard += A * psign * tard / (2.0 * TPM::norm(K_tp,k_a,k_c) * TPM::norm(K_tp,k_d_,k_z_) );

               }

               //tp(6)
               if(k_a == k_e_){

                  K_tp = (k_b + k_c)%Tools::gL();

                  tard = 0.0;

                  for(int pi = 0;pi < 2;++pi)
                     tard += tpm(1,K_tp,pi,k_b,k_c,k_d_,k_z_);

                  kard -= A * psign * tard / (2.0 * TPM::norm(K_tp,k_b,k_c) * TPM::norm(K_tp,k_d_,k_z_) );

               }

               //tp(7)
               if(k_c == k_d_){

                  K_tp = (k_a + k_b)%Tools::gL();

                  tard = 0.0;

                  for(int pi = 0;pi < 2;++pi)
                     tard += tpm(1,K_tp,pi,k_a,k_b,k_e_,k_z_);

                  kard += A * psign * tard / (2.0 * TPM::norm(K_tp,k_a,k_b) * TPM::norm(K_tp,k_e_,k_z_) );

               }

               //tp(8)
               if(k_b == k_d_){

                  K_tp = (k_a + k_c)%Tools::gL();

                  tard = 0.0;

                  for(int pi = 0;pi < 2;++pi)
                     tard += tpm(1,K_tp,pi,k_a,k_c,k_e_,k_z_);

                  kard -= A * psign * tard / (2.0 * TPM::norm(K_tp,k_a,k_c) * TPM::norm(K_tp,k_e_,k_z_) );

               }

               //tp(9)
               if(k_a == k_d_){

                  K_tp = (k_b + k_c)%Tools::gL();

                  tard = 0.0;

                  for(int pi = 0;pi < 2;++pi)
                     tard += tpm(1,K_tp,pi,k_b,k_c,k_e_,k_z_);

                  kard += A * psign * tard / (2.0 * TPM::norm(K_tp,k_b,k_c) * TPM::norm(K_tp,k_e_,k_z_) );

               }

            }

            //tp(1)
            if(k_c == k_z){

               K_tp = (k_a + k_b)%Tools::gL();

               tard = 0.0;

               for(int pi = 0;pi < 2;++pi)
                  tard += tpm(1,K_tp,pi,k_a,k_b,k_d,k_e);

               kard += A * tard / ( 2.0 * TPM::norm(K_tp,k_a,k_b) * TPM::norm(K_tp,k_d,k_e) ) ;

            }


            //tp(2)
            if(k_b == k_z){

               K_tp = (k_a + k_c)%Tools::gL();

               tard = 0.0;

               for(int pi = 0;pi < 2;++pi)
                  tard += tpm(1,K_tp,pi,k_a,k_c,k_d,k_e);

               kard -= A * tard / (2.0 * TPM::norm(K_tp,k_a,k_c) * TPM::norm(K_tp,k_d,k_e) );

            }

            //tp(4)
            if(k_c == k_e){

               K_tp = (k_a + k_b)%Tools::gL();

               tard = 0.0;

               for(int pi = 0;pi < 2;++pi)
                  tard += tpm(1,K_tp,pi,k_a,k_b,k_d,k_z);

               kard -= A * tard / (2.0 * TPM::norm(K_tp,k_a,k_b) * TPM::norm(K_tp,k_d,k_z) ) ;

            }

            //tp(5)
            if(k_b == k_e){

               K_tp = (k_a + k_c)%Tools::gL();

               tard = 0.0;

               for(int pi = 0;pi < 2;++pi)
                  tard += tpm(1,K_tp,pi,k_a,k_c,k_d,k_z);

               kard += A * tard / (2.0 * TPM::norm(K_tp,k_a,k_c) * TPM::norm(K_tp,k_d,k_z) );

            }

            //tp(7)
            if(k_c == k_d){

               K_tp = (k_a + k_b)%Tools::gL();

               tard = 0.0;

               for(int pi = 0;pi < 2;++pi)
                  tard += tpm(1,K_tp,pi,k_a,k_b,k_e,k_z);

               kard += A * tard / (2.0 * TPM::norm(K_tp,k_a,k_b) * TPM::norm(K_tp,k_e,k_z) );

            }

            //tp(8)
            if(k_b == k_d){

               K_tp = (k_a + k_c)%Tools::gL();

               tard = 0.0;

               for(int pi = 0;pi < 2;++pi)
                  tard += tpm(1,K_tp,pi,k_a,k_c,k_e,k_z);

               kard -= A * tard / (2.0 * TPM::norm(K_tp,k_a,k_c) * TPM::norm(K_tp,k_e,k_z) );

            }

            //tp(9)
            if(k_a == k_d){

               K_tp = (k_b + k_c)%Tools::gL();

               tard = 0.0;

               for(int pi = 0;pi < 2;++pi)
                  tard += tpm(1,K_tp,pi,k_b,k_c,k_e,k_z);

               kard += A * tard / (2.0 * TPM::norm(K_tp,k_b,k_c) * TPM::norm(K_tp,k_e,k_z) );

            }

            //finally the DPM norm
            kard *= DPM::norm(1,K,p,1,k_a,k_b,k_c) * DPM::norm(1,K,p,1,k_d,k_e,k_z);

            (*this)(B,i,j) += kard;

         }
      }

   }

   this->symmetrize();

}

/**
 * The T1-map: maps a TPM object (tpm) on a DPM object (*this). 
 * @param tpm input TPM
 */
void DPM::T(const TPM &tpm){

   double a = 1.0;
   double b = 1.0/(Tools::gN()*(Tools::gN() - 1.0));
   double c = 1.0/(Tools::gN() - 1.0);

   this->T(a,b,c,tpm);

}

/** 
 * The hat function maps a TPM object tpm to a DPM object (*this) so that bar(this) = tpm,
 * The inverse of the TPM::bar function. It is a T1-like map.
 * @param tpm input TPM
 */
void DPM::hat(const TPM &tpm){

   double a = 1.0/(2*Tools::gL() - 4.0);
   double b = 1.0/((2*Tools::gL() - 4.0)*(2*Tools::gL() - 3.0)*(2*Tools::gL() - 2.0));
   double c = 1.0/((2*Tools::gL() - 4.0)*(2*Tools::gL() - 3.0));

   this->T(a,b,c,tpm);

}

/**
 * @param S dp spin
 * @param K dp momentum
 * @param p dp parity
 * @param S_ab intermediate spin
 * @param k_a first momentum index
 * @param k_b second momentum index
 * @param k_c third momentum index
 * @return the norm of the basisstate
 */
double DPM::norm(int S,int K,int p,int S_ab,int k_a,int k_b,int k_c){

   if(S == 0){//S = 1/2

   if(K == 0){

      if(k_a == 0 || k_b == 0 || k_c == 0){

         if( (k_a == k_b) || (k_b == k_c) || (k_c == k_a))//(Tools::gL()/2 0 Tools::gL()/2) is mapped onto itsself
            return 0.5;

         if(k_c == 0)
            return 0.5;
         else{//recoupling needed!

            if(p == 0){//positive parity

               if(S_ab == 0)
                  return 1;
               else
                  return 1.0/std::sqrt(3.0);

            }
            else{//negative parity

               if(S_ab == 0)
                  return 1.0/std::sqrt(3.0);
               else
                  return 1;

            }

         }

      }
      else
         return 1.0/std::sqrt(2.0);

   }
   else if(K == Tools::gL()/2){

      if(k_a == Tools::gL()/2 || k_b == Tools::gL()/2 || k_c == Tools::gL()/2){

         if( (k_a == k_b) || (k_b == k_c) || (k_c == k_a))//(0 0 Tools::gL()/2) is mapped onto itsself
            return 0.5;

         if(k_c == Tools::gL()/2)//normal: mapped onto itsself
            return 0.5;
         else{//recoupling needed!

            if(p == 0){//positive parity

               if(S_ab == 0)
                  return 1;
               else
                  return 1.0/std::sqrt(3.0);

            }
            else{//negative parity

               if(S_ab == 0)
                  return 1.0/std::sqrt(3.0);
               else
                  return 1;

            }

         }

      }
      else
         return 1.0/std::sqrt(2.0);

   }
   else
      return 1.0/std::sqrt(2.0);

   }
   else{//S = 3/2

      if(K == 0){

         if(k_a == 0 || k_b == 0 || k_c == 0)
            return 0.5;
         else
            return 1.0/std::sqrt(2.0);

      }
      else if (K == Tools::gL()/2){

         if(k_a == Tools::gL()/2 || k_b == Tools::gL()/2 || k_c == Tools::gL()/2)
            return 0.5;
         else
            return 1.0/std::sqrt(2.0);

      }
      else
         return 1.0/std::sqrt(2.0);

   }

}
