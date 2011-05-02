#ifndef PHM_H
#define PHM_H

#include <iostream>

using std::ostream;

#include "BlockMatrix.h"
#include "TPM.h"
#include "PPHM.h"

/**
 * @author Brecht Verstichel
 * @date 23-03-2011\n\n
 * This class, PHM, is a class written for spinsymmetrical, translationally invariant particle-hole matrices with parity taken into account.
 * It inherits all the functions from its mother class BlockMatrix, some special member functions 
 * and some lists that give the relationship between the sp and the ph basis.
 */
class PHM : public BlockMatrix {

    /**
    * Output stream operator overloaded,
    * @param output The stream to which you are writing (e.g. cout)
    * @param phm_p the PHM you want to print
    */
   friend ostream &operator<<(ostream &output,const PHM &phm_p);

   public:
      
      //constructor
      PHM();

      //copy constructor
      PHM(const PHM &);

      //destructor
      virtual ~PHM();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      double operator()(int B,int k_a,int k_b,int k_c,int k_d) const;

      double operator()(int S,int K,int p,int k_a,int k_b,int k_c,int k_d) const;

      static int get_phase_order(int S,int &K,int p,int &a,int &b);

      int gN() const;

      int gM() const;

      int gL() const;

      void G(const TPM &);

      static void init(int,int);

      static void clear();

      static double norm(int,int,int);

      void bar(const PPHM &);

   private:

      //!static list that takes in a blockindex B, a ph index i and returns two sp indices
      static vector< vector<int> > *ph2s;

      //!static list that takes in a blockindex B and two sp indices a,b and returns a ph index i
      static int ***s2ph;

      //!list of 6j symbols needed.
      static double **_6j;

      //!static list that takes a blockindex B and returns the ph spin S, the ph momentum K and ph parity p.
      static int **block_char;

      //!static list that takes a the block characteristics S,K,p and returns a blockindex B
      static int ***char_block;

      //!number of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

      //!nr of sites
      static int L;

};

#endif
