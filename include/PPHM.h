#ifndef PPHM_H
#define PPHM_H

#include <iostream>

using std::ostream;

#include "BlockMatrix.h"
#include "TPM.h"

/**
 * @author Brecht Verstichel
 * @date 03-06-2010\n\n
 * This class, PPHM, is a class written for spinsymmetrical, translationally invaraint two-particle-one-hole matrices. It is written specially for the T_2 condition. 
 * It inherits all the functions from its mother class BlockMatrix, some special member functions and two lists that give
 * the relationship between the pph (two-particle one hole) and the sp basis. This matrix has M blocks, M/2 for S = 1/2 block with degeneracy 2
 * and M/2 for S = 3/2 block with degeneracy 4.
 */
class PPHM : public BlockMatrix {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param pphm_p the PPHM you want to print
    */
   friend ostream &operator<<(ostream &output,const PPHM &pphm_p);

   public:
      
      //constructor
      PPHM();

      //copy constructor
      PPHM(const PPHM &);

      //destructor
      virtual ~PPHM();

      void construct_lists();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      double operator()(int S,int K,int p,int S_ab,int k_a,int k_b,int k_c,int S_de,int k_d,int k_e,int k_z) const;

      double operator()(int B,int S_ab,int k_a,int k_b,int k_c,int S_de,int k_d,int k_e,int k_z) const;

      int get_inco(int B,int S,int S_ab,int k_a,int k_b,int k_c,int &i) const;

      static int get_phase_order(int S,int &K,int p,int &S_ab,int &a,int &b,int &c);

      //geef N terug
      int gN() const;

      //geef M terug
      int gM() const;

      int gL() const;

      //get spin of block
      int gS(int block) const;

      //get momentum of block
      int gK(int block) const;

      //get parity of block
      int gp(int block) const;

      //maak een PPHM van een TPM via de T2 conditie
      void T(const TPM &);

      void print_basis();

      static void init(int L_in,int N_in);
       
      static void clear();

      int total_dim();

      static double norm(int,int,int,int);

   private:

      //!static list that takes in a pph index i and a blockindex B and returns three momentum sp indices and an intermediate spin S_ab.
      static vector< vector<int> > *pph2s;

      //!static list that takes three momentum sp indices k_a, k_b and k_c, a blockindex B and an intermediate spinindex S_ab, and returns a pph index i:
      static int *****s2pph;

      //!list of 6j symbols needed.
      static double **_6j;

      //!list of block characteristics: takes in a block and return pph spin, momentum and parity
      static int **block_char;

      //!list that takes in pph spin, momentum and parity and returns the blockindex
      static int ***char_block;

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

      //!dimension of the lattice
      static int L;

};

#endif
