#ifndef TPM_H
#define TPM_H

#include <iostream>
#include <vector>
#include <fstream>

using std::ostream;
using std::ifstream;
using std::vector;

#include "BlockMatrix.h"

class SUP;

/**
 * @author Brecht Verstichel
 * @date 10-05-2010\n\n
 * This class TPM is a class written for two particle matrices with spinsymmetry, translational symmetry and parity included, it inherits all the functions from its mother 
 * BlockMatrix, some special member functions and two lists that give the relationship between the sp and the tp basis.
 */
class TPM : public BlockMatrix {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << tpm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << tpm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param tpm_p the TPM you want to print
    */
   friend ostream &operator<<(ostream &output,const TPM &tpm_p);

   public:
      
      //constructor
      TPM();

      //copy constructor
      TPM(const TPM &);

      //destructor
      virtual ~TPM();

      void constr_lists();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      //easy to access the numbers, in sp mode and blockindex
      double operator()(int B,int a,int b,int c,int d) const;

      //easy to access the numbers, in sp mode and with tp spin and momentum quantumnumber
      double operator()(int S,int K,int p,int a,int b,int c,int d) const;

      static int get_phase_order(int S,int &K,int p,int &a,int &b);

      //geef N terug
      int gN() const;

      //geef M terug
      int gM() const;

      int gL() const;

      void hubbard(double U);

      //Q afbeelding en zijn inverse
      void Q(int option,const TPM &);

      //Q like afbeelding Q(A,B,C,tpm_d)
      void Q(int option,double A,double B,double C,const TPM &);

      //overlapmatrix afbeelding en zijn inverse
      void S(int option,const TPM &);

      void unit();

      void proj_Tr();

      //de hessiaan afbeelding:
      void H(const TPM &b,const SUP &D);

      //los het stelsel op
      int solve(TPM &b,const SUP &D);

      void min_unit(double scale);

      void collaps(int option,const SUP &);

      //return the spin
      double spin() const;

      //output to file
      void out_sp(const char *) const;

      //input from file
      void in(ifstream &);

      static void init(int,int);

      static void clear();

      static void init_overlap();

      static double norm(int,int,int);

   private:

      //!static list that takes in a tp index i and a blockindex B, and returns two sp indices.
      static vector< vector<int> > *t2s;

      //!static that takes two sp momentum indices k_a,k_b and a blockindex B, and returns a tp index i:
      static int ***s2t;

      //!static list that takes a blockindex B and returns the tp spin S, the tp pseudo-momentum K~ and the tp parity p.
      static int **block_char;

      //!static list that returns the blockindex when given the S, K~ and p.
      static int ***char_block;

      //!overlapmatrix parameters
      static double Sa,Sc;

      //!list of 6j symbols needed.
      static double **_6j;

      //!nr of particles
      static int N;

      //!nr of sites
      static int L;

      //!dimension of sp hilbert space
      static int M;

};

#endif
