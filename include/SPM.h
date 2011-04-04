#ifndef SPM_H
#define SPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "Matrix.h"
#include "TPM.h"
#include "PHM.h"

/**
 * @author Brecht Verstichel
 * @date 21-03-2011\n\n
 * This class SPM was written for single particle matrices in a spinsymmetrical and translationally invariant system with parity included.
 */
class SPM {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * @param output The stream to which you are writing (e.g. cout)
    * @param spm_p de SPM you want to print
    */
   friend ostream &operator<<(ostream &output,const SPM &spm_p);

   public:
      
      //constructor
      SPM();

      //copy constructor
      SPM(const SPM &);

      //TPM constructor
      SPM(double ,const TPM &);

      //PHM constructor
      SPM(double ,const PHM &);

      //destructor
      virtual ~SPM();

      double operator[](int) const;

      double &operator[](int);

      int gN() const;

      int gM() const;

      int gL() const;

      const double *gspm() const;

      void bar(double,const TPM &);

      void bar(double,const PHM &);

      static void init(int,int);

      static void clear();

   private:

      //!dimension of single particle space
      static int M;

      //!nr of particles
      static int N;

      //!nr of sites
      static int L;

      //!dimension of the vector
      static int dim;

      //!list of degeneracies
      static int *deg;

      //!the object containing the numbers
      double *spm;

};

#endif
