#ifndef SUP_H
#define SUP_H

#include <iostream> 
#include <fstream> 

using std::ostream;
using std::ofstream;
using std::ifstream;

#include "TPM.h"
#include "PHM.h"
#include "DPM.h"
#include "PPHM.h"

#ifdef PQG

#define __G_CON

#endif

#ifdef PQGT1

#define __G_CON
#define __T1_CON

#endif

#ifdef PQGT2

#define __G_CON
#define __T2_CON

#endif

#ifdef PQGT

#define __G_CON
#define __T1_CON
#define __T2_CON

#endif

class EIG;

/**
 * @author Brecht Verstichel
 * @date 09-03-2010\n\n
 * This class, SUP is a blockmatrix over the carrierspace's of active N-representability conditions. 
 * This class contains two TPM objects, and if compiled with the right option a PHM or DPM object, 
 * You have to remember that these matrices are independent of each other (by which I mean that TPM::Q(SUP_PQ::tpm (0))
 * is not neccesarily equal to SUP_PQ::tpm (1)) etc. .
 */
class SUP{
  
   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param SZ_p the SUP you want to print
    */
   friend ostream &operator<<(ostream &output,const SUP &SZ_p);

   public:

      //constructor
      SUP();

      //copy constructor
      SUP(const SUP &);

      //destructor
      ~SUP();

      //overload += operator
      SUP &operator+=(const SUP &);

      //overload -= operator
      SUP &operator-=(const SUP &);

      //overload equality operator
      SUP &operator=(const SUP &);

      //overload equality operator
      SUP &operator=(double );

      TPM &tpm(int i);

      const TPM &tpm(int i) const;

      double ddot(const SUP &) const;

      void invert();

      void dscal(double alpha);

      //positieve of negatieve vierkantswortel uit een supermatrix
      void sqrt(int option);

      void L_map(const SUP &,const SUP &);

      void daxpy(double alpha,const SUP &);

      SUP &mprod(const SUP &,const SUP &);

      void fill(const TPM &);

      void fill();

      void fill_Random();

      void out(ofstream &) const;

      void in(ifstream &);

#ifdef __G_CON
      PHM &phm();

      const PHM &phm() const;
#endif

#ifdef __T1_CON
      DPM &dpm();

      const DPM &dpm() const;
#endif

#ifdef __T2_CON
      PPHM &pphm();

      const PPHM &pphm() const;
#endif

      void sep_pm(SUP &,SUP &);

   private:

      //!double pointer of TPM's, will contain the P and Q block of the SUP in the first and second block.
      TPM **SZ_tp;

#ifdef __G_CON
      //!pointer to the particle hole matrix
      PHM *SZ_ph;
#endif

#ifdef __T1_CON
      //!pointer tot he three particles matrix DPM
      DPM *SZ_dp;
#endif

#ifdef __T2_CON
      //!pointer tot he three particles matrix DPM
      PPHM *SZ_pph;
#endif

};

#endif
