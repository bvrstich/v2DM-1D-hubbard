#ifndef EIG_H
#define EIG_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "BlockVector.h"
#include "SUP.h"

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

/**
 * @author Brecht Verstichel
 * @date 06-05-2010\n\n
 * This class, EIG is a "block"-vector over the carrierspace's of the active condtions. It contains room
 * to store the eigenvalues and special member function that work with these eigenvalues.
 * This class should only be used when a SUP matrix has been diagonalized, some functions could give strange results when the EIG object is filled
 * with random numbers.\n\n
 */
class EIG{

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param eig_p the EIG you want to print
    */
   friend ostream &operator<<(ostream &output,const EIG &eig_p);

   public:

   //constructor met initialisatie op 
   EIG(SUP &);

   //copy constructor
   EIG(const EIG &);

   //destructor
   ~EIG();

   void diagonalize(SUP &);

   //overload equality operator
   EIG &operator=(const EIG &);

   BlockVector<TPM> &tpv(int);

   const BlockVector<TPM> &tpv(int) const;

#ifdef __G_CON

   BlockVector<PHM> &phv();

   const BlockVector<PHM> &phv() const;

#endif

#ifdef __T1_CON

   BlockVector<DPM> &dpv();

   const BlockVector<DPM> &dpv() const;

#endif

#ifdef __T2_CON

   BlockVector<PPHM> &pphv();

   const BlockVector<PPHM> &pphv() const;

#endif

   double min() const;

   double max() const;

   static void init(int,int);

   private:

   //!double pointer to a BlockVector<TPM> object, the eigenvalues of the P and Q part of a SUP matrix will be stored here.
   BlockVector<TPM> **v_tp;

#ifdef __G_CON

   //!single pointer to a BlockVector<PHM> object, the eigenvalues of G part of a SUP matrix will be stored here.
   BlockVector<PHM> *v_ph;

#endif

#ifdef __T1_CON

   //!single pointer to a BlockVector<DPM> object, the eigenvalues of T1 part of a SUP matrix will be stored here.
   BlockVector<DPM> *v_dp;

#endif

#ifdef __T2_CON

   //!single pointer to a BlockVector<PPHM> object, the eigenvalues of T2 part of a SUP matrix will be stored here.
   BlockVector<PPHM> *v_pph;

#endif

};

#endif
