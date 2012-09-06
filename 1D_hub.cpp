/**
 * @mainpage 
 * This is an implementation of the dual-only potential reduction method
 * for optimizing the second order density matrix using the P Q G T1 and T2 N-representability conditions for the special
 * case of the 1 dimensional Hubbard model with periodic boundary condition, which means that all the symmetries for this case have been 
 * implemented.
 * At compile time you can decide which condtions will be active compile with make PQ, PQG, PQGT1, PQGT2 or PQGT=(for all conditions).
 * @author Brecht Verstichel
 * @date 10-05-2010
 */
#include <iostream>
#include <fstream>
#include <time.h>
#include <cmath>

using std::cout;
using std::endl;
using std::ofstream;

#include "include.h"

/**
 * In the main the actual program is run.\n 
 * We start from the unity density matrix normed on the particle number and minimize the 
 * ojective function:\n\n
 * Tr (Gamma H) - t * ln(det P(Gamma)) \n\n
 * Once the minimum is found the parameter t is reduced and a new search is initiated,
 * this goes on until convergence is reached.\n
 * The potential is minimized using the Newton-Raphson method and the resulting linear system
 * is solved via the linear conjugate gradient method.
 */
int main(int argc,char *argv[]){

   //initialize the random nr generator
   srand(time(NULL));

   cout.precision(10);

   int L = atoi(argv[1]);//dimension of the lattice, nr of sites
   int N = atoi(argv[2]);//nr of particles

   double U = atof(argv[3]);//onsite repulsion

   Tools::init(L,N);

   TPM::init();
   SPM::init(L,N);

#ifdef __G_CON
   PHM::init(L,N);
#endif

#ifdef __T1_CON
   DPM::init(L,N);
#endif

#ifdef __T2_CON
   PPHM::init(L,N);
#endif

   SUP::init(L,N);
   EIG::init(L,N);

   TPM ham;
   ham.hubbard(U);

   TPM rdm;
   rdm.unit();

   TPM backup_rdm(rdm);

   double t = 1.0;
   double tolerance = 1.0e-15;

   int n = 2*L*(2*L  - 1) + 4*L*L;

   //outer iteration: scaling of the potential barrier
   while(t * n > 1.0e-5){

      cout << n*t << "\t" << rdm.trace() << "\t" << rdm.ddot(ham) << "\t";

      int nr_cg_iter = 0;
      int nr_newton_iter = 0;

      double convergence = 1.0;

      //inner iteration: 
      //Newton's method for finding the minimum of the current potential
      while(convergence > tolerance){

         ++nr_newton_iter;

         SUP P;

         P.fill(rdm);

         P.invert();

         //eerst -gradient aanmaken:
         TPM grad;

         grad.constr_grad(t,ham,P);

         //dit wordt de stap:
         TPM delta;

         //los het hessiaan stelsel op:
         nr_cg_iter += delta.solve(t,P,grad);

         //line search
         double a = delta.line_search(t,P,ham);

         //rdm += a*delta;
         rdm.daxpy(a,delta);

         convergence = a*a*delta.ddot(delta);

      }

      cout << nr_newton_iter << "\t" << nr_cg_iter << endl;

      t /= 1.5;

   }

   cout << endl;

   cout << "Final Energy:\t" << ham.ddot(rdm) << endl;
   cout << endl;
   cout << "Final Spin:\t" << rdm.spin() << endl;

#ifdef __T2_CON
   PPHM::clear();
#endif

#ifdef __T1_CON
   DPM::clear();
#endif

#ifdef __G_CON
   PHM::clear();
#endif

   SPM::clear();
   TPM::clear();

   Tools::clear();

   return 0;

}
