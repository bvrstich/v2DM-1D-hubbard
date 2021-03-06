/**
 * @mainpage 
 * This is an implementation of a primal dual interior point method
 * for optimizing the second order density matrix using the P Q G T1 and T2 N-representability conditions for the special
 * case of the 1 dimensional Hubbard model with periodic boundary condition, which means that all the symmetries for this case have been 
 * implemented.
 * The method used is a path following algorithm with predictor corrector steps.
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
 * 
 * In the main the actual program is run.\n 
 * Part 1: An easy initial point is taken and then centered to the required precision (flag == 0)\n
 * Part 2: When the primal dual point is sufficiently centered steps are taken to reduce the
 * primal dual gap and take a large step in that direction (predictor) (flag == 1)\n
 * After each step a correcting step (flag == 2) is taken that brings the primal dual point closer to
 * the central path.\n
 * Part 3: When the primal dual gap is smaller that the required accuracy exit the while. (flag == 3)\n
 * For more information on the actual method, see primal_dual.pdf
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
   SPM::init();

#ifdef __G_CON
   PHM::init();
#endif

#ifdef __T1_CON
   DPM::init();
#endif

#ifdef __T2_CON
   PPHM::init();
#endif

   //hamiltoniaan
   TPM ham;
   ham.hubbard(U);

   TPM ham_copy(ham);

   //only traceless hamiltonian needed in program.
   ham.proj_Tr();

   //primal
   SUP X;

   //dual
   SUP Z;

   //Lagrange multiplier
   SUP V;

   //just dubya
   SUP W;

   SUP u_0;

   //little help
   TPM hulp;

   u_0.tpm(0).unit();

   u_0.fill();

   X = 0.0;
   Z = 0.0;

   //what does this do?
   double sigma = 1.0;

   double tolerance = 1.0e-5;

   double D_conv(1.0),P_conv(1.0),convergence(1.0);

   double mazzy = 1.6;

   int iter;
   int max_iter = 1;

   int tot_iter;

   while(D_conv > tolerance || P_conv > tolerance || fabs(convergence) > tolerance){

      D_conv = 1.0;

      iter = 0;

      while(D_conv > tolerance && iter <= max_iter){

         tot_iter++;

         ++iter;

         //solve system
         SUP B(Z);

         B -= u_0;

         B.daxpy(mazzy/sigma,X);

         TPM b;

         b.collaps(1,B);

         b.daxpy(-mazzy/sigma,ham);

         hulp.S(-1,b);

         //hulp is the matrix containing the gamma_i's
         hulp.proj_Tr();

         //construct W
         W.fill(hulp);

         W += u_0;

         W.daxpy(-1.0/sigma,X);

         //update Z and V with eigenvalue decomposition:
         W.sep_pm(Z,V);

         V.dscal(-sigma);

         //check infeasibility of the dual problem:
         TPM v;

         v.collaps(1,V);

         v -= ham;

         D_conv = sqrt(v.ddot(v));

     }

      //update primal:
      X = V;

      //check primal feasibility (W is a helping variable now)
      W.fill(hulp);

      W += u_0;

      W -= Z;

      P_conv = sqrt(W.ddot(W));

      convergence = Z.tpm(0).ddot(ham) + X.ddot(u_0);

      cout << P_conv << "\t" << D_conv << "\t" << sigma << "\t" << convergence << "\t" << Z.tpm(0).ddot(ham_copy) << endl;

      if(D_conv < P_conv)
         sigma *= 1.01;
      else
         sigma /= 1.01;

   }

   cout << endl;
   cout << "Energy: " << ham_copy.ddot(Z.tpm(0)) << endl;
   cout << "pd gap: " << Z.ddot(X) << endl;
   cout << "dual conv: " << D_conv << endl;
   cout << "primal conv: " << P_conv << endl;

   cout << endl;
   cout << tot_iter << endl;

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
