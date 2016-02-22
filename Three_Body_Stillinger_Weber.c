/*
*
* CDDL HEADER START
*
* The contents of this file are subject to the terms of the Common Development
* and Distribution License Version 1.0 (the "License").
*
* You can obtain a copy of the license at
* http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
* specific language governing permissions and limitations under the License.
*
* When distributing Covered Code, include this CDDL HEADER in each file and
* include the License file in a prominent location with the name LICENSE.CDDL.
* If applicable, add the following below this CDDL HEADER, with the fields
* enclosed by brackets "[]" replaced with your own identifying information:
*
* Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
*
* CDDL HEADER END
*

*
* Copyright (c) 2012, Regents of the University of Minnesota.  All rights reserved.
*
* Contributors:
*    Ryan S. Elliott
*    Ellad B. Tadmor
*    Valeriu Smirichinski
*    Amit Singh
*/

/*******************************************************************************
*
*  Three_Body_Stillinger_Weber
*
*  Stillinger-Weber type Three-Body potential KIM Model Driver
*
*  Language: C
*
*  Release:
*
*******************************************************************************/
/*******************************************************************************

* For three-body potential any potential energy function for N-particle system can be written in the
* following form of two-body and three-body terms:
     F(1, 2, ...., N)   = sum_{i; 1 <= i < j <= N} phi_two(r; r = r_ij)
                          + sum_{i; 1 <= i \neq j < k <= N} phi_three(r_ij, r_ik, r_jk);

* For Stillinger-Weber, two-body term phi_two is

     phi_two(r) = epsilon * f_2(r_cap; r_cap = r/sigma); where f_2 is

     f_2(r_cap) = A * ( B*r_cap^(-p) - r_cap^-q ) * exp(1/(r_cap - a)); when  r_cap <  a
                = 0   when r_cap >= a

* and three-body term phi_three is

     phi_three(r_ij, r_ik, r_jk) = epsilon * f_3(r1_cap, r2_cap, r3_cap); where

     r1_cap = r_ij/sigma, r2_cap = r_ik/sigma, r3_cap = r_jk/sigma,      and

     f_3(r1_cap, r2_cap, r3_cap) = lambda * (exp(gamma((1/(r1_cap - a)) + (1/(r2_cap - a))))) * (costheta_jik - costheta_0)^2;  when r1_cap < a && r2_cap < a
                                 = 0;                                                                                           otherwise

     costheta_jik = (r_ij^2 + r_ik^2 - r_jk^2) / (2*r_ij*r_ik).

*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "KIM_API_C.h"
#include "KIM_API_status.h"

#define THIRD 1.0/3.0
/******************************************************************************
* Below are the definitions and values of all Model parameters
*******************************************************************************/
#define DIM 3       /* dimensionality of space */
#define SPECCODE 1  /* internal species code */

/* Define prototypes for Model Driver init */
/* must be all lowercase to be compatible with the KIM API (to support Fortran Tests) */
/**/
int model_driver_init(void* km, char* paramfile_names, int* nmstrlen, int* numparamfiles);

/* Define prototypes for Model (Driver) reinit, compute, and destroy */
/* defined as static to avoid namespace clashes with other Models    */
/**/
static int reinit(void* km);
static int destroy(void* km);
static int compute(void* km);
/**/
static void calc_phi_two(double* A, double* B, double* p, double* q, double* a, double* sigma, double* epsilon,
                         double r, double* phi);
static void calc_phi_dphi_two(double* A, double* B, double* p, double* q, double* a, double* sigma, double* epsilon,
                              double r, double* phi, double* dphi);
static void calc_phi_d2phi_two(double* A, double* B, double* p, double* q, double* a, double* sigma, double* epsilon,
                               double r, double* phi, double* dphi, double* d2phi);

static void calc_phi_three(double* a, double* lambda, double* gamma, double* sigma, double* epsilon, double* costhetat,
                           double rij, double rik, double rjk, double* phi);
static void calc_phi_dphi_three(double* a, double* lambda, double* gamma, double* sigma, double* epsilon, double* costhetat,
                                double rij, double rik, double rjk, double* phi, double* dphi);
static void calc_phi_d2phi_three(double* a, double* lambda, double* gamma, double* sigma, double* epsilon, double* costhetat,
                                 double rij, double rik, double rjk, double* phi, double* dphi, double* d2phi);

/* Define model_buffer structure */
struct model_buffer {
   int NBC;
   int HalfOrFull;
   int IterOrLoca;
   int energy_ind;
   int forces_ind;
   int particleEnergy_ind;
   int process_dEdr_ind;
   int process_d2Edr2_ind;
   int model_index_shift;
   int numberOfParticles_ind;
   int particleSpecies_ind;
   int coordinates_ind;
   int numberContributingParticles_ind;
   int boxSideLengths_ind;
   int get_neigh_ind;
   int cutoff_ind;


   double* Pcutoff;
   double* cutsq;
   double* A;
   double* B;
   double* p;
   double* q;
   double* a;
   double* lambda;
   double* gamma;
   double* sigma;
   double* epsilon;
   double* costhetat;
};


/* Calculate pair potential phi_two(r) */
static void calc_phi_two(double* A, double* B, double* p, double* q, double* a, double* sigma, double* epsilon,
                         double r, double* phi)
{
   /* Local variables */
   double r_cap;

   r_cap = r/(*sigma);

   if (r_cap >=  *a)
   {
      *phi = 0.0;
   }
   else
   {
      *phi = (*epsilon) * (*A) * ( (*B) * pow(r_cap,-(*p)) - pow(r_cap,-(*q)) ) * exp(1/(r_cap - *a));
   }

   return;
}

/* Calculate pair potential phi_two(r) and its derivative dphi_two(r) */
static void calc_phi_dphi_two(double* A, double* B, double* p, double* q, double* a, double* sigma, double* epsilon,
                              double r, double* phi, double* dphi)
{
   /* Local variables */
   double r_cap;

   r_cap = r/(*sigma);

   if (r_cap >= *a)
   {
      *phi = 0.0;
      *dphi = 0.0;
   }
   else
   {
      *phi = (*epsilon) * (*A) * ( (*B) * pow(r_cap,-(*p)) - pow(r_cap,-(*q)) ) * exp(1/(r_cap - *a));

      *dphi =  ( (*q) * pow(r_cap,-((*q)+1)) - (*p * (*B)) * pow(r_cap,-((*p)+1)) )
                  - ( (*B) * pow(r_cap,-(*p)) - pow(r_cap,-(*q)) ) * pow((r_cap - *a),-2);
      *dphi *= (*epsilon / *sigma) * (*A) * exp(1/(r_cap - *a));
   }

   return;
}

/* Calculate pair potential phi_two(r) and its 1st & 2nd derivatives dphi_two(r), d2phi_two(r) */
static void calc_phi_d2phi_two(double* A, double* B, double* p, double* q, double* a, double* sigma, double* epsilon,
                               double r, double* phi, double* dphi, double* d2phi)
{
   /* Local variables */
   double r_cap;

   r_cap = r/(*sigma);

   if (r_cap >= *a)
   {
      *phi = 0.0;
      *dphi = 0.0;
      *d2phi = 0.0;
   }
   else
   {
      *phi = (*epsilon) * (*A) * ( (*B) * pow(r_cap,-(*p)) - pow(r_cap,-(*q)) ) * exp(1/(r_cap - *a));

      *dphi =  ( (*q) * pow(r_cap,-((*q)+1)) - (*p * *B) * pow(r_cap,-((*p)+1)) )
                  - ( (*B) * pow(r_cap,-(*p)) - pow(r_cap,-(*q)) ) * pow((r_cap - *a),-2);
      *dphi *= (*epsilon / *sigma) * (*A) * exp(1/(r_cap - *a));

      *d2phi = ( *B * pow(r_cap,-(*p)) - pow(r_cap,-(*q)) ) * ( pow((r_cap - *a),-4) + 2*pow((r_cap - *a),-3) )
                  + 2 * ( *p * *B *pow(r_cap,-(*p + 1)) - *q * pow(r_cap,-(*q + 1)) ) * pow((r_cap - *a),-2)
                  + ( *p * (*p + 1) *  *B * pow(r_cap,-(*p + 2)) - *q * (*q + 1) * pow(r_cap,-(*q + 2)) );
      *d2phi *= (*epsilon / (*sigma * *sigma)) * (*A) * exp(1/(r_cap - *a));
   }

   return;

}

/* Calculate pair potential phi_three(rij, rik, rjk) */
static void calc_phi_three(double* a, double* lambda, double* gamma, double* sigma, double* epsilon, double* costhetat,
                           double rij, double rik, double rjk, double* phi)
{
   /* local variables */
   double c1;

   double rij_cap;
   double rik_cap;
   double rjk_cap;

   double costhetajik;

   double diff_costhetajik;

   double exp_ij_ik;

   c1 = *lambda * *epsilon;

   rij_cap = rij/(*sigma);
   rik_cap = rik/(*sigma);
   rjk_cap = rjk/(*sigma);

   costhetajik = (pow(rij,2) + pow(rik,2) - pow(rjk,2))/(2*rij*rik);

   diff_costhetajik = costhetajik - *costhetat;

   exp_ij_ik = exp((*gamma) * (1/(rij_cap - *a) + 1/(rik_cap - *a)));

   if ((rij_cap < *a) && (rik_cap < *a))
   {
     *phi  = c1 * exp_ij_ik * diff_costhetajik *  diff_costhetajik;
   }

   return;
}

/* Calculate pair potential phi_three(rij, rik, rjk) and its 1st derivative
  dphi_three(rij, rik, rjk)

  dphi has three components as derivatives of phi w.r.t. rij, rik, rjk
*/
static void calc_phi_dphi_three(double* a, double* lambda, double* gamma, double* sigma, double* epsilon, double* costhetat,
                                double rij, double rik, double rjk, double* phi, double* dphi)
{
   /* local variables */
   double c1;
   double c2;

   double rij_cap;
   double rik_cap;
   double rjk_cap;

   double costhetajik;

   double diff_costhetajik;

   double costhetajik_ij;
   double costhetajik_ik;
   double costhetajik_jk;

   double exp_ij_ik;

   double d_ij;
   double d_ik;
   double d_jk;

   c1 = *lambda * *epsilon;
   c2 = c1 / *sigma;

   rij_cap = rij/(*sigma);
   rik_cap = rik/(*sigma);
   rjk_cap = rjk/(*sigma);

   costhetajik = (pow(rij,2) + pow(rik,2) - pow(rjk,2))/(2*rij*rik);

   diff_costhetajik = costhetajik - *costhetat;

   /* Derivatives of cosines w.r.t rij, rik, rjk */
   costhetajik_ij = (*sigma)*(pow(rij,2) - pow(rik,2) + pow(rjk,2))/(2*rij*rij*rik);
   costhetajik_ik = (*sigma)*(pow(rik,2) - pow(rij,2) + pow(rjk,2))/(2*rij*rik*rik);
   costhetajik_jk = -(*sigma)*rjk/(rij*rik);

   /* Variables for simplifying terms */
   exp_ij_ik = exp((*gamma) * (1/(rij_cap - *a) + 1/(rik_cap - *a)));

   d_ij = -(*gamma)*pow((rij_cap - *a),-2);
   d_ik = -(*gamma)*pow((rik_cap - *a),-2);
   d_jk = -(*gamma)*pow((rjk_cap - *a),-2);


   if ((rij_cap < *a) && (rik_cap < *a))
   {
     *phi  = c1 * exp_ij_ik * diff_costhetajik *  diff_costhetajik;
     dphi[0] = c2 * diff_costhetajik * exp_ij_ik *  (d_ij * diff_costhetajik   +  2 * costhetajik_ij);
     dphi[1] = c2 * diff_costhetajik * exp_ij_ik *  (d_ik * diff_costhetajik   +  2 * costhetajik_ik);
     dphi[2] = c2 * diff_costhetajik * exp_ij_ik *  2 * costhetajik_jk;
   }

   return;
}

/* Calculate pair potential phi_three(rij, rik, rjk) and its 1st & 2nd derivatives
  dphi_three(rij, rik, rjk), d2phi_three(rij, rik, rjk)

  dphi has three components as derivatives of phi w.r.t. rij, rik, rjk

  d2phi as symmetric Hessian matrix of phi has six components: [0]=(ij,ij), [3]=(ij,ik), [4]=(ij,jk)
                                                                            [1]=(ik,ik), [5]=(ik,jk)
                                                                                         [2]=(jk,jk)
*/
static void calc_phi_d2phi_three(double* a, double* lambda, double* gamma, double* sigma, double* epsilon, double* costhetat,
                                 double rij, double rik, double rjk, double* phi, double* dphi, double* d2phi)
{
    /* local variables */
   double c1;
   double c2;
   double c3;

   double rij_cap;
   double rik_cap;
   double rjk_cap;

   double costhetajik;

   double diff_costhetajik;

   double diff_costhetajik_2;

   double costhetajik_ij;
   double costhetajik_ik;
   double costhetajik_jk;

   double costhetajik_ij_ij;
   double costhetajik_ik_ik;
   double costhetajik_jk_jk;
   double costhetajik_ij_ik;
   double costhetajik_ij_jk;
   double costhetajik_ik_jk;

   double exp_ij_ik;

   double d_ij;
   double d_ik;
   double d_jk;

   double d_ij_2;
   double d_ik_2;
   double d_jk_2;

   double dd_ij;
   double dd_ik;
   double dd_jk;

   c1 = *lambda * *epsilon;
   c2 = c1 / *sigma;
   c3 = c2 / *sigma;

   rij_cap = rij/(*sigma);
   rik_cap = rik/(*sigma);
   rjk_cap = rjk/(*sigma);

   costhetajik = (pow(rij,2) + pow(rik,2) - pow(rjk,2))/(2*rij*rik);

   diff_costhetajik = costhetajik - *costhetat;

   diff_costhetajik_2 = diff_costhetajik *  diff_costhetajik;

   /* Derivatives of cosines w.r.t. r_ij, r_ik, r_jk */
   costhetajik_ij = (*sigma)*(pow(rij,2) - pow(rik,2) + pow(rjk,2))/(2*rij*rij*rik);
   costhetajik_ik = (*sigma)*(pow(rik,2) - pow(rij,2) + pow(rjk,2))/(2*rij*rik*rik);
   costhetajik_jk = -(*sigma)*rjk/(rij*rik);

   /* Hessian matrix of cosine */
   costhetajik_ij_ij = (*sigma * *sigma)*(pow(rik,2) - pow(rjk,2))/(rij*rij*rij*rik);
   costhetajik_ik_ik = (*sigma * *sigma)*(pow(rij,2) - pow(rjk,2))/(rij*rik*rik*rik);
   costhetajik_jk_jk = -(*sigma * *sigma)/(rij*rik);
   costhetajik_ij_ik = -(*sigma * *sigma)*(pow(rij,2) + pow(rik,2) + pow(rjk,2))/(2*rij*rij*rik*rik);
   costhetajik_ij_jk = (*sigma * *sigma)*rjk/(rij*rij*rik);
   costhetajik_ik_jk = (*sigma * *sigma)*rjk/(rik*rik*rij);

   /* Variables for simplifying terms */
   exp_ij_ik = exp((*gamma) * (1/(rij_cap - *a) + 1/(rik_cap - *a)));

   d_ij = -(*gamma)*pow((rij_cap - *a),-2);
   d_ik = -(*gamma)*pow((rik_cap - *a),-2);
   d_jk = -(*gamma)*pow((rjk_cap - *a),-2);

   d_ij_2 = d_ij * d_ij;
   d_ik_2 = d_ik * d_ik;
   d_jk_2 = d_jk * d_jk;

   dd_ij = (2 * *gamma)*pow((rij_cap - *a),-3);
   dd_ik = (2 * *gamma)*pow((rik_cap - *a),-3);
   dd_jk = (2 * *gamma)*pow((rjk_cap - *a),-3);

   *phi  = 0.0;
   dphi[0] = dphi[1] = dphi[2] = 0.0;
   d2phi[0] = d2phi[1] = d2phi[2] = 0.0;
   d2phi[3] = d2phi[4] = d2phi[5] = 0.0;

   if ((rij_cap < *a) && (rik_cap < *a))
   {
     *phi  += exp_ij_ik * diff_costhetajik_2;
     dphi[0] += diff_costhetajik * exp_ij_ik *  (d_ij * diff_costhetajik   +  2 * costhetajik_ij);
     dphi[1] += diff_costhetajik * exp_ij_ik *  (d_ik * diff_costhetajik   +  2 * costhetajik_ik);
     dphi[2] += diff_costhetajik * exp_ij_ik *  2 * costhetajik_jk;

     d2phi[0] += exp_ij_ik * ((d_ij_2 + dd_ij) * diff_costhetajik_2 +
                  (4 * d_ij * costhetajik_ij + 2 * costhetajik_ij_ij) * diff_costhetajik + 2 * costhetajik_ij * costhetajik_ij);
     d2phi[1] += exp_ij_ik * ((d_ik_2 + dd_ik) * diff_costhetajik_2 +
                  (4 * d_ik * costhetajik_ik + 2 * costhetajik_ik_ik) * diff_costhetajik + 2 * costhetajik_ik * costhetajik_ik);
     d2phi[2] += 2 * exp_ij_ik * ( costhetajik_jk_jk * diff_costhetajik +  costhetajik_jk * costhetajik_jk);
     d2phi[3] += exp_ij_ik * (d_ij * d_ik * diff_costhetajik_2 +
                  (d_ij * costhetajik_ik + d_ik * costhetajik_ij + costhetajik_ij_ik) * 2 * diff_costhetajik + 2 * costhetajik_ij * costhetajik_ik);
     d2phi[4] += exp_ij_ik * ((d_ij * costhetajik_jk + costhetajik_ij_jk) * 2 * diff_costhetajik + 2 * costhetajik_ij * costhetajik_jk);
     d2phi[5] += exp_ij_ik * ((d_ik * costhetajik_jk + costhetajik_ik_jk) * 2 * diff_costhetajik + 2 * costhetajik_ik * costhetajik_jk);
   }

   *phi  *= c1;
   dphi[0] *= c2; /*  derivative is w.r.t. rij */
   dphi[1] *= c2; /*  derivative is w.r.t. rik */
   dphi[2] *= c2; /*  derivative is w.r.t. rjk */

   d2phi[0]  *= c3; /*  derivative is w.r.t. rij, rij */
   d2phi[1]  *= c3; /*  derivative is w.r.t. rik, rik */
   d2phi[2]  *= c3; /*  derivative is w.r.t. rjk, rjk */
   d2phi[3]  *= c3; /*  derivative is w.r.t. rij, rik */
   d2phi[4]  *= c3; /*  derivative is w.r.t. rij, rjk */
   d2phi[5]  *= c3; /*  derivative is w.r.t. rik, rjk */

   return;
}

/* compute function */
static int compute(void* km)
{
   /* local variables */
   intptr_t* pkim = *((intptr_t**) km);
   double R1;
   double R2;
   double R3;
   double R_pairs[2];
   double *pR_pairs = &(R_pairs[0]);
   double Rsqij;
   double Rsqjk;
   double Rsqik;
   double phi_two;
   double dphi_two;
   double d2phi_two;
   double dEidr_two;
   double d2Eidr_two;
   double phi_three;
   double* dphi_three;
   double* d2phi_three;
   double* dEidr_three;
   double* d2Eidr_three;
   double Rij[DIM];
   double Rjk[DIM];
   double Rik[DIM];
   double *pRij = &(Rij[0]);
   double *pRjk = &(Rjk[0]);
   double *pRik = &(Rik[0]);
   double Rij_pairs[2][3];
   double *pRij_pairs = &(Rij_pairs[0][0]);
   int ier;
   int i;
   int i_pairs[2];
   int *pi_pairs = &(i_pairs[0]);
   int j;
   int j_pairs[2];
   int *pj_pairs = &(j_pairs[0]);
   int k;
   int jj;
   int kk;
   int kdim;
   int currentAtom;
   int* neighListOfCurrentAtom;
   struct model_buffer* buffer;
   int comp_energy;
   int comp_force;
   int comp_particleEnergy;
   int comp_process_dEdr;
   int comp_process_d2Edr2;
   int NBC;
   int HalfOrFull;
   int IterOrLoca;
   int model_index_shift;
   int zero = 0;
   int one = 1;
   int request;

   int* nAtoms;
   int* particleSpecies;
   double* cutoff;
   double* cutsq;
   double* A;
   double* B;
   double* p;
   double* q;
   double* a;
   double* lambda;
   double* gamma;
   double* sigma;
   double* epsilon;
   double* costhetat;
   double* Rij_list;
   double* coords;
   double* energy;
   double* force;
   double* particleEnergy;
   double* boxSideLengths;
   int* numContrib;
   int numberContrib;
   int numOfAtomNeigh;
   typedef int (*get_neigh_ptr)(void *,int *,int *,int *, int *, int **, double **);
   get_neigh_ptr get_neigh;

   dphi_three = (double *) malloc(3*sizeof(double));
   d2phi_three = (double *) malloc(6*sizeof(double));
   dEidr_three = (double *) malloc(3*sizeof(double));
   d2Eidr_three = (double *) malloc(6*sizeof(double));

   /* get buffer from KIM object */
   buffer = (struct model_buffer*) KIM_API_get_model_buffer(pkim, &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
      return ier;
   }

   /* unpack info from the buffer */
   NBC = buffer->NBC;
   HalfOrFull = buffer->HalfOrFull;
   IterOrLoca = buffer->IterOrLoca;
   model_index_shift = buffer->model_index_shift;
   /* unpack the Model's parameters stored in the KIM API object */
   cutsq  = buffer->cutsq;
   A      = buffer->A;
   B      = buffer->B;
   p      = buffer->p;
   q      = buffer->q;
   a      = buffer->a;
   lambda = buffer->lambda;
   gamma  = buffer->gamma;
   sigma  = buffer->sigma;
   epsilon = buffer->epsilon;
   costhetat = buffer->costhetat;

   /* check to see if we have been asked to compute the forces, particleEnergy, and dEdr */
   KIM_API_getm_compute_by_index(pkim, &ier, 5*3,
                                 buffer->energy_ind,         &comp_energy,         1,
                                 buffer->forces_ind,         &comp_force,          1,
                                 buffer->particleEnergy_ind, &comp_particleEnergy, 1,
                                 buffer->process_dEdr_ind,   &comp_process_dEdr,   1,
                                 buffer->process_d2Edr2_ind, &comp_process_d2Edr2, 1);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_compute_by_index", ier);
      return ier;
   }

   KIM_API_getm_data_by_index(pkim, &ier, 9*3,
                              buffer->cutoff_ind,                      &cutoff,         1,
                              buffer->numberOfParticles_ind,           &nAtoms,         1,
                              buffer->particleSpecies_ind,               &particleSpecies,  1,
                              buffer->coordinates_ind,                 &coords,         1,
                              buffer->numberContributingParticles_ind, &numContrib,     (HalfOrFull==1),
                              buffer->boxSideLengths_ind,              &boxSideLengths, (NBC==1),
                              buffer->energy_ind,                      &energy,         comp_energy,
                              buffer->forces_ind,                      &force,          comp_force,
                              buffer->particleEnergy_ind,              &particleEnergy, comp_particleEnergy);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data_by_index", ier);
      return ier;
   }
   if (NBC!=3)
   {
      get_neigh = (get_neigh_ptr) KIM_API_get_method_by_index(pkim, buffer->get_neigh_ind, &ier);
      if (KIM_STATUS_OK > ier)
      {
         KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_method_by_index", ier);
         return ier;
      }
   }

   if (HalfOrFull == 1)
   {
      if (3 != NBC) /*    non-CLUSTER cases */
      {
         numberContrib = *numContrib;
      }
      else /*   CLUSTER cases */
      {
         numberContrib = *nAtoms;
      }
   }
   else
   {  /* provide initialization even if not used */
      numberContrib = *nAtoms;
   }

   /* Check to be sure that the atom types are correct */
   /**/
   ier = KIM_STATUS_FAIL; /* assume an error */
   for (i = 0; i < *nAtoms; ++i)
   {
      if ( SPECCODE != particleSpecies[i])
      {
         KIM_API_report_error(__LINE__, __FILE__, "Unexpected species type detected", ier);
         return ier;
      }
   }
   ier = KIM_STATUS_OK; /* everything is ok */

   /* initialize potential energies, forces, and virial term */
   if (comp_particleEnergy)
   {
      for (i = 0; i < *nAtoms; ++i)
      {
         particleEnergy[i] = 0.0;
      }
   }
   else if (comp_energy)
   {
      *energy = 0.0;
   }

   if (comp_force)
   {
      for (i = 0; i < *nAtoms; ++i)
      {
         for (kdim = 0; kdim < DIM; ++kdim)
         {
            force[i*DIM + kdim] = 0.0;
         }
      }
   }

   /* Initialize neighbor handling for CLUSTER NBC */
   if (3 == NBC) /* CLUSTER */
   {
      neighListOfCurrentAtom = (int *) malloc((*nAtoms)*sizeof(int));
   }

   /* Initialize neighbor handling for Iterator mode */

   if (1 == IterOrLoca)
   {
      ier = (*get_neigh)(&pkim, &zero, &zero, &currentAtom, &numOfAtomNeigh,
                          &neighListOfCurrentAtom, &Rij_list);
      /* check for successful initialization */
      if (KIM_STATUS_NEIGH_ITER_INIT_OK != ier)
      {
         KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh", ier);
         ier = KIM_STATUS_FAIL;
         return ier;
      }
   }

   /* Compute enery and forces */

   /* loop over particles and compute enregy and forces */
   i = -1;
   while( 1 )
   {

      /* Set up neighbor list for next atom for all NBC methods */
      if (1 == IterOrLoca) /* ITERATOR mode */
      {
         ier = (*get_neigh)(&pkim, &zero, &one, &currentAtom, &numOfAtomNeigh,
                             &neighListOfCurrentAtom, &Rij_list);
         if (KIM_STATUS_NEIGH_ITER_PAST_END == ier) /* the end of the list, terminate loop */
         {
            break;
         }
         if (KIM_STATUS_OK > ier) /* some sort of problem, exit */
         {
            KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh", ier);
            return ier;
         }

         i = currentAtom + model_index_shift;
      }
      else
      {
         i++;
         if (*nAtoms <= i) /* incremented past end of list, terminate loop */
         {
            break;
         }

         if (3 == NBC) /* CLUSTER NBC method */
         {
            numOfAtomNeigh = *nAtoms - (i + 1);
            for (kdim = 0; kdim < numOfAtomNeigh; ++kdim)
            {
               neighListOfCurrentAtom[kdim] = i + kdim + 1 - model_index_shift;
            }
            ier = KIM_STATUS_OK;
         }
         else
         {
            request = i - model_index_shift;
            ier = (*get_neigh)(&pkim, &one, &request,
                                &currentAtom, &numOfAtomNeigh,
                                &neighListOfCurrentAtom, &Rij_list);
            if (KIM_STATUS_OK != ier) /* some sort of problem, exit */
            {
               KIM_API_report_error(__LINE__, __FILE__, "get_neigh", ier);
               ier = KIM_STATUS_FAIL;
               return ier;
            }
         }
      }

      /* loop over the neighbors of atom i */
      for (jj = 0; jj < numOfAtomNeigh; ++ jj)
      {

         j = neighListOfCurrentAtom[jj] + model_index_shift; /* get neighbor ID */

         /* compute relative position vector and squared distance */
         Rsqij = 0.0;
         for (kdim = 0; kdim < DIM; ++kdim)
         {
            if (0 != NBC) /* all methods except NEIGH_RVEC */
            {
               Rij[kdim] = coords[j*DIM + kdim] - coords[i*DIM + kdim];
            }
            else          /* NEIGH_RVEC method */
            {
               Rij[kdim] = Rij_list[jj*DIM + kdim];
            }

            /* apply periodic boundary conditions if required */
            if (2 == NBC)
            {
               if (abs(Rij[kdim]) > 0.5*boxSideLengths[kdim])
               {
                  Rij[kdim] -= (Rij[kdim]/fabs(Rij[kdim]))*boxSideLengths[kdim];
               }
            }

            /* compute squared distance */
            Rsqij += Rij[kdim]*Rij[kdim];
         }

         /* compute energy and force */
         if (Rsqij > *cutsq) continue; /* particles are not interacting  */
         R1 = sqrt(Rsqij);
         if (comp_process_d2Edr2)
         {
             /* compute pair potential and its derivatives */
             calc_phi_d2phi_two(A, B, p, q, a, sigma, epsilon,
                                R1, &phi_two, &dphi_two, &d2phi_two);

             /* compute dEidr */
             if ((1 == HalfOrFull) && (j < numberContrib))
             {
                  /* Half mode -- double contribution */
                   dEidr_two = dphi_two;
                   d2Eidr_two = d2phi_two;
             }
             else
             {
                  /* Full mode -- regular contribution */
                  dEidr_two  = 0.5*dphi_two;
                  d2Eidr_two = 0.5*d2phi_two;
             }
         }
         else if (comp_force || comp_process_dEdr)
         {
             /* compute pair potential and its derivative */
             calc_phi_dphi_two(A, B, p, q, a, sigma, epsilon,
                               R1, &phi_two, &dphi_two);

             /* compute dEidr */
             if ((1 == HalfOrFull) && (j < numberContrib))
             {
                  /* Half mode -- double contribution */
                  dEidr_two = dphi_two;
             }
             else
             {

                  /* Full mode -- regular contribution */
                  dEidr_two = 0.5*dphi_two;
             }
         }
         else
         {
            /* compute just pair potential */
             calc_phi_two(A, B, p, q, a, sigma, epsilon,
                          R1, &phi_two);
         }

         /* contribution to energy */
         if (comp_particleEnergy)
         {
            particleEnergy[i] += 0.5*phi_two;
            /* if half list add energy for the other atom in the pair */
            if ((1 == HalfOrFull) && (j < numberContrib)) particleEnergy[j] += 0.5*phi_two;
         }
         if (comp_energy)
         {
             if ((1 == HalfOrFull) && (j < numberContrib))
             {
                 /* Half mode -- add v to total energy */
                 *energy += phi_two;
             }
             else
             {
                 /* Full mode -- add half v to total energy */
                 *energy += 0.5*phi_two;
             }
         }

         /* contribution to process_dEdr */
         if (comp_process_dEdr)
         {
            ier = KIM_API_process_dEdr(km, &dEidr_two, &R1, &pRij, &i, &j);
         }

         /* contribution to process_d2Edr2 */
         if (comp_process_d2Edr2)
         {
             R_pairs[0] = R_pairs[1] = R1;
             Rij_pairs[0][0] = Rij_pairs[1][0] = Rij[0];
             Rij_pairs[0][1] = Rij_pairs[1][1] = Rij[1];
             Rij_pairs[0][2] = Rij_pairs[1][2] = Rij[2];
             i_pairs[0] = i_pairs[1] = i;
             j_pairs[0] = j_pairs[1] = j;

             ier = KIM_API_process_d2Edr2(km, &d2Eidr_two, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
         }

         /* contribution to forces */
         if (comp_force)
         {
           for (kdim = 0; kdim < DIM; ++kdim)
           {
               force[i*DIM + kdim] += dEidr_two*Rij[kdim]/R1; /* accumulate force on atom i */
               force[j*DIM + kdim] -= dEidr_two*Rij[kdim]/R1; /* accumulate force on atom j */
           }
         }


         /* Start adding three body terms */
         /*********************************/
         if( jj == numOfAtomNeigh-1) continue;

         for (kk = jj+1; kk < numOfAtomNeigh; ++kk)
         {

           k = neighListOfCurrentAtom[kk] + model_index_shift; /* get neighbor ID */

           /* compute relative position vector and squared distance */
           Rsqik = 0.0;
           Rsqjk = 0.0;
           for (kdim = 0; kdim < DIM; ++kdim)
           {
              if (0 != NBC) /* all methods except NEIGH_RVEC */
              {
                 Rik[kdim] = coords[k*DIM + kdim] - coords[i*DIM + kdim];
                 Rjk[kdim] = Rik[kdim] - Rij[kdim];
              }
              else          /* NEIGH_RVEC method */
              {
                 Rik[kdim] = Rij_list[kk*DIM + kdim];
                 Rjk[kdim] = Rik[kdim] - Rij[kdim];
              }

              /* apply periodic boundary conditions if required */
              if (2 == NBC)
              {
                 if (abs(Rik[kdim]) > 0.5*boxSideLengths[kdim])
                 {
                    Rik[kdim] -= (Rik[kdim]/fabs(Rik[kdim]))*boxSideLengths[kdim];
                 }
                 if (abs(Rjk[kdim]) > 0.5*boxSideLengths[kdim])
                 {
                    Rjk[kdim] -= (Rjk[kdim]/fabs(Rjk[kdim]))*boxSideLengths[kdim];
                 }
              }

              /* compute squared distance */
              Rsqik += Rik[kdim]*Rik[kdim];
              Rsqjk += Rjk[kdim]*Rjk[kdim];
           }

           /* compute energy and force */

           if (Rsqik > *cutsq) continue; /* particles are interacting ? */

            R2 = sqrt(Rsqik);
            R3 = sqrt(Rsqjk);

            if (comp_process_d2Edr2)
            {
               /* compute three-body potential and its derivatives */
               calc_phi_d2phi_three(a, lambda, gamma, sigma, epsilon, costhetat,
                                    R1, R2, R3, &phi_three, dphi_three, d2phi_three);

               /* compute dEidr */
               if ((1 == HalfOrFull))
               {
                  /*  Half mode */
                  ier = KIM_STATUS_FAIL;
                  KIM_API_report_error(__LINE__, __FILE__, "The driver does not support half list mode", ier);
                  return ier;
                }
                else
                {
                  /* Full mode -- regular contribution */
                  dEidr_three[0]  = dphi_three[0];
                  dEidr_three[1]  = dphi_three[1];
                  dEidr_three[2]  = dphi_three[2];

                  d2Eidr_three[0] = d2phi_three[0];
                  d2Eidr_three[1] = d2phi_three[1];
                  d2Eidr_three[2] = d2phi_three[2];
                  d2Eidr_three[3] = d2phi_three[3];
                  d2Eidr_three[4] = d2phi_three[4];
                  d2Eidr_three[5] = d2phi_three[5];
                }
            }
            else if (comp_force || comp_process_dEdr)
            {
               /* compute three-body potential and its derivative */
               calc_phi_dphi_three(a, lambda, gamma, sigma, epsilon, costhetat,
                                   R1, R2, R3, &phi_three, dphi_three);

               /* compute dEidr */
               if ((1 == HalfOrFull))
               {
                  /* Half mode */
                  ier = KIM_STATUS_FAIL;
                  KIM_API_report_error(__LINE__, __FILE__, "The driver does not support half list mode", ier);
                  return ier;
               }
               else
               {
                  /* Full mode -- regular contribution */
                  dEidr_three[0]  =  dphi_three[0];
                  dEidr_three[1]  =  dphi_three[1];
                  dEidr_three[2]  =  dphi_three[2];
               }

            }
            else
            {
               /* compute just three-body potential */
               calc_phi_three(a, lambda, gamma, sigma, epsilon, costhetat,
                             R1, R2, R3, &phi_three);
            }

            /* contribution to energy */
            if (comp_particleEnergy)
            {
               particleEnergy[i] += phi_three;
               /* if half list add energy for other atoms in the triplet */
               if ((1 == HalfOrFull)) {
                  ier = KIM_STATUS_FAIL;
                  KIM_API_report_error(__LINE__, __FILE__, "The driver does not support half list mode", ier);
                  return ier;
               }
            }
            if (comp_energy)
            {
               if ((1 == HalfOrFull))
               {
                  /* Half mode */
                  ier = KIM_STATUS_FAIL;
                  KIM_API_report_error(__LINE__, __FILE__, "The driver does not support half list mode", ier);
                  return ier;
               }
               else
               {
                  /* Full mode -- add half v to total energy */
                  *energy +=  phi_three;
               }
            }

            /* contribution to process_dEdr */
            if (comp_process_dEdr)
            {
               ier = KIM_API_process_dEdr(km, dEidr_three, &R1, &pRij, &i, &j);
               ier = KIM_API_process_dEdr(km, dEidr_three+1, &R2, &pRik, &i, &k);
               ier = KIM_API_process_dEdr(km, dEidr_three+2, &R3, &pRjk, &j, &k);
            }

            /* contribution to process_d2Edr2 */
            if (comp_process_d2Edr2)
            {
               R_pairs[0] = R_pairs[1] = R1;
               Rij_pairs[0][0] = Rij_pairs[1][0] = Rij[0];
               Rij_pairs[0][1] = Rij_pairs[1][1] = Rij[1];
               Rij_pairs[0][2] = Rij_pairs[1][2] = Rij[2];
               i_pairs[0] = i_pairs[1] = i;
               j_pairs[0] = j_pairs[1] = j;

               ier = KIM_API_process_d2Edr2(km, d2Eidr_three, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);

               R_pairs[0] = R_pairs[1] = R2;
               Rij_pairs[0][0] = Rij_pairs[1][0] = Rik[0];
               Rij_pairs[0][1] = Rij_pairs[1][1] = Rik[1];
               Rij_pairs[0][2] = Rij_pairs[1][2] = Rik[2];
               i_pairs[0] = i_pairs[1] = i;
               j_pairs[0] = j_pairs[1] = k;

               ier = KIM_API_process_d2Edr2(km, d2Eidr_three+1, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);

               R_pairs[0] = R_pairs[1] = R3;
               Rij_pairs[0][0] = Rij_pairs[1][0] = Rjk[0];
               Rij_pairs[0][1] = Rij_pairs[1][1] = Rjk[1];
               Rij_pairs[0][2] = Rij_pairs[1][2] = Rjk[2];
               i_pairs[0] = i_pairs[1] = j;
               j_pairs[0] = j_pairs[1] = k;

               ier = KIM_API_process_d2Edr2(km, d2Eidr_three+2, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);


               R_pairs[0] = R1;
               R_pairs[1] = R2;
               Rij_pairs[0][0] =  Rij[0];
               Rij_pairs[0][1] =  Rij[1];
               Rij_pairs[0][2] =  Rij[2];
               Rij_pairs[1][0] =  Rik[0];
               Rij_pairs[1][1] =  Rik[1];
               Rij_pairs[1][2] =  Rik[2];
               i_pairs[0]  = i;
               j_pairs[0]  = j;
               i_pairs[1]  = i;
               j_pairs[1]  = k;

               ier = KIM_API_process_d2Edr2(km, d2Eidr_three+3, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);

               R_pairs[0] = R2;
               R_pairs[1] = R1;
               Rij_pairs[0][0] =  Rik[0];
               Rij_pairs[0][1] =  Rik[1];
               Rij_pairs[0][2] =  Rik[2];
               Rij_pairs[1][0] =  Rij[0];
               Rij_pairs[1][1] =  Rij[1];
               Rij_pairs[1][2] =  Rij[2];
               i_pairs[0]  = i;
               j_pairs[0]  = k;
               i_pairs[1]  = i;
               j_pairs[1]  = j;

               ier = KIM_API_process_d2Edr2(km, d2Eidr_three+3, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);

               R_pairs[0] = R1;
               R_pairs[1] = R3;
               Rij_pairs[0][0] =  Rij[0];
               Rij_pairs[0][1] =  Rij[1];
               Rij_pairs[0][2] =  Rij[2];
               Rij_pairs[1][0] =  Rjk[0];
               Rij_pairs[1][1] =  Rjk[1];
               Rij_pairs[1][2] =  Rjk[2];
               i_pairs[0]  = i;
               j_pairs[0]  = j;
               i_pairs[1]  = j;
               j_pairs[1]  = k;

               ier = KIM_API_process_d2Edr2(km, d2Eidr_three+4, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);

               R_pairs[0] = R3;
               R_pairs[1] = R1;
               Rij_pairs[0][0] =  Rjk[0];
               Rij_pairs[0][1] =  Rjk[1];
               Rij_pairs[0][2] =  Rjk[2];
               Rij_pairs[1][0] =  Rij[0];
               Rij_pairs[1][1] =  Rij[1];
               Rij_pairs[1][2] =  Rij[2];
               i_pairs[0]  = j;
               j_pairs[0]  = k;
               i_pairs[1]  = i;
               j_pairs[1]  = j;

               ier = KIM_API_process_d2Edr2(km, d2Eidr_three+4, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);

               R_pairs[0] = R2;
               R_pairs[1] = R3;
               Rij_pairs[0][0] =  Rik[0];
               Rij_pairs[0][1] =  Rik[1];
               Rij_pairs[0][2] =  Rik[2];
               Rij_pairs[1][0] =  Rjk[0];
               Rij_pairs[1][1] =  Rjk[1];
               Rij_pairs[1][2] =  Rjk[2];
               i_pairs[0]  = i;
               j_pairs[0]  = k;
               i_pairs[1]  = j;
               j_pairs[1]  = k;

               ier = KIM_API_process_d2Edr2(km, d2Eidr_three+5, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);

               R_pairs[0] = R3;
               R_pairs[1] = R2;
               Rij_pairs[0][0] =  Rjk[0];
               Rij_pairs[0][1] =  Rjk[1];
               Rij_pairs[0][2] =  Rjk[2];
               Rij_pairs[1][0] =  Rik[0];
               Rij_pairs[1][1] =  Rik[1];
               Rij_pairs[1][2] =  Rik[2];
               i_pairs[0]  = j;
               j_pairs[0]  = k;
               i_pairs[1]  = i;
               j_pairs[1]  = k;

               ier = KIM_API_process_d2Edr2(km, d2Eidr_three+5, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
            }
            /* contribution to forces */
            if (comp_force)
            {
               for (kdim = 0; kdim < DIM; ++kdim)
               {
                  force[i*DIM + kdim] += dEidr_three[0]*Rij[kdim]/R1 + dEidr_three[1]*Rik[kdim]/R2; /* accumulate force on atom i */
                  force[j*DIM + kdim] -= dEidr_three[0]*Rij[kdim]/R1 - dEidr_three[2]*Rjk[kdim]/R3;
                  force[k*DIM + kdim] -= dEidr_three[2]*Rjk[kdim]/R3 + dEidr_three[1]*Rik[kdim]/R2;
               }
            }
         } /* loop on kk */

      /* End adding three body terms */
      /*******************************/

      } /* loop on jj */

   }    /* infinite while loop terminated by break statements above */

   /* perform final tasks */

   if (comp_particleEnergy && comp_energy)
   {
      *energy = 0.0;
      for (k = 0; k < *nAtoms; ++k)
      {
         *energy += particleEnergy[k];
      }
   }

   /* Free temporary storage */
   if (3 == NBC)
   {
      free(neighListOfCurrentAtom);
   }

   free(dphi_three);
   free(d2phi_three);
   free(dEidr_three);
   free(d2Eidr_three);

   /* everything is great */
   ier = KIM_STATUS_OK;

   return ier;
}

/* Initialization function */
int model_driver_init(void *km, char* paramfile_names, int* nmstrlen, int* numparamfiles)
{
   /* KIM variables */
   intptr_t* pkim = *((intptr_t**) km);
   char* paramfile1name;

   /* Local variables */
   FILE* fid;
   double cutoff;
   double A;
   double B;
   double p;
   double q;
   double a;
   double lambda;
   double gamma;
   double sigma;
   double epsilon;
   double costhetat;
   /* double <FILL parameter 1>; */
   /* double <FILL parameter 2>; */
   /* FILL as many parameters as needed */
   double* model_Pcutoff;
   double* model_cutoff;
   double* model_cutsq;
   double* model_A;
   double* model_B;
   double* model_p;
   double* model_q;
   double* model_a;
   double* model_lambda;
   double* model_gamma;
   double* model_sigma;
   double* model_epsilon;
   double* model_costhetat;
   /* double* model_<FILL parameter 1>; */
   /* double* model_<FILL parameter 2>; */
   /* FILL as many parameters as needed */
   int ier;
   struct model_buffer* buffer;
   const char* NBCstr;

   /* set paramfile1name */
   if (*numparamfiles != 1)
   {
       ier = KIM_STATUS_FAIL;
       KIM_API_report_error(__LINE__, __FILE__, "Incorrect number of parameter files.", ier);
       return ier;
   }
   paramfile1name = paramfile_names;

   /* store pointer to functions in KIM object */
   KIM_API_setm_method(pkim, &ier, 3*4,
                     "compute", 1, &compute, 1,
                     "reinit",  1, &reinit,  1,
                     "destroy", 1, &destroy, 1);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_setm_method", ier);
      return ier;
   }

   /* Read in model parameters from parameter file */
   fid = fopen(paramfile1name, "r");
   if (fid == NULL)
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "Unable to open parameter file for Three_Body_Stillinger_Weber parameters", ier);
      return ier;
   }

   ier = fscanf(fid, "%lf\n %lf\n %lf\n %lf\n %lf\n %lf\n %lf\n %lf\n %lf\n %lf\n",
                        &A,
                        &B,
                        &p,
                        &q,
                        &a,
                        &lambda,
                        &gamma,
                        &sigma,
                        &epsilon,
                        &costhetat
                        /* &<FILL parameter 1>, */
                        /* &<FILL parameter 2>, */
                        /* FILL as many parameters as needed */
                        );
   fclose(fid);

   /* check that we read the right number of parameters */
   if (10 != ier)
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "Unable to read all SW parameters", ier);
      return ier;
   }

   /* convert to appropriate units */
   A *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "fs",
                                               0.0, 0.0,  0.0, 0.0, 0.0, &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
      return ier;
   }

   B *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "fs",
                                               0.0, 0.0,  0.0, 0.0, 0.0, &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
      return ier;
   }

   p *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "fs",
                                               0.0, 0.0,  0.0, 0.0, 0.0, &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
      return ier;
   }

   q *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "fs",
                                               0.0, 0.0,  0.0, 0.0, 0.0, &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
      return ier;
   }

   a *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "fs",
                                               0.0, 0.0,  0.0, 0.0, 0.0, &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
      return ier;
   }

   lambda *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "fs",
                                               0.0, 0.0,  0.0, 0.0, 0.0, &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
      return ier;
   }

   gamma *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "fs",
                                               0.0, 0.0,  0.0, 0.0, 0.0, &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
      return ier;
   }

   sigma *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "fs",
                                               1.0, 0.0,  0.0, 0.0, 0.0, &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
      return ier;
   }

   epsilon *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "fs",
                                               0.0, 1.0,  0.0, 0.0, 0.0, &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
      return ier;
   }

   costhetat *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "fs",
                                               0.0, 0.0,  0.0, 0.0, 0.0, &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
      return ier;
   }
   /* FILL as many parameters as necessary */

   /* store model cutoff in KIM object */
   model_cutoff = (double*) KIM_API_get_data(pkim, "cutoff", &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_data", ier);
      return ier;
   }
   *model_cutoff = a*sigma;

   /* allocate memory for parameters */
   model_Pcutoff = (double*) malloc(1*sizeof(double));
   if (NULL == model_Pcutoff)
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "malloc", ier);
      return ier;
   }
   model_cutsq = (double*) malloc(1*sizeof(double));
   if (NULL == model_cutsq)
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "malloc", ier);
      return ier;
   }
   model_A = (double*) malloc(1*sizeof(double));
   if (NULL == model_A)
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "malloc", ier);
      return ier;
   }
   model_B = (double*) malloc(1*sizeof(double));
   if (NULL == model_B)
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "malloc", ier);
      return ier;
   }
   model_p = (double*) malloc(1*sizeof(double));
   if (NULL == model_p)
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "malloc", ier);
      return ier;
   }
   model_q = (double*) malloc(1*sizeof(double));
   if (NULL == model_q)
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "malloc", ier);
      return ier;
   }
   model_a = (double*) malloc(1*sizeof(double));
   if (NULL == model_a)
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "malloc", ier);
      return ier;
   }
   model_lambda = (double*) malloc(1*sizeof(double));
   if (NULL == model_lambda)
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "malloc", ier);
      return ier;
   }
   model_gamma = (double*) malloc(1*sizeof(double));
   if (NULL == model_gamma)
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "malloc", ier);
      return ier;
   }
   model_sigma = (double*) malloc(1*sizeof(double));
   if (NULL == model_sigma)
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "malloc", ier);
      return ier;
   }
   model_epsilon = (double*) malloc(1*sizeof(double));
   if (NULL == model_epsilon)
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "malloc", ier);
      return ier;
   }
   model_costhetat = (double*) malloc(1*sizeof(double));
   if (NULL == model_costhetat)
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "malloc", ier);
      return ier;
   }
   /* FILL: repeat above statements as many times as necessary for all parameters.
      Use "FREE" and "FIXED" as appropriate. (Recall FREE parameters can be modified by
      the calling routine. FIXED parameters depend on the FREE parameters and must be
      appropriately adjusted in the reinit() method.)  */

   /* store parameters in KIM object */
   KIM_API_setm_data(pkim, &ier, 12*4,
                             "PARAM_FREE_cutoff",   1,  model_Pcutoff,              1,
                             "PARAM_FIXED_cutsq",   1,  model_cutsq,                1,
                             "PARAM_FREE_A",        1,  model_A,                    1,
                             "PARAM_FREE_B",        1,  model_B,                    1,
                             "PARAM_FREE_p",        1,  model_p,                    1,
                             "PARAM_FREE_q",        1,  model_q,                    1,
                             "PARAM_FREE_a",        1,  model_a,                    1,
                             "PARAM_FREE_lambda",   1,  model_lambda,               1,
                             "PARAM_FREE_gamma",    1,  model_gamma,                1,
                             "PARAM_FREE_sigma",    1,  model_sigma,                1,
                             "PARAM_FREE_epsilon",  1,  model_epsilon,              1,
                             "PARAM_FREE_costhetat", 1,  model_costhetat,           1
                             /* FILL as many parameters as needed */
                            );

   /* set value of parameters */
   *model_Pcutoff = *model_cutoff;
   *model_cutsq = (*model_cutoff)*(*model_cutoff);
   *model_A = A;
   *model_B = B;
   *model_p = p;
   *model_q = q;
   *model_a = a;
   *model_lambda = lambda;
   *model_gamma = gamma;
   *model_sigma = sigma;
   *model_epsilon = epsilon;
   *model_costhetat = costhetat;
   /* *model_<FILL parameter 1> = <FILL parameter 1>; */
   /* *model_<FILL parameter 2> = <FILL parameter 2>; */
   /* FILL as many parameters as needed */

   /* allocate buffer */
   buffer = (struct model_buffer*) malloc(sizeof(struct model_buffer));
   if (NULL == buffer)
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "malloc", ier);
      return ier;
   }

   /* setup buffer */
   /* Determine neighbor list boundary condition (NBC) */
   ier = KIM_API_get_NBC_method(pkim, &NBCstr);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_NBC_method", ier);
      return ier;
   }
   if ((!strcmp("NEIGH_RVEC_H",NBCstr)) || (!strcmp("NEIGH_RVEC_F",NBCstr)))
   {
      buffer->NBC = 0;
   }
   else if ((!strcmp("NEIGH_PURE_H",NBCstr)) || (!strcmp("NEIGH_PURE_F",NBCstr)))
   {
      buffer->NBC = 1;
   }
   else if ((!strcmp("MI_OPBC_H",NBCstr)) || (!strcmp("MI_OPBC_F",NBCstr)))
   {
      buffer->NBC = 2;
   }
   else if (!strcmp("CLUSTER",NBCstr))
   {
      buffer->NBC = 3;
   }
   else
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "Unknown NBC method", ier);
      return ier;
   }

   /* Determine if Half or Full neighbor lists are being used */
   /*****************************
    * HalfOrFull = 1 -- Half
    *            = 2 -- Full
    *****************************/

   if (KIM_API_is_half_neighbors(pkim, &ier))
   {
      buffer->HalfOrFull = 1;
   }
   else
   {
     buffer->HalfOrFull = 2;
   }


   /* determine neighbor list handling mode */
   if (buffer->NBC != 3)
   {
      /*****************************
       * IterOrLoca = 1 -- Iterator
       *            = 2 -- Locator
       *****************************/
      buffer->IterOrLoca = KIM_API_get_neigh_mode(pkim, &ier);
      if (KIM_STATUS_OK > ier)
      {
         KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh_mode", ier);
         return ier;
      }
      if ((buffer->IterOrLoca != 1) && (buffer->IterOrLoca != 2))
      {
         printf("* ERROR: Unsupported IterOrLoca mode = %i\n", buffer->IterOrLoca);
         return ier;
      }
   }
   else
   {
      buffer->IterOrLoca = 2;   /* for CLUSTER NBC */
   }

   buffer->model_index_shift = KIM_API_get_model_index_shift(pkim);

   KIM_API_getm_index(pkim, &ier, 12*3,
                      "cutoff",                      &(buffer->cutoff_ind),                      1,
                      "numberOfParticles",           &(buffer->numberOfParticles_ind),           1,
                      "particleSpecies",               &(buffer->particleSpecies_ind),               1,
                      "numberContributingParticles", &(buffer->numberContributingParticles_ind), 1,
                      "coordinates",                 &(buffer->coordinates_ind),                 1,
                      "get_neigh",                   &(buffer->get_neigh_ind),                   1,
                      "boxSideLengths",              &(buffer->boxSideLengths_ind),              1,
                      "energy",                      &(buffer->energy_ind),                      1,
                      "forces",                      &(buffer->forces_ind),                      1,
                      "particleEnergy",              &(buffer->particleEnergy_ind),              1,
                      "process_dEdr",                &(buffer->process_dEdr_ind),                1,
                      "process_d2Edr2",              &(buffer->process_d2Edr2_ind),              1
                     );
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_index", ier);
      return ier;
   }
   /* store in model buffer */
   KIM_API_set_model_buffer(pkim, (void*) buffer, &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_set_model_buffer", ier);
      return ier;
   }

   /* store parameters in buffer */
   buffer->Pcutoff    = model_Pcutoff;
   buffer->cutsq      = model_cutsq;
   buffer->A          = model_A;
   buffer->B          = model_B;
   buffer->p          = model_p;
   buffer->q          = model_q;
   buffer->a          = model_a;
   buffer->lambda     = model_lambda;
   buffer->gamma      = model_gamma;
   buffer->sigma      = model_sigma;
   buffer->epsilon    = model_epsilon;
   buffer->costhetat  = model_costhetat;
   /* buffer-><FILL parameter 1> = <FILL parameter 1>; */
   /* buffer-><FILL parameter 2> = <FILL parameter 2>; */
   /* FILL as many parameters as needed */

   ier = KIM_STATUS_OK;
   return ier;
}

/* Reinitialization function */
static int reinit(void *km)
{
    /* Local variables */
   intptr_t* pkim = *((intptr_t**) km);
   int ier;
   double *cutoff;
   struct model_buffer* buffer;

   /* get buffer from KIM object */
   buffer = (struct model_buffer*) KIM_API_get_model_buffer(pkim, &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
      return ier;
   }

   /* set new values in KIM object     */
   /*                                  */
   /* store model cutoff in KIM object */
   cutoff = KIM_API_get_data(pkim, "cutoff", &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_data", ier);
      return ier;
   }
   *cutoff = *buffer->Pcutoff;

   /* set value of parameter cutsq */
   *buffer->cutsq = (*cutoff)*(*cutoff);

   /* FILL: store any other FIXED parameters whose values depend on FREE parameters */

   ier = KIM_STATUS_OK;
   return ier;
}

/* destroy function */
static int destroy(void *km)
{
   /* Local variables */
   intptr_t* pkim = *((intptr_t**) km);
   struct model_buffer* buffer;
   int ier;

   /* get model buffer from KIM object */
   buffer = (struct model_buffer*) KIM_API_get_model_buffer(pkim, &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
      return ier;
   }

   /* free parameters */
   free(buffer->Pcutoff);
   free(buffer->cutsq);
   free(buffer->A);
   free(buffer->B);
   free(buffer->p);
   free(buffer->q);
   free(buffer->a);
   free(buffer->lambda);
   free(buffer->gamma);
   free(buffer->sigma);
   free(buffer->epsilon);
   free(buffer->costhetat);
   /* FILL: repeat above statements as many times as necessary for all FREE and FIXED parameters. */

   /* destroy the buffer */
   free(buffer);

   ier = KIM_STATUS_OK;
   return ier;
}
