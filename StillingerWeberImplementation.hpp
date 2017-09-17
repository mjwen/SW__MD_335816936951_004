//
// CDDL HEADER START
//
// The contents of this file are subject to the terms of the Common Development
// and Distribution License Version 1.0 (the "License").
//
// You can obtain a copy of the license at
// http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
// specific language governing permissions and limitations under the License.
//
// When distributing Covered Code, include this CDDL HEADER in each file and
// include the License file in a prominent location with the name LICENSE.CDDL.
// If applicable, add the following below this CDDL HEADER, with the fields
// enclosed by brackets "[]" replaced with your own identifying information:
//
// Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
//
// CDDL HEADER END
//

//
// Copyright (c) 2017, Regents of the University of Minnesota.
// All rights reserved.
//
// Contributors:
//    Mingjian Wen
//


#ifndef STILLINGER_WEBER_IMPLEMENTATION_HPP_
#define STILLINGER_WEBER_IMPLEMENTATION_HPP_

#include <vector>
#include <unordered_map>
#include "KIM_LogVerbosity.hpp"
#include "StillingerWeber.hpp"

#define DIMENSION 3
#define ONE 1.0
#define HALF 0.5

#define MAX_PARAMETER_FILES 1


//==============================================================================
//
// Type definitions, enumerations, and helper function prototypes
//
//==============================================================================

// type declaration for get neighbor functions
typedef int (GetNeighborFunction)(void const * const, int const,
                                  int * const, int const ** const);

// type declaration for vector of constant dimension
typedef double VectorOfSizeDIM[DIMENSION];

// helper routine declarations
void AllocateAndInitialize2DArray(double**& arrayPtr, int const extentZero,
                                  int const extentOne);
void Deallocate2DArray(double**& arrayPtr);

//==============================================================================
//
// Declaration of StillingerWeberImplementation class
//
//==============================================================================

//******************************************************************************
class StillingerWeberImplementation
{
 public:
  StillingerWeberImplementation(
      KIM::ModelDriverCreate * const modelDriverCreate,
      KIM::LengthUnit const requestedLengthUnit,
      KIM::EnergyUnit const requestedEnergyUnit,
      KIM::ChargeUnit const requestedChargeUnit,
      KIM::TemperatureUnit const requestedTemperatureUnit,
      KIM::TimeUnit const requestedTimeUnit,
      int * const ier);
  ~StillingerWeberImplementation();  // no explicit Destroy() needed here

  int Refresh(KIM::ModelRefresh * const modelRefresh);
  int Compute(KIM::ModelCompute const * const modelCompute);

 private:
  // Constant values that never change
  //   Set in constructor (via SetConstantValues)
  //
  //
  // StillingerWeberImplementation: constants
  int numberModelSpecies_;
  std::vector<int> modelSpeciesCodeList_;
  int numberUniqueSpeciesPairs_;


  // Constant values that are read from the input files and never change
  //   Set in constructor (via functions listed below)
  //
  //
  // KIM API: Model Fixed Parameters
  //   Memory allocated in   AllocateFixedParameterMemory()
  //   Memory deallocated in destructor
  //   Data set in ReadParameterFile routines
  // none
  //
  // KIM API: Model Free Parameters whose (pointer) values never change
  //   Memory allocated in   AllocateFreeParameterMemory() (from constructor)
  //   Memory deallocated in destructor
  //   Data set in ReadParameterFile routines OR by KIM Simulator
  double* cutoff_;
  double* sigma_;


  // Mutable values that only change when Refresh() executes
  //   Set in Refresh (via SetRefreshMutableValues)
  //
  //
  // KIM API: Model Fixed Parameters
  // none
  //
  // KIM API: Model Free Parameters
  // none
  //
  // StillingerWeberImplementation: values
  double influenceDistance_;
  double** cutoffSq_2D_;
  double** sigma_2D_;


  // Mutable values that can change with each call to Refresh() and Compute()
  //   Memory may be reallocated on each call
  //
  //
  // StillingerWeberImplementation: values that change
  int cachedNumberOfParticles_;


  // Helper methods
  //
  //
  // Related to constructor
  void AllocateFreeParameterMemory();
  static int OpenParameterFiles(
      KIM::ModelDriverCreate * const modelDriverCreate,
      int const numberParameterFiles,
      FILE* parameterFilePointers[MAX_PARAMETER_FILES]);
  static void CloseParameterFiles(
      int const numberParameterFiles,
      FILE* const parameterFilePointers[MAX_PARAMETER_FILES]);
  int ProcessParameterFiles(
      KIM::ModelDriverCreate * const modelDriverCreate,
      int const numberParameterFiles,
      FILE* const parameterFilePointers[MAX_PARAMETER_FILES]);
  void getNextDataLine(FILE* const filePtr, char* const nextLine,
                       int const maxSize, int* endOfFileFlag);
  int ConvertUnits(
      KIM::ModelDriverCreate * const modelDriverCreate,
      KIM::LengthUnit const requestedLengthUnit,
      KIM::EnergyUnit const requestedEnergyUnit,
      KIM::ChargeUnit const requestedChargeUnit,
      KIM::TemperatureUnit const requestedTemperatureUnit,
      KIM::TimeUnit const requestedTimeUnit);
  int RegisterKIMModelSettings(KIM::ModelDriverCreate * const modelDriverCreate);
  int RegisterKIMParameters(KIM::ModelDriverCreate * const modelDriverCreate);
  int RegisterKIMFunctions(KIM::ModelDriverCreate * const modelDriverCreate) const;
  //
  // Related to Refresh()
  template<class ModelObj>
  int SetRefreshMutableValues(ModelObj * const modelObj);
  //
  // Related to Compute()
  int SetComputeMutableValues(KIM::ModelCompute const * const modelCompute,
                              bool& isComputeProcess_dEdr,
                              bool& isComputeProcess_d2Edr2,
                              bool& isComputeEnergy,
                              bool& isComputeForces,
                              bool& isComputeParticleEnergy,
                              int const*& particleSpecies,
                              int const*& particleContributing,
                              VectorOfSizeDIM const*& coordinates,
                              double*& energy,
                              double*& particleEnergy,
                              VectorOfSizeDIM*& forces);
  int CheckParticleSpecies(KIM::ModelCompute const * const modelCompute,
                           int const* const particleSpecies) const;
  int GetComputeIndex(const bool& isComputeProcess_dEdr,
                      const bool& isComputeProcess_d2Edr2,
                      const bool& isComputeEnergy,
                      const bool& isComputeForces,
                      const bool& isComputeParticleEnergy) const;

  // compute functions
  template< bool isComputeProcess_dEdr, bool isComputeProcess_d2Edr2,
            bool isComputeEnergy, bool isComputeForces,
            bool isComputeParticleEnergy>
  int Compute(KIM::ModelCompute const * const modelCompute,
              const int* const particleSpecies,
              const int* const particleContributing,
              const VectorOfSizeDIM* const coordinates,
              double* const energy,
              VectorOfSizeDIM* const forces,
              double* const particleEnergy);
};

//==============================================================================
//
// Definition of StillingerWeberImplementation::Compute functions
//
// NOTE: Here we rely on the compiler optimizations to prune dead code
//       after the template expansions.  This provides high efficiency
//       and easy maintenance.
//
//==============================================================================

#include "KIM_ModelComputeLogMacros.hpp"
template< bool isComputeProcess_dEdr, bool isComputeProcess_d2Edr2,
          bool isComputeEnergy, bool isComputeForces,
          bool isComputeParticleEnergy>
int StillingerWeberImplementation::Compute(
    KIM::ModelCompute const * const modelCompute,
    const int* const particleSpecies,
    const int* const particleContributing,
    const VectorOfSizeDIM* const coordinates,
    double* const energy,
    VectorOfSizeDIM* const forces,
    double* const particleEnergy)
{
  int ier = false;

  if ((isComputeEnergy == false) &&
      (isComputeParticleEnergy == false) &&
      (isComputeForces == false) &&
      (isComputeProcess_dEdr == false) &&
      (isComputeProcess_d2Edr2 == false))
    return ier;

  // initialize energy and forces
  if (isComputeEnergy == true) {
    *energy = 0.0;
  }

  if (isComputeParticleEnergy == true) {
    for (int i = 0; i < cachedNumberOfParticles_; ++i) {
      particleEnergy[i] = 0.0;
    }
  }

  if (isComputeForces == true) {
    for (int i = 0; i < cachedNumberOfParticles_; ++i) {
      for (int j = 0; j < DIMENSION; ++j)
        forces[i][j] = 0.0;
    }
  }

  // calculate contribution from pair function
  //
  // Setup loop over contributing particles
  int i = 0;
  int numnei = 0;
  int const * n1atom = 0;

  for (i = 0; i < cachedNumberOfParticles_; ++i) {
    if (particleContributing[i]) {
      modelCompute->GetNeighborList(0, i, &numnei, &n1atom);
      int const iSpecies = particleSpecies[i];

      // Setup loop over neighbors of current particle
      for (int jj = 0; jj < numnei; ++jj) {
        int const j = n1atom[jj];
        int const jSpecies = particleSpecies[j];
        double rij[DIMENSION];

        // Compute rij
        for (int k = 0; k < DIMENSION; ++k)
          rij[k] = coordinates[j][k] - coordinates[i][k];

        // compute distance squared
        double const rij_sq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];

        if (rij_sq <= cutoffSq_2D_[iSpecies][jSpecies]) {
          // compute contribution to energy, force, etc.
          double phi = 0.0;
          double dphiByR = 0.0;
          double d2phi = 0.0;
          double dEidrByR = 0.0;
          double d2Eidr2 = 0.0;

          // Compute pair potential and its derivatives
          if (isComputeProcess_d2Edr2 == true) {
            // Compute d2phi
            d2phi = 0.0;
            d2Eidr2 = 0.5*d2phi;
          }

          if ((isComputeProcess_dEdr == true) || (isComputeForces == true)) {
            // Compute dphi
            dphiByR = 0.0;
            dEidrByR = 0.5*dphiByR;
          }

          if ((isComputeEnergy == true) || (isComputeParticleEnergy == true)) {
            // Compute phi
            phi = 0.0;
          }

          // Contribution to energy
          if (isComputeEnergy == true) {
            *energy += 0.5*phi;
          }

          // Contribution to particleEnergy
          if (isComputeParticleEnergy == true) {
            double const halfPhi = 0.5*phi;
            particleEnergy[i] += halfPhi;
          }

          // Contribution to forces
          if (isComputeForces == true) {
            for (int k = 0; k < DIMENSION; ++k) {
              double const contrib = dEidrByR * rij[k];
              forces[i][k] += contrib;
              forces[j][k] -= contrib;
            }
          }

          // Call process_dEdr
          if (isComputeProcess_dEdr == true) {
            double const rij_mag = sqrt(rij_sq);
            double const dEidr = dEidrByR*rij_mag;
            ier = modelCompute->ProcessDEDrTerm(dEidr, rij_mag, rij, i, j);
            if (ier) {
              LOG_ERROR("process_dEdr");
              return ier;
            }
          }

          // Call process_d2Edr2
          if (isComputeProcess_d2Edr2 == true) {
            double const rij_mag = sqrt(rij_sq);
            double const R_pairs[2] = {rij_mag, rij_mag};
            double const* const pRs = &R_pairs[0];
            double const Rij_pairs[6]
                = {rij[0], rij[1], rij[2],
                   rij[0], rij[1], rij[2]};
            double const* const pRijConsts = &Rij_pairs[0];
            int const i_pairs[2] = {i, i};
            int const j_pairs[2] = {j, j};
            int const* const pis = &i_pairs[0];
            int const* const pjs = &j_pairs[0];

            ier = modelCompute->ProcessD2EDr2Term(d2Eidr2, pRs, pRijConsts, pis,
                                                  pjs);
            if (ier) {
              LOG_ERROR("process_d2Edr2");
              return ier;
            }
          }
        }  // if particleContributing
      }  // if particles i and j interact
    }  // end of first neighbor loop
  }  // end of loop over contributing particles

  // everything is good
  ier = false;
  return ier;
}

#endif  // STILLINGER_WEBER_IMPLEMENTATION_HPP_
