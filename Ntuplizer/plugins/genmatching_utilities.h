#ifndef GENMATCHING_UTILITIES_H
#define GENMATCHING_UTILITIES_H

#include <limits>

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TMatrixF.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TrajectoryParametrization/interface/CartesianTrajectoryError.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "propagate_utilities.h"
#include "propagation_definitions.h"

/// 5 parameter covariance matrix
typedef math::Error<5>::type CovarianceMatrix;

bool hasMotherWithPdgId(const reco::Candidate* particle, int pdgId);

typedef std::pair<TrajectoryStateOnSurface, double> TsosPath;
/**
 * @brief matrix contains the chi2-like values. boolMatrix is filled to represent the
 * best values while avoiding duplicates.
 *
 * boolMatrix is filled with 0.0, 1.0s are used to mark the best minimal values.
 * There can be at most one true value per row and one true value per column. In case
 * of duplicates either in a column or in a row, the one with the lowest value is kept
 * and the others are set to 0.0.
 *
 * @param matrix the input matrix with chi2-like values
 * @param boolMatrix the output boolean matrix with 1.0 for the minimal values (modified by this
 * function)
 */
void markUniqueBestMatches(const TMatrixF& matrix, TMatrixF& boolMatrix);

/**
 * @brief Matches a reconstructed track to a generated surface.
 *
 * Try to propagate the recoTrack to the genSurface. If the propagation is successful
 * and the deltaR between the gen and reco final states is below deltaR_thr, the
 * matching is successful.
 *
 * @param genSurface The generated surface to match to.
 * @param recoTrack The reconstructed track to be matched.
 * @param magField The magnetic field used for propagation.
 * @param propagatorAlong The propagator for along direction.
 * @param propagatorOpposite The propagator for opposite direction.
 * @param genFinalParams The final parameters of the generated trajectory.
 * @param recoFinalParams The final parameters of the reconstructed trajectory (modified by this
 * function).
 * @param finalRecoError The error of the reconstructed trajectory (modified by this function).
 * @param deltaR_thr reco and gen deltaR must be below this to be considered a match.
 * @return The matching result as a GenMatchResults enum value.
 */
GenMatchResults matchRecoTrackToGenSurface(
    const PropagationSurface genSurface, const reco::Track* recoTrack,
    const MagneticField* magField, const Propagator* propagatorAlong,
    const Propagator* propagatorOpposite, const GlobalTrajectoryParameters& genFinalParams,
    GlobalTrajectoryParameters& recoFinalParams, CartesianTrajectoryError& finalRecoError,
    Float_t deltaR_thr);
/**
 * @brief Calculates the chi-square between reco and gen given reco errors
 *
 * The chi-square is computed as (gen-reco)/(reco_error*reco_error) where the
 * elements are 3-component position and 3-component momentum. The final result is
 * reduced chi-square (chi-square divided by 6).
 */
AlgebraicVector6 calculateChi2Vector(GlobalTrajectoryParameters genParams,
                                     GlobalTrajectoryParameters recoParams,
                                     CartesianTrajectoryError recoError);

#endif