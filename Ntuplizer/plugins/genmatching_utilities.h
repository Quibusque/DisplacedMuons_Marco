#ifndef GENMATCHING_UTILITIES_H
#define GENMATCHING_UTILITIES_H

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TrajectoryParametrization/interface/CartesianTrajectoryError.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "propagate_utilities.h"

/// 5 parameter covariance matrix
typedef math::Error<5>::type CovarianceMatrix;

bool hasMotherWithPdgId(const reco::Candidate* particle, int pdgId);

kPropagationSurface matchGenParticleToRecoTrack(
    const reco::GenParticle& genParticle, const reco::Track* recoTrack,
    const MagneticField* magField, const Propagator* propagatorAlong,
    const Propagator* propagatorOpposite, GlobalTrajectoryParameters& genFinalParams,
    GlobalTrajectoryParameters& recoFinalParams, CartesianTrajectoryError& finalRecoError);

AlgebraicVector6 calculateChi2Vector(GlobalTrajectoryParameters genParams,
                                     GlobalTrajectoryParameters recoParams,
                                     CartesianTrajectoryError recoError);

#endif