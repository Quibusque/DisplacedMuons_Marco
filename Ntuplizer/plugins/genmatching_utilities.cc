#include "genmatching_utilities.h"

#include "propagate_utilities.h"
#include "propagation_definitions.h"

bool hasMotherWithPdgId(const reco::Candidate* particle, int pdgId) {
    // Loop on mothers, if any, and return true if a mother with the specified pdgId is found
    for (size_t i = 0; i < particle->numberOfMothers(); i++) {
        const reco::Candidate* mother = particle->mother(i);
        if (mother->pdgId() == pdgId || hasMotherWithPdgId(mother, pdgId)) {
            return true;
        }
    }
    return false;
}


GenMatchResults matchRecoTrackToGenSurface(
    const PropagationSurface genSurface, const reco::Track* recoTrack, const MagneticField* magField,
    const Propagator* propagatorAlong, const Propagator* propagatorOpposite,
    const GlobalTrajectoryParameters& genFinalParams, GlobalTrajectoryParameters& recoFinalParams,
    CartesianTrajectoryError& finalRecoError, Float_t deltaR_thr) {
    if (genSurface.genMatchResult == GenMatchResults::NONE) {
        return GenMatchResults::NONE;
    }

    // Build reco initial state, including error from track covariance matrix
    GlobalPoint recoVertex(recoTrack->vx(), recoTrack->vy(), recoTrack->vz());
    GlobalVector recoMomentum(recoTrack->px(), recoTrack->py(), recoTrack->pz());
    int recoCharge = recoTrack->charge();
    FreeTrajectoryState recoFTS(recoVertex, recoMomentum, recoCharge, magField);
    recoFTS.setCurvilinearError(CurvilinearTrajectoryError(recoTrack->covariance()));
    TsosPath recoTsosPath = TsosPath();

    bool recoPropagationGood = propagateToSurface(recoFTS, recoTsosPath, genSurface, magField,
                                                  propagatorAlong, propagatorOpposite);
    if (recoPropagationGood) {
        recoFinalParams =
            GlobalTrajectoryParameters(recoTsosPath.first.globalPosition(),
                                       recoTsosPath.first.globalMomentum(), recoCharge, magField);

        Float_t deltaR_value = reco::deltaR(genFinalParams.momentum(), recoFinalParams.momentum());
        if (deltaR_value > deltaR_thr) {
            return GenMatchResults::DELTA_R_FAIL;
        }

        finalRecoError = CartesianTrajectoryError(recoTsosPath.first.cartesianError());
        return genSurface.genMatchResult;
    }

    return GenMatchResults::NONE;
}

AlgebraicVector6 calculateChi2Vector(GlobalTrajectoryParameters genParams,
                                     GlobalTrajectoryParameters recoParams,
                                     CartesianTrajectoryError recoError) {
    AlgebraicVector6 genVector = genParams.vector();
    AlgebraicVector6 recoVector = recoParams.vector();

    AlgebraicVector6 delta = genVector - recoVector;

    const AlgebraicSymMatrix66& recoCovMatrix = recoError.matrix();
    AlgebraicVector6 errors = recoCovMatrix.Diagonal();
    AlgebraicVector6 chi2_vec;
    for (int i = 0; i < 6; i++) {
        chi2_vec(i) = delta(i) * delta(i) / errors(i);
    }

    return chi2_vec/6;
}