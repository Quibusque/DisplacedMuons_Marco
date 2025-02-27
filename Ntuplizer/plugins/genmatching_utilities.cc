#include "genmatching_utilities.h"
#include "propagate_utilities.h"

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

// FOR REFERENCE, > 0 are valid surfaces
// enum kPropagationSurface {
//     RECO_OUTSIDE_CMS = -4, GEN_OUTSIDE_CMS = -3, OUTSIDE_DELTAR = -2,
//     NONE = -1, CYLINDER = 1, MAX_Z = 2, MIN_Z = 3
// }
kPropagationSurface matchGenParticleToRecoTrack(
    const reco::GenParticle& genParticle, const reco::Track* recoTrack,
    const MagneticField* magField, const Propagator* propagatorAlong,
    const Propagator* propagatorOpposite, GlobalTrajectoryParameters& genFinalParams,
    GlobalTrajectoryParameters& recoFinalParams, CartesianTrajectoryError& finalRecoError) {
    // Very loose initial cut based on deltaR between gen and reco
    Float_t deltaRValue = reco::deltaR(genParticle, *recoTrack);
    if (!(deltaRValue < 0.5)) {
        return OUTSIDE_DELTAR;
    }

    // Initial states
    // Gen vertices are in mm, while reco vertices are in cm, make gen vertices in cm
    GlobalPoint genVertex(genParticle.vx() * 0.1, genParticle.vy() * 0.1, genParticle.vz() * 0.1);
    GlobalVector genMomentum(genParticle.px(), genParticle.py(), genParticle.pz());
    int genCharge = genParticle.charge();
    FreeTrajectoryState genFTS(genVertex, genMomentum, genCharge, magField);
    GlobalPoint recoVertex(recoTrack->vx(), recoTrack->vy(), recoTrack->vz());
    GlobalVector recoMomentum(recoTrack->px(), recoTrack->py(), recoTrack->pz());
    int recoCharge = recoTrack->charge();
    FreeTrajectoryState recoFTS(recoVertex, recoMomentum, recoCharge, magField);

    // recoFTS error from track covariance matrix
    CovarianceMatrix initialRecoCovMatrix = recoTrack->covariance();
    CurvilinearTrajectoryError initialRecoError(initialRecoCovMatrix);
    recoFTS.setCurvilinearError(initialRecoError);

    // Propagate the particles to the target surface
    TsosPath genTsosPath, recoTsosPath;
    kPropagationSurface propagationSurface =
        propagateToCommonSurface(genFTS, recoFTS, genTsosPath, recoTsosPath, magField,
                                 propagatorAlong, propagatorOpposite);

    // This is necessary not to crash because of invalid access to the TsosPath
    if (propagationSurface > 0) {
        finalRecoError = CartesianTrajectoryError(recoTsosPath.first.cartesianError());

        // recoTsosPath.first
        genFinalParams =
            GlobalTrajectoryParameters(genTsosPath.first.globalPosition(),
                                       genTsosPath.first.globalMomentum(), genCharge, magField);
        recoFinalParams =
            GlobalTrajectoryParameters(recoTsosPath.first.globalPosition(),
                                       recoTsosPath.first.globalMomentum(), recoCharge, magField);
    }

    return propagationSurface;
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

    return chi2_vec;
}