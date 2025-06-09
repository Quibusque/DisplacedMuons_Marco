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

GenMatchResults matchRecoTrackToGenSurface(const PropagationSurface genSurface,
                                           const reco::Track* recoTrack,
                                           const MagneticField* magField,
                                           const Propagator* propagatorAlong,
                                           const Propagator* propagatorOpposite,
                                           GlobalTrajectoryParameters& recoFinalParams,
                                           CartesianTrajectoryError& finalRecoError) {
    if (static_cast<int>(genSurface.genMatchResult) < 0) {
        return genSurface.genMatchResult;
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
        chi2_vec(i) = delta(i) * delta(i) / (errors(i) * errors(i));
    }

    return chi2_vec / 6;
}

bool isGenMatch(const GlobalTrajectoryParameters& genFinalParams,
                const GlobalTrajectoryParameters& recoFinalParams,
                const PropagationSurface& surface) {
    // Proper angle difference in [0,pi] range
    Float_t deltaPhi = std::abs(recoFinalParams.momentum().phi() - genFinalParams.momentum().phi());
    deltaPhi = (deltaPhi < M_PI) ? deltaPhi : (2 * M_PI - deltaPhi);
    deltaPhi = std::abs(deltaPhi);
    switch (surface.genMatchResult) {
        case GenMatchResults::CYLINDER: {
            Float_t deltaZ =
                std::abs(recoFinalParams.position().z() - genFinalParams.position().z());
            return (deltaPhi < 0.025 && deltaZ < 150.0);
        }
        case GenMatchResults::POS_ENDCAP:
        case GenMatchResults::NEG_ENDCAP: {
            Float_t deltaR =
                std::abs(recoFinalParams.position().perp() - genFinalParams.position().perp());
            // M_PI/36 == 5 degrees
            return (deltaPhi < M_PI / 36. && deltaR < 50.0);
        }

        default:
            return false;
    }
}

Float_t genMatchDistance(const GlobalTrajectoryParameters& genFinalParams,
                         const GlobalTrajectoryParameters& recoFinalParams,
                         const PropagationSurface& surface) {
    // Proper angle difference in [0,pi] range
    Float_t deltaPhi = std::abs(recoFinalParams.momentum().phi() - genFinalParams.momentum().phi());
    deltaPhi = (deltaPhi < M_PI) ? deltaPhi : (2 * M_PI - deltaPhi);
    deltaPhi = std::abs(deltaPhi);
    switch (surface.genMatchResult) {
        case GenMatchResults::CYLINDER: {
            Float_t deltaZ =
                std::abs(recoFinalParams.position().z() - genFinalParams.position().z());
            return std::sqrt(deltaPhi * deltaPhi + deltaZ * deltaZ);
        }
        case GenMatchResults::POS_ENDCAP:
        case GenMatchResults::NEG_ENDCAP: {
            Float_t deltaR =
                std::abs(recoFinalParams.position().perp() - genFinalParams.position().perp());
            return std::sqrt(deltaPhi * deltaPhi + deltaR * deltaR);
        }
        default:
            return std::numeric_limits<float>::infinity();
    }
}
