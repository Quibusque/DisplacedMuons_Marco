#ifndef PROPAGATION_DEFINITIONS_H
#define PROPAGATION_DEFINITIONS_H

enum class GenMatchResults {
    GEN_OUTSIDE_CMS = -3,
    DELTA_R_FAIL = -2,
    NONE = -1,
    CYLINDER = 1,
    POS_ENDCAP = 2,
    NEG_ENDCAP = 3
};

struct PropagationSurface {
    Float_t radius;
    Float_t maxZ;
    Float_t minZ;
    GenMatchResults genMatchResult;

    bool operator==(const PropagationSurface& other) const {
        return radius == other.radius && maxZ == other.maxZ && minZ == other.minZ &&
               genMatchResult == other.genMatchResult;
    }
    bool operator!=(const PropagationSurface& other) const { return !(*this == other); }
};

namespace PropagationConstants {
const PropagationSurface CYLINDER = {420.0, 700.0, -700.0, GenMatchResults::CYLINDER};
const PropagationSurface POS_ENDCAP = {700.0, 700.0, 0.0, GenMatchResults::POS_ENDCAP};
const PropagationSurface NEG_ENDCAP = {700.0, 0.0, -700.0, GenMatchResults::NEG_ENDCAP};
const PropagationSurface NONE = {0.0, 0.0, 0.0, GenMatchResults::NONE};
const PropagationSurface GEN_OUTSIDE_CMS = {0.0, 0.0, 0.0, GenMatchResults::GEN_OUTSIDE_CMS};
const Float_t MAX_CMS_CYLINDER_RADIUS = 800.0;
const Float_t MAX_CMS_Z = 1200.0;
const Float_t MIN_CMS_Z = -1200.0;
}  // namespace PropagationConstants

#endif  // PROPAGATION_DEFINITIONS_H