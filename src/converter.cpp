#include "converter.h"

// Include all converter model headers
#include "converter_models/Flyback.h"
#include "converter_models/Buck.h"
#include "converter_models/Boost.h"
#include "converter_models/SingleSwitchForward.h"
#include "converter_models/TwoSwitchForward.h"
#include "converter_models/ActiveClampForward.h"
#include "converter_models/PushPull.h"
#include "converter_models/Llc.h"
#include "converter_models/Cllc.h"
#include "converter_models/Dab.h"
#include "converter_models/PhaseShiftedFullBridge.h"
#include "converter_models/PhaseShiftedHalfBridge.h"
#include "converter_models/IsolatedBuck.h"
#include "converter_models/IsolatedBuckBoost.h"
#include "converter_models/CurrentTransformer.h"
#include "converter_models/PowerFactorCorrection.h"
#include "converter_models/Cuk.h"
#include "converter_models/Sepic.h"
#include "converter_models/Zeta.h"
#include "converter_models/FourSwitchBuckBoost.h"
#include "converter_models/AsymmetricHalfBridge.h"
#include "converter_models/Weinberg.h"
#include "converter_models/Vienna.h"
#include "converter_models/Clllc.h"
#include "converter_models/Src.h"
#include "converter_models/CommonModeChoke.h"
#include "converter_models/DifferentialModeChoke.h"
#include "constructive_models/MasMigration.h"

#include <set>
#include <unordered_map>
#include <cmath>

namespace PyMKF {

// ------------------------------------------------------------------
// Topology name normalization.
//
// dispatch_converter() below uses legacy short forms ("flyback", "buck",
// ...). External callers may pass:
//   * the same short form (canonical internal),
//   * the MAS 1.0 camelCase enum value ("flybackConverter", ...),
//   * or the pre-1.0 Title Case label ("Flyback Converter", ...).
//
// Normalize all three to the internal short form so process_converter()
// and design_magnetics_from_converter() transparently accept any of them.
// ------------------------------------------------------------------
static std::string normalize_topology_name(const std::string& s) {
    static const std::set<std::string> short_forms = {
        "flyback", "advanced_flyback",
        "buck", "advanced_buck",
        "boost", "advanced_boost",
        "single_switch_forward",
        "two_switch_forward",
        "active_clamp_forward",
        "push_pull",
        "llc", "advanced_llc",
        "cllc", "advanced_cllc",
        "dab", "advanced_dab",
        "phase_shifted_full_bridge", "psfb",
        "phase_shifted_half_bridge", "pshb",
        "isolated_buck",
        "isolated_buck_boost",
        "current_transformer",
        "power_factor_correction", "pfc",
        "cuk", "advanced_cuk",
        "sepic", "advanced_sepic",
        "zeta", "advanced_zeta",
        "four_switch_buck_boost", "advanced_four_switch_buck_boost",
        "asymmetric_half_bridge", "advanced_asymmetric_half_bridge",
        "weinberg", "advanced_weinberg",
        "vienna", "advanced_vienna",
        "clllc", "advanced_clllc",
        "src", "advanced_src"
    };
    if (short_forms.count(s)) return s;

    static const std::unordered_map<std::string, std::string> aliases = {
        // MAS 1.0 camelCase (canonical for designRequirements.topology)
        {"flybackConverter",                 "flyback"},
        {"buckConverter",                    "buck"},
        {"boostConverter",                   "boost"},
        {"singleSwitchForwardConverter",     "single_switch_forward"},
        {"twoSwitchForwardConverter",        "two_switch_forward"},
        {"activeClampForwardConverter",      "active_clamp_forward"},
        {"pushPullConverter",                "push_pull"},
        {"llcResonantConverter",             "llc"},
        {"llcConverter",                     "llc"},
        {"cllcResonantConverter",            "cllc"},
        {"cllcConverter",                    "cllc"},
        {"dualActiveBridgeConverter",        "dab"},
        {"dual_active_bridge",               "dab"},
        {"phaseShiftedFullBridgeConverter",  "phase_shifted_full_bridge"},
        {"phaseShiftedHalfBridgeConverter",  "phase_shifted_half_bridge"},
        {"isolatedBuckConverter",            "isolated_buck"},
        {"isolatedBuckBoostConverter",       "isolated_buck_boost"},
        {"currentTransformer",               "current_transformer"},
        {"powerFactorCorrection",            "power_factor_correction"},
        {"cukConverter",                      "cuk"},
        {"sepicConverter",                    "sepic"},
        {"zetaConverter",                     "zeta"},
        {"fourSwitchBuckBoostConverter",      "four_switch_buck_boost"},
        {"asymmetricHalfBridgeConverter",     "asymmetric_half_bridge"},
        {"weinbergConverter",                 "weinberg"},
        {"viennaRectifierConverter",          "vienna"},
        {"clllcResonantConverter",            "clllc"},
        {"seriesResonantConverter",           "src"},

        // Pre-1.0 Title Case labels
        {"Flyback Converter",                 "flyback"},
        {"Buck Converter",                    "buck"},
        {"Boost Converter",                   "boost"},
        {"Single-Switch Forward Converter",   "single_switch_forward"},
        {"Two-Switch Forward Converter",      "two_switch_forward"},
        {"Active-Clamp Forward Converter",    "active_clamp_forward"},
        {"Push-Pull Converter",               "push_pull"},
        {"LLC Resonant Converter",            "llc"},
        {"LLC Converter",                     "llc"},
        {"CLLC Resonant Converter",           "cllc"},
        {"Dual Active Bridge Converter",      "dab"},
        {"Phase-Shifted Full Bridge Converter", "phase_shifted_full_bridge"},
        {"Phase-Shifted Half Bridge Converter", "phase_shifted_half_bridge"},
        {"Isolated Buck Converter",           "isolated_buck"},
        {"Isolated Buck-Boost Converter",     "isolated_buck_boost"},
        {"Current Transformer",               "current_transformer"},
        {"Power Factor Correction",           "power_factor_correction"},
        {"Cuk Converter",                     "cuk"},
        {"SEPIC Converter",                   "sepic"},
        {"Sepic Converter",                   "sepic"},
        {"Zeta Converter",                    "zeta"},
        {"Four-Switch Buck-Boost Converter",  "four_switch_buck_boost"},
        {"Asymmetric Half-Bridge Converter",  "asymmetric_half_bridge"},
        {"Weinberg Converter",                "weinberg"},
        {"Vienna Rectifier Converter",        "vienna"},
        {"Vienna Rectifier",                  "vienna"},
        {"CLLLC Resonant Converter",          "clllc"},
        {"Series Resonant Converter",         "src"},
        {"SRC",                               "src"}
    };

    auto it = aliases.find(s);
    if (it != aliases.end()) return it->second;

    // Pass through unchanged; dispatch_converter() will throw if it is
    // truly unknown. We never silently substitute an unrelated topology.
    return s;
}

OpenMagnetics::Inputs process_flyback_internal(const json& converterJson, bool useNgspice) {
    OpenMagnetics::AdvancedFlyback converter(converterJson);
    converter._assertErrors = true;
    
    if (useNgspice) {
        // Get design parameters directly from AdvancedFlyback
        std::vector<double> turnsRatios = converter.get_desired_turns_ratios();
        double inductance = converter.get_desired_inductance();
        auto operatingPoints = converter.simulate_and_extract_operating_points(turnsRatios, inductance);
        
        // Build DesignRequirements manually like AdvancedFlyback::process() does
        OpenMagnetics::DesignRequirements designReqs;
        designReqs.get_mutable_turns_ratios().clear();
        for (auto turnsRatio : turnsRatios) {
            OpenMagnetics::DimensionWithTolerance turnsRatioWithTolerance;
            turnsRatioWithTolerance.set_nominal(turnsRatio);
            designReqs.get_mutable_turns_ratios().push_back(turnsRatioWithTolerance);
        }
        OpenMagnetics::DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_nominal(inductance);
        designReqs.set_magnetizing_inductance(inductanceWithTolerance);
        designReqs.set_topology(MAS::Topologies::FLYBACK_CONVERTER);
        
        OpenMagnetics::Inputs result;
        result.set_design_requirements(designReqs);
        result.set_operating_points(operatingPoints);
        return result;
    } else {
        OpenMagnetics::Inputs result = converter.process();
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::FLYBACK_CONVERTER);
        return result;
    }
}

OpenMagnetics::Inputs process_buck_internal(const json& converterJson, bool useNgspice) {
    OpenMagnetics::AdvancedBuck converter(converterJson);
    converter._assertErrors = true;
    
    if (useNgspice) {
        auto designReqs = converter.process_design_requirements();
        double inductance = OpenMagnetics::resolve_dimensional_values(designReqs.get_magnetizing_inductance());
        auto operatingPoints = converter.simulate_and_extract_operating_points(inductance);
        OpenMagnetics::Inputs result;
        result.set_design_requirements(designReqs);
        result.set_operating_points(operatingPoints);
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::BUCK_CONVERTER);
        return result;
    } else {
        OpenMagnetics::Inputs result = converter.process();
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::BUCK_CONVERTER);
        return result;
    }
}

OpenMagnetics::Inputs process_boost_internal(const json& converterJson, bool useNgspice) {
    OpenMagnetics::AdvancedBoost converter(converterJson);
    converter._assertErrors = true;
    
    if (useNgspice) {
        auto designReqs = converter.process_design_requirements();
        double inductance = OpenMagnetics::resolve_dimensional_values(designReqs.get_magnetizing_inductance());
        auto operatingPoints = converter.simulate_and_extract_operating_points(inductance);
        OpenMagnetics::Inputs result;
        result.set_design_requirements(designReqs);
        result.set_operating_points(operatingPoints);
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::BOOST_CONVERTER);
        return result;
    } else {
        OpenMagnetics::Inputs result = converter.process();
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::BOOST_CONVERTER);
        return result;
    }
}

OpenMagnetics::Inputs process_single_switch_forward_internal(const json& converterJson, bool useNgspice) {
    OpenMagnetics::AdvancedSingleSwitchForward converter(converterJson);
    converter._assertErrors = true;
    
    if (useNgspice) {
        auto designReqs = converter.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(OpenMagnetics::resolve_dimensional_values(tr));
        }
        double inductance = OpenMagnetics::resolve_dimensional_values(designReqs.get_magnetizing_inductance());
        auto operatingPoints = converter.simulate_and_extract_operating_points(turnsRatios, inductance);
        OpenMagnetics::Inputs result;
        result.set_design_requirements(designReqs);
        result.set_operating_points(operatingPoints);
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::SINGLE_SWITCH_FORWARD_CONVERTER);
        return result;
    } else {
        OpenMagnetics::Inputs result = converter.process();
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::SINGLE_SWITCH_FORWARD_CONVERTER);
        return result;
    }
}

OpenMagnetics::Inputs process_two_switch_forward_internal(const json& converterJson, bool useNgspice) {
    OpenMagnetics::AdvancedTwoSwitchForward converter(converterJson);
    converter._assertErrors = true;
    
    if (useNgspice) {
        auto designReqs = converter.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(OpenMagnetics::resolve_dimensional_values(tr));
        }
        double inductance = OpenMagnetics::resolve_dimensional_values(designReqs.get_magnetizing_inductance());
        auto operatingPoints = converter.simulate_and_extract_operating_points(turnsRatios, inductance);
        OpenMagnetics::Inputs result;
        result.set_design_requirements(designReqs);
        result.set_operating_points(operatingPoints);
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::TWO_SWITCH_FORWARD_CONVERTER);
        return result;
    } else {
        OpenMagnetics::Inputs result = converter.process();
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::TWO_SWITCH_FORWARD_CONVERTER);
        return result;
    }
}

OpenMagnetics::Inputs process_active_clamp_forward_internal(const json& converterJson, bool useNgspice) {
    OpenMagnetics::AdvancedActiveClampForward converter(converterJson);
    converter._assertErrors = true;
    
    if (useNgspice) {
        auto designReqs = converter.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(OpenMagnetics::resolve_dimensional_values(tr));
        }
        double inductance = OpenMagnetics::resolve_dimensional_values(designReqs.get_magnetizing_inductance());
        auto operatingPoints = converter.simulate_and_extract_operating_points(turnsRatios, inductance);
        OpenMagnetics::Inputs result;
        result.set_design_requirements(designReqs);
        result.set_operating_points(operatingPoints);
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::ACTIVE_CLAMP_FORWARD_CONVERTER);
        return result;
    } else {
        OpenMagnetics::Inputs result = converter.process();
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::ACTIVE_CLAMP_FORWARD_CONVERTER);
        return result;
    }
}

OpenMagnetics::Inputs process_push_pull_internal(const json& converterJson, bool useNgspice) {
    OpenMagnetics::AdvancedPushPull converter(converterJson);
    converter._assertErrors = true;
    
    if (useNgspice) {
        auto designReqs = converter.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(OpenMagnetics::resolve_dimensional_values(tr));
        }
        double inductance = OpenMagnetics::resolve_dimensional_values(designReqs.get_magnetizing_inductance());
        auto operatingPoints = converter.simulate_and_extract_operating_points(turnsRatios, inductance);
        OpenMagnetics::Inputs result;
        result.set_design_requirements(designReqs);
        result.set_operating_points(operatingPoints);
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::PUSH_PULL_CONVERTER);
        return result;
    } else {
        OpenMagnetics::Inputs result = converter.process();
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::PUSH_PULL_CONVERTER);
        return result;
    }
}

OpenMagnetics::Inputs process_llc_internal(const json& converterJson, bool useNgspice) {
    OpenMagnetics::AdvancedLlc converter(converterJson);
    converter._assertErrors = true;
    
    if (useNgspice) {
        auto designReqs = converter.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(OpenMagnetics::resolve_dimensional_values(tr));
        }
        double inductance = OpenMagnetics::resolve_dimensional_values(designReqs.get_magnetizing_inductance());
        auto operatingPoints = converter.simulate_and_extract_operating_points(turnsRatios, inductance);
        OpenMagnetics::Inputs result;
        result.set_design_requirements(designReqs);
        result.set_operating_points(operatingPoints);
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::LLC_RESONANT_CONVERTER);
        return result;
    } else {
        OpenMagnetics::Inputs result = converter.process();
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::LLC_RESONANT_CONVERTER);
        return result;
    }
}

OpenMagnetics::Inputs process_cllc_internal(const json& converterJson, bool useNgspice) {
    OpenMagnetics::AdvancedCllcConverter converter(converterJson);
    converter._assertErrors = true;
    
    if (useNgspice) {
        auto designReqs = converter.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(OpenMagnetics::resolve_dimensional_values(tr));
        }
        double inductance = OpenMagnetics::resolve_dimensional_values(designReqs.get_magnetizing_inductance());
        auto operatingPoints = converter.simulate_and_extract_operating_points(turnsRatios, inductance);
        OpenMagnetics::Inputs result;
        result.set_design_requirements(designReqs);
        result.set_operating_points(operatingPoints);
        return result;
    } else {
        return converter.process();
    }
}

OpenMagnetics::Inputs process_dab_internal(const json& converterJson, bool useNgspice) {
    OpenMagnetics::AdvancedDab converter(converterJson);
    converter._assertErrors = true;
    
    if (useNgspice) {
        auto designReqs = converter.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(OpenMagnetics::resolve_dimensional_values(tr));
        }
        double inductance = OpenMagnetics::resolve_dimensional_values(designReqs.get_magnetizing_inductance());
        auto operatingPoints = converter.simulate_and_extract_operating_points(turnsRatios, inductance);
        OpenMagnetics::Inputs result;
        result.set_design_requirements(designReqs);
        result.set_operating_points(operatingPoints);
        return result;
    } else {
        return converter.process();
    }
}

OpenMagnetics::Inputs process_psfb_internal(const json& converterJson, bool useNgspice) {
    OpenMagnetics::AdvancedPsfb converter(converterJson);
    converter._assertErrors = true;
    
    if (useNgspice) {
        auto designReqs = converter.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(OpenMagnetics::resolve_dimensional_values(tr));
        }
        double inductance = OpenMagnetics::resolve_dimensional_values(designReqs.get_magnetizing_inductance());
        auto operatingPoints = converter.simulate_and_extract_operating_points(turnsRatios, inductance);
        OpenMagnetics::Inputs result;
        result.set_design_requirements(designReqs);
        result.set_operating_points(operatingPoints);
        return result;
    } else {
        return converter.process();
    }
}

OpenMagnetics::Inputs process_pshb_internal(const json& converterJson, bool useNgspice) {
    OpenMagnetics::AdvancedPshb converter(converterJson);
    converter._assertErrors = true;
    
    if (useNgspice) {
        auto designReqs = converter.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(OpenMagnetics::resolve_dimensional_values(tr));
        }
        double inductance = OpenMagnetics::resolve_dimensional_values(designReqs.get_magnetizing_inductance());
        auto operatingPoints = converter.simulate_and_extract_operating_points(turnsRatios, inductance);
        OpenMagnetics::Inputs result;
        result.set_design_requirements(designReqs);
        result.set_operating_points(operatingPoints);
        return result;
    } else {
        return converter.process();
    }
}

OpenMagnetics::Inputs process_isolated_buck_internal(const json& converterJson, bool useNgspice) {
    OpenMagnetics::AdvancedIsolatedBuck converter(converterJson);
    converter._assertErrors = true;
    
    if (useNgspice) {
        auto designReqs = converter.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(OpenMagnetics::resolve_dimensional_values(tr));
        }
        double inductance = OpenMagnetics::resolve_dimensional_values(designReqs.get_magnetizing_inductance());
        auto operatingPoints = converter.simulate_and_extract_operating_points(turnsRatios, inductance);
        OpenMagnetics::Inputs result;
        result.set_design_requirements(designReqs);
        result.set_operating_points(operatingPoints);
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::ISOLATED_BUCK_CONVERTER);
        return result;
    } else {
        OpenMagnetics::Inputs result = converter.process();
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::ISOLATED_BUCK_CONVERTER);
        return result;
    }
}

OpenMagnetics::Inputs process_isolated_buck_boost_internal(const json& converterJson, bool useNgspice) {
    OpenMagnetics::AdvancedIsolatedBuckBoost converter(converterJson);
    converter._assertErrors = true;
    
    if (useNgspice) {
        auto designReqs = converter.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(OpenMagnetics::resolve_dimensional_values(tr));
        }
        double inductance = OpenMagnetics::resolve_dimensional_values(designReqs.get_magnetizing_inductance());
        auto operatingPoints = converter.simulate_and_extract_operating_points(turnsRatios, inductance);
        OpenMagnetics::Inputs result;
        result.set_design_requirements(designReqs);
        result.set_operating_points(operatingPoints);
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::ISOLATED_BUCK_BOOST_CONVERTER);
        return result;
    } else {
        OpenMagnetics::Inputs result = converter.process();
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::ISOLATED_BUCK_BOOST_CONVERTER);
        return result;
    }
}

OpenMagnetics::Inputs process_current_transformer_internal(const json& converterJson) {
    OpenMagnetics::CurrentTransformer converter(converterJson);
    converter._assertErrors = true;
    
    double turnsRatio = converterJson.contains("turnsRatio") ? converterJson["turnsRatio"].get<double>() : 100.0;
    double secondaryResistance = converterJson.contains("secondaryResistance") ? converterJson["secondaryResistance"].get<double>() : 0.0;
    OpenMagnetics::Inputs result = converter.process(turnsRatio, secondaryResistance);
    result.get_mutable_design_requirements().set_topology(MAS::Topologies::CURRENT_TRANSFORMER);
    return result;
}

OpenMagnetics::Inputs process_pfc_internal(const json& converterJson) {
    OpenMagnetics::PowerFactorCorrection converter(converterJson);
    converter._assertErrors = true;
    return converter.process();
}

// ─────────────────────────────────────────────────────────────────────
// 9 new topology bindings (2026-05): Cuk, Sepic, Zeta,
// FourSwitchBuckBoost, AsymmetricHalfBridge, Weinberg, Vienna, Clllc,
// Src. All have first-class implementations in MKF's converter_models/
// — they were just not wired through PyMKF until now.
//
// Pattern groups:
//   * Single-inductor non-isolated (Cuk, Sepic, Zeta, FSBB) follow the
//     Buck/Boost pattern: AdvancedXxx::process() + simulate via inductance.
//   * Isolated with vector turns ratios (AsymHB, Clllc) follow the LLC
//     pattern.
//   * Weinberg has scalar turns ratio (single-secondary topology) — its
//     ngspice dispatch lives separately below.
//   * Vienna and Src have AdvancedXxx classes that only override
//     process_design_requirements(); the inherited Topology::process()
//     handles them. They follow the PFC pattern (no useNgspice branch).
// ─────────────────────────────────────────────────────────────────────

OpenMagnetics::Inputs process_cuk_internal(const json& converterJson, bool useNgspice) {
    OpenMagnetics::AdvancedCuk converter(converterJson);
    converter._assertErrors = true;
    if (useNgspice) {
        auto designReqs = converter.process_design_requirements();
        double inductance = OpenMagnetics::resolve_dimensional_values(designReqs.get_magnetizing_inductance());
        auto operatingPoints = converter.simulate_and_extract_operating_points(inductance);
        OpenMagnetics::Inputs result;
        result.set_design_requirements(designReqs);
        result.set_operating_points(operatingPoints);
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::CUK_CONVERTER);
        return result;
    } else {
        OpenMagnetics::Inputs result = converter.process();
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::CUK_CONVERTER);
        return result;
    }
}

OpenMagnetics::Inputs process_sepic_internal(const json& converterJson, bool useNgspice) {
    OpenMagnetics::AdvancedSepic converter(converterJson);
    converter._assertErrors = true;
    if (useNgspice) {
        auto designReqs = converter.process_design_requirements();
        double inductance = OpenMagnetics::resolve_dimensional_values(designReqs.get_magnetizing_inductance());
        auto operatingPoints = converter.simulate_and_extract_operating_points(inductance);
        OpenMagnetics::Inputs result;
        result.set_design_requirements(designReqs);
        result.set_operating_points(operatingPoints);
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::SEPIC_CONVERTER);
        return result;
    } else {
        OpenMagnetics::Inputs result = converter.process();
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::SEPIC_CONVERTER);
        return result;
    }
}

OpenMagnetics::Inputs process_zeta_internal(const json& converterJson, bool useNgspice) {
    OpenMagnetics::AdvancedZeta converter(converterJson);
    converter._assertErrors = true;
    if (useNgspice) {
        auto designReqs = converter.process_design_requirements();
        double inductance = OpenMagnetics::resolve_dimensional_values(designReqs.get_magnetizing_inductance());
        auto operatingPoints = converter.simulate_and_extract_operating_points(inductance);
        OpenMagnetics::Inputs result;
        result.set_design_requirements(designReqs);
        result.set_operating_points(operatingPoints);
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::ZETA_CONVERTER);
        return result;
    } else {
        OpenMagnetics::Inputs result = converter.process();
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::ZETA_CONVERTER);
        return result;
    }
}

OpenMagnetics::Inputs process_four_switch_buck_boost_internal(const json& converterJson, bool useNgspice) {
    OpenMagnetics::AdvancedFourSwitchBuckBoost converter(converterJson);
    converter._assertErrors = true;
    if (useNgspice) {
        auto designReqs = converter.process_design_requirements();
        double inductance = OpenMagnetics::resolve_dimensional_values(designReqs.get_magnetizing_inductance());
        auto operatingPoints = converter.simulate_and_extract_operating_points(inductance);
        OpenMagnetics::Inputs result;
        result.set_design_requirements(designReqs);
        result.set_operating_points(operatingPoints);
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::FOUR_SWITCH_BUCK_BOOST_CONVERTER);
        return result;
    } else {
        OpenMagnetics::Inputs result = converter.process();
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::FOUR_SWITCH_BUCK_BOOST_CONVERTER);
        return result;
    }
}

OpenMagnetics::Inputs process_asymmetric_half_bridge_internal(const json& converterJson, bool useNgspice) {
    OpenMagnetics::AdvancedAsymmetricHalfBridge converter(converterJson);
    converter._assertErrors = true;
    if (useNgspice) {
        auto designReqs = converter.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(OpenMagnetics::resolve_dimensional_values(tr));
        }
        double inductance = OpenMagnetics::resolve_dimensional_values(designReqs.get_magnetizing_inductance());
        auto operatingPoints = converter.simulate_and_extract_operating_points(turnsRatios, inductance);
        OpenMagnetics::Inputs result;
        result.set_design_requirements(designReqs);
        result.set_operating_points(operatingPoints);
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::ASYMMETRIC_HALF_BRIDGE_CONVERTER);
        return result;
    } else {
        OpenMagnetics::Inputs result = converter.process();
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::ASYMMETRIC_HALF_BRIDGE_CONVERTER);
        return result;
    }
}

OpenMagnetics::Inputs process_weinberg_internal(const json& converterJson, bool useNgspice) {
    OpenMagnetics::AdvancedWeinberg converter(converterJson);
    converter._assertErrors = true;
    if (useNgspice) {
        auto designReqs = converter.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(OpenMagnetics::resolve_dimensional_values(tr));
        }
        if (turnsRatios.empty()) {
            throw std::runtime_error("Weinberg converter requires at least one turns ratio");
        }
        double inductance = OpenMagnetics::resolve_dimensional_values(designReqs.get_magnetizing_inductance());
        auto operatingPoints = converter.simulate_and_extract_operating_points(turnsRatios[0], inductance);
        OpenMagnetics::Inputs result;
        result.set_design_requirements(designReqs);
        result.set_operating_points(operatingPoints);
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::WEINBERG_CONVERTER);
        return result;
    } else {
        OpenMagnetics::Inputs result = converter.process();
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::WEINBERG_CONVERTER);
        return result;
    }
}

OpenMagnetics::Inputs process_clllc_internal(const json& converterJson, bool useNgspice) {
    OpenMagnetics::AdvancedClllc converter(converterJson);
    converter._assertErrors = true;
    if (useNgspice) {
        auto designReqs = converter.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(OpenMagnetics::resolve_dimensional_values(tr));
        }
        double inductance = OpenMagnetics::resolve_dimensional_values(designReqs.get_magnetizing_inductance());
        auto operatingPoints = converter.simulate_and_extract_operating_points(turnsRatios, inductance);
        OpenMagnetics::Inputs result;
        result.set_design_requirements(designReqs);
        result.set_operating_points(operatingPoints);
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::CLLLC_RESONANT_CONVERTER);
        return result;
    } else {
        OpenMagnetics::Inputs result = converter.process();
        result.get_mutable_design_requirements().set_topology(MAS::Topologies::CLLLC_RESONANT_CONVERTER);
        return result;
    }
}

// Vienna and Src follow the PFC pattern: AdvancedXxx overrides only
// process_design_requirements(); we invoke the inherited Topology::process()
// via the simple base class on the AdvancedXxx instance. useNgspice has
// no separate path here — process() already handles operating-point
// simulation through the topology pipeline.
OpenMagnetics::Inputs process_vienna_internal(const json& converterJson, bool useNgspice) {
    (void)useNgspice;
    OpenMagnetics::AdvancedVienna converter(converterJson);
    converter._assertErrors = true;
    OpenMagnetics::Inputs result = converter.process();
    result.get_mutable_design_requirements().set_topology(MAS::Topologies::VIENNA_RECTIFIER_CONVERTER);
    return result;
}

OpenMagnetics::Inputs process_src_internal(const json& converterJson, bool useNgspice) {
    (void)useNgspice;
    OpenMagnetics::AdvancedSrc converter(converterJson);
    converter._assertErrors = true;
    OpenMagnetics::Inputs result = converter.process();
    result.get_mutable_design_requirements().set_topology(MAS::Topologies::SERIES_RESONANT_CONVERTER);
    return result;
}

OpenMagnetics::Inputs dispatch_converter(const std::string& topologyName, const json& converterJson, bool useNgspice) {
    if (topologyName == "flyback" || topologyName == "advanced_flyback") {
        return process_flyback_internal(converterJson, useNgspice);
    }
    else if (topologyName == "buck" || topologyName == "advanced_buck") {
        return process_buck_internal(converterJson, useNgspice);
    }
    else if (topologyName == "boost" || topologyName == "advanced_boost") {
        return process_boost_internal(converterJson, useNgspice);
    }
    else if (topologyName == "single_switch_forward") {
        return process_single_switch_forward_internal(converterJson, useNgspice);
    }
    else if (topologyName == "two_switch_forward") {
        return process_two_switch_forward_internal(converterJson, useNgspice);
    }
    else if (topologyName == "active_clamp_forward") {
        return process_active_clamp_forward_internal(converterJson, useNgspice);
    }
    else if (topologyName == "push_pull") {
        return process_push_pull_internal(converterJson, useNgspice);
    }
    else if (topologyName == "llc" || topologyName == "advanced_llc") {
        return process_llc_internal(converterJson, useNgspice);
    }
    else if (topologyName == "cllc" || topologyName == "advanced_cllc") {
        return process_cllc_internal(converterJson, useNgspice);
    }
    else if (topologyName == "dab" || topologyName == "advanced_dab") {
        return process_dab_internal(converterJson, useNgspice);
    }
    else if (topologyName == "phase_shifted_full_bridge" || topologyName == "psfb" ||
             topologyName == "advanced_phase_shifted_full_bridge" || topologyName == "advanced_psfb") {
        return process_psfb_internal(converterJson, useNgspice);
    }
    else if (topologyName == "phase_shifted_half_bridge" || topologyName == "pshb" ||
             topologyName == "advanced_phase_shifted_half_bridge" || topologyName == "advanced_pshb") {
        return process_pshb_internal(converterJson, useNgspice);
    }
    else if (topologyName == "isolated_buck") {
        return process_isolated_buck_internal(converterJson, useNgspice);
    }
    else if (topologyName == "isolated_buck_boost") {
        return process_isolated_buck_boost_internal(converterJson, useNgspice);
    }
    else if (topologyName == "current_transformer") {
        return process_current_transformer_internal(converterJson);
    }
    else if (topologyName == "power_factor_correction" || topologyName == "pfc") {
        return process_pfc_internal(converterJson);
    }
    else if (topologyName == "cuk" || topologyName == "advanced_cuk") {
        return process_cuk_internal(converterJson, useNgspice);
    }
    else if (topologyName == "sepic" || topologyName == "advanced_sepic") {
        return process_sepic_internal(converterJson, useNgspice);
    }
    else if (topologyName == "zeta" || topologyName == "advanced_zeta") {
        return process_zeta_internal(converterJson, useNgspice);
    }
    else if (topologyName == "four_switch_buck_boost" || topologyName == "advanced_four_switch_buck_boost") {
        return process_four_switch_buck_boost_internal(converterJson, useNgspice);
    }
    else if (topologyName == "asymmetric_half_bridge" || topologyName == "advanced_asymmetric_half_bridge") {
        return process_asymmetric_half_bridge_internal(converterJson, useNgspice);
    }
    else if (topologyName == "weinberg" || topologyName == "advanced_weinberg") {
        return process_weinberg_internal(converterJson, useNgspice);
    }
    else if (topologyName == "vienna" || topologyName == "advanced_vienna") {
        return process_vienna_internal(converterJson, useNgspice);
    }
    else if (topologyName == "clllc" || topologyName == "advanced_clllc") {
        return process_clllc_internal(converterJson, useNgspice);
    }
    else if (topologyName == "src" || topologyName == "advanced_src") {
        return process_src_internal(converterJson, useNgspice);
    }
    else {
        throw std::invalid_argument("Unknown topology: " + topologyName);
    }
}

json process_converter_internal(const std::string& topologyName, const json& converterJson, bool useNgspice) {
    try {
        OpenMagnetics::Inputs inputs = dispatch_converter(topologyName, converterJson, useNgspice);
        json result;
        to_json(result, inputs);
        return result;
    }
    catch (const std::exception& e) {
        json error;
        error["error"] = "Exception: " + std::string(e.what());
        return error;
    }
}

json process_converter(const std::string& topologyName, json converterJson, bool useNgspice) {
    // Normalize topology name (accept MAS 1.0 camelCase, pre-1.0 Title Case,
    // and the internal short form) and migrate any pre-1.0 enum strings
    // embedded in the converter JSON (e.g. flyback "mode": "Continuous
    // Conduction Mode" -> "continuousConductionMode") in-place.
    std::string normalized = normalize_topology_name(topologyName);
    OpenMagnetics::compat::migrate_pre_1_0(converterJson);
    return process_converter_internal(normalized, converterJson, useNgspice);
}

// Simplified design_magnetics_from_converter using MKF template method
json design_magnetics_from_converter(
    const std::string& topologyNameRaw,
    json converterJson, 
    int maxResults, 
    json coreModeJson, 
    bool useNgspice, 
    json weightsJson) {
    
    (void)useNgspice;  // Template method always uses ngspice

    // Accept MAS 1.0 camelCase, pre-1.0 Title Case, or internal short form.
    const std::string topologyName = normalize_topology_name(topologyNameRaw);
    OpenMagnetics::compat::migrate_pre_1_0(converterJson);

    try {
        OpenMagnetics::CoreAdviser::CoreAdviserModes coreMode;
        from_json(coreModeJson, coreMode);
        
        std::map<OpenMagnetics::MagneticFilters, double> weights;
        if (!weightsJson.is_null()) {
            for (auto& [key, value] : weightsJson.items()) {
                OpenMagnetics::MagneticFilters filter;
                OpenMagnetics::from_json(key, filter);
                weights[filter] = value;
            }
        }
        
        OpenMagnetics::MagneticAdviser magneticAdviser;
        magneticAdviser.set_core_mode(coreMode);
        
        std::vector<std::pair<OpenMagnetics::Mas, double>> masMagnetics;
        
        // Use MKF template method - handles all converters automatically.
        // For Buck / Boost the basic process_design_requirements() sets
        // only `magnetizingInductance.minimum` (a ripple-floor), and
        // the CoreAdviser picks the smallest core satisfying it
        // regardless of the L the user asked for. Dispatch to the
        // Advanced* subclass so spec.desiredInductance becomes the
        // nominal target (AdvancedBuck/Boost::process_design_requirements
        // override does this, MKF commit 8acd72c7).
        //
        // Flyback stays on the basic ctor: AdvancedFlyback's from_json
        // requires extra fields (desiredDutyCycle, desiredTurnsRatios)
        // that aren't in a normal converter spec. The basic
        // Flyback::process_design_requirements already sets nominal
        // from its V·s + duty-cycle derivation, so the magnetic gets
        // sized correctly without the Advanced* path.
        if (topologyName == "flyback" || topologyName == "advanced_flyback") {
            OpenMagnetics::Flyback converter(converterJson);
            converter._assertErrors = true;
            masMagnetics = weights.empty()
                ? magneticAdviser.get_advised_magnetic_from_converter(converter, maxResults)
                : magneticAdviser.get_advised_magnetic_from_converter(converter, weights, maxResults);
        }
        else if (topologyName == "buck" || topologyName == "advanced_buck") {
            OpenMagnetics::AdvancedBuck converter(converterJson);
            converter._assertErrors = true;
            masMagnetics = weights.empty()
                ? magneticAdviser.get_advised_magnetic_from_converter(converter, maxResults)
                : magneticAdviser.get_advised_magnetic_from_converter(converter, weights, maxResults);
        }
        else if (topologyName == "boost" || topologyName == "advanced_boost") {
            OpenMagnetics::AdvancedBoost converter(converterJson);
            converter._assertErrors = true;
            masMagnetics = weights.empty()
                ? magneticAdviser.get_advised_magnetic_from_converter(converter, maxResults)
                : magneticAdviser.get_advised_magnetic_from_converter(converter, weights, maxResults);
        }
        else if (topologyName == "single_switch_forward") {
            OpenMagnetics::SingleSwitchForward converter(converterJson);
            converter._assertErrors = true;
            masMagnetics = weights.empty()
                ? magneticAdviser.get_advised_magnetic_from_converter(converter, maxResults)
                : magneticAdviser.get_advised_magnetic_from_converter(converter, weights, maxResults);
        }
        else if (topologyName == "two_switch_forward") {
            OpenMagnetics::TwoSwitchForward converter(converterJson);
            converter._assertErrors = true;
            masMagnetics = weights.empty()
                ? magneticAdviser.get_advised_magnetic_from_converter(converter, maxResults)
                : magneticAdviser.get_advised_magnetic_from_converter(converter, weights, maxResults);
        }
        else if (topologyName == "active_clamp_forward") {
            OpenMagnetics::ActiveClampForward converter(converterJson);
            converter._assertErrors = true;
            masMagnetics = weights.empty()
                ? magneticAdviser.get_advised_magnetic_from_converter(converter, maxResults)
                : magneticAdviser.get_advised_magnetic_from_converter(converter, weights, maxResults);
        }
        else if (topologyName == "push_pull") {
            OpenMagnetics::PushPull converter(converterJson);
            converter._assertErrors = true;
            masMagnetics = weights.empty()
                ? magneticAdviser.get_advised_magnetic_from_converter(converter, maxResults)
                : magneticAdviser.get_advised_magnetic_from_converter(converter, weights, maxResults);
        }
        else if (topologyName == "isolated_buck") {
            OpenMagnetics::IsolatedBuck converter(converterJson);
            converter._assertErrors = true;
            masMagnetics = weights.empty()
                ? magneticAdviser.get_advised_magnetic_from_converter(converter, maxResults)
                : magneticAdviser.get_advised_magnetic_from_converter(converter, weights, maxResults);
        }
        else if (topologyName == "isolated_buck_boost") {
            OpenMagnetics::IsolatedBuckBoost converter(converterJson);
            converter._assertErrors = true;
            masMagnetics = weights.empty()
                ? magneticAdviser.get_advised_magnetic_from_converter(converter, maxResults)
                : magneticAdviser.get_advised_magnetic_from_converter(converter, weights, maxResults);
        }
        else if (topologyName == "llc" || topologyName == "advanced_llc") {
            OpenMagnetics::Llc converter(converterJson);
            converter._assertErrors = true;
            masMagnetics = weights.empty()
                ? magneticAdviser.get_advised_magnetic_from_converter(converter, maxResults)
                : magneticAdviser.get_advised_magnetic_from_converter(converter, weights, maxResults);
        }
        else if (topologyName == "cuk" || topologyName == "advanced_cuk") {
            OpenMagnetics::Cuk converter(converterJson);
            converter._assertErrors = true;
            masMagnetics = weights.empty()
                ? magneticAdviser.get_advised_magnetic_from_converter(converter, maxResults)
                : magneticAdviser.get_advised_magnetic_from_converter(converter, weights, maxResults);
        }
        else if (topologyName == "sepic" || topologyName == "advanced_sepic") {
            OpenMagnetics::Sepic converter(converterJson);
            converter._assertErrors = true;
            masMagnetics = weights.empty()
                ? magneticAdviser.get_advised_magnetic_from_converter(converter, maxResults)
                : magneticAdviser.get_advised_magnetic_from_converter(converter, weights, maxResults);
        }
        else if (topologyName == "zeta" || topologyName == "advanced_zeta") {
            OpenMagnetics::Zeta converter(converterJson);
            converter._assertErrors = true;
            masMagnetics = weights.empty()
                ? magneticAdviser.get_advised_magnetic_from_converter(converter, maxResults)
                : magneticAdviser.get_advised_magnetic_from_converter(converter, weights, maxResults);
        }
        else if (topologyName == "four_switch_buck_boost" || topologyName == "advanced_four_switch_buck_boost") {
            OpenMagnetics::FourSwitchBuckBoost converter(converterJson);
            converter._assertErrors = true;
            masMagnetics = weights.empty()
                ? magneticAdviser.get_advised_magnetic_from_converter(converter, maxResults)
                : magneticAdviser.get_advised_magnetic_from_converter(converter, weights, maxResults);
        }
        else if (topologyName == "asymmetric_half_bridge" || topologyName == "advanced_asymmetric_half_bridge") {
            OpenMagnetics::AsymmetricHalfBridge converter(converterJson);
            converter._assertErrors = true;
            masMagnetics = weights.empty()
                ? magneticAdviser.get_advised_magnetic_from_converter(converter, maxResults)
                : magneticAdviser.get_advised_magnetic_from_converter(converter, weights, maxResults);
        }
        else if (topologyName == "phase_shifted_full_bridge" || topologyName == "psfb" ||
                 topologyName == "advanced_phase_shifted_full_bridge" || topologyName == "advanced_psfb") {
            // Use the basic Psfb converter (mirrors the AHB / LLC branches):
            // its process_design_requirements() derives turns ratios and
            // magnetizing inductance from the converter spec. Falling through
            // to the generic else branch would dispatch to AdvancedPsfb, whose
            // from_json requires desiredTurnsRatios / desiredMagnetizingInductance
            // (the design *outputs*), which a plain converter spec does not carry.
            OpenMagnetics::Psfb converter(converterJson);
            converter._assertErrors = true;
            masMagnetics = weights.empty()
                ? magneticAdviser.get_advised_magnetic_from_converter(converter, maxResults)
                : magneticAdviser.get_advised_magnetic_from_converter(converter, weights, maxResults);
        }
        else if (topologyName == "weinberg" || topologyName == "advanced_weinberg") {
            OpenMagnetics::Weinberg converter(converterJson);
            converter._assertErrors = true;
            masMagnetics = weights.empty()
                ? magneticAdviser.get_advised_magnetic_from_converter(converter, maxResults)
                : magneticAdviser.get_advised_magnetic_from_converter(converter, weights, maxResults);
        }
        else if (topologyName == "vienna" || topologyName == "advanced_vienna") {
            OpenMagnetics::Vienna converter(converterJson);
            converter._assertErrors = true;
            masMagnetics = weights.empty()
                ? magneticAdviser.get_advised_magnetic_from_converter(converter, maxResults)
                : magneticAdviser.get_advised_magnetic_from_converter(converter, weights, maxResults);
        }
        else if (topologyName == "clllc" || topologyName == "advanced_clllc") {
            OpenMagnetics::Clllc converter(converterJson);
            converter._assertErrors = true;
            masMagnetics = weights.empty()
                ? magneticAdviser.get_advised_magnetic_from_converter(converter, maxResults)
                : magneticAdviser.get_advised_magnetic_from_converter(converter, weights, maxResults);
        }
        else if (topologyName == "src" || topologyName == "advanced_src") {
            OpenMagnetics::Src converter(converterJson);
            converter._assertErrors = true;
            masMagnetics = weights.empty()
                ? magneticAdviser.get_advised_magnetic_from_converter(converter, maxResults)
                : magneticAdviser.get_advised_magnetic_from_converter(converter, weights, maxResults);
        }
        else if (topologyName == "dab" || topologyName == "advanced_dab" ||
                 topologyName == "dual_active_bridge") {
            // Use the base Dab model (mirrors llc / src / clllc above): it
            // derives turns ratios + magnetizing inductance itself from the
            // V1/V2 windows, so it does not require the AdvancedDab
            // desiredTurnsRatios / desiredMagnetizingInductance inputs.
            OpenMagnetics::Dab converter(converterJson);
            converter._assertErrors = true;
            masMagnetics = weights.empty()
                ? magneticAdviser.get_advised_magnetic_from_converter(converter, maxResults)
                : magneticAdviser.get_advised_magnetic_from_converter(converter, weights, maxResults);
        }
        else {
            // Fall back to old approach for other topologies
            json inputsJson = process_converter_internal(topologyName, converterJson, useNgspice);
            if (inputsJson.contains("error")) return inputsJson;
            
            OpenMagnetics::Inputs inputs(inputsJson);
            masMagnetics = weights.empty()
                ? magneticAdviser.get_advised_magnetic(inputs, maxResults)
                : magneticAdviser.get_advised_magnetic(inputs, weights, maxResults);
        }
        
        // Build results (same as before)
        auto scoringsPerFilter = magneticAdviser.get_scorings();
        
        json results = json();
        results["data"] = json::array();
        for (size_t i = 0; i < masMagnetics.size(); ++i) {
            auto& masMagnetic = masMagnetics[i].first;
            double scoring = masMagnetics[i].second;
            std::string name = masMagnetic.get_magnetic().get_manufacturer_info().value().get_reference().value();
            json result;
            json masJson;
            to_json(masJson, masMagnetic);
            result["mas"] = masJson;
            result["scoring"] = scoring;
            if (scoringsPerFilter.count(name)) {
                json filterScorings;
                for (auto& filterPair : scoringsPerFilter[name]) {
                    auto filter = filterPair.first;
                    double filterScore = filterPair.second;
                    filterScorings[std::string(magic_enum::enum_name(filter))] = filterScore;
                }
                result["scoringPerFilter"] = filterScorings;
            }
            results["data"].push_back(result);
        }
        
        sort(results["data"].begin(), results["data"].end(), [](json& b1, json& b2) {
            return b1["scoring"] > b2["scoring"];
        });
        
        OpenMagnetics::settings.reset();
        return results;
    }
    catch (const std::exception& e) {
        json error;
        error["error"] = "Exception: " + std::string(e.what());
        return error;
    }
}


json process_flyback(json flybackJson) {
    return process_converter("flyback", flybackJson, true);
}

json process_buck(json buckJson) {
    return process_converter("buck", buckJson, true);
}

json process_boost(json boostJson) {
    return process_converter("boost", boostJson, true);
}

json process_single_switch_forward(json forwardJson) {
    return process_converter("single_switch_forward", forwardJson, true);
}

json process_two_switch_forward(json forwardJson) {
    return process_converter("two_switch_forward", forwardJson, true);
}

json process_active_clamp_forward(json forwardJson) {
    return process_converter("active_clamp_forward", forwardJson, true);
}

json process_push_pull(json pushPullJson) {
    return process_converter("push_pull", pushPullJson, true);
}

json process_isolated_buck(json isolatedBuckJson) {
    return process_converter("isolated_buck", isolatedBuckJson, true);
}

json process_isolated_buck_boost(json isolatedBuckBoostJson) {
    return process_converter("isolated_buck_boost", isolatedBuckBoostJson, true);
}

json process_current_transformer(json ctJson, double turnsRatio, double secondaryResistance) {
    ctJson["turnsRatio"] = turnsRatio;
    ctJson["secondaryResistance"] = secondaryResistance;
    return process_converter("current_transformer", ctJson, true);
}

json process_cuk(json cukJson) { return process_converter("cuk", cukJson, true); }
json process_sepic(json sepicJson) { return process_converter("sepic", sepicJson, true); }
json process_zeta(json zetaJson) { return process_converter("zeta", zetaJson, true); }
json process_four_switch_buck_boost(json j) { return process_converter("four_switch_buck_boost", j, true); }
json process_asymmetric_half_bridge(json j) { return process_converter("asymmetric_half_bridge", j, true); }
json process_weinberg(json j) { return process_converter("weinberg", j, true); }
json process_vienna(json j) { return process_converter("vienna", j, true); }
json process_clllc(json j) { return process_converter("clllc", j, true); }
json process_src(json j) { return process_converter("src", j, true); }

// ─────────────────────────────────────────────────────────────────────
// get_extra_components_inputs (MKF 2026-04-29 API)
//
// Returns a JSON list describing the design requirements / operating
// points of any extra components a topology brings along besides its
// main magnetic — resonant tank caps for LLC, snubber nets, auxiliary
// inductors, etc. The downstream caller can pipe each entry into a
// component-search (TAS) or a fresh PyOM design call.
//
// Each list entry is one of:
//   { "kind": "magnetic",  "inputs": { ...MAS::Inputs JSON... } }
//   { "kind": "capacitor", "inputs": { ...CAS::Inputs JSON... } }
//
// Arguments:
//   topology_name      — same key set as dispatch_converter (e.g. "llc",
//                        "psfb", "active_clamp_forward").
//   converter_json     — the converter spec used to construct the simple
//                        topology class (NOT the AdvancedXxx version).
//   mode               — "IDEAL" or "REAL". IDEAL returns archetypal
//                        requirements without parasitics; REAL refines
//                        with the main magnetic's actual leakage,
//                        winding resistance, etc.
//   magnetic_json      — optional MAS::Magnetic JSON, required by some
//                        topologies in REAL mode. Pass null in IDEAL.
// ─────────────────────────────────────────────────────────────────────

namespace {

// Construct the simple topology subclass (NOT AdvancedXxx), run the
// design pipeline (process() / process_design_requirements()) so the
// internal resonant-tank / extra-component values are populated, and
// then invoke get_extra_components_inputs on it. The simple classes
// own the design physics, so they're what the new API was designed
// against. Some topologies require process() to have run first; we
// invoke it via the SFINAE-detected member that exists on each.
template <typename TopologyT>
std::vector<std::variant<OpenMagnetics::Inputs, CAS::Inputs>>
collect_extra_components(const json& converterJson,
                         OpenMagnetics::ExtraComponentsMode mode,
                         std::optional<OpenMagnetics::Magnetic> magnetic) {
    TopologyT topology(converterJson);
    topology._assertErrors = true;
    // Run process() — the Topology base class entry point that wires
    // process_design_requirements() + process_operating_points() and
    // populates the internal Lr / Cr / inductance ratio that
    // get_extra_components_inputs() reads. The Inputs return value is
    // discarded; we only care about the side effects on `topology`.
    topology.process();
    return topology.get_extra_components_inputs(mode, magnetic);
}

std::vector<std::variant<OpenMagnetics::Inputs, CAS::Inputs>>
dispatch_extra_components(const std::string& topologyName,
                          const json& converterJson,
                          OpenMagnetics::ExtraComponentsMode mode,
                          std::optional<OpenMagnetics::Magnetic> magnetic) {
    if (topologyName == "flyback")
        return collect_extra_components<OpenMagnetics::Flyback>(converterJson, mode, magnetic);
    if (topologyName == "buck")
        return collect_extra_components<OpenMagnetics::Buck>(converterJson, mode, magnetic);
    if (topologyName == "boost")
        return collect_extra_components<OpenMagnetics::Boost>(converterJson, mode, magnetic);
    if (topologyName == "single_switch_forward")
        return collect_extra_components<OpenMagnetics::SingleSwitchForward>(converterJson, mode, magnetic);
    if (topologyName == "two_switch_forward")
        return collect_extra_components<OpenMagnetics::TwoSwitchForward>(converterJson, mode, magnetic);
    if (topologyName == "active_clamp_forward")
        return collect_extra_components<OpenMagnetics::ActiveClampForward>(converterJson, mode, magnetic);
    if (topologyName == "push_pull")
        return collect_extra_components<OpenMagnetics::PushPull>(converterJson, mode, magnetic);
    if (topologyName == "llc")
        return collect_extra_components<OpenMagnetics::Llc>(converterJson, mode, magnetic);
    if (topologyName == "dab")
        return collect_extra_components<OpenMagnetics::Dab>(converterJson, mode, magnetic);
    if (topologyName == "phase_shifted_full_bridge" || topologyName == "psfb")
        return collect_extra_components<OpenMagnetics::Psfb>(converterJson, mode, magnetic);
    if (topologyName == "phase_shifted_half_bridge" || topologyName == "pshb")
        return collect_extra_components<OpenMagnetics::Pshb>(converterJson, mode, magnetic);
    if (topologyName == "cllc")
        return collect_extra_components<OpenMagnetics::CllcConverter>(converterJson, mode, magnetic);
    if (topologyName == "isolated_buck")
        return collect_extra_components<OpenMagnetics::IsolatedBuck>(converterJson, mode, magnetic);
    if (topologyName == "isolated_buck_boost")
        return collect_extra_components<OpenMagnetics::IsolatedBuckBoost>(converterJson, mode, magnetic);
    if (topologyName == "cuk" || topologyName == "advanced_cuk")
        return collect_extra_components<OpenMagnetics::Cuk>(converterJson, mode, magnetic);
    if (topologyName == "sepic" || topologyName == "advanced_sepic")
        return collect_extra_components<OpenMagnetics::Sepic>(converterJson, mode, magnetic);
    if (topologyName == "zeta" || topologyName == "advanced_zeta")
        return collect_extra_components<OpenMagnetics::Zeta>(converterJson, mode, magnetic);
    if (topologyName == "four_switch_buck_boost" || topologyName == "advanced_four_switch_buck_boost")
        return collect_extra_components<OpenMagnetics::FourSwitchBuckBoost>(converterJson, mode, magnetic);
    if (topologyName == "asymmetric_half_bridge" || topologyName == "advanced_asymmetric_half_bridge")
        return collect_extra_components<OpenMagnetics::AsymmetricHalfBridge>(converterJson, mode, magnetic);
    if (topologyName == "weinberg" || topologyName == "advanced_weinberg")
        return collect_extra_components<OpenMagnetics::Weinberg>(converterJson, mode, magnetic);
    if (topologyName == "clllc" || topologyName == "advanced_clllc")
        return collect_extra_components<OpenMagnetics::Clllc>(converterJson, mode, magnetic);
    if (topologyName == "src" || topologyName == "advanced_src")
        return collect_extra_components<OpenMagnetics::Src>(converterJson, mode, magnetic);
    if (topologyName == "vienna" || topologyName == "advanced_vienna") {
        // Vienna 3-φ rectifier doesn't carry detached extras beyond
        // the boost inductor (the main magnetic itself). Return an
        // empty list so the pipeline can proceed instead of throwing.
        return {};
    }
    throw std::invalid_argument(
        "get_extra_components_inputs: topology '" + topologyName +
        "' has no extra-components dispatch (or hasn't been wired in PyMKF)."
    );
}

}  // namespace

// ─────────────────────────────────────────────────────────────────────
// generate_ngspice_circuit (MKF native SPICE generator)
//
// Returns the canonical ngspice deck for a converter at a given operating
// point, sized from the magnetic that was designed in Phase 3. This
// replaces hand-rolled netlist templates on the Proteus side: MKF owns
// the topology wiring (rectifier polarity, gate-drive timing, dead-time,
// snubbers, K-coupling) so the simulation always validates the actual
// designed magnetic.
//
// Args mirror the C++ method: turns_ratios (Np:Ns vector), magnetizing
// inductance, optional vin / operating-point indices. The simple
// topology class is constructed and process()'d first so internal state
// (Lr / Cr for LLC, etc.) is populated before deck generation.
// ─────────────────────────────────────────────────────────────────────
namespace {
// Apply optional bridge simulation mode override on a topology instance.
// Must be called BEFORE topology.process() because process() consults the
// mode when computing dead time and other timing parameters.
//
// modeStr accepts:
//   "" | "default" | "pulse" | "behavioral_pulse" → no-op (BEHAVIORAL_PULSE
//       stays as default; safe for non-bridge topologies because
//       set_bridge_simulation_mode throws if you try to force PULSE on a
//       non-bridge topology).
//   "switch" | "voltage_controlled_switch" → set VOLTAGE_CONTROLLED_SWITCH,
//       producing a real SW1/body-diode/snubber bridge instead of a single
//       dependent V-source. Required when downstream tooling needs to size
//       the bridge MOSFETs (e.g. Heaviside's TAS decomposer).
//   anything else → returns false so the caller can raise a clear error.
template <typename TopologyT>
bool apply_bridge_simulation_mode(TopologyT& topology, const std::string& modeStr) {
    if (modeStr.empty() || modeStr == "default" ||
        modeStr == "pulse" || modeStr == "behavioral_pulse") {
        return true;
    }
    if (modeStr == "switch" || modeStr == "voltage_controlled_switch") {
        topology.set_bridge_simulation_mode(
            OpenMagnetics::BridgeSimulationMode::VOLTAGE_CONTROLLED_SWITCH);
        return true;
    }
    return false;
}

// Apply a partial SpiceSimulationConfig override from a JSON dict. Missing
// keys keep the per-topology default. Unknown keys are silently ignored so
// that callers built against an older binding still work after new fields
// are added.
//
// Must be called BEFORE topology.process() because the config is read by
// the netlist emitters in generate_ngspice_circuit (which runs after
// process()).
template <typename TopologyT>
void apply_spice_simulation_config(TopologyT& topology, const json& cfgJson) {
    if (cfgJson.is_null()) return;
    if (!cfgJson.is_object()) {
        throw std::runtime_error(
            "generate_ngspice_circuit: spice_config must be a JSON object");
    }
    if (cfgJson.empty()) return;

    OpenMagnetics::SpiceSimulationConfig cfg = topology.spice_config();
    auto take_d = [&](const char* k, double& dst) {
        auto it = cfgJson.find(k);
        if (it != cfgJson.end() && it->is_number()) dst = it->get<double>();
    };
    // Overload for fields that recently switched to std::optional<double> in
    // MAS (e.g. snubRReal). If JSON provides a value we set it; otherwise the
    // optional keeps its existing state (Topology default).
    auto take_od = [&](const char* k, std::optional<double>& dst) {
        auto it = cfgJson.find(k);
        if (it != cfgJson.end() && it->is_number()) dst = it->get<double>();
    };
    auto take_i = [&](const char* k, int& dst) {
        auto it = cfgJson.find(k);
        if (it != cfgJson.end() && it->is_number_integer()) dst = it->get<int>();
    };
    auto take_s = [&](const char* k, std::string& dst) {
        auto it = cfgJson.find(k);
        if (it != cfgJson.end() && it->is_string()) dst = it->get<std::string>();
    };
    take_d("pwmHigh",                        cfg.pwmHigh);
    take_d("pwmRise",                        cfg.pwmRise);
    take_d("pwmFall",                        cfg.pwmFall);
    take_d("swModelVT",                      cfg.swModelVT);
    take_d("swModelVH",                      cfg.swModelVH);
    take_d("swModelRON",                     cfg.swModelRON);
    take_d("swModelROFF",                    cfg.swModelROFF);
    take_d("snubR",                          cfg.snubR);
    take_d("snubC",                          cfg.snubC);
    take_od("snubRReal",                     cfg.snubRReal);
    take_d("diodeIS",                        cfg.diodeIS);
    take_d("diodeRS",                        cfg.diodeRS);
    take_d("outputCapacitance",              cfg.outputCapacitance);
    take_d("outputCapInitialChargeFraction", cfg.outputCapInitialChargeFraction);
    take_i("samplesPerPeriod",               cfg.samplesPerPeriod);
    take_d("relTol",                         cfg.relTol);
    take_d("absTol",                         cfg.absTol);
    take_d("vnTol",                          cfg.vnTol);
    take_i("itl1",                           cfg.itl1);
    take_i("itl4",                           cfg.itl4);
    take_s("method",                         cfg.method);
    take_d("trTol",                          cfg.trTol);
    topology.set_spice_config(std::move(cfg));
}

// Isolated topologies (transformer): pass turns ratios + magnetizing inductance.
template <typename TopologyT>
std::string generate_spice_isolated(const json& converterJson,
                                     const std::vector<double>& turnsRatios,
                                     double magnetizingInductance,
                                     size_t vinIdx, size_t opIdx,
                                     const std::string& bridgeMode,
                                     const json& spiceConfig) {
    TopologyT topology(converterJson);
    topology._assertErrors = true;
    if (!apply_bridge_simulation_mode(topology, bridgeMode)) {
        throw std::runtime_error(
            "generate_ngspice_circuit: unknown bridge_simulation_mode '" +
            bridgeMode + "' — expected '', 'pulse', or 'switch'");
    }
    apply_spice_simulation_config(topology, spiceConfig);
    topology.process();
    return topology.generate_ngspice_circuit(turnsRatios, magnetizingInductance, vinIdx, opIdx);
}

// Non-isolated topologies (single inductor): just pass the inductance.
template <typename TopologyT>
std::string generate_spice_inductor(const json& converterJson,
                                     double inductance,
                                     size_t vinIdx, size_t opIdx,
                                     const std::string& bridgeMode,
                                     const json& spiceConfig) {
    TopologyT topology(converterJson);
    topology._assertErrors = true;
    if (!apply_bridge_simulation_mode(topology, bridgeMode)) {
        throw std::runtime_error(
            "generate_ngspice_circuit: unknown bridge_simulation_mode '" +
            bridgeMode + "' — expected '', 'pulse', or 'switch'");
    }
    apply_spice_simulation_config(topology, spiceConfig);
    topology.process();
    return topology.generate_ngspice_circuit(inductance, vinIdx, opIdx);
}

// Isolated topologies with a single secondary winding (Weinberg): the
// generate_ngspice_circuit signature takes a scalar turns ratio rather
// than a vector. We pass turnsRatios[0] if available, otherwise 1.0.
template <typename TopologyT>
std::string generate_spice_isolated_scalar(const json& converterJson,
                                            const std::vector<double>& turnsRatios,
                                            double magnetizingInductance,
                                            size_t vinIdx, size_t opIdx,
                                            const std::string& bridgeMode,
                                            const json& spiceConfig) {
    TopologyT topology(converterJson);
    topology._assertErrors = true;
    if (!apply_bridge_simulation_mode(topology, bridgeMode)) {
        throw std::runtime_error(
            "generate_ngspice_circuit: unknown bridge_simulation_mode '" +
            bridgeMode + "' — expected '', 'pulse', or 'switch'");
    }
    apply_spice_simulation_config(topology, spiceConfig);
    topology.process();
    double turnsRatio = turnsRatios.empty() ? 1.0 : turnsRatios[0];
    return topology.generate_ngspice_circuit(turnsRatio, magnetizingInductance, vinIdx, opIdx);
}

}  // namespace

json generate_ngspice_circuit(const std::string& topologyName,
                              json converterJson,
                              std::vector<double> turnsRatios,
                              double magnetizingInductance,
                              size_t vinIdx,
                              size_t opIdx,
                              std::string bridgeMode,
                              json spiceConfig) {
    try {
        std::string spice;
        // Non-isolated single-inductor topologies — magnetizingInductance
        // arg is interpreted as the main inductor value.
        if (topologyName == "buck")
            spice = generate_spice_inductor<OpenMagnetics::Buck>(converterJson, magnetizingInductance, vinIdx, opIdx, bridgeMode, spiceConfig);
        else if (topologyName == "boost")
            spice = generate_spice_inductor<OpenMagnetics::Boost>(converterJson, magnetizingInductance, vinIdx, opIdx, bridgeMode, spiceConfig);
        // Isolated topologies — turnsRatios + magnetizing inductance.
        else if (topologyName == "flyback") spice = generate_spice_isolated<OpenMagnetics::Flyback>(converterJson, turnsRatios, magnetizingInductance, vinIdx, opIdx, bridgeMode, spiceConfig);
        else if (topologyName == "single_switch_forward") spice = generate_spice_isolated<OpenMagnetics::SingleSwitchForward>(converterJson, turnsRatios, magnetizingInductance, vinIdx, opIdx, bridgeMode, spiceConfig);
        else if (topologyName == "two_switch_forward")    spice = generate_spice_isolated<OpenMagnetics::TwoSwitchForward>(converterJson, turnsRatios, magnetizingInductance, vinIdx, opIdx, bridgeMode, spiceConfig);
        else if (topologyName == "active_clamp_forward")  spice = generate_spice_isolated<OpenMagnetics::ActiveClampForward>(converterJson, turnsRatios, magnetizingInductance, vinIdx, opIdx, bridgeMode, spiceConfig);
        else if (topologyName == "push_pull") spice = generate_spice_isolated<OpenMagnetics::PushPull>(converterJson, turnsRatios, magnetizingInductance, vinIdx, opIdx, bridgeMode, spiceConfig);
        else if (topologyName == "llc")       spice = generate_spice_isolated<OpenMagnetics::Llc>(converterJson, turnsRatios, magnetizingInductance, vinIdx, opIdx, bridgeMode, spiceConfig);
        // CLLC has a distinct signature: generate_ngspice_circuit takes
        // (double turnsRatio, const CllcResonantParameters&, ...) instead
        // of the uniform (vector<double> turnsRatios, double Lm, ...)
        // shape, because CLLC needs the full resonant-parameter set
        // (Lr1, Cr1, Lr2, Cr2, Lm) computed by
        // CllcConverter::calculate_resonant_parameters() — a single
        // magnetizing inductance scalar would not determine the tank.
        else if (topologyName == "cllc" || topologyName == "advanced_cllc") {
            OpenMagnetics::AdvancedCllcConverter topology(converterJson);
            topology._assertErrors = true;
            if (!apply_bridge_simulation_mode(topology, bridgeMode)) {
                throw std::runtime_error(
                    "generate_ngspice_circuit: unknown bridge_simulation_mode '" +
                    bridgeMode + "' — expected '', 'pulse', or 'switch'");
            }
            apply_spice_simulation_config(topology, spiceConfig);
            topology.process();
            auto params = topology.calculate_resonant_parameters();
            double turnsRatio = turnsRatios.empty() ? 1.0 : turnsRatios[0];
            spice = topology.generate_ngspice_circuit(turnsRatio, params, vinIdx, opIdx);
        }
        else if (topologyName == "dab" || topologyName == "dual_active_bridge" || topologyName == "advanced_dab")
            spice = generate_spice_isolated<OpenMagnetics::Dab>(converterJson, turnsRatios, magnetizingInductance, vinIdx, opIdx, bridgeMode, spiceConfig);
        else if (topologyName == "phase_shifted_full_bridge" || topologyName == "psfb") spice = generate_spice_isolated<OpenMagnetics::Psfb>(converterJson, turnsRatios, magnetizingInductance, vinIdx, opIdx, bridgeMode, spiceConfig);
        else if (topologyName == "phase_shifted_half_bridge" || topologyName == "pshb") spice = generate_spice_isolated<OpenMagnetics::Pshb>(converterJson, turnsRatios, magnetizingInductance, vinIdx, opIdx, bridgeMode, spiceConfig);
        else if (topologyName == "isolated_buck")       spice = generate_spice_isolated<OpenMagnetics::IsolatedBuck>(converterJson, turnsRatios, magnetizingInductance, vinIdx, opIdx, bridgeMode, spiceConfig);
        else if (topologyName == "isolated_buck_boost") spice = generate_spice_isolated<OpenMagnetics::IsolatedBuckBoost>(converterJson, turnsRatios, magnetizingInductance, vinIdx, opIdx, bridgeMode, spiceConfig);
        // Single-inductor non-isolated (2026-05): Cuk, Sepic, Zeta, FSBB.
        else if (topologyName == "cuk" || topologyName == "advanced_cuk")
            spice = generate_spice_inductor<OpenMagnetics::Cuk>(converterJson, magnetizingInductance, vinIdx, opIdx, bridgeMode, spiceConfig);
        else if (topologyName == "sepic" || topologyName == "advanced_sepic")
            spice = generate_spice_inductor<OpenMagnetics::Sepic>(converterJson, magnetizingInductance, vinIdx, opIdx, bridgeMode, spiceConfig);
        else if (topologyName == "zeta" || topologyName == "advanced_zeta")
            spice = generate_spice_inductor<OpenMagnetics::Zeta>(converterJson, magnetizingInductance, vinIdx, opIdx, bridgeMode, spiceConfig);
        else if (topologyName == "four_switch_buck_boost" || topologyName == "advanced_four_switch_buck_boost")
            spice = generate_spice_inductor<OpenMagnetics::FourSwitchBuckBoost>(converterJson, magnetizingInductance, vinIdx, opIdx, bridgeMode, spiceConfig);
        // Isolated with vector turns ratios (2026-05): AsymHB, Clllc, Src, Vienna.
        else if (topologyName == "asymmetric_half_bridge" || topologyName == "advanced_asymmetric_half_bridge")
            spice = generate_spice_isolated<OpenMagnetics::AsymmetricHalfBridge>(converterJson, turnsRatios, magnetizingInductance, vinIdx, opIdx, bridgeMode, spiceConfig);
        else if (topologyName == "clllc" || topologyName == "advanced_clllc")
            spice = generate_spice_isolated<OpenMagnetics::Clllc>(converterJson, turnsRatios, magnetizingInductance, vinIdx, opIdx, bridgeMode, spiceConfig);
        else if (topologyName == "src" || topologyName == "advanced_src")
            spice = generate_spice_isolated<OpenMagnetics::Src>(converterJson, turnsRatios, magnetizingInductance, vinIdx, opIdx, bridgeMode, spiceConfig);
        else if (topologyName == "vienna" || topologyName == "advanced_vienna")
            spice = generate_spice_isolated<OpenMagnetics::Vienna>(converterJson, turnsRatios, magnetizingInductance, vinIdx, opIdx, bridgeMode, spiceConfig);
        // Isolated with scalar turns ratio (single secondary): Weinberg.
        else if (topologyName == "weinberg" || topologyName == "advanced_weinberg")
            spice = generate_spice_isolated_scalar<OpenMagnetics::Weinberg>(converterJson, turnsRatios, magnetizingInductance, vinIdx, opIdx, bridgeMode, spiceConfig);
        // PFC has a distinct method signature — (inductance, dcResistance,
        // simulationTime, timeStep) — because its deck simulates over a
        // line cycle rather than a switching cycle. Use the method's
        // defaults for dcResistance/simulationTime/timeStep; vinIdx/opIdx
        // are not meaningful for PFC (line voltage drives the operating
        // point).
        else if (topologyName == "power_factor_correction"
              || topologyName == "advanced_power_factor_correction"
              || topologyName == "pfc") {
            OpenMagnetics::PowerFactorCorrection topology(converterJson);
            topology._assertErrors = true;
            apply_spice_simulation_config(topology, spiceConfig);
            topology.process();
            spice = topology.generate_ngspice_circuit(magnetizingInductance);
        }
        else return json{{"error", "generate_ngspice_circuit: unknown topology '" + topologyName + "'"}};
        return json{{"netlist", spice}};
    } catch (const std::exception& exc) {
        return json{{"error", std::string("generate_ngspice_circuit: ") + exc.what()}};
    }
}

json get_extra_components_inputs(const std::string& topologyName,
                                 json converterJson,
                                 const std::string& modeStr,
                                 json magneticJson) {
    try {
        OpenMagnetics::ExtraComponentsMode mode;
        if (modeStr == "IDEAL" || modeStr == "ideal")
            mode = OpenMagnetics::ExtraComponentsMode::IDEAL;
        else if (modeStr == "REAL" || modeStr == "real")
            mode = OpenMagnetics::ExtraComponentsMode::REAL;
        else
            return json{{"error", "mode must be 'IDEAL' or 'REAL', got: " + modeStr}};

        std::optional<OpenMagnetics::Magnetic> magnetic;
        if (!magneticJson.is_null()) {
            magnetic = OpenMagnetics::Magnetic(magneticJson);
        }

        auto components = dispatch_extra_components(topologyName, converterJson, mode, magnetic);

        json out = json::array();
        for (const auto& comp : components) {
            json entry;
            std::visit([&entry](const auto& inputs) {
                using T = std::decay_t<decltype(inputs)>;
                json inputsJson;
                to_json(inputsJson, inputs);
                if constexpr (std::is_same_v<T, OpenMagnetics::Inputs>) {
                    entry["kind"] = "magnetic";
                } else {
                    // CAS::Inputs
                    entry["kind"] = "capacitor";
                }
                entry["inputs"] = inputsJson;
            }, comp);
            out.push_back(entry);
        }
        return out;
    } catch (const std::exception& exc) {
        return json{{"error", std::string("get_extra_components_inputs: ") + exc.what()}};
    }
}

// ============================================================================
// Converter input parity aliases (Phase B port from WebLibMKF 2026-05).
// WASM exposes `calculate_<topology>_inputs(json)` and
// `calculate_advanced_<topology>_inputs(json)` per topology; PyMKF's
// `process_converter("<topology>", json, true)` already covers the same
// dispatch, so each alias is a one-line forward. Diagnostic-only fields
// added by WASM (per-OP duty cycle / peak current / etc.) are intentionally
// omitted — they're for wizard display, not for the magnetics adviser.
// ============================================================================

#define DEFINE_CONVERTER_ALIAS(fn_name, topology) \
    json fn_name(json converterJson) { return process_converter(topology, converterJson, true); }

DEFINE_CONVERTER_ALIAS(calculate_flyback_inputs, "flyback")
DEFINE_CONVERTER_ALIAS(calculate_advanced_flyback_inputs, "advanced_flyback")
DEFINE_CONVERTER_ALIAS(calculate_buck_inputs, "buck")
DEFINE_CONVERTER_ALIAS(calculate_advanced_buck_inputs, "advanced_buck")
DEFINE_CONVERTER_ALIAS(calculate_boost_inputs, "boost")
DEFINE_CONVERTER_ALIAS(calculate_advanced_boost_inputs, "advanced_boost")
DEFINE_CONVERTER_ALIAS(calculate_single_switch_forward_inputs, "single_switch_forward")
DEFINE_CONVERTER_ALIAS(calculate_advanced_single_switch_forward_inputs, "single_switch_forward")
DEFINE_CONVERTER_ALIAS(calculate_two_switch_forward_inputs, "two_switch_forward")
DEFINE_CONVERTER_ALIAS(calculate_advanced_two_switch_forward_inputs, "two_switch_forward")
DEFINE_CONVERTER_ALIAS(calculate_active_clamp_forward_inputs, "active_clamp_forward")
DEFINE_CONVERTER_ALIAS(calculate_advanced_active_clamp_forward_inputs, "active_clamp_forward")
DEFINE_CONVERTER_ALIAS(calculate_push_pull_inputs, "push_pull")
DEFINE_CONVERTER_ALIAS(calculate_advanced_push_pull_inputs, "push_pull")
DEFINE_CONVERTER_ALIAS(calculate_isolated_buck_inputs, "isolated_buck")
DEFINE_CONVERTER_ALIAS(calculate_advanced_isolated_buck_inputs, "isolated_buck")
DEFINE_CONVERTER_ALIAS(calculate_isolated_buck_boost_inputs, "isolated_buck_boost")
DEFINE_CONVERTER_ALIAS(calculate_advanced_isolated_buck_boost_inputs, "isolated_buck_boost")
DEFINE_CONVERTER_ALIAS(calculate_cuk_inputs, "cuk")
DEFINE_CONVERTER_ALIAS(calculate_advanced_cuk_inputs, "advanced_cuk")
DEFINE_CONVERTER_ALIAS(calculate_sepic_inputs, "sepic")
DEFINE_CONVERTER_ALIAS(calculate_advanced_sepic_inputs, "advanced_sepic")
DEFINE_CONVERTER_ALIAS(calculate_zeta_inputs, "zeta")
DEFINE_CONVERTER_ALIAS(calculate_advanced_zeta_inputs, "advanced_zeta")
DEFINE_CONVERTER_ALIAS(calculate_four_switch_buck_boost_inputs, "four_switch_buck_boost")
DEFINE_CONVERTER_ALIAS(calculate_advanced_four_switch_buck_boost_inputs, "advanced_four_switch_buck_boost")
DEFINE_CONVERTER_ALIAS(calculate_weinberg_inputs, "weinberg")
DEFINE_CONVERTER_ALIAS(calculate_advanced_weinberg_inputs, "advanced_weinberg")
DEFINE_CONVERTER_ALIAS(calculate_clllc_inputs, "clllc")
DEFINE_CONVERTER_ALIAS(calculate_advanced_clllc_inputs, "advanced_clllc")
DEFINE_CONVERTER_ALIAS(calculate_vienna_inputs, "vienna")
DEFINE_CONVERTER_ALIAS(calculate_advanced_vienna_inputs, "advanced_vienna")
DEFINE_CONVERTER_ALIAS(calculate_src_inputs, "src")
DEFINE_CONVERTER_ALIAS(calculate_advanced_src_inputs, "advanced_src")
DEFINE_CONVERTER_ALIAS(calculate_llc_inputs, "llc")
DEFINE_CONVERTER_ALIAS(calculate_advanced_llc_inputs, "advanced_llc")
DEFINE_CONVERTER_ALIAS(calculate_cllc_inputs, "cllc")
DEFINE_CONVERTER_ALIAS(calculate_advanced_cllc_inputs, "advanced_cllc")
DEFINE_CONVERTER_ALIAS(calculate_dab_inputs, "dab")
DEFINE_CONVERTER_ALIAS(calculate_advanced_dab_inputs, "advanced_dab")
// PFC has no AdvancedPowerFactorCorrection class in MKF — basic-only by design.
DEFINE_CONVERTER_ALIAS(calculate_pfc_inputs, "power_factor_correction")
DEFINE_CONVERTER_ALIAS(calculate_psfb_inputs, "phase_shifted_full_bridge")
DEFINE_CONVERTER_ALIAS(calculate_advanced_psfb_inputs, "advanced_phase_shifted_full_bridge")
DEFINE_CONVERTER_ALIAS(calculate_pshb_inputs, "phase_shifted_half_bridge")
DEFINE_CONVERTER_ALIAS(calculate_advanced_pshb_inputs, "advanced_phase_shifted_half_bridge")
DEFINE_CONVERTER_ALIAS(calculate_ahb_inputs, "asymmetric_half_bridge")
DEFINE_CONVERTER_ALIAS(calculate_advanced_ahb_inputs, "advanced_asymmetric_half_bridge")

#undef DEFINE_CONVERTER_ALIAS

// ============================================================================
// Common-Mode / Differential-Mode Choke wrappers (parity port from WebLibMKF
// 2026-05). Each wrapper mirrors the same-named function in
// `/home/alf/OpenMagnetics/WebLibMKF/src/libMKF.cpp`. Per-period waveform
// repetition (only used by the display layer) is intentionally dropped —
// the asgard adviser only needs a single period.
// ============================================================================

json calculate_cmc_inputs(json cmcInputsJson) {
    try {
        OpenMagnetics::CommonModeChoke cmcInputs(cmcInputsJson);

        size_t numberOfPeriods = 1;
        if (cmcInputsJson.contains("numberOfPeriods")) {
            numberOfPeriods = cmcInputsJson["numberOfPeriods"].get<size_t>();
        }
        cmcInputs.set_num_periods_to_extract(numberOfPeriods);

        auto inputs = cmcInputs.process();
        json result;
        to_json(result, inputs);

        json diag;
        diag["computedInductance"] = cmcInputs.get_computed_inductance();
        result["cmcDiagnostics"] = diag;
        return result;
    }
    catch (const std::exception& exc) {
        return json{{"error", std::string("calculate_cmc_inputs: ") + exc.what()}};
    }
}

json calculate_advanced_cmc_inputs(json cmcInputsJson) {
    try {
        OpenMagnetics::AdvancedCommonModeChoke cmcInputs(cmcInputsJson);

        size_t numberOfPeriods = 1;
        if (cmcInputsJson.contains("numberOfPeriods")) {
            numberOfPeriods = cmcInputsJson["numberOfPeriods"].get<size_t>();
        }
        cmcInputs.set_num_periods_to_extract(numberOfPeriods);

        auto inputs = cmcInputs.process();
        json result;
        to_json(result, inputs);

        json diag;
        diag["computedInductance"] = cmcInputs.get_computed_inductance();
        result["cmcDiagnostics"] = diag;
        return result;
    }
    catch (const std::exception& exc) {
        return json{{"error", std::string("calculate_advanced_cmc_inputs: ") + exc.what()}};
    }
}

std::string generate_cmc_ngspice_circuit(json cmcInputsJson) {
    try {
        bool isAdvancedCmc = cmcInputsJson.contains("desiredInductance");

        std::unique_ptr<OpenMagnetics::CommonModeChoke> cmcPtr;
        double inductance;
        double frequency = 150000;

        if (isAdvancedCmc) {
            auto advanced = std::make_unique<OpenMagnetics::AdvancedCommonModeChoke>(cmcInputsJson);
            inductance = advanced->get_desired_inductance();
            cmcPtr = std::move(advanced);
        }
        else {
            cmcPtr = std::make_unique<OpenMagnetics::CommonModeChoke>(cmcInputsJson);
            auto designRequirements = cmcPtr->process_design_requirements();
            inductance = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
            if (!(inductance > 0)) {
                throw std::runtime_error("Unable to calculate CMC inductance");
            }
        }

        return cmcPtr->generate_ngspice_circuit(inductance, frequency);
    }
    catch (const std::exception& exc) {
        return "Exception: " + std::string{exc.what()};
    }
}

json simulate_cmc_lisn_waveforms(json cmcInputsJson, double inductance) {
    try {
        OpenMagnetics::CommonModeChoke cmc(cmcInputsJson);
        auto designRequirements = cmc.process_design_requirements();

        std::vector<double> frequencies;
        if (cmcInputsJson.contains("impedancePoints") && cmcInputsJson["impedancePoints"].is_array()) {
            for (const auto& point : cmcInputsJson["impedancePoints"]) {
                if (point.contains("frequency")) {
                    frequencies.push_back(point["frequency"].get<double>());
                }
            }
        }
        if (frequencies.empty()) {
            frequencies.push_back(150000);
        }

        auto waveforms = cmc.simulate_and_extract_waveforms(inductance, frequencies);
        auto operatingPoints = cmc.simulate_and_extract_operating_points(inductance);

        json result;
        json inputs;
        inputs["designRequirements"] = json();
        to_json(inputs["designRequirements"], designRequirements);
        inputs["operatingPoints"] = json::array();
        for (const auto& op : operatingPoints) {
            json opJson;
            to_json(opJson, op);
            inputs["operatingPoints"].push_back(opJson);
        }
        result["inputs"] = inputs;

        result["converterWaveforms"] = json::array();
        for (const auto& wf : waveforms) {
            json cwJson;
            cwJson["frequency"] = wf.frequency;
            cwJson["time"] = wf.time;
            cwJson["inputVoltage"] = wf.inputVoltage;
            cwJson["windingCurrents"] = wf.windingCurrents;
            cwJson["lisnVoltage"] = wf.lisnVoltage;
            cwJson["operatingPointName"] = wf.operatingPointName;
            cwJson["commonModeAttenuation"] = wf.commonModeAttenuation;
            cwJson["commonModeImpedance"] = wf.commonModeImpedance;
            cwJson["theoreticalImpedance"] = wf.theoreticalImpedance;
            result["converterWaveforms"].push_back(cwJson);
        }
        return result;
    }
    catch (const std::exception& exc) {
        return json{{"error", std::string("simulate_cmc_lisn_waveforms: ") + exc.what()}};
    }
}

json simulate_cmc_ideal_waveforms(json cmcInputsJson, double inductance, double parasiticCap_pF, double dvdt_V_ns) {
    try {
        OpenMagnetics::CommonModeChoke cmc(cmcInputsJson);
        auto designRequirements = cmc.process_design_requirements();

        int numberOfPeriods            = cmcInputsJson.value("numberOfPeriods", 2);
        int numberOfSteadyStatePeriods = cmcInputsJson.value("numberOfSteadyStatePeriods", 10);

        auto operatingPoints = cmc.simulate_realistic_cmc(
            inductance, parasiticCap_pF, dvdt_V_ns,
            numberOfPeriods, numberOfSteadyStatePeriods);

        json result;
        json inputs;
        inputs["designRequirements"] = json();
        to_json(inputs["designRequirements"], designRequirements);
        inputs["operatingPoints"] = json::array();
        for (const auto& op : operatingPoints) {
            json opJson;
            to_json(opJson, op);
            inputs["operatingPoints"].push_back(opJson);
        }
        result["inputs"] = inputs;
        result["converterWaveforms"] = json::array();

        cmc.process();
        json diag;
        diag["computedInductance"] = cmc.get_computed_inductance();
        result["cmcDiagnostics"] = diag;
        return result;
    }
    catch (const std::exception& exc) {
        return json{{"error", std::string("simulate_cmc_ideal_waveforms: ") + exc.what()}};
    }
}

json calculate_dmc_inputs(json dmcInputsJson) {
    try {
        OpenMagnetics::DifferentialModeChoke dmcInputs(dmcInputsJson);
        auto inputs = dmcInputs.process();
        json result;
        to_json(result, inputs);

        json diag;
        diag["computedInductance"]      = dmcInputs.get_computed_inductance();
        diag["computedMinFrequency"]    = dmcInputs.get_computed_min_frequency();
        diag["computedMaxFrequency"]    = dmcInputs.get_computed_max_frequency();
        diag["impedanceAtMinFrequency"] = dmcInputs.get_computed_impedance_at_min_freq();
        diag["numberWindings"]          = dmcInputs.get_computed_number_windings();
        result["dmcDiagnostics"] = diag;
        return result;
    }
    catch (const std::exception& exc) {
        return json{{"error", std::string("calculate_dmc_inputs: ") + exc.what()}};
    }
}

json verify_dmc_attenuation(json dmcInputsJson, double inductance, double capacitance) {
    try {
        OpenMagnetics::DifferentialModeChoke dmc(dmcInputsJson);

        std::optional<double> cap = (capacitance > 0) ? std::optional<double>(capacitance) : std::nullopt;
        auto results = dmc.verify_attenuation(inductance, cap);

        json result = json::array();
        for (const auto& r : results) {
            json rJson;
            rJson["frequency"] = r.frequency;
            rJson["requiredAttenuation"] = r.requiredAttenuation;
            rJson["measuredAttenuation"] = r.measuredAttenuation;
            rJson["theoreticalAttenuation"] = r.theoreticalAttenuation;
            rJson["passed"] = r.passed;
            rJson["message"] = r.message;
            result.push_back(rJson);
        }
        return result;
    }
    catch (const std::exception& exc) {
        return json{{"error", std::string("verify_dmc_attenuation: ") + exc.what()}};
    }
}

json propose_dmc_design(json dmcInputsJson) {
    try {
        OpenMagnetics::DifferentialModeChoke dmc(dmcInputsJson);
        return dmc.propose_design();
    }
    catch (const std::exception& exc) {
        return json{{"error", std::string("propose_dmc_design: ") + exc.what()}};
    }
}

json simulate_dmc_waveforms(json dmcInputsJson, double inductance) {
    try {
        OpenMagnetics::DifferentialModeChoke dmc(dmcInputsJson);

        std::vector<double> frequencies;
        auto minimumImpedance = dmc.get_minimum_impedance();
        if (minimumImpedance) {
            for (const auto& imp : *minimumImpedance) {
                frequencies.push_back(imp.get_frequency());
            }
        }
        if (frequencies.empty()) {
            frequencies = {150000, 500000, 1000000, 10000000, 30000000};
        }

        size_t numberOfPeriods = 1;
        if (dmcInputsJson.contains("numberOfPeriods")) {
            numberOfPeriods = dmcInputsJson["numberOfPeriods"].get<size_t>();
        }

        auto waveforms = dmc.simulate_and_extract_waveforms(inductance, frequencies, numberOfPeriods);

        json result = json::array();
        for (const auto& wf : waveforms) {
            json wfJson;
            wfJson["time"] = wf.time;
            wfJson["frequency"] = wf.frequency;
            wfJson["inputVoltage"] = wf.inputVoltage;
            wfJson["outputVoltage"] = wf.outputVoltage;
            wfJson["inductorCurrent"] = wf.inductorCurrent;
            wfJson["operatingPointName"] = wf.operatingPointName;
            wfJson["dmAttenuation"] = wf.dmAttenuation;
            result.push_back(wfJson);
        }
        return result;
    }
    catch (const std::exception& exc) {
        return json{{"error", std::string("simulate_dmc_waveforms: ") + exc.what()}};
    }
}

std::string generate_dmc_ngspice_circuit(json dmcInputsJson) {
    try {
        OpenMagnetics::DifferentialModeChoke dmc(dmcInputsJson);

        // Inductance: prefer explicit `minimumInductance` from the wizard,
        // fall back to propose_design()'s proposal. Mirrors WASM
        // generate_dmc_ngspice_circuit at libMKF.cpp:8113-8136.
        double inductance;
        if (dmcInputsJson.contains("minimumInductance") && dmcInputsJson["minimumInductance"].is_number()) {
            inductance = dmcInputsJson["minimumInductance"].get<double>();
        }
        else {
            auto proposal = dmc.propose_design();
            if (proposal.contains("inductance") && proposal["inductance"].is_number()) {
                inductance = proposal["inductance"].get<double>();
            }
            else if (proposal.contains("minimumInductance") && proposal["minimumInductance"].is_number()) {
                inductance = proposal["minimumInductance"].get<double>();
            }
            else {
                throw std::runtime_error("Unable to determine DMC inductance for SPICE generation");
            }
        }

        constexpr double kDmcSpiceTestFrequency = 150000.0;
        return dmc.generate_ngspice_circuit(inductance, kDmcSpiceTestFrequency);
    }
    catch (const std::exception& exc) {
        return "Exception: " + std::string{exc.what()};
    }
}

// ============================================================================
// Phase B simulate_*_ideal_waveforms functions (parity port from WebLibMKF).
//
// Each function runs ngspice simulation to produce per-period waveforms and
// operating points, plus topology-specific diagnostics. Adapted for PyMKF:
//   - json in / json out (not std::string)
//   - json{{"error", ...}} for errors (not "Exception: ")
//   - No json::parse() / .dump(4)
//   - No WASM timing debug prints
// ============================================================================

// ------ Flyback (isolated, transformer-based) ------
json simulate_flyback_ideal_waveforms(json inputsJson) {
    try {
        bool isAdvanced = inputsJson.contains("desiredInductance");

        OpenMagnetics::DesignRequirements designRequirements;
        std::vector<double> turnsRatios;
        double magnetizingInductance;

        std::unique_ptr<OpenMagnetics::Flyback> flybackPtr;

        if (isAdvanced) {
            auto advPtr = std::make_unique<OpenMagnetics::AdvancedFlyback>(inputsJson);
            magnetizingInductance = advPtr->get_desired_inductance();
            turnsRatios = advPtr->get_desired_turns_ratios();

            designRequirements.get_mutable_turns_ratios().clear();
            for (auto tr : turnsRatios) {
                OpenMagnetics::DimensionWithTolerance trTol;
                trTol.set_nominal(tr);
                designRequirements.get_mutable_turns_ratios().push_back(trTol);
            }
            OpenMagnetics::DimensionWithTolerance inductanceTol;
            inductanceTol.set_nominal(magnetizingInductance);
            designRequirements.set_magnetizing_inductance(inductanceTol);
            std::vector<OpenMagnetics::IsolationSide> isolationSides;
            for (size_t w = 0; w < turnsRatios.size() + 1; ++w) {
                isolationSides.push_back(OpenMagnetics::get_isolation_side_from_index(w));
            }
            designRequirements.set_isolation_sides(isolationSides);
            designRequirements.set_topology(MAS::Topologies::FLYBACK_CONVERTER);

            flybackPtr = std::move(advPtr);
        } else {
            flybackPtr = std::make_unique<OpenMagnetics::Flyback>(inputsJson);
            designRequirements = flybackPtr->process_design_requirements();
            for (const auto& tr : designRequirements.get_turns_ratios()) {
                turnsRatios.push_back(tr.get_nominal().value());
            }
            magnetizingInductance = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
            if (!(magnetizingInductance > 0)) {
                throw std::runtime_error("Unable to calculate magnetizing inductance");
            }
        }

        OpenMagnetics::NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("ngspice simulation is required but ngspice is not available");
        }

        size_t numberOfPeriods = inputsJson.value("numberOfPeriods", 2);
        size_t numberOfSteadyStatePeriods = inputsJson.value("numberOfSteadyStatePeriods", 5);
        flybackPtr->set_num_periods_to_extract(numberOfPeriods);
        flybackPtr->set_num_steady_state_periods(numberOfSteadyStatePeriods);

        auto topologyWaveforms = flybackPtr->simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance, numberOfPeriods);
        auto operatingPoints = flybackPtr->simulate_and_extract_operating_points(turnsRatios, magnetizingInductance, numberOfPeriods);

        json result;
        json inputs;
        inputs["designRequirements"] = json();
        to_json(inputs["designRequirements"], designRequirements);
        inputs["operatingPoints"] = json::array();
        for (const auto& op : operatingPoints) { json j; to_json(j, op); inputs["operatingPoints"].push_back(j); }
        result["inputs"] = inputs;

        result["converterWaveforms"] = json::array();
        for (const auto& tw : topologyWaveforms) { json j; to_json(j, tw); result["converterWaveforms"].push_back(j); }

        flybackPtr->process_operating_points(turnsRatios, magnetizingInductance);
        {
            json diag;
            const auto& names = flybackPtr->get_per_op_name();
            const auto& dC   = flybackPtr->get_per_op_duty_cycle();
            const auto& fsw  = flybackPtr->get_per_op_switching_frequency();
            const auto& iAvg = flybackPtr->get_per_op_primary_average_current();
            const auto& iPP  = flybackPtr->get_per_op_primary_peak_to_peak();
            const auto& iPk  = flybackPtr->get_per_op_primary_peak_current();
            const auto& iSec = flybackPtr->get_per_op_secondary_peak_current();
            const auto& ccm  = flybackPtr->get_per_op_is_ccm();
            diag["dutyCycle"]             = dC.empty()   ? flybackPtr->get_last_duty_cycle()              : dC.front();
            diag["switchingFrequency"]    = fsw.empty()  ? 0.0                                            : fsw.front();
            diag["primaryAverageCurrent"] = iAvg.empty() ? 0.0                                            : iAvg.front();
            diag["primaryPeakToPeak"]     = iPP.empty()  ? 0.0                                            : iPP.front();
            diag["primaryPeakCurrent"]    = iPk.empty()  ? 0.0                                            : iPk.front();
            diag["secondaryPeakCurrent"]  = iSec.empty() ? 0.0                                            : iSec.front();
            diag["isCcm"]                 = ccm.empty()  ? false                                          : (bool)ccm.front();
            json perOp = json::array();
            for (size_t i = 0; i < dC.size(); ++i) {
                json row;
                row["operatingPointName"]    = (i < names.size()) ? names[i] : ("OP " + std::to_string(i));
                row["dutyCycle"]             = dC[i];
                row["switchingFrequency"]    = fsw[i];
                row["primaryAverageCurrent"] = iAvg[i];
                row["primaryPeakToPeak"]     = iPP[i];
                row["primaryPeakCurrent"]    = iPk[i];
                row["secondaryPeakCurrent"]  = iSec[i];
                row["isCcm"]                 = (bool)ccm[i];
                perOp.push_back(row);
            }
            diag["perOp"] = perOp;
            result["flybackDiagnostics"] = diag;
        }

        return result;
    }
    catch (const std::exception& exc) {
        return json{{"error", std::string("simulate_flyback_ideal_waveforms: ") + exc.what()}};
    }
}

json simulate_flyback_with_magnetic(json inputsJson, json magneticJson) {
    try {
        OpenMagnetics::AdvancedFlyback converter(inputsJson);

        double magnetizingInductance = converter.get_desired_inductance();
        std::vector<double> turnsRatios = converter.get_desired_turns_ratios();

        OpenMagnetics::NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("ngspice simulation is required but ngspice is not available");
        }

        size_t numberOfPeriods = inputsJson.value("numberOfPeriods", 2);
        size_t numberOfSteadyStatePeriods = inputsJson.value("numberOfSteadyStatePeriods", 5);
        converter.set_num_periods_to_extract(numberOfPeriods);
        converter.set_num_steady_state_periods(numberOfSteadyStatePeriods);

        auto topologyWaveforms = converter.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance, numberOfPeriods);
        auto operatingPoints = converter.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance, numberOfPeriods);

        json result;
        json inputs;
        auto designRequirements = converter.process_design_requirements();
        inputs["designRequirements"] = json();
        to_json(inputs["designRequirements"], designRequirements);
        inputs["operatingPoints"] = json::array();
        for (const auto& op : operatingPoints) { json j; to_json(j, op); inputs["operatingPoints"].push_back(j); }
        result["inputs"] = inputs;

        result["converterWaveforms"] = json::array();
        for (const auto& tw : topologyWaveforms) { json j; to_json(j, tw); result["converterWaveforms"].push_back(j); }

        return result;
    }
    catch (const std::exception& exc) {
        return json{{"error", std::string("simulate_flyback_with_magnetic: ") + exc.what()}};
    }
}

// ------ Helper macro for non-isolated single-inductor topologies (Buck/Boost pattern) ------
#define DEFINE_SIMULATE_NON_ISOLATED(func_name, BaseType, AdvancedType, TopologyEnum, diagKey, \
    DIAG_FLAT, DIAG_PEROP) \
json func_name(json inputsJson) { \
    try { \
        bool isAdvanced = inputsJson.contains("desiredInductance"); \
        OpenMagnetics::DesignRequirements designRequirements; \
        double inductance; \
        std::unique_ptr<BaseType> ptr; \
        if (isAdvanced) { \
            auto advPtr = std::make_unique<AdvancedType>(inputsJson); \
            inductance = advPtr->get_desired_inductance(); \
            designRequirements.get_mutable_turns_ratios().clear(); \
            OpenMagnetics::DimensionWithTolerance inductanceTol; \
            inductanceTol.set_nominal(inductance); \
            designRequirements.set_magnetizing_inductance(inductanceTol); \
            std::vector<OpenMagnetics::IsolationSide> isolationSides; \
            isolationSides.push_back(OpenMagnetics::get_isolation_side_from_index(0)); \
            designRequirements.set_isolation_sides(isolationSides); \
            designRequirements.set_topology(TopologyEnum); \
            ptr = std::move(advPtr); \
        } else { \
            ptr = std::make_unique<BaseType>(inputsJson); \
            designRequirements = ptr->process_design_requirements(); \
            inductance = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance()); \
            if (!(inductance > 0)) { \
                throw std::runtime_error("Unable to calculate inductance"); \
            } \
        } \
        OpenMagnetics::NgspiceRunner runner; \
        if (!runner.is_available()) { \
            throw std::runtime_error("ngspice simulation is required but ngspice is not available"); \
        } \
        size_t numberOfPeriods = inputsJson.value("numberOfPeriods", 2); \
        size_t numberOfSteadyStatePeriods = inputsJson.value("numberOfSteadyStatePeriods", 5); \
        ptr->set_num_periods_to_extract(numberOfPeriods); \
        ptr->set_num_steady_state_periods(numberOfSteadyStatePeriods); \
        auto topologyWaveforms = ptr->simulate_and_extract_topology_waveforms(inductance); \
        auto operatingPoints = ptr->simulate_and_extract_operating_points(inductance); \
        json result; \
        json inputs; \
        inputs["designRequirements"] = json(); \
        to_json(inputs["designRequirements"], designRequirements); \
        inputs["operatingPoints"] = json::array(); \
        for (const auto& op : operatingPoints) { json j; to_json(j, op); inputs["operatingPoints"].push_back(j); } \
        result["inputs"] = inputs; \
        result["converterWaveforms"] = json::array(); \
        for (const auto& tw : topologyWaveforms) { json j; to_json(j, tw); result["converterWaveforms"].push_back(j); } \
        ptr->process(); \
        { \
            json diag; \
            DIAG_FLAT \
            json perOp = json::array(); \
            DIAG_PEROP \
            diag["perOp"] = perOp; \
            result[diagKey] = diag; \
        } \
        return result; \
    } \
    catch (const std::exception& exc) { \
        return json{{"error", std::string(#func_name ": ") + exc.what()}}; \
    } \
}

// Buck diagnostics
#define BUCK_BOOST_DIAG_FLAT \
    const auto& names = ptr->get_per_op_name(); \
    const auto& dC    = ptr->get_per_op_duty_cycle(); \
    const auto& iAvg  = ptr->get_per_op_inductor_average_current(); \
    const auto& iPP   = ptr->get_per_op_inductor_peak_to_peak(); \
    const auto& iPk   = ptr->get_per_op_peak_inductor_current(); \
    const auto& ccm   = ptr->get_per_op_is_ccm(); \
    const auto& cRat  = ptr->get_per_op_conduction_ratio(); \
    diag["dutyCycle"]              = dC.empty()   ? ptr->get_last_duty_cycle()              : dC.front(); \
    diag["inductorAverageCurrent"] = iAvg.empty() ? ptr->get_last_inductor_average_current() : iAvg.front(); \
    diag["inductorPeakToPeak"]     = iPP.empty()  ? ptr->get_last_inductor_peak_to_peak()    : iPP.front(); \
    diag["peakInductorCurrent"]    = iPk.empty()  ? ptr->get_last_peak_inductor_current()    : iPk.front(); \
    diag["conductionRatio"]        = cRat.empty() ? ptr->get_last_conduction_ratio()         : cRat.front(); \
    diag["isCcm"]                  = ccm.empty()  ? ptr->get_last_is_ccm()                   : (bool)ccm.front();

#define BUCK_BOOST_DIAG_PEROP \
    for (size_t i = 0; i < dC.size(); ++i) { \
        json row; \
        row["operatingPointName"]    = (i < names.size()) ? names[i] : ("OP " + std::to_string(i)); \
        row["dutyCycle"]             = dC[i]; \
        row["inductorAverageCurrent"] = iAvg[i]; \
        row["inductorPeakToPeak"]    = iPP[i]; \
        row["peakInductorCurrent"]   = iPk[i]; \
        row["conductionRatio"]       = cRat[i]; \
        row["isCcm"]                 = (bool)ccm[i]; \
        perOp.push_back(row); \
    }

DEFINE_SIMULATE_NON_ISOLATED(simulate_buck_ideal_waveforms, OpenMagnetics::Buck, OpenMagnetics::AdvancedBuck,
    MAS::Topologies::BUCK_CONVERTER, "buckDiagnostics", BUCK_BOOST_DIAG_FLAT, BUCK_BOOST_DIAG_PEROP)

DEFINE_SIMULATE_NON_ISOLATED(simulate_boost_ideal_waveforms, OpenMagnetics::Boost, OpenMagnetics::AdvancedBoost,
    MAS::Topologies::BOOST_CONVERTER, "boostDiagnostics", BUCK_BOOST_DIAG_FLAT, BUCK_BOOST_DIAG_PEROP)

#undef BUCK_BOOST_DIAG_FLAT
#undef BUCK_BOOST_DIAG_PEROP

// ------ Coupling topologies (Sepic/Cuk/Zeta): single inductance L1, coupling diagnostics ------

#define COUPLING_DIAG_FLAT \
    const auto& names = ptr->get_per_op_name(); \
    const auto& v_duty_cycle = ptr->get_per_op_duty_cycle(); \
    const auto& v_conversion_ratio = ptr->get_per_op_conversion_ratio(); \
    const auto& v_coupling_cap_voltage = ptr->get_per_op_coupling_cap_voltage(); \
    const auto& v_input_inductor_average = ptr->get_per_op_input_inductor_average(); \
    const auto& v_output_inductor_average = ptr->get_per_op_output_inductor_average(); \
    const auto& v_input_inductor_ripple = ptr->get_per_op_input_inductor_ripple(); \
    const auto& v_output_inductor_ripple = ptr->get_per_op_output_inductor_ripple(); \
    const auto& v_switch_peak_voltage = ptr->get_per_op_switch_peak_voltage(); \
    const auto& v_switch_peak_current = ptr->get_per_op_switch_peak_current(); \
    const auto& v_diode_peak_reverse_voltage = ptr->get_per_op_diode_peak_reverse_voltage(); \
    const auto& v_diode_peak_current = ptr->get_per_op_diode_peak_current(); \
    const auto& v_coupling_cap_rms_current = ptr->get_per_op_coupling_cap_rms_current(); \
    const auto& v_is_ccm = ptr->get_per_op_is_ccm(); \
    diag["dutyCycle"] = v_duty_cycle.empty() ? ptr->get_last_duty_cycle() : v_duty_cycle.front(); \
    diag["conversionRatio"] = v_conversion_ratio.empty() ? ptr->get_last_conversion_ratio() : v_conversion_ratio.front(); \
    diag["couplingCapVoltage"] = v_coupling_cap_voltage.empty() ? ptr->get_last_coupling_cap_voltage() : v_coupling_cap_voltage.front(); \
    diag["inputInductorAverage"] = v_input_inductor_average.empty() ? ptr->get_last_input_inductor_average() : v_input_inductor_average.front(); \
    diag["outputInductorAverage"] = v_output_inductor_average.empty() ? ptr->get_last_output_inductor_average() : v_output_inductor_average.front(); \
    diag["inputInductorRipple"] = v_input_inductor_ripple.empty() ? ptr->get_last_input_inductor_ripple() : v_input_inductor_ripple.front(); \
    diag["outputInductorRipple"] = v_output_inductor_ripple.empty() ? ptr->get_last_output_inductor_ripple() : v_output_inductor_ripple.front(); \
    diag["switchPeakVoltage"] = v_switch_peak_voltage.empty() ? ptr->get_last_switch_peak_voltage() : v_switch_peak_voltage.front(); \
    diag["switchPeakCurrent"] = v_switch_peak_current.empty() ? ptr->get_last_switch_peak_current() : v_switch_peak_current.front(); \
    diag["diodePeakReverseVoltage"] = v_diode_peak_reverse_voltage.empty() ? ptr->get_last_diode_peak_reverse_voltage() : v_diode_peak_reverse_voltage.front(); \
    diag["diodePeakCurrent"] = v_diode_peak_current.empty() ? ptr->get_last_diode_peak_current() : v_diode_peak_current.front(); \
    diag["couplingCapRmsCurrent"] = v_coupling_cap_rms_current.empty() ? ptr->get_last_coupling_cap_rms_current() : v_coupling_cap_rms_current.front(); \
    diag["isCcm"] = v_is_ccm.empty() ? ptr->get_last_is_ccm() : (bool)v_is_ccm.front();

#define COUPLING_DIAG_PEROP_BASE \
    for (size_t i = 0; i < v_duty_cycle.size(); ++i) { \
        json row; \
        row["operatingPointName"] = (i < names.size()) ? names[i] : ("OP " + std::to_string(i)); \
        row["dutyCycle"] = v_duty_cycle[i]; \
        row["conversionRatio"] = v_conversion_ratio[i]; \
        row["couplingCapVoltage"] = v_coupling_cap_voltage[i]; \
        row["inputInductorAverage"] = v_input_inductor_average[i]; \
        row["outputInductorAverage"] = v_output_inductor_average[i]; \
        row["inputInductorRipple"] = v_input_inductor_ripple[i]; \
        row["outputInductorRipple"] = v_output_inductor_ripple[i]; \
        row["switchPeakVoltage"] = v_switch_peak_voltage[i]; \
        row["switchPeakCurrent"] = v_switch_peak_current[i]; \
        row["diodePeakReverseVoltage"] = v_diode_peak_reverse_voltage[i]; \
        row["diodePeakCurrent"] = v_diode_peak_current[i]; \
        row["couplingCapRmsCurrent"] = v_coupling_cap_rms_current[i]; \
        row["isCcm"] = (bool)v_is_ccm[i];

// Sepic: sizedCs, sizedCo
json simulate_sepic_ideal_waveforms(json inputsJson) {
    try {
        bool isAdvanced = inputsJson.contains("desiredInductance");
        OpenMagnetics::DesignRequirements designRequirements;
        double inductanceL1;
        std::unique_ptr<OpenMagnetics::Sepic> ptr;
        if (isAdvanced) {
            auto advPtr = std::make_unique<OpenMagnetics::AdvancedSepic>(inputsJson);
            inductanceL1 = advPtr->get_desired_inductance();
            designRequirements.get_mutable_turns_ratios().clear();
            OpenMagnetics::DimensionWithTolerance inductanceTol;
            inductanceTol.set_nominal(inductanceL1);
            designRequirements.set_magnetizing_inductance(inductanceTol);
            std::vector<OpenMagnetics::IsolationSide> isolationSides;
            isolationSides.push_back(OpenMagnetics::get_isolation_side_from_index(0));
            designRequirements.set_isolation_sides(isolationSides);
            designRequirements.set_topology(MAS::Topologies::SEPIC_CONVERTER);
            ptr = std::move(advPtr);
        } else {
            ptr = std::make_unique<OpenMagnetics::Sepic>(inputsJson);
            designRequirements = ptr->process_design_requirements();
            inductanceL1 = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
            if (!(inductanceL1 > 0)) throw std::runtime_error("Unable to calculate inductance");
        }
        OpenMagnetics::NgspiceRunner runner;
        if (!runner.is_available()) throw std::runtime_error("ngspice is not available");
        size_t numberOfPeriods = inputsJson.value("numberOfPeriods", 2);
        size_t numberOfSteadyStatePeriods = inputsJson.value("numberOfSteadyStatePeriods", 5);
        ptr->set_num_periods_to_extract(numberOfPeriods);
        ptr->set_num_steady_state_periods(numberOfSteadyStatePeriods);
        auto topologyWaveforms = ptr->simulate_and_extract_topology_waveforms(inductanceL1, numberOfPeriods);
        auto operatingPoints = ptr->simulate_and_extract_operating_points(inductanceL1);
        json result;
        json inputs;
        inputs["designRequirements"] = json();
        to_json(inputs["designRequirements"], designRequirements);
        inputs["operatingPoints"] = json::array();
        for (const auto& op : operatingPoints) { json j; to_json(j, op); inputs["operatingPoints"].push_back(j); }
        result["inputs"] = inputs;
        result["converterWaveforms"] = json::array();
        for (const auto& tw : topologyWaveforms) { json j; to_json(j, tw); result["converterWaveforms"].push_back(j); }
        ptr->process();
        {
            json diag;
            COUPLING_DIAG_FLAT
            const auto& v_sized_cs = ptr->get_per_op_sized_cs();
            const auto& v_sized_co = ptr->get_per_op_sized_co();
            diag["sizedCs"] = v_sized_cs.empty() ? ptr->get_last_sized_cs() : v_sized_cs.front();
            diag["sizedCo"] = v_sized_co.empty() ? ptr->get_last_sized_co() : v_sized_co.front();
            json perOp = json::array();
            COUPLING_DIAG_PEROP_BASE
                row["sizedCs"] = v_sized_cs[i];
                row["sizedCo"] = v_sized_co[i];
                perOp.push_back(row);
            }
            diag["perOp"] = perOp;
            result["sepicDiagnostics"] = diag;
        }
        return result;
    }
    catch (const std::exception& exc) {
        return json{{"error", std::string("simulate_sepic_ideal_waveforms: ") + exc.what()}};
    }
}

// Cuk: sizedCa, sizedCb, sizedCo, rhpZeroFrequency
json simulate_cuk_ideal_waveforms(json inputsJson) {
    try {
        bool isAdvanced = inputsJson.contains("desiredInductance");
        OpenMagnetics::DesignRequirements designRequirements;
        double inductanceL1;
        std::unique_ptr<OpenMagnetics::Cuk> ptr;
        if (isAdvanced) {
            auto advPtr = std::make_unique<OpenMagnetics::AdvancedCuk>(inputsJson);
            inductanceL1 = advPtr->get_desired_inductance();
            designRequirements.get_mutable_turns_ratios().clear();
            OpenMagnetics::DimensionWithTolerance inductanceTol;
            inductanceTol.set_nominal(inductanceL1);
            designRequirements.set_magnetizing_inductance(inductanceTol);
            std::vector<OpenMagnetics::IsolationSide> isolationSides;
            isolationSides.push_back(OpenMagnetics::get_isolation_side_from_index(0));
            designRequirements.set_isolation_sides(isolationSides);
            designRequirements.set_topology(MAS::Topologies::CUK_CONVERTER);
            ptr = std::move(advPtr);
        } else {
            ptr = std::make_unique<OpenMagnetics::Cuk>(inputsJson);
            designRequirements = ptr->process_design_requirements();
            inductanceL1 = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
            if (!(inductanceL1 > 0)) throw std::runtime_error("Unable to calculate inductance");
        }
        OpenMagnetics::NgspiceRunner runner;
        if (!runner.is_available()) throw std::runtime_error("ngspice is not available");
        size_t numberOfPeriods = inputsJson.value("numberOfPeriods", 2);
        size_t numberOfSteadyStatePeriods = inputsJson.value("numberOfSteadyStatePeriods", 5);
        ptr->set_num_periods_to_extract(numberOfPeriods);
        ptr->set_num_steady_state_periods(numberOfSteadyStatePeriods);
        auto operatingPoints = ptr->simulate_and_extract_operating_points(inductanceL1);
        auto topologyWaveforms = ptr->simulate_and_extract_topology_waveforms(inductanceL1, numberOfPeriods);
        json result;
        json inputs;
        inputs["designRequirements"] = json();
        to_json(inputs["designRequirements"], designRequirements);
        inputs["operatingPoints"] = json::array();
        for (const auto& op : operatingPoints) { json j; to_json(j, op); inputs["operatingPoints"].push_back(j); }
        result["inputs"] = inputs;
        result["converterWaveforms"] = json::array();
        for (const auto& tw : topologyWaveforms) { json j; to_json(j, tw); result["converterWaveforms"].push_back(j); }
        ptr->process();
        {
            json diag;
            COUPLING_DIAG_FLAT
            const auto& v_sized_ca = ptr->get_per_op_sized_ca();
            const auto& v_sized_cb = ptr->get_per_op_sized_cb();
            const auto& v_sized_co = ptr->get_per_op_sized_co();
            const auto& v_rhp_zero_frequency = ptr->get_per_op_rhp_zero_frequency();
            diag["sizedCa"] = v_sized_ca.empty() ? ptr->get_last_sized_ca() : v_sized_ca.front();
            diag["sizedCb"] = v_sized_cb.empty() ? ptr->get_last_sized_cb() : v_sized_cb.front();
            diag["sizedCo"] = v_sized_co.empty() ? ptr->get_last_sized_co() : v_sized_co.front();
            diag["rhpZeroFrequency"] = v_rhp_zero_frequency.empty() ? ptr->get_last_rhp_zero_frequency() : v_rhp_zero_frequency.front();
            json perOp = json::array();
            COUPLING_DIAG_PEROP_BASE
                row["sizedCa"] = v_sized_ca[i];
                row["sizedCb"] = v_sized_cb[i];
                row["sizedCo"] = v_sized_co[i];
                row["rhpZeroFrequency"] = v_rhp_zero_frequency[i];
                perOp.push_back(row);
            }
            diag["perOp"] = perOp;
            result["cukDiagnostics"] = diag;
        }
        return result;
    }
    catch (const std::exception& exc) {
        return json{{"error", std::string("simulate_cuk_ideal_waveforms: ") + exc.what()}};
    }
}

// Zeta: sizedCc, sizedCo, outputVoltageRipple, inputCurrentRipple
json simulate_zeta_ideal_waveforms(json inputsJson) {
    try {
        bool isAdvanced = inputsJson.contains("desiredInductance");
        OpenMagnetics::DesignRequirements designRequirements;
        double inductanceL1;
        std::unique_ptr<OpenMagnetics::Zeta> ptr;
        if (isAdvanced) {
            auto advPtr = std::make_unique<OpenMagnetics::AdvancedZeta>(inputsJson);
            inductanceL1 = advPtr->get_desired_inductance();
            designRequirements.get_mutable_turns_ratios().clear();
            OpenMagnetics::DimensionWithTolerance inductanceTol;
            inductanceTol.set_nominal(inductanceL1);
            designRequirements.set_magnetizing_inductance(inductanceTol);
            std::vector<OpenMagnetics::IsolationSide> isolationSides;
            isolationSides.push_back(OpenMagnetics::get_isolation_side_from_index(0));
            designRequirements.set_isolation_sides(isolationSides);
            designRequirements.set_topology(MAS::Topologies::ZETA_CONVERTER);
            ptr = std::move(advPtr);
        } else {
            ptr = std::make_unique<OpenMagnetics::Zeta>(inputsJson);
            designRequirements = ptr->process_design_requirements();
            inductanceL1 = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
            if (!(inductanceL1 > 0)) throw std::runtime_error("Unable to calculate inductance");
        }
        OpenMagnetics::NgspiceRunner runner;
        if (!runner.is_available()) throw std::runtime_error("ngspice is not available");
        size_t numberOfPeriods = inputsJson.value("numberOfPeriods", 2);
        size_t numberOfSteadyStatePeriods = inputsJson.value("numberOfSteadyStatePeriods", 5);
        ptr->set_num_periods_to_extract(numberOfPeriods);
        ptr->set_num_steady_state_periods(numberOfSteadyStatePeriods);
        auto operatingPoints = ptr->simulate_and_extract_operating_points(inductanceL1);
        auto topologyWaveforms = ptr->simulate_and_extract_topology_waveforms(inductanceL1, numberOfPeriods);
        json result;
        json inputs;
        inputs["designRequirements"] = json();
        to_json(inputs["designRequirements"], designRequirements);
        inputs["operatingPoints"] = json::array();
        for (const auto& op : operatingPoints) { json j; to_json(j, op); inputs["operatingPoints"].push_back(j); }
        result["inputs"] = inputs;
        result["converterWaveforms"] = json::array();
        for (const auto& tw : topologyWaveforms) { json j; to_json(j, tw); result["converterWaveforms"].push_back(j); }
        ptr->process();
        {
            json diag;
            COUPLING_DIAG_FLAT
            const auto& v_sized_cc = ptr->get_per_op_sized_cc();
            const auto& v_sized_co = ptr->get_per_op_sized_co();
            const auto& v_output_voltage_ripple = ptr->get_per_op_output_voltage_ripple();
            const auto& v_input_current_ripple = ptr->get_per_op_input_current_ripple();
            diag["sizedCc"] = v_sized_cc.empty() ? ptr->get_last_sized_cc() : v_sized_cc.front();
            diag["sizedCo"] = v_sized_co.empty() ? ptr->get_last_sized_co() : v_sized_co.front();
            diag["outputVoltageRipple"] = v_output_voltage_ripple.empty() ? ptr->get_last_output_voltage_ripple() : v_output_voltage_ripple.front();
            diag["inputCurrentRipple"] = v_input_current_ripple.empty() ? ptr->get_last_input_current_ripple() : v_input_current_ripple.front();
            json perOp = json::array();
            COUPLING_DIAG_PEROP_BASE
                row["sizedCc"] = v_sized_cc[i];
                row["sizedCo"] = v_sized_co[i];
                row["outputVoltageRipple"] = v_output_voltage_ripple[i];
                row["inputCurrentRipple"] = v_input_current_ripple[i];
                perOp.push_back(row);
            }
            diag["perOp"] = perOp;
            result["zetaDiagnostics"] = diag;
        }
        return result;
    }
    catch (const std::exception& exc) {
        return json{{"error", std::string("simulate_zeta_ideal_waveforms: ") + exc.what()}};
    }
}

#undef COUPLING_DIAG_FLAT
#undef COUPLING_DIAG_PEROP_BASE

// ------ FourSwitchBuckBoost (single inductance, FSBB diagnostics) ------
json simulate_four_switch_buck_boost_ideal_waveforms(json inputsJson) {
    try {
        bool isAdvanced = inputsJson.contains("desiredInductance");
        OpenMagnetics::DesignRequirements designRequirements;
        double inductance;
        std::unique_ptr<OpenMagnetics::FourSwitchBuckBoost> ptr;
        if (isAdvanced) {
            auto advPtr = std::make_unique<OpenMagnetics::AdvancedFourSwitchBuckBoost>(inputsJson);
            inductance = advPtr->get_desired_inductance();
            designRequirements.get_mutable_turns_ratios().clear();
            OpenMagnetics::DimensionWithTolerance inductanceTol;
            inductanceTol.set_nominal(inductance);
            designRequirements.set_magnetizing_inductance(inductanceTol);
            std::vector<OpenMagnetics::IsolationSide> isolationSides;
            isolationSides.push_back(OpenMagnetics::get_isolation_side_from_index(0));
            designRequirements.set_isolation_sides(isolationSides);
            designRequirements.set_topology(MAS::Topologies::FOUR_SWITCH_BUCK_BOOST_CONVERTER);
            ptr = std::move(advPtr);
        } else {
            ptr = std::make_unique<OpenMagnetics::FourSwitchBuckBoost>(inputsJson);
            designRequirements = ptr->process_design_requirements();
            inductance = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
            if (!(inductance > 0)) throw std::runtime_error("Unable to calculate inductance");
        }
        OpenMagnetics::NgspiceRunner runner;
        if (!runner.is_available()) throw std::runtime_error("ngspice is not available");
        size_t numberOfPeriods = inputsJson.value("numberOfPeriods", 2);
        size_t numberOfSteadyStatePeriods = inputsJson.value("numberOfSteadyStatePeriods", 5);
        ptr->set_num_periods_to_extract(numberOfPeriods);
        ptr->set_num_steady_state_periods(numberOfSteadyStatePeriods);
        auto operatingPoints = ptr->simulate_and_extract_operating_points(inductance);
        auto topologyWaveforms = ptr->simulate_and_extract_topology_waveforms(inductance, numberOfPeriods);
        json result;
        json inputs;
        inputs["designRequirements"] = json();
        to_json(inputs["designRequirements"], designRequirements);
        inputs["operatingPoints"] = json::array();
        for (const auto& op : operatingPoints) { json j; to_json(j, op); inputs["operatingPoints"].push_back(j); }
        result["inputs"] = inputs;
        result["converterWaveforms"] = json::array();
        for (const auto& tw : topologyWaveforms) { json j; to_json(j, tw); result["converterWaveforms"].push_back(j); }
        ptr->process();
        {
            json diag;
            const auto& names = ptr->get_per_op_name();
            const auto& v_inductor_average_current = ptr->get_per_op_inductor_average_current();
            const auto& v_sized_output_capacitance = ptr->get_per_op_sized_output_capacitance();
            diag["inductorAverageCurrent"] = v_inductor_average_current.empty() ? ptr->get_last_inductor_average_current() : v_inductor_average_current.front();
            diag["sizedOutputCapacitance"] = v_sized_output_capacitance.empty() ? ptr->get_last_sized_output_capacitance() : v_sized_output_capacitance.front();
            json perOp = json::array();
            for (size_t i = 0; i < v_inductor_average_current.size(); ++i) {
                json row;
                row["operatingPointName"] = (i < names.size()) ? names[i] : ("OP " + std::to_string(i));
                row["inductorAverageCurrent"] = v_inductor_average_current[i];
                row["sizedOutputCapacitance"] = v_sized_output_capacitance[i];
                perOp.push_back(row);
            }
            diag["perOp"] = perOp;
            result["fourSwitchBuckBoostDiagnostics"] = diag;
        }
        return result;
    }
    catch (const std::exception& exc) {
        return json{{"error", std::string("simulate_four_switch_buck_boost_ideal_waveforms: ") + exc.what()}};
    }
}

// ------ Helper macro for isolated topologies with turnsRatios + magnetizingInductance ------
// Used by Forward, TwoSwitchForward, ActiveClampForward, PushPull, IsolatedBuck, IsolatedBuckBoost, Weinberg
#define DEFINE_SIMULATE_ISOLATED(func_name, BaseType, AdvancedType, TopologyEnum, diagKey, \
    SIM_WAVEFORM_CALL, SIM_OP_CALL, DIAG_BLOCK) \
json func_name(json inputsJson) { \
    try { \
        bool isAdvanced = inputsJson.contains("desiredInductance"); \
        size_t numberOfPeriods = inputsJson.value("numberOfPeriods", 2); \
        size_t numberOfSteadyStatePeriods = inputsJson.value("numberOfSteadyStatePeriods", 5); \
        OpenMagnetics::DesignRequirements designRequirements; \
        double magnetizingInductance; \
        std::vector<double> turnsRatios; \
        std::unique_ptr<BaseType> ptr; \
        if (isAdvanced) { \
            auto advPtr = std::make_unique<AdvancedType>(inputsJson); \
            magnetizingInductance = advPtr->get_desired_inductance(); \
            turnsRatios = advPtr->get_desired_turns_ratios(); \
            designRequirements.get_mutable_turns_ratios().clear(); \
            for (auto tr : turnsRatios) { \
                OpenMagnetics::DimensionWithTolerance trTol; \
                trTol.set_nominal(tr); \
                designRequirements.get_mutable_turns_ratios().push_back(trTol); \
            } \
            OpenMagnetics::DimensionWithTolerance inductanceTol; \
            inductanceTol.set_nominal(magnetizingInductance); \
            designRequirements.set_magnetizing_inductance(inductanceTol); \
            designRequirements.set_topology(TopologyEnum); \
            ptr = std::move(advPtr); \
        } else { \
            ptr = std::make_unique<BaseType>(inputsJson); \
            designRequirements = ptr->process_design_requirements(); \
            for (const auto& tr : designRequirements.get_turns_ratios()) { \
                if (tr.get_nominal()) turnsRatios.push_back(tr.get_nominal().value()); \
            } \
            magnetizingInductance = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance()); \
            if (!(magnetizingInductance > 0)) throw std::runtime_error("Unable to calculate inductance"); \
        } \
        ptr->set_num_periods_to_extract(numberOfPeriods); \
        ptr->set_num_steady_state_periods(numberOfSteadyStatePeriods); \
        OpenMagnetics::NgspiceRunner runner; \
        if (!runner.is_available()) throw std::runtime_error("ngspice is not available"); \
        auto topologyWaveforms = SIM_WAVEFORM_CALL; \
        auto operatingPoints = SIM_OP_CALL; \
        json result; \
        json inputs; \
        inputs["designRequirements"] = json(); \
        to_json(inputs["designRequirements"], designRequirements); \
        inputs["operatingPoints"] = json::array(); \
        for (const auto& op : operatingPoints) { json j; to_json(j, op); inputs["operatingPoints"].push_back(j); } \
        result["inputs"] = inputs; \
        result["converterWaveforms"] = json::array(); \
        for (const auto& tw : topologyWaveforms) { json j; to_json(j, tw); result["converterWaveforms"].push_back(j); } \
        ptr->process(); \
        DIAG_BLOCK \
        return result; \
    } \
    catch (const std::exception& exc) { \
        return json{{"error", std::string(#func_name ": ") + exc.what()}}; \
    } \
}

// SingleSwitchForward diagnostics (has primaryTurnsRatio + resetVoltage)
#define SINGLE_SWITCH_FORWARD_DIAG_BLOCK \
    { \
        json diag; \
        const auto& names = ptr->get_per_op_name(); \
        const auto& v_maximum_duty_cycle = ptr->get_per_op_maximum_duty_cycle(); \
        const auto& v_computed_magnetizing_inductance = ptr->get_per_op_computed_magnetizing_inductance(); \
        const auto& v_computed_secondary_turns_ratio = ptr->get_per_op_computed_secondary_turns_ratio(); \
        const auto& v_primary_peak_current = ptr->get_per_op_primary_peak_current(); \
        const auto& v_secondary_peak_current = ptr->get_per_op_secondary_peak_current(); \
        const auto& v_magnetizing_peak_current = ptr->get_per_op_magnetizing_peak_current(); \
        const auto& v_is_ccm = ptr->get_per_op_is_ccm(); \
        const auto& v_computed_primary_turns_ratio = ptr->get_per_op_computed_primary_turns_ratio(); \
        const auto& v_reset_voltage = ptr->get_per_op_reset_voltage(); \
        diag["maximumDutyCycle"] = v_maximum_duty_cycle.empty() ? ptr->get_last_maximum_duty_cycle() : v_maximum_duty_cycle.front(); \
        diag["magnetizingInductance"] = v_computed_magnetizing_inductance.empty() ? ptr->get_last_computed_magnetizing_inductance() : v_computed_magnetizing_inductance.front(); \
        diag["secondaryTurnsRatio"] = v_computed_secondary_turns_ratio.empty() ? ptr->get_last_computed_secondary_turns_ratio() : v_computed_secondary_turns_ratio.front(); \
        diag["primaryPeakCurrent"] = v_primary_peak_current.empty() ? ptr->get_last_primary_peak_current() : v_primary_peak_current.front(); \
        diag["secondaryPeakCurrent"] = v_secondary_peak_current.empty() ? ptr->get_last_secondary_peak_current() : v_secondary_peak_current.front(); \
        diag["magnetizingPeakCurrent"] = v_magnetizing_peak_current.empty() ? ptr->get_last_magnetizing_peak_current() : v_magnetizing_peak_current.front(); \
        diag["isCcm"] = v_is_ccm.empty() ? ptr->get_last_is_ccm() : (bool)v_is_ccm.front(); \
        diag["primaryTurnsRatio"] = v_computed_primary_turns_ratio.empty() ? ptr->get_last_computed_primary_turns_ratio() : v_computed_primary_turns_ratio.front(); \
        diag["resetVoltage"] = v_reset_voltage.empty() ? ptr->get_last_reset_voltage() : v_reset_voltage.front(); \
        json perOp = json::array(); \
        for (size_t i = 0; i < v_maximum_duty_cycle.size(); ++i) { \
            json row; \
            row["operatingPointName"] = (i < names.size()) ? names[i] : ("OP " + std::to_string(i)); \
            row["maximumDutyCycle"] = v_maximum_duty_cycle[i]; \
            row["magnetizingInductance"] = v_computed_magnetizing_inductance[i]; \
            row["secondaryTurnsRatio"] = v_computed_secondary_turns_ratio[i]; \
            row["primaryPeakCurrent"] = v_primary_peak_current[i]; \
            row["secondaryPeakCurrent"] = v_secondary_peak_current[i]; \
            row["magnetizingPeakCurrent"] = v_magnetizing_peak_current[i]; \
            row["isCcm"] = (bool)v_is_ccm[i]; \
            row["primaryTurnsRatio"] = v_computed_primary_turns_ratio[i]; \
            row["resetVoltage"] = v_reset_voltage[i]; \
            perOp.push_back(row); \
        } \
        diag["perOp"] = perOp; \
        result["singleSwitchForwardDiagnostics"] = diag; \
    }

// Common forward diagnostics (TwoSwitchForward, ActiveClampForward, PushPull — no primaryTurnsRatio/resetVoltage)
#define FORWARD_COMMON_DIAG_BLOCK(diagKeyStr) \
    { \
        json diag; \
        const auto& names = ptr->get_per_op_name(); \
        const auto& v_maximum_duty_cycle = ptr->get_per_op_maximum_duty_cycle(); \
        const auto& v_computed_magnetizing_inductance = ptr->get_per_op_computed_magnetizing_inductance(); \
        const auto& v_computed_secondary_turns_ratio = ptr->get_per_op_computed_secondary_turns_ratio(); \
        const auto& v_primary_peak_current = ptr->get_per_op_primary_peak_current(); \
        const auto& v_secondary_peak_current = ptr->get_per_op_secondary_peak_current(); \
        const auto& v_magnetizing_peak_current = ptr->get_per_op_magnetizing_peak_current(); \
        const auto& v_is_ccm = ptr->get_per_op_is_ccm(); \
        diag["maximumDutyCycle"] = v_maximum_duty_cycle.empty() ? ptr->get_last_maximum_duty_cycle() : v_maximum_duty_cycle.front(); \
        diag["magnetizingInductance"] = v_computed_magnetizing_inductance.empty() ? ptr->get_last_computed_magnetizing_inductance() : v_computed_magnetizing_inductance.front(); \
        diag["secondaryTurnsRatio"] = v_computed_secondary_turns_ratio.empty() ? ptr->get_last_computed_secondary_turns_ratio() : v_computed_secondary_turns_ratio.front(); \
        diag["primaryPeakCurrent"] = v_primary_peak_current.empty() ? ptr->get_last_primary_peak_current() : v_primary_peak_current.front(); \
        diag["secondaryPeakCurrent"] = v_secondary_peak_current.empty() ? ptr->get_last_secondary_peak_current() : v_secondary_peak_current.front(); \
        diag["magnetizingPeakCurrent"] = v_magnetizing_peak_current.empty() ? ptr->get_last_magnetizing_peak_current() : v_magnetizing_peak_current.front(); \
        diag["isCcm"] = v_is_ccm.empty() ? ptr->get_last_is_ccm() : (bool)v_is_ccm.front(); \
        json perOp = json::array(); \
        for (size_t i = 0; i < v_maximum_duty_cycle.size(); ++i) { \
            json row; \
            row["operatingPointName"] = (i < names.size()) ? names[i] : ("OP " + std::to_string(i)); \
            row["maximumDutyCycle"] = v_maximum_duty_cycle[i]; \
            row["magnetizingInductance"] = v_computed_magnetizing_inductance[i]; \
            row["secondaryTurnsRatio"] = v_computed_secondary_turns_ratio[i]; \
            row["primaryPeakCurrent"] = v_primary_peak_current[i]; \
            row["secondaryPeakCurrent"] = v_secondary_peak_current[i]; \
            row["magnetizingPeakCurrent"] = v_magnetizing_peak_current[i]; \
            row["isCcm"] = (bool)v_is_ccm[i]; \
            perOp.push_back(row); \
        } \
        diag["perOp"] = perOp; \
        result[diagKeyStr] = diag; \
    }

// IsolatedBuck/IsolatedBuckBoost diagnostics — takes diagKeyStr so the
// expansion works inside DEFINE_SIMULATE_ISOLATED (diagKey is a macro
// parameter there, not a C++ variable).
#define ISOLATED_BUCK_DIAG_BLOCK(diagKeyStr) \
    { \
        json diag; \
        const auto& names = ptr->get_per_op_name(); \
        const auto& v_duty_cycle = ptr->get_per_op_duty_cycle(); \
        const auto& v_magnetizing_current_ripple = ptr->get_per_op_magnetizing_current_ripple(); \
        const auto& v_primary_average_current = ptr->get_per_op_primary_average_current(); \
        const auto& v_primary_peak_current = ptr->get_per_op_primary_peak_current(); \
        const auto& v_secondary_peak_current = ptr->get_per_op_secondary_peak_current(); \
        const auto& v_is_ccm = ptr->get_per_op_is_ccm(); \
        diag["dutyCycle"] = v_duty_cycle.empty() ? ptr->get_last_duty_cycle() : v_duty_cycle.front(); \
        diag["magnetizingCurrentRipple"] = v_magnetizing_current_ripple.empty() ? ptr->get_last_magnetizing_current_ripple() : v_magnetizing_current_ripple.front(); \
        diag["primaryAverageCurrent"] = v_primary_average_current.empty() ? ptr->get_last_primary_average_current() : v_primary_average_current.front(); \
        diag["primaryPeakCurrent"] = v_primary_peak_current.empty() ? ptr->get_last_primary_peak_current() : v_primary_peak_current.front(); \
        diag["secondaryPeakCurrent"] = v_secondary_peak_current.empty() ? ptr->get_last_secondary_peak_current() : v_secondary_peak_current.front(); \
        diag["isCcm"] = v_is_ccm.empty() ? ptr->get_last_is_ccm() : (bool)v_is_ccm.front(); \
        json perOp = json::array(); \
        for (size_t i = 0; i < v_duty_cycle.size(); ++i) { \
            json row; \
            row["operatingPointName"] = (i < names.size()) ? names[i] : ("OP " + std::to_string(i)); \
            row["dutyCycle"] = v_duty_cycle[i]; \
            row["magnetizingCurrentRipple"] = v_magnetizing_current_ripple[i]; \
            row["primaryAverageCurrent"] = v_primary_average_current[i]; \
            row["primaryPeakCurrent"] = v_primary_peak_current[i]; \
            row["secondaryPeakCurrent"] = v_secondary_peak_current[i]; \
            row["isCcm"] = (bool)v_is_ccm[i]; \
            perOp.push_back(row); \
        } \
        diag["perOp"] = perOp; \
        result[diagKeyStr] = diag; \
    }

DEFINE_SIMULATE_ISOLATED(simulate_forward_ideal_waveforms, OpenMagnetics::SingleSwitchForward, OpenMagnetics::AdvancedSingleSwitchForward,
    MAS::Topologies::SINGLE_SWITCH_FORWARD_CONVERTER, "singleSwitchForwardDiagnostics",
    ptr->simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance),
    ptr->simulate_and_extract_operating_points(turnsRatios, magnetizingInductance),
    SINGLE_SWITCH_FORWARD_DIAG_BLOCK)

DEFINE_SIMULATE_ISOLATED(simulate_two_switch_forward_ideal_waveforms, OpenMagnetics::TwoSwitchForward, OpenMagnetics::AdvancedTwoSwitchForward,
    MAS::Topologies::TWO_SWITCH_FORWARD_CONVERTER, "twoSwitchForwardDiagnostics",
    ptr->simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance),
    ptr->simulate_and_extract_operating_points(turnsRatios, magnetizingInductance),
    FORWARD_COMMON_DIAG_BLOCK("twoSwitchForwardDiagnostics"))

DEFINE_SIMULATE_ISOLATED(simulate_active_clamp_forward_ideal_waveforms, OpenMagnetics::ActiveClampForward, OpenMagnetics::AdvancedActiveClampForward,
    MAS::Topologies::ACTIVE_CLAMP_FORWARD_CONVERTER, "activeClampForwardDiagnostics",
    ptr->simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance),
    ptr->simulate_and_extract_operating_points(turnsRatios, magnetizingInductance),
    FORWARD_COMMON_DIAG_BLOCK("activeClampForwardDiagnostics"))

DEFINE_SIMULATE_ISOLATED(simulate_push_pull_ideal_waveforms, OpenMagnetics::PushPull, OpenMagnetics::AdvancedPushPull,
    MAS::Topologies::PUSH_PULL_CONVERTER, "pushPullDiagnostics",
    ptr->simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance),
    ptr->simulate_and_extract_operating_points(turnsRatios, magnetizingInductance),
    { \
        json diag; \
        const auto& names = ptr->get_per_op_name(); \
        const auto& v_duty_cycle = ptr->get_per_op_duty_cycle(); \
        const auto& v_switching_frequency = ptr->get_per_op_switching_frequency(); \
        const auto& v_primary_average_current = ptr->get_per_op_primary_average_current(); \
        const auto& v_primary_peak_current = ptr->get_per_op_primary_peak_current(); \
        const auto& v_magnetizing_peak_current = ptr->get_per_op_magnetizing_peak_current(); \
        const auto& v_is_ccm = ptr->get_per_op_is_ccm(); \
        diag["dutyCycle"] = v_duty_cycle.empty() ? ptr->get_last_duty_cycle() : v_duty_cycle.front(); \
        diag["switchingFrequency"] = v_switching_frequency.empty() ? ptr->get_last_switching_frequency() : v_switching_frequency.front(); \
        diag["primaryAverageCurrent"] = v_primary_average_current.empty() ? ptr->get_last_primary_average_current() : v_primary_average_current.front(); \
        diag["primaryPeakCurrent"] = v_primary_peak_current.empty() ? ptr->get_last_primary_peak_current() : v_primary_peak_current.front(); \
        diag["magnetizingPeakCurrent"] = v_magnetizing_peak_current.empty() ? ptr->get_last_magnetizing_peak_current() : v_magnetizing_peak_current.front(); \
        diag["isCcm"] = v_is_ccm.empty() ? ptr->get_last_is_ccm() : (bool)v_is_ccm.front(); \
        json perOp = json::array(); \
        for (size_t i = 0; i < v_duty_cycle.size(); ++i) { \
            json row; \
            row["operatingPointName"] = (i < names.size()) ? names[i] : ("OP " + std::to_string(i)); \
            row["dutyCycle"] = v_duty_cycle[i]; \
            row["switchingFrequency"] = v_switching_frequency[i]; \
            row["primaryAverageCurrent"] = v_primary_average_current[i]; \
            row["primaryPeakCurrent"] = v_primary_peak_current[i]; \
            row["magnetizingPeakCurrent"] = v_magnetizing_peak_current[i]; \
            row["isCcm"] = (bool)v_is_ccm[i]; \
            perOp.push_back(row); \
        } \
        diag["perOp"] = perOp; \
        result["pushPullDiagnostics"] = diag; \
    })

DEFINE_SIMULATE_ISOLATED(simulate_isolated_buck_ideal_waveforms, OpenMagnetics::IsolatedBuck, OpenMagnetics::AdvancedIsolatedBuck,
    MAS::Topologies::ISOLATED_BUCK_CONVERTER, "isolatedBuckDiagnostics",
    ptr->simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance),
    ptr->simulate_and_extract_operating_points(turnsRatios, magnetizingInductance),
    ISOLATED_BUCK_DIAG_BLOCK("isolatedBuckDiagnostics"))

DEFINE_SIMULATE_ISOLATED(simulate_isolated_buck_boost_ideal_waveforms, OpenMagnetics::IsolatedBuckBoost, OpenMagnetics::AdvancedIsolatedBuckBoost,
    MAS::Topologies::ISOLATED_BUCK_BOOST_CONVERTER, "isolatedBuckBoostDiagnostics",
    ptr->simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance),
    ptr->simulate_and_extract_operating_points(turnsRatios, magnetizingInductance),
    ISOLATED_BUCK_DIAG_BLOCK("isolatedBuckBoostDiagnostics"))

#undef SINGLE_SWITCH_FORWARD_DIAG_BLOCK
#undef FORWARD_COMMON_DIAG_BLOCK
#undef ISOLATED_BUCK_DIAG_BLOCK

// ------ Weinberg (scalar turnsRatio + magnetizingInductance) ------
json simulate_weinberg_ideal_waveforms(json inputsJson) {
    try {
        bool isAdvanced = inputsJson.contains("desiredInductance");
        size_t numberOfPeriods = inputsJson.value("numberOfPeriods", 2);
        size_t numberOfSteadyStatePeriods = inputsJson.value("numberOfSteadyStatePeriods", 5);
        OpenMagnetics::DesignRequirements designRequirements;
        double magnetizingInductance;
        double turnsRatio;
        std::unique_ptr<OpenMagnetics::Weinberg> ptr;
        if (isAdvanced) {
            auto advPtr = std::make_unique<OpenMagnetics::AdvancedWeinberg>(inputsJson);
            magnetizingInductance = advPtr->get_desired_inductance();
            turnsRatio = advPtr->get_desired_turns_ratio();
            designRequirements = advPtr->process_design_requirements();
            ptr = std::move(advPtr);
        } else {
            ptr = std::make_unique<OpenMagnetics::Weinberg>(inputsJson);
            designRequirements = ptr->process_design_requirements();
            if (designRequirements.get_turns_ratios().empty() || !designRequirements.get_turns_ratios()[0].get_nominal())
                throw std::runtime_error("Weinberg: no turns ratio available");
            turnsRatio = designRequirements.get_turns_ratios()[0].get_nominal().value();
            magnetizingInductance = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
            if (!(magnetizingInductance > 0)) throw std::runtime_error("Unable to calculate inductance");
        }
        ptr->set_num_periods_to_extract(numberOfPeriods);
        ptr->set_num_steady_state_periods(numberOfSteadyStatePeriods);
        OpenMagnetics::NgspiceRunner runner;
        if (!runner.is_available()) throw std::runtime_error("ngspice is not available");
        auto operatingPoints = ptr->simulate_and_extract_operating_points(turnsRatio, magnetizingInductance);
        auto topologyWaveforms = ptr->simulate_and_extract_topology_waveforms(turnsRatio, magnetizingInductance, numberOfPeriods);
        json result;
        json inputs;
        inputs["designRequirements"] = json();
        to_json(inputs["designRequirements"], designRequirements);
        inputs["operatingPoints"] = json::array();
        for (const auto& op : operatingPoints) { json j; to_json(j, op); inputs["operatingPoints"].push_back(j); }
        result["inputs"] = inputs;
        result["converterWaveforms"] = json::array();
        for (const auto& tw : topologyWaveforms) { json j; to_json(j, tw); result["converterWaveforms"].push_back(j); }
        ptr->process();
        {
            json diag;
            const auto& names = ptr->get_per_op_name();
            const auto& v_duty_cycle = ptr->get_per_op_duty_cycle();
            const auto& v_conversion_ratio = ptr->get_per_op_conversion_ratio();
            const auto& v_operating_regime = ptr->get_per_op_operating_regime();
            const auto& v_overlap_fraction = ptr->get_per_op_overlap_fraction();
            const auto& v_switch_peak_voltage = ptr->get_per_op_switch_peak_voltage();
            const auto& v_switch_peak_current = ptr->get_per_op_switch_peak_current();
            const auto& v_diode_peak_reverse_voltage = ptr->get_per_op_diode_peak_reverse_voltage();
            const auto& v_diode_peak_current = ptr->get_per_op_diode_peak_current();
            const auto& v_energy_recovery_avg_current = ptr->get_per_op_energy_recovery_avg_current();
            const auto& v_input_inductor_average = ptr->get_per_op_input_inductor_average();
            const auto& v_input_inductor_ripple = ptr->get_per_op_input_inductor_ripple();
            const auto& v_magnetizing_ripple = ptr->get_per_op_magnetizing_ripple();
            const auto& v_flux_imbalance_margin = ptr->get_per_op_flux_imbalance_margin();
            const auto& v_rhp_zero_frequency = ptr->get_per_op_rhp_zero_frequency();
            const auto& v_is_ccm = ptr->get_per_op_is_ccm();
            const auto& v_sized_co = ptr->get_per_op_sized_co();
            const auto& v_output_voltage_ripple = ptr->get_per_op_output_voltage_ripple();
            diag["dutyCycle"] = v_duty_cycle.empty() ? ptr->get_last_duty_cycle() : v_duty_cycle.front();
            diag["conversionRatio"] = v_conversion_ratio.empty() ? ptr->get_last_conversion_ratio() : v_conversion_ratio.front();
            diag["operatingRegime"] = v_operating_regime.empty() ? ptr->get_last_operating_regime() : v_operating_regime.front();
            diag["overlapFraction"] = v_overlap_fraction.empty() ? ptr->get_last_overlap_fraction() : v_overlap_fraction.front();
            diag["switchPeakVoltage"] = v_switch_peak_voltage.empty() ? ptr->get_last_switch_peak_voltage() : v_switch_peak_voltage.front();
            diag["switchPeakCurrent"] = v_switch_peak_current.empty() ? ptr->get_last_switch_peak_current() : v_switch_peak_current.front();
            diag["diodePeakReverseVoltage"] = v_diode_peak_reverse_voltage.empty() ? ptr->get_last_diode_peak_reverse_voltage() : v_diode_peak_reverse_voltage.front();
            diag["diodePeakCurrent"] = v_diode_peak_current.empty() ? ptr->get_last_diode_peak_current() : v_diode_peak_current.front();
            diag["energyRecoveryAvgCurrent"] = v_energy_recovery_avg_current.empty() ? ptr->get_last_energy_recovery_avg_current() : v_energy_recovery_avg_current.front();
            diag["inputInductorAverage"] = v_input_inductor_average.empty() ? ptr->get_last_input_inductor_average() : v_input_inductor_average.front();
            diag["inputInductorRipple"] = v_input_inductor_ripple.empty() ? ptr->get_last_input_inductor_ripple() : v_input_inductor_ripple.front();
            diag["magnetizingRipple"] = v_magnetizing_ripple.empty() ? ptr->get_last_magnetizing_ripple() : v_magnetizing_ripple.front();
            diag["fluxImbalanceMargin"] = v_flux_imbalance_margin.empty() ? ptr->get_last_flux_imbalance_margin() : v_flux_imbalance_margin.front();
            diag["rhpZeroFrequency"] = v_rhp_zero_frequency.empty() ? ptr->get_last_rhp_zero_frequency() : v_rhp_zero_frequency.front();
            diag["isCcm"] = v_is_ccm.empty() ? ptr->get_last_is_ccm() : (bool)v_is_ccm.front();
            diag["sizedCo"] = v_sized_co.empty() ? ptr->get_last_sized_co() : v_sized_co.front();
            diag["outputVoltageRipple"] = v_output_voltage_ripple.empty() ? ptr->get_last_output_voltage_ripple() : v_output_voltage_ripple.front();
            json perOp = json::array();
            for (size_t i = 0; i < v_duty_cycle.size(); ++i) {
                json row;
                row["operatingPointName"] = (i < names.size()) ? names[i] : ("OP " + std::to_string(i));
                row["dutyCycle"] = v_duty_cycle[i]; row["conversionRatio"] = v_conversion_ratio[i];
                row["operatingRegime"] = v_operating_regime[i]; row["overlapFraction"] = v_overlap_fraction[i];
                row["switchPeakVoltage"] = v_switch_peak_voltage[i]; row["switchPeakCurrent"] = v_switch_peak_current[i];
                row["diodePeakReverseVoltage"] = v_diode_peak_reverse_voltage[i]; row["diodePeakCurrent"] = v_diode_peak_current[i];
                row["energyRecoveryAvgCurrent"] = v_energy_recovery_avg_current[i];
                row["inputInductorAverage"] = v_input_inductor_average[i]; row["inputInductorRipple"] = v_input_inductor_ripple[i];
                row["magnetizingRipple"] = v_magnetizing_ripple[i]; row["fluxImbalanceMargin"] = v_flux_imbalance_margin[i];
                row["rhpZeroFrequency"] = v_rhp_zero_frequency[i]; row["isCcm"] = (bool)v_is_ccm[i];
                row["sizedCo"] = v_sized_co[i]; row["outputVoltageRipple"] = v_output_voltage_ripple[i];
                perOp.push_back(row);
            }
            diag["perOp"] = perOp;
            result["weinbergDiagnostics"] = diag;
        }
        return result;
    }
    catch (const std::exception& exc) {
        return json{{"error", std::string("simulate_weinberg_ideal_waveforms: ") + exc.what()}};
    }
}

// ------ LLC (special: not Advanced, just Llc with optional user overrides) ------
json simulate_llc_ideal_waveforms(json inputsJson) {
    try {
        OpenMagnetics::Llc llcInputs(inputsJson);
        if (inputsJson.contains("desiredResonantInductance") && inputsJson["desiredResonantInductance"].is_number())
            llcInputs.set_user_resonant_inductance(inputsJson["desiredResonantInductance"].get<double>());
        if (inputsJson.contains("desiredResonantCapacitance") && inputsJson["desiredResonantCapacitance"].is_number())
            llcInputs.set_user_resonant_capacitance(inputsJson["desiredResonantCapacitance"].get<double>());

        auto designRequirements = llcInputs.process_design_requirements();
        double magnetizingInductance = inputsJson.value("magnetizingInductance", 200e-6);
        std::vector<double> turnsRatios;
        for (const auto& tr : designRequirements.get_turns_ratios()) {
            if (tr.get_nominal()) turnsRatios.push_back(tr.get_nominal().value());
        }
        OpenMagnetics::NgspiceRunner runner;
        if (!runner.is_available()) throw std::runtime_error("ngspice is not available");
        size_t numberOfPeriods = inputsJson.value("numberOfPeriods", 2);
        size_t numberOfSteadyStatePeriods = inputsJson.value("numberOfSteadyStatePeriods", 3);
        llcInputs.set_num_periods_to_extract(static_cast<int>(numberOfPeriods));
        llcInputs.set_num_steady_state_periods(static_cast<int>(numberOfSteadyStatePeriods));
        auto topologyWaveforms = llcInputs.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance, numberOfPeriods);
        auto operatingPoints = llcInputs.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance, numberOfPeriods);
        json result;
        json inputsJ;
        inputsJ["designRequirements"] = json();
        to_json(inputsJ["designRequirements"], designRequirements);
        inputsJ["operatingPoints"] = json::array();
        for (const auto& op : operatingPoints) { json j; to_json(j, op); inputsJ["operatingPoints"].push_back(j); }
        result["inputs"] = inputsJ;
        result["converterWaveforms"] = json::array();
        for (const auto& tw : topologyWaveforms) { json j; to_json(j, tw); result["converterWaveforms"].push_back(j); }
        json diag;
        diag["computedResonantInductance"]  = llcInputs.get_computed_resonant_inductance();
        diag["computedResonantCapacitance"] = llcInputs.get_computed_resonant_capacitance();
        diag["computedInductanceRatio"]     = llcInputs.get_computed_inductance_ratio();
        diag["lipFrequency"]                = llcInputs.get_lip_frequency();
        diag["lipInputVoltage"]             = llcInputs.get_lip_input_voltage();
        diag["lastSubStateSequence"]        = llcInputs.get_last_sub_state_sequence();
        {
            const auto& names = llcInputs.get_per_op_name();
            const auto& mode  = llcInputs.get_per_op_mode();
            const auto& res   = llcInputs.get_per_op_steady_state_residual();
            const auto& zvs   = llcInputs.get_per_op_zvs_margin_lagging();
            const auto& ipk   = llcInputs.get_per_op_primary_peak_current();
            diag["lastMode"]                = mode.empty() ? llcInputs.get_last_mode()                  : mode.front();
            diag["lastSteadyStateResidual"] = res.empty()  ? llcInputs.get_last_steady_state_residual() : res.front();
            diag["lastZvsMarginLagging"]    = zvs.empty()  ? llcInputs.get_last_zvs_margin_lagging()    : zvs.front();
            diag["lastPrimaryPeakCurrent"]  = ipk.empty()  ? llcInputs.get_last_primary_peak_current()  : ipk.front();
            json perOp = json::array();
            for (size_t i = 0; i < mode.size(); ++i) {
                json row;
                row["operatingPointName"]    = (i < names.size()) ? names[i] : ("OP " + std::to_string(i));
                row["lastMode"]              = mode[i];
                row["steadyStateResidual"]   = res[i];
                row["zvsMarginLagging"]      = zvs[i];
                row["primaryPeakCurrent"]    = ipk[i];
                perOp.push_back(row);
            }
            diag["perOp"] = perOp;
        }
        result["llcDiagnostics"] = diag;
        return result;
    }
    catch (const std::exception& exc) {
        return json{{"error", std::string("simulate_llc_ideal_waveforms: ") + exc.what()}};
    }
}

// ------ CLLC ------
json simulate_cllc_ideal_waveforms(json inputsJson) {
    try {
        OpenMagnetics::CllcConverter cllc(inputsJson);
        auto designRequirements = cllc.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designRequirements.get_turns_ratios()) {
            if (tr.get_nominal()) turnsRatios.push_back(tr.get_nominal().value());
        }
        if (turnsRatios.empty()) throw std::runtime_error("CLLC: no turns ratios");
        double magnetizingInductance;
        if (inputsJson.contains("magnetizingInductance") && inputsJson["magnetizingInductance"].is_number())
            magnetizingInductance = inputsJson["magnetizingInductance"].get<double>();
        else {
            magnetizingInductance = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
            if (!(magnetizingInductance > 0)) throw std::runtime_error("CLLC: no magnetizing inductance");
        }
        OpenMagnetics::NgspiceRunner runner;
        if (!runner.is_available()) throw std::runtime_error("ngspice is not available");
        size_t numberOfPeriods = inputsJson.value("numberOfPeriods", 2);
        size_t numberOfSteadyStatePeriods = inputsJson.value("numberOfSteadyStatePeriods", 3);
        cllc.set_num_periods_to_extract(static_cast<int>(numberOfPeriods));
        cllc.set_num_steady_state_periods(static_cast<int>(numberOfSteadyStatePeriods));
        auto topologyWaveforms = cllc.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance, numberOfPeriods);
        auto operatingPoints = cllc.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
        cllc.process_operating_points(turnsRatios, magnetizingInductance);
        json result;
        json inputsJ;
        inputsJ["designRequirements"] = json();
        to_json(inputsJ["designRequirements"], designRequirements);
        inputsJ["operatingPoints"] = json::array();
        for (const auto& op : operatingPoints) { json j; to_json(j, op); inputsJ["operatingPoints"].push_back(j); }
        result["inputs"] = inputsJ;
        result["converterWaveforms"] = json::array();
        for (const auto& tw : topologyWaveforms) { json j; to_json(j, tw); result["converterWaveforms"].push_back(j); }
        json diag;
        diag["lipFrequency"]              = cllc.get_lip_frequency();
        diag["lastSubStateSequence"]      = cllc.get_last_sub_state_sequence();
        diag["bridgeVoltageFactor"]       = cllc.get_bridge_voltage_factor();
        {
            const auto& names = cllc.get_per_op_name();
            const auto& mode  = cllc.get_per_op_mode();
            const auto& res   = cllc.get_per_op_steady_state_residual();
            const auto& zp    = cllc.get_per_op_zvs_margin_primary();
            const auto& zs    = cllc.get_per_op_zvs_margin_secondary();
            const auto& rtt   = cllc.get_per_op_resonant_transition_time();
            const auto& ipk   = cllc.get_per_op_primary_peak_current();
            const auto& vcr   = cllc.get_per_op_resonant_cap_peak_voltage();
            diag["lastMode"]                  = mode.empty() ? cllc.get_last_mode()                      : mode.front();
            diag["lastSteadyStateResidual"]   = res.empty()  ? cllc.get_last_steady_state_residual()     : res.front();
            diag["lastZvsMarginPrimary"]      = zp.empty()   ? cllc.get_last_zvs_margin_primary()        : zp.front();
            diag["lastZvsMarginSecondary"]    = zs.empty()   ? cllc.get_last_zvs_margin_secondary()      : zs.front();
            diag["lastResonantTransitionTime"]= rtt.empty()  ? cllc.get_last_resonant_transition_time()  : rtt.front();
            diag["lastPrimaryPeakCurrent"]    = ipk.empty()  ? cllc.get_last_primary_peak_current()      : ipk.front();
            diag["lastResonantCapPeakVoltage"]= vcr.empty()  ? cllc.get_last_resonant_cap_peak_voltage() : vcr.front();
            json perOp = json::array();
            for (size_t i = 0; i < mode.size(); ++i) {
                json row;
                row["operatingPointName"]     = (i < names.size()) ? names[i] : ("OP " + std::to_string(i));
                row["lastMode"]               = mode[i]; row["steadyStateResidual"] = res[i];
                row["zvsMarginPrimary"]       = zp[i]; row["zvsMarginSecondary"] = zs[i];
                row["resonantTransitionTime"] = rtt[i]; row["primaryPeakCurrent"] = ipk[i];
                row["resonantCapPeakVoltage"] = vcr[i];
                perOp.push_back(row);
            }
            diag["perOp"] = perOp;
        }
        result["cllcDiagnostics"] = diag;
        return result;
    }
    catch (const std::exception& exc) {
        return json{{"error", std::string("simulate_cllc_ideal_waveforms: ") + exc.what()}};
    }
}

// ------ CLLLC ------
json simulate_clllc_ideal_waveforms(json inputsJson) {
    try {
        OpenMagnetics::Clllc model(inputsJson);
        auto designRequirements = model.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designRequirements.get_turns_ratios()) {
            if (tr.get_nominal()) turnsRatios.push_back(tr.get_nominal().value());
        }
        if (turnsRatios.empty()) throw std::runtime_error("Clllc: no turns ratios");
        double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
        if (!(magnetizingInductance > 0)) throw std::runtime_error("Clllc: no magnetizing inductance");
        OpenMagnetics::NgspiceRunner runner;
        if (!runner.is_available()) throw std::runtime_error("ngspice is not available");
        size_t numberOfPeriods = inputsJson.value("numberOfPeriods", 2);
        size_t numberOfSteadyStatePeriods = inputsJson.value("numberOfSteadyStatePeriods", 5);
        model.set_num_periods_to_extract(static_cast<int>(numberOfPeriods));
        model.set_num_steady_state_periods(static_cast<int>(numberOfSteadyStatePeriods));
        auto operatingPoints = model.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
        auto topologyWaveforms = model.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance, numberOfPeriods);
        json result;
        json inputs;
        inputs["designRequirements"] = json();
        to_json(inputs["designRequirements"], designRequirements);
        inputs["operatingPoints"] = json::array();
        for (const auto& op : operatingPoints) { json j; to_json(j, op); inputs["operatingPoints"].push_back(j); }
        result["inputs"] = inputs;
        result["converterWaveforms"] = json::array();
        for (const auto& tw : topologyWaveforms) { json j; to_json(j, tw); result["converterWaveforms"].push_back(j); }
        model.process();
        {
            json diag;
            diag["computedPrimarySeriesInductance"]      = model.get_computed_primary_series_inductance();
            diag["computedSecondarySeriesInductance"]    = model.get_computed_secondary_series_inductance();
            diag["computedPrimaryResonantCapacitance"]   = model.get_computed_primary_resonant_capacitance();
            diag["computedSecondaryResonantCapacitance"] = model.get_computed_secondary_resonant_capacitance();
            diag["computedMagnetizingInductance"]        = model.get_computed_magnetizing_inductance();
            diag["computedTurnsRatio"]                   = model.get_computed_turns_ratio();
            diag["computedDeadTime"]                     = model.get_computed_dead_time();
            diag["computedInductanceRatioK"]             = model.get_computed_inductance_ratio_k();
            diag["computedQualityFactor"]                = model.get_computed_quality_factor();
            diag["computedPrimaryResonantFrequency"]     = model.get_computed_primary_resonant_frequency();
            const auto& names = model.get_per_op_name();
            const auto& mf    = model.get_per_op_mode_forward();
            const auto& mr    = model.get_per_op_mode_reverse();
            const auto& zpl   = model.get_per_op_zvs_margin_primary_lagging();
            const auto& zsl   = model.get_per_op_zvs_margin_secondary_lagging();
            const auto& zltp  = model.get_per_op_zvs_load_threshold_primary();
            const auto& zlts  = model.get_per_op_zvs_load_threshold_secondary();
            const auto& rtt   = model.get_per_op_resonant_transition_time();
            const auto& ipk   = model.get_per_op_primary_peak_current();
            const auto& spk   = model.get_per_op_secondary_peak_current();
            const auto& irms  = model.get_per_op_primary_rms_current();
            const auto& srms  = model.get_per_op_secondary_rms_current();
            const auto& mpk   = model.get_per_op_magnetizing_peak_current();
            const auto& vc1   = model.get_per_op_cr1_peak_voltage();
            const auto& vc2   = model.get_per_op_cr2_peak_voltage();
            const auto& shr   = model.get_per_op_current_sharing_ratio();
            const auto& res   = model.get_per_op_steady_state_residual();
            diag["lastModeForward"]                = mf.empty()   ? model.get_last_mode_forward()                  : mf.front();
            diag["lastModeReverse"]                = mr.empty()   ? model.get_last_mode_reverse()                  : mr.front();
            diag["lastZvsMarginPrimaryLagging"]    = zpl.empty()  ? model.get_last_zvs_margin_primary_lagging()    : zpl.front();
            diag["lastZvsMarginSecondaryLagging"]  = zsl.empty()  ? model.get_last_zvs_margin_secondary_lagging()  : zsl.front();
            diag["lastZvsLoadThresholdPrimary"]    = zltp.empty() ? model.get_last_zvs_load_threshold_primary()    : zltp.front();
            diag["lastZvsLoadThresholdSecondary"]  = zlts.empty() ? model.get_last_zvs_load_threshold_secondary()  : zlts.front();
            diag["lastResonantTransitionTime"]     = rtt.empty()  ? model.get_last_resonant_transition_time()      : rtt.front();
            diag["lastPrimaryPeakCurrent"]         = ipk.empty()  ? model.get_last_primary_peak_current()          : ipk.front();
            diag["lastSecondaryPeakCurrent"]       = spk.empty()  ? model.get_last_secondary_peak_current()        : spk.front();
            diag["lastPrimaryRmsCurrent"]          = irms.empty() ? model.get_last_primary_rms_current()           : irms.front();
            diag["lastSecondaryRmsCurrent"]        = srms.empty() ? model.get_last_secondary_rms_current()         : srms.front();
            diag["lastMagnetizingPeakCurrent"]     = mpk.empty()  ? model.get_last_magnetizing_peak_current()      : mpk.front();
            diag["lastCr1PeakVoltage"]             = vc1.empty()  ? model.get_last_cr1_peak_voltage()              : vc1.front();
            diag["lastCr2PeakVoltage"]             = vc2.empty()  ? model.get_last_cr2_peak_voltage()              : vc2.front();
            diag["lastCurrentSharingRatio"]        = shr.empty()  ? model.get_last_current_sharing_ratio()         : shr.front();
            diag["lastSteadyStateResidual"]        = res.empty()  ? model.get_last_steady_state_residual()         : res.front();
            json perOp = json::array();
            for (size_t i = 0; i < mf.size(); ++i) {
                json row;
                row["operatingPointName"]              = (i < names.size()) ? names[i] : ("OP " + std::to_string(i));
                row["modeForward"] = mf[i]; row["modeReverse"] = mr[i];
                row["zvsMarginPrimaryLagging"] = zpl[i]; row["zvsMarginSecondaryLagging"] = zsl[i];
                row["zvsLoadThresholdPrimary"] = zltp[i]; row["zvsLoadThresholdSecondary"] = zlts[i];
                row["resonantTransitionTime"] = rtt[i]; row["primaryPeakCurrent"] = ipk[i];
                row["secondaryPeakCurrent"] = spk[i]; row["primaryRmsCurrent"] = irms[i];
                row["secondaryRmsCurrent"] = srms[i]; row["magnetizingPeakCurrent"] = mpk[i];
                row["cr1PeakVoltage"] = vc1[i]; row["cr2PeakVoltage"] = vc2[i];
                row["currentSharingRatio"] = shr[i]; row["steadyStateResidual"] = res[i];
                perOp.push_back(row);
            }
            diag["perOp"] = perOp;
            result["clllcDiagnostics"] = diag;
        }
        return result;
    }
    catch (const std::exception& exc) {
        return json{{"error", std::string("simulate_clllc_ideal_waveforms: ") + exc.what()}};
    }
}

// ------ SRC ------
json simulate_src_ideal_waveforms(json inputsJson) {
    try {
        bool isAdvanced = inputsJson.contains("desiredTurnsRatios") ||
                          inputsJson.contains("desiredResonantInductance") ||
                          inputsJson.contains("desiredResonantCapacitance");
        std::unique_ptr<OpenMagnetics::Src> model;
        if (isAdvanced) model = std::make_unique<OpenMagnetics::AdvancedSrc>(inputsJson);
        else model = std::make_unique<OpenMagnetics::Src>(inputsJson);
        size_t numberOfPeriods = inputsJson.value("numberOfPeriods", 2);
        size_t numberOfSteadyStatePeriods = inputsJson.value("numberOfSteadyStatePeriods", 5);
        model->set_num_periods_to_extract(static_cast<int>(numberOfPeriods));
        model->set_num_steady_state_periods(static_cast<int>(numberOfSteadyStatePeriods));
        auto designRequirements = model->process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designRequirements.get_turns_ratios()) {
            if (tr.get_nominal()) turnsRatios.push_back(tr.get_nominal().value());
        }
        if (turnsRatios.empty()) throw std::runtime_error("SRC: no turns ratios");
        double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
        if (!(magnetizingInductance > 0)) throw std::runtime_error("SRC: no magnetizing inductance");
        OpenMagnetics::NgspiceRunner runner;
        if (!runner.is_available()) throw std::runtime_error("ngspice is not available");
        auto topologyWaveforms = model->simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance, numberOfPeriods);
        auto operatingPoints = model->simulate_and_extract_operating_points(turnsRatios, magnetizingInductance, numberOfPeriods);
        json result;
        json inputsJ;
        inputsJ["designRequirements"] = json();
        to_json(inputsJ["designRequirements"], designRequirements);
        inputsJ["operatingPoints"] = json::array();
        for (const auto& op : operatingPoints) { json j; to_json(j, op); inputsJ["operatingPoints"].push_back(j); }
        result["inputs"] = inputsJ;
        result["converterWaveforms"] = json::array();
        for (const auto& tw : topologyWaveforms) { json j; to_json(j, tw); result["converterWaveforms"].push_back(j); }
        json diag;
        diag["computedResonantInductance"]  = model->get_computed_resonant_inductance();
        diag["computedResonantCapacitance"] = model->get_computed_resonant_capacitance();
        diag["computedResonantFrequency"]   = model->get_computed_resonant_frequency();
        {
            const auto& names = model->get_per_op_name();
            const auto& gm    = model->get_per_op_gain_m();
            const auto& nfsw  = model->get_per_op_normalized_fsw();
            const auto& ir    = model->get_per_op_ir_peak();
            const auto& vc    = model->get_per_op_vcr_peak();
            const auto& abv   = model->get_per_op_is_above_resonance();
            diag["lastGainM"]            = gm.empty()   ? model->get_last_gain()               : gm.front();
            diag["lastNormalizedFsw"]    = nfsw.empty() ? model->get_last_normalized_fsw()     : nfsw.front();
            diag["lastIrPeak"]           = ir.empty()   ? model->get_last_ir_peak()            : ir.front();
            diag["lastVcrPeak"]          = vc.empty()   ? model->get_last_vcr_peak()           : vc.front();
            diag["lastIsAboveResonance"] = abv.empty()  ? model->get_last_is_above_resonance() : (abv.front() != 0);
            json perOp = json::array();
            for (size_t i = 0; i < gm.size(); ++i) {
                json row;
                row["operatingPointName"] = (i < names.size()) ? names[i] : ("OP " + std::to_string(i));
                row["gainM"] = gm[i]; row["normalizedFsw"] = nfsw[i]; row["irPeak"] = ir[i];
                row["vcrPeak"] = vc[i]; row["isAboveResonance"] = (abv[i] != 0);
                perOp.push_back(row);
            }
            diag["perOp"] = perOp;
        }
        result["srcDiagnostics"] = diag;
        return result;
    }
    catch (const std::exception& exc) {
        return json{{"error", std::string("simulate_src_ideal_waveforms: ") + exc.what()}};
    }
}

// ------ DAB (always AdvancedDab) ------
json simulate_dab_ideal_waveforms(json inputsJson) {
    try {
        OpenMagnetics::AdvancedDab dabInputs(inputsJson);
        auto inputs = dabInputs.process();
        auto designRequirements = inputs.get_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designRequirements.get_turns_ratios()) {
            if (tr.get_nominal()) turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
        if (!(magnetizingInductance > 0)) throw std::runtime_error("DAB: no magnetizing inductance");
        OpenMagnetics::NgspiceRunner runner;
        if (!runner.is_available()) throw std::runtime_error("ngspice is not available");
        size_t numberOfPeriods = inputsJson.value("numberOfPeriods", 2);
        size_t numberOfSteadyStatePeriods = inputsJson.value("numberOfSteadyStatePeriods", 3);
        dabInputs.set_num_periods_to_extract(static_cast<int>(numberOfPeriods));
        dabInputs.set_num_steady_state_periods(static_cast<int>(numberOfSteadyStatePeriods));
        auto topologyWaveforms = dabInputs.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance, numberOfPeriods);
        auto operatingPoints = dabInputs.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
        json result;
        json inputsJ;
        inputsJ["designRequirements"] = json();
        to_json(inputsJ["designRequirements"], designRequirements);
        inputsJ["operatingPoints"] = json::array();
        for (const auto& op : operatingPoints) { json j; to_json(j, op); inputsJ["operatingPoints"].push_back(j); }
        result["inputs"] = inputsJ;
        result["converterWaveforms"] = json::array();
        for (const auto& tw : topologyWaveforms) { json j; to_json(j, tw); result["converterWaveforms"].push_back(j); }
        {
            json dabDiag;
            const auto& names = dabInputs.get_per_op_name();
            const auto& mt    = dabInputs.get_per_op_modulation_type();
            const auto& zp    = dabInputs.get_per_op_zvs_margin_primary();
            const auto& zs    = dabInputs.get_per_op_zvs_margin_secondary();
            const auto& d3    = dabInputs.get_per_op_d3_rad();
            const auto& vr    = dabInputs.get_per_op_voltage_conversion_ratio();
            dabDiag["modulationType"]           = mt.empty() ? dabInputs.get_last_modulation_type()       : mt.front();
            dabDiag["computedD3Deg"]            = (d3.empty() ? dabInputs.get_last_d3_rad()               : d3.front()) * 180.0 / M_PI;
            dabDiag["zvsMarginPrimaryDeg"]      = (zp.empty() ? dabInputs.get_last_zvs_margin_primary()   : zp.front()) * 180.0 / M_PI;
            dabDiag["zvsMarginSecondaryDeg"]    = (zs.empty() ? dabInputs.get_last_zvs_margin_secondary() : zs.front()) * 180.0 / M_PI;
            dabDiag["computedSeriesInductance"] = dabInputs.get_computed_series_inductance();
            dabDiag["voltageConversionRatio"]   = vr.empty() ? dabInputs.get_last_voltage_conversion_ratio() : vr.front();
            json perOp = json::array();
            for (size_t i = 0; i < mt.size(); ++i) {
                json row;
                row["operatingPointName"]    = (i < names.size()) ? names[i] : ("OP " + std::to_string(i));
                row["modulationType"]        = mt[i]; row["zvsMarginPrimaryDeg"] = zp[i] * 180.0 / M_PI;
                row["zvsMarginSecondaryDeg"] = zs[i] * 180.0 / M_PI; row["computedD3Deg"] = d3[i] * 180.0 / M_PI;
                row["voltageConversionRatio"]= vr[i];
                perOp.push_back(row);
            }
            dabDiag["perOp"] = perOp;
            result["dabDiagnostics"] = dabDiag;
        }
        return result;
    }
    catch (const std::exception& exc) {
        return json{{"error", std::string("simulate_dab_ideal_waveforms: ") + exc.what()}};
    }
}

// ------ PSFB / PSHB (always AdvancedPsfb/AdvancedPshb, same diagnostic schema) ------
#define DEFINE_SIMULATE_PHASE_SHIFTED(func_name, AdvancedType, diagKey) \
json func_name(json inputsJson) { \
    try { \
        AdvancedType model(inputsJson); \
        auto inputs = model.process(); \
        auto designRequirements = inputs.get_design_requirements(); \
        std::vector<double> turnsRatios; \
        for (const auto& tr : designRequirements.get_turns_ratios()) { \
            if (tr.get_nominal()) turnsRatios.push_back(tr.get_nominal().value()); \
        } \
        double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance()); \
        if (!(magnetizingInductance > 0)) throw std::runtime_error(std::string(#func_name) + ": no magnetizing inductance"); \
        OpenMagnetics::NgspiceRunner runner; \
        if (!runner.is_available()) throw std::runtime_error("ngspice is not available"); \
        size_t numberOfPeriods = inputsJson.value("numberOfPeriods", 2); \
        size_t numberOfSteadyStatePeriods = inputsJson.value("numberOfSteadyStatePeriods", 3); \
        model.set_num_periods_to_extract(static_cast<int>(numberOfPeriods)); \
        model.set_num_steady_state_periods(static_cast<int>(numberOfSteadyStatePeriods)); \
        auto topologyWaveforms = model.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance, numberOfPeriods); \
        auto operatingPoints = model.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance); \
        json result; \
        json inputsJ; \
        inputsJ["designRequirements"] = json(); \
        to_json(inputsJ["designRequirements"], designRequirements); \
        inputsJ["operatingPoints"] = json::array(); \
        for (const auto& op : operatingPoints) { json j; to_json(j, op); inputsJ["operatingPoints"].push_back(j); } \
        result["inputs"] = inputsJ; \
        result["converterWaveforms"] = json::array(); \
        for (const auto& tw : topologyWaveforms) { json j; to_json(j, tw); result["converterWaveforms"].push_back(j); } \
        { \
            json diag; \
            const auto& names = model.get_per_op_name(); \
            const auto& dcl   = model.get_per_op_duty_cycle_loss(); \
            const auto& deff  = model.get_per_op_effective_duty_cycle(); \
            const auto& zvs   = model.get_per_op_zvs_margin_lagging(); \
            const auto& zlt   = model.get_per_op_zvs_load_threshold(); \
            const auto& rtt   = model.get_per_op_resonant_transition_time(); \
            const auto& ipk   = model.get_per_op_primary_peak_current(); \
            diag["effectiveDutyCycle"]            = deff.empty() ? model.get_last_effective_duty_cycle()     : deff.front(); \
            diag["dutyCycleLoss"]                 = dcl.empty()  ? model.get_last_duty_cycle_loss()          : dcl.front(); \
            diag["zvsMarginLagging"]              = zvs.empty()  ? model.get_last_zvs_margin_lagging()       : zvs.front(); \
            diag["zvsLoadThreshold"]              = zlt.empty()  ? model.get_last_zvs_load_threshold()       : zlt.front(); \
            diag["resonantTransitionTime"]        = rtt.empty()  ? model.get_last_resonant_transition_time() : rtt.front(); \
            diag["primaryPeakCurrent"]            = ipk.empty()  ? model.get_last_primary_peak_current()     : ipk.front(); \
            diag["computedSeriesInductance"]      = model.get_computed_series_inductance(); \
            diag["computedOutputInductance"]      = model.get_computed_output_inductance(); \
            diag["computedMagnetizingInductance"] = model.get_computed_magnetizing_inductance(); \
            diag["computedDeadTime"]              = model.get_computed_dead_time(); \
            json perOp = json::array(); \
            for (size_t i = 0; i < dcl.size(); ++i) { \
                json row; \
                row["operatingPointName"]    = (i < names.size()) ? names[i] : ("OP " + std::to_string(i)); \
                row["dutyCycleLoss"] = dcl[i]; row["effectiveDutyCycle"] = deff[i]; \
                row["zvsMarginLagging"] = zvs[i]; row["zvsLoadThreshold"] = zlt[i]; \
                row["resonantTransitionTime"] = rtt[i]; row["primaryPeakCurrent"] = ipk[i]; \
                perOp.push_back(row); \
            } \
            diag["perOp"] = perOp; \
            result[diagKey] = diag; \
        } \
        return result; \
    } \
    catch (const std::exception& exc) { \
        return json{{"error", std::string(#func_name ": ") + exc.what()}}; \
    } \
}

DEFINE_SIMULATE_PHASE_SHIFTED(simulate_psfb_ideal_waveforms, OpenMagnetics::AdvancedPsfb, "psfbDiagnostics")
DEFINE_SIMULATE_PHASE_SHIFTED(simulate_pshb_ideal_waveforms, OpenMagnetics::AdvancedPshb, "pshbDiagnostics")
#undef DEFINE_SIMULATE_PHASE_SHIFTED

// ------ AHB (always AdvancedAsymmetricHalfBridge) ------
json simulate_ahb_ideal_waveforms(json inputsJson) {
    try {
        OpenMagnetics::AdvancedAsymmetricHalfBridge ahbInputs(inputsJson);
        auto inputs = ahbInputs.process();
        auto designRequirements = inputs.get_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designRequirements.get_turns_ratios()) {
            if (tr.get_nominal()) turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
        if (!(magnetizingInductance > 0)) throw std::runtime_error("AHB: no magnetizing inductance");
        OpenMagnetics::NgspiceRunner runner;
        if (!runner.is_available()) throw std::runtime_error("ngspice is not available");
        size_t numberOfPeriods = inputsJson.value("numberOfPeriods", 2);
        size_t numberOfSteadyStatePeriods = inputsJson.value("numberOfSteadyStatePeriods", 3);
        ahbInputs.set_num_periods_to_extract(static_cast<int>(numberOfPeriods));
        ahbInputs.set_num_steady_state_periods(static_cast<int>(numberOfSteadyStatePeriods));
        auto topologyWaveforms = ahbInputs.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance, numberOfPeriods);
        auto operatingPoints = ahbInputs.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
        json result;
        json inputsJ;
        inputsJ["designRequirements"] = json();
        to_json(inputsJ["designRequirements"], designRequirements);
        inputsJ["operatingPoints"] = json::array();
        for (const auto& op : operatingPoints) { json j; to_json(j, op); inputsJ["operatingPoints"].push_back(j); }
        result["inputs"] = inputsJ;
        result["converterWaveforms"] = json::array();
        for (const auto& tw : topologyWaveforms) { json j; to_json(j, tw); result["converterWaveforms"].push_back(j); }
        {
            json diag;
            const auto& names = ahbInputs.get_per_op_name();
            const auto& dc    = ahbInputs.get_per_op_duty_cycle();
            const auto& cr    = ahbInputs.get_per_op_conversion_ratio();
            const auto& cbv   = ahbInputs.get_per_op_dc_blocking_cap_voltage();
            const auto& cbr   = ahbInputs.get_per_op_dc_blocking_cap_ripple();
            const auto& ppvp  = ahbInputs.get_per_op_primary_peak_voltage_positive();
            const auto& ppvn  = ahbInputs.get_per_op_primary_peak_voltage_negative();
            const auto& spvq1 = ahbInputs.get_per_op_switch_peak_voltage_q1();
            const auto& spvq2 = ahbInputs.get_per_op_switch_peak_voltage_q2();
            const auto& sirq1 = ahbInputs.get_per_op_switch_rms_current_q1();
            const auto& sirq2 = ahbInputs.get_per_op_switch_rms_current_q2();
            const auto& zm    = ahbInputs.get_per_op_zvs_margin();
            const auto& rtt   = ahbInputs.get_per_op_resonant_transition_time();
            const auto& ssfe  = ahbInputs.get_per_op_steady_state_flux_excursion();
            const auto& tfee  = ahbInputs.get_per_op_transient_flux_excursion_estimate();
            const auto& mcr   = ahbInputs.get_per_op_magnetizing_current_ripple();
            const auto& oir   = ahbInputs.get_per_op_output_inductor_ripple();
            const auto& om    = ahbInputs.get_per_op_operating_mode();
            const auto& rt    = ahbInputs.get_per_op_rectifier_type();
            diag["operatingMode"]                = om.empty()    ? ahbInputs.get_last_operating_mode()                     : om.front();
            diag["rectifierType"]                = rt.empty()    ? ahbInputs.get_last_rectifier_type()                     : rt.front();
            diag["dutyCycle"]                    = dc.empty()    ? ahbInputs.get_last_duty_cycle()                         : dc.front();
            diag["conversionRatio"]              = cr.empty()    ? ahbInputs.get_last_conversion_ratio()                   : cr.front();
            diag["dcBlockingCapVoltage"]         = cbv.empty()   ? ahbInputs.get_last_dc_blocking_cap_voltage()             : cbv.front();
            diag["dcBlockingCapRipple"]          = cbr.empty()   ? ahbInputs.get_last_dc_blocking_cap_ripple()              : cbr.front();
            diag["primaryPeakVoltagePositive"]   = ppvp.empty()  ? ahbInputs.get_last_primary_peak_voltage_positive()       : ppvp.front();
            diag["primaryPeakVoltageNegative"]   = ppvn.empty()  ? ahbInputs.get_last_primary_peak_voltage_negative()       : ppvn.front();
            diag["switchPeakVoltageQ1"]          = spvq1.empty() ? ahbInputs.get_last_switch_peak_voltage_q1()              : spvq1.front();
            diag["switchPeakVoltageQ2"]          = spvq2.empty() ? ahbInputs.get_last_switch_peak_voltage_q2()              : spvq2.front();
            diag["switchRmsCurrentQ1"]           = sirq1.empty() ? ahbInputs.get_last_switch_rms_current_q1()               : sirq1.front();
            diag["switchRmsCurrentQ2"]           = sirq2.empty() ? ahbInputs.get_last_switch_rms_current_q2()               : sirq2.front();
            diag["zvsMargin"]                    = zm.empty()    ? ahbInputs.get_last_zvs_margin()                          : zm.front();
            diag["resonantTransitionTime"]       = rtt.empty()   ? ahbInputs.get_last_resonant_transition_time()            : rtt.front();
            diag["steadyStateFluxExcursion"]     = ssfe.empty()  ? ahbInputs.get_last_steady_state_flux_excursion()         : ssfe.front();
            diag["transientFluxExcursionEstimate"] = tfee.empty() ? ahbInputs.get_last_transient_flux_excursion_estimate()  : tfee.front();
            diag["magnetizingCurrentRipple"]     = mcr.empty()   ? ahbInputs.get_last_magnetizing_current_ripple()          : mcr.front();
            diag["outputInductorRipple"]         = oir.empty()   ? ahbInputs.get_last_output_inductor_ripple()              : oir.front();
            json perOp = json::array();
            for (size_t i = 0; i < dc.size(); ++i) {
                json row;
                row["operatingPointName"]            = (i < names.size()) ? names[i] : ("OP " + std::to_string(i));
                row["dutyCycle"] = dc[i]; row["conversionRatio"] = cr[i];
                row["dcBlockingCapVoltage"] = cbv[i]; row["dcBlockingCapRipple"] = cbr[i];
                row["primaryPeakVoltagePositive"] = ppvp[i]; row["primaryPeakVoltageNegative"] = ppvn[i];
                row["switchPeakVoltageQ1"] = spvq1[i]; row["switchPeakVoltageQ2"] = spvq2[i];
                row["switchRmsCurrentQ1"] = sirq1[i]; row["switchRmsCurrentQ2"] = sirq2[i];
                row["zvsMargin"] = zm[i]; row["resonantTransitionTime"] = rtt[i];
                row["steadyStateFluxExcursion"] = ssfe[i]; row["transientFluxExcursionEstimate"] = tfee[i];
                row["magnetizingCurrentRipple"] = mcr[i]; row["outputInductorRipple"] = oir[i];
                row["operatingMode"] = om[i]; row["rectifierType"] = rt[i];
                perOp.push_back(row);
            }
            diag["perOp"] = perOp;
            result["ahbDiagnostics"] = diag;
        }
        return result;
    }
    catch (const std::exception& exc) {
        return json{{"error", std::string("simulate_ahb_ideal_waveforms: ") + exc.what()}};
    }
}

// ------ Vienna ------
json simulate_vienna_ideal_waveforms(json inputsJson) {
    try {
        bool isAdvanced = inputsJson.contains("desiredBoostInductance");
        std::unique_ptr<OpenMagnetics::Vienna> model;
        if (isAdvanced) model = std::make_unique<OpenMagnetics::AdvancedVienna>(inputsJson);
        else model = std::make_unique<OpenMagnetics::Vienna>(inputsJson);
        size_t numberOfPeriods = inputsJson.value("numberOfPeriods", 1);
        size_t numberOfSteadyStatePeriods = inputsJson.value("numberOfSteadyStatePeriods", 3);
        model->set_num_periods_to_extract(static_cast<int>(numberOfPeriods));
        model->set_num_steady_state_periods(static_cast<int>(numberOfSteadyStatePeriods));
        auto designRequirements = model->process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designRequirements.get_turns_ratios()) {
            if (tr.get_nominal()) turnsRatios.push_back(tr.get_nominal().value());
        }
        double boostInductance = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
        if (!(boostInductance > 0)) throw std::runtime_error("Vienna: no boost inductance");
        OpenMagnetics::NgspiceRunner runner;
        if (!runner.is_available()) throw std::runtime_error("ngspice is not available");
        auto topologyWaveforms = model->simulate_and_extract_topology_waveforms(turnsRatios, boostInductance, numberOfPeriods);
        auto operatingPoints = model->simulate_and_extract_operating_points(turnsRatios, boostInductance, numberOfPeriods);
        json result;
        json inputsJ;
        inputsJ["designRequirements"] = json();
        to_json(inputsJ["designRequirements"], designRequirements);
        inputsJ["operatingPoints"] = json::array();
        for (const auto& op : operatingPoints) { json j; to_json(j, op); inputsJ["operatingPoints"].push_back(j); }
        result["inputs"] = inputsJ;
        result["converterWaveforms"] = json::array();
        for (const auto& tw : topologyWaveforms) { json j; to_json(j, tw); result["converterWaveforms"].push_back(j); }
        json diag;
        diag["computedBoostInductance"]      = model->get_computed_boost_inductance();
        diag["computedModulationIndex"]      = model->get_computed_modulation_index();
        diag["computedLinePeakCurrent"]      = model->get_computed_line_peak_current();
        diag["computedSwitchVoltageStress"]  = model->get_computed_switch_voltage_stress();
        result["viennaDiagnostics"] = diag;
        return result;
    }
    catch (const std::exception& exc) {
        return json{{"error", std::string("simulate_vienna_ideal_waveforms: ") + exc.what()}};
    }
}

// ------ PFC ------
json simulate_pfc_waveforms(json inputsJson) {
    try {
        OpenMagnetics::PowerFactorCorrection pfcInputs(inputsJson);
        double inductance;
        if (inputsJson.contains("inductance")) {
            inductance = inputsJson["inductance"].get<double>();
        } else {
            std::string mode = inputsJson.value("mode", "Continuous Conduction Mode");
            if (mode == "Continuous Conduction Mode") inductance = pfcInputs.calculate_inductance_ccm();
            else if (mode == "Critical Conduction Mode") inductance = pfcInputs.calculate_inductance_crcm();
            else inductance = pfcInputs.calculate_inductance_dcm();
        }
        auto designRequirements = pfcInputs.process_design_requirements();
        double dcResistance = inputsJson.value("dcResistance", 0.1);
        int numberOfCycles = inputsJson.value("numberOfPeriods", 2);
        pfcInputs.set_num_periods_to_extract(numberOfCycles);
        auto simWaveforms = pfcInputs.simulate_and_extract_waveforms(inductance, dcResistance, numberOfCycles);
        auto operatingPoints = pfcInputs.simulate_and_extract_operating_points(inductance, dcResistance);
        json result;
        result["inductance"] = inductance;
        json inputs;
        inputs["designRequirements"] = json();
        to_json(inputs["designRequirements"], designRequirements);
        inputs["operatingPoints"] = json::array();
        for (const auto& op : operatingPoints) { json j; to_json(j, op); inputs["operatingPoints"].push_back(j); }
        result["inputs"] = inputs;
        result["converterWaveforms"] = json::array();
        json convOp;
        convOp["operatingPointName"] = simWaveforms.operatingPointName;
        convOp["switchingFrequency"] = simWaveforms.switchingFrequency;
        if (!simWaveforms.inputVoltage.empty()) {
            json iv; iv["time"] = simWaveforms.time; iv["data"] = simWaveforms.inputVoltage;
            convOp["inputVoltage"] = iv;
        }
        if (!simWaveforms.inputCurrent.empty()) {
            json ic; ic["time"] = simWaveforms.time; ic["data"] = simWaveforms.inputCurrent;
            convOp["inputCurrent"] = ic;
        }
        if (!simWaveforms.outputCurrent.empty()) {
            json oc; oc["time"] = simWaveforms.time; oc["data"] = simWaveforms.outputCurrent;
            convOp["outputCurrents"] = json::array({oc});
        }
        result["converterWaveforms"].push_back(convOp);
        result["powerFactor"] = simWaveforms.powerFactor;
        result["efficiency"] = simWaveforms.efficiency;
        result["currentThd"] = simWaveforms.currentThd;
        pfcInputs.process();
        {
            json diag;
            diag["computedInductance"] = pfcInputs.get_computed_inductance();
            diag["actualMode"]         = pfcInputs.get_computed_actual_mode();
            const auto& names = pfcInputs.get_per_op_name();
            const auto& dc    = pfcInputs.get_per_op_duty_cycle_peak();
            const auto& ipk   = pfcInputs.get_per_op_peak_inductor_current();
            const auto& ir    = pfcInputs.get_per_op_inductor_ripple();
            const auto& irms  = pfcInputs.get_per_op_line_rms_current();
            const auto& pin   = pfcInputs.get_per_op_input_power();
            diag["dutyCyclePeak"]       = dc.empty()   ? pfcInputs.get_last_duty_cycle_peak()       : dc.front();
            diag["peakInductorCurrent"] = ipk.empty()  ? pfcInputs.get_last_peak_inductor_current() : ipk.front();
            diag["inductorRipple"]      = ir.empty()   ? pfcInputs.get_last_inductor_ripple()       : ir.front();
            diag["lineRmsCurrent"]      = irms.empty() ? pfcInputs.get_last_line_rms_current()      : irms.front();
            diag["inputPower"]          = pin.empty()  ? pfcInputs.get_last_input_power()           : pin.front();
            json perOp = json::array();
            for (size_t i = 0; i < dc.size(); ++i) {
                json row;
                row["operatingPointName"]  = (i < names.size()) ? names[i] : ("OP " + std::to_string(i));
                row["dutyCyclePeak"] = dc[i]; row["peakInductorCurrent"] = ipk[i];
                row["inductorRipple"] = ir[i]; row["lineRmsCurrent"] = irms[i]; row["inputPower"] = pin[i];
                perOp.push_back(row);
            }
            diag["perOp"] = perOp;
            result["pfcDiagnostics"] = diag;
        }
        return result;
    }
    catch (const std::exception& exc) {
        return json{{"error", std::string("simulate_pfc_waveforms: ") + exc.what()}};
    }
}

#undef DEFINE_SIMULATE_NON_ISOLATED
#undef DEFINE_SIMULATE_ISOLATED

// ============================================================================
// Phase B generate_*_ngspice_circuit functions (per-topology SPICE generators).
//
// These are the per-topology convenience wrappers that match the WASM API.
// They take a converter JSON + optional vin/op indices and return a raw
// std::string netlist. The unified generate_ngspice_circuit() above already
// dispatches all topologies; these are individual entry points for parity.
// ============================================================================

template<typename ConverterType, typename AdvancedConverterType>
std::string generate_converter_ngspice_circuit_helper_py(json inputsJson, size_t inputVoltageIndex, size_t operatingPointIndex, const std::string& desiredFieldName) {
    try {
        bool isAdvanced = inputsJson.contains(desiredFieldName);
        std::vector<double> turnsRatios;
        double magnetizingInductance;
        std::string netlist;
        if (isAdvanced) {
            AdvancedConverterType converter(inputsJson);
            magnetizingInductance = converter.get_desired_inductance();
            turnsRatios = converter.get_desired_turns_ratios();
            netlist = converter.generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, operatingPointIndex);
        } else {
            ConverterType converter(inputsJson);
            auto designRequirements = converter.process_design_requirements();
            for (const auto& tr : designRequirements.get_turns_ratios()) {
                turnsRatios.push_back(tr.get_nominal().value());
            }
            magnetizingInductance = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
            if (!(magnetizingInductance > 0)) throw std::runtime_error("Unable to calculate magnetizing inductance");
            netlist = converter.generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, operatingPointIndex);
        }
        return netlist;
    } catch (const std::exception& exc) {
        return "Exception: " + std::string{exc.what()};
    }
}

std::string generate_flyback_ngspice_circuit(json inputsJson, size_t inputVoltageIndex, size_t operatingPointIndex) {
    return generate_converter_ngspice_circuit_helper_py<OpenMagnetics::Flyback, OpenMagnetics::AdvancedFlyback>(
        inputsJson, inputVoltageIndex, operatingPointIndex, "desiredInductance");
}

std::string generate_buck_ngspice_circuit(json inputsJson, size_t inputVoltageIndex, size_t operatingPointIndex) {
    try {
        bool isAdvanced = inputsJson.contains("desiredInductance");
        double inductance;
        std::string netlist;
        if (isAdvanced) {
            OpenMagnetics::AdvancedBuck converter(inputsJson);
            inductance = converter.get_desired_inductance();
            netlist = converter.generate_ngspice_circuit(inductance, inputVoltageIndex, operatingPointIndex);
        } else {
            OpenMagnetics::Buck converter(inputsJson);
            auto designRequirements = converter.process_design_requirements();
            inductance = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
            if (!(inductance > 0)) throw std::runtime_error("Unable to calculate inductance");
            netlist = converter.generate_ngspice_circuit(inductance, inputVoltageIndex, operatingPointIndex);
        }
        return netlist;
    } catch (const std::exception& exc) {
        return "Exception: " + std::string{exc.what()};
    }
}

std::string generate_boost_ngspice_circuit(json inputsJson, size_t inputVoltageIndex, size_t operatingPointIndex) {
    try {
        bool isAdvanced = inputsJson.contains("desiredInductance");
        double inductance;
        std::string netlist;
        if (isAdvanced) {
            OpenMagnetics::AdvancedBoost converter(inputsJson);
            inductance = converter.get_desired_inductance();
            netlist = converter.generate_ngspice_circuit(inductance, inputVoltageIndex, operatingPointIndex);
        } else {
            OpenMagnetics::Boost converter(inputsJson);
            auto designRequirements = converter.process_design_requirements();
            inductance = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
            if (!(inductance > 0)) throw std::runtime_error("Unable to calculate inductance");
            netlist = converter.generate_ngspice_circuit(inductance, inputVoltageIndex, operatingPointIndex);
        }
        return netlist;
    } catch (const std::exception& exc) {
        return "Exception: " + std::string{exc.what()};
    }
}

std::string generate_sepic_ngspice_circuit(json inputsJson, size_t inputVoltageIndex, size_t operatingPointIndex) {
    try {
        bool isAdvanced = inputsJson.contains("desiredInductance");
        double inductanceL1;
        std::string netlist;
        if (isAdvanced) {
            OpenMagnetics::AdvancedSepic converter(inputsJson);
            inductanceL1 = converter.get_desired_inductance();
            netlist = converter.generate_ngspice_circuit(inductanceL1, inputVoltageIndex, operatingPointIndex);
        } else {
            OpenMagnetics::Sepic converter(inputsJson);
            auto designRequirements = converter.process_design_requirements();
            inductanceL1 = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
            if (!(inductanceL1 > 0)) throw std::runtime_error("Unable to calculate inductance");
            netlist = converter.generate_ngspice_circuit(inductanceL1, inputVoltageIndex, operatingPointIndex);
        }
        return netlist;
    } catch (const std::exception& exc) {
        return "Exception: " + std::string{exc.what()};
    }
}

std::string generate_forward_ngspice_circuit(json inputsJson, size_t inputVoltageIndex, size_t operatingPointIndex) {
    return generate_converter_ngspice_circuit_helper_py<OpenMagnetics::SingleSwitchForward, OpenMagnetics::AdvancedSingleSwitchForward>(
        inputsJson, inputVoltageIndex, operatingPointIndex, "desiredInductance");
}

std::string generate_two_switch_forward_ngspice_circuit(json inputsJson, size_t inputVoltageIndex, size_t operatingPointIndex) {
    return generate_converter_ngspice_circuit_helper_py<OpenMagnetics::TwoSwitchForward, OpenMagnetics::AdvancedTwoSwitchForward>(
        inputsJson, inputVoltageIndex, operatingPointIndex, "desiredInductance");
}

std::string generate_active_clamp_forward_ngspice_circuit(json inputsJson, size_t inputVoltageIndex, size_t operatingPointIndex) {
    return generate_converter_ngspice_circuit_helper_py<OpenMagnetics::ActiveClampForward, OpenMagnetics::AdvancedActiveClampForward>(
        inputsJson, inputVoltageIndex, operatingPointIndex, "desiredInductance");
}

std::string generate_push_pull_ngspice_circuit(json inputsJson, size_t inputVoltageIndex, size_t operatingPointIndex) {
    return generate_converter_ngspice_circuit_helper_py<OpenMagnetics::PushPull, OpenMagnetics::AdvancedPushPull>(
        inputsJson, inputVoltageIndex, operatingPointIndex, "desiredInductance");
}

std::string generate_isolated_buck_ngspice_circuit(json inputsJson, size_t inputVoltageIndex, size_t operatingPointIndex) {
    return generate_converter_ngspice_circuit_helper_py<OpenMagnetics::IsolatedBuck, OpenMagnetics::AdvancedIsolatedBuck>(
        inputsJson, inputVoltageIndex, operatingPointIndex, "desiredInductance");
}

std::string generate_isolated_buck_boost_ngspice_circuit(json inputsJson, size_t inputVoltageIndex, size_t operatingPointIndex) {
    return generate_converter_ngspice_circuit_helper_py<OpenMagnetics::IsolatedBuckBoost, OpenMagnetics::AdvancedIsolatedBuckBoost>(
        inputsJson, inputVoltageIndex, operatingPointIndex, "desiredInductance");
}

std::string generate_llc_ngspice_circuit(json inputsJson, size_t inputVoltageIndex, size_t operatingPointIndex) {
    try {
        OpenMagnetics::Llc converter(inputsJson);
        if (inputsJson.contains("desiredResonantInductance") && inputsJson["desiredResonantInductance"].is_number())
            converter.set_user_resonant_inductance(inputsJson["desiredResonantInductance"].get<double>());
        if (inputsJson.contains("desiredResonantCapacitance") && inputsJson["desiredResonantCapacitance"].is_number())
            converter.set_user_resonant_capacitance(inputsJson["desiredResonantCapacitance"].get<double>());
        auto designRequirements = converter.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designRequirements.get_turns_ratios()) {
            turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = inputsJson.value("magnetizingInductance", 200e-6);
        return converter.generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, operatingPointIndex);
    } catch (const std::exception& exc) {
        return "Exception: " + std::string{exc.what()};
    }
}

std::string generate_cllc_ngspice_circuit(json inputsJson, size_t inputVoltageIndex, size_t operatingPointIndex) {
    try {
        OpenMagnetics::CllcConverter converter(inputsJson);
        auto designRequirements = converter.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designRequirements.get_turns_ratios()) {
            if (tr.get_nominal()) turnsRatios.push_back(tr.get_nominal().value());
        }
        if (turnsRatios.empty()) throw std::runtime_error("CLLC: no turns ratios");
        auto params = converter.calculate_resonant_parameters();
        double turnsRatio = turnsRatios[0];
        return converter.generate_ngspice_circuit(turnsRatio, params, inputVoltageIndex, operatingPointIndex);
    } catch (const std::exception& exc) {
        return "Exception: " + std::string{exc.what()};
    }
}

std::string generate_src_ngspice_circuit(json inputsJson, size_t inputVoltageIndex, size_t operatingPointIndex) {
    try {
        bool isAdvanced = inputsJson.contains("desiredTurnsRatios") ||
                          inputsJson.contains("desiredResonantInductance") ||
                          inputsJson.contains("desiredResonantCapacitance");
        std::unique_ptr<OpenMagnetics::Src> model;
        if (isAdvanced) model = std::make_unique<OpenMagnetics::AdvancedSrc>(inputsJson);
        else model = std::make_unique<OpenMagnetics::Src>(inputsJson);
        auto designRequirements = model->process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designRequirements.get_turns_ratios()) {
            if (tr.get_nominal()) turnsRatios.push_back(tr.get_nominal().value());
        }
        if (turnsRatios.empty()) throw std::runtime_error("SRC: no turns ratios");
        double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
        if (!(magnetizingInductance > 0)) throw std::runtime_error("SRC: no magnetizing inductance");
        return model->generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, operatingPointIndex);
    } catch (const std::exception& exc) {
        return "Exception: " + std::string{exc.what()};
    }
}

std::string generate_dab_ngspice_circuit(json inputsJson, size_t inputVoltageIndex, size_t operatingPointIndex) {
    try {
        bool isAdvanced = inputsJson.contains("desiredMagnetizingInductance");
        std::vector<double> turnsRatios;
        double magnetizingInductance;
        std::string netlist;
        if (isAdvanced) {
            OpenMagnetics::AdvancedDab converter(inputsJson);
            magnetizingInductance = converter.get_desired_magnetizing_inductance();
            turnsRatios = converter.get_desired_turns_ratios();
            netlist = converter.generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, operatingPointIndex);
        } else {
            OpenMagnetics::Dab converter(inputsJson);
            auto designRequirements = converter.process_design_requirements();
            for (const auto& tr : designRequirements.get_turns_ratios()) {
                turnsRatios.push_back(tr.get_nominal().value());
            }
            magnetizingInductance = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
            if (!(magnetizingInductance > 0)) throw std::runtime_error("Unable to calculate magnetizing inductance");
            netlist = converter.generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, operatingPointIndex);
        }
        return netlist;
    } catch (const std::exception& exc) {
        return "Exception: " + std::string{exc.what()};
    }
}

std::string generate_psfb_ngspice_circuit(json inputsJson, size_t inputVoltageIndex, size_t operatingPointIndex) {
    try {
        bool isAdvanced = inputsJson.contains("desiredMagnetizingInductance");
        std::vector<double> turnsRatios;
        double magnetizingInductance;
        std::string netlist;
        if (isAdvanced) {
            OpenMagnetics::AdvancedPsfb converter(inputsJson);
            magnetizingInductance = converter.get_desired_magnetizing_inductance();
            turnsRatios = converter.get_desired_turns_ratios();
            netlist = converter.generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, operatingPointIndex);
        } else {
            OpenMagnetics::Psfb converter(inputsJson);
            auto designRequirements = converter.process_design_requirements();
            for (const auto& tr : designRequirements.get_turns_ratios()) {
                turnsRatios.push_back(tr.get_nominal().value());
            }
            magnetizingInductance = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
            if (!(magnetizingInductance > 0)) throw std::runtime_error("Unable to calculate magnetizing inductance");
            netlist = converter.generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, operatingPointIndex);
        }
        return netlist;
    } catch (const std::exception& exc) {
        return "Exception: " + std::string{exc.what()};
    }
}

// ============================================================================
// 2026-05 parity completion: per-topology ngspice generators for the remaining
// topologies whose MKF class already exposes generate_ngspice_circuit but had
// no PyOpenMagnetics wrapper. Each mirrors the parameter extraction of the
// matching simulate_*_ideal_waveforms function so the two stay consistent.
// ============================================================================

// --- Single-inductor topologies (Cuk / Zeta / FourSwitchBuckBoost) ---------
std::string generate_cuk_ngspice_circuit(json inputsJson, size_t inputVoltageIndex, size_t operatingPointIndex) {
    try {
        bool isAdvanced = inputsJson.contains("desiredInductance");
        double inductanceL1;
        std::string netlist;
        if (isAdvanced) {
            OpenMagnetics::AdvancedCuk converter(inputsJson);
            inductanceL1 = converter.get_desired_inductance();
            netlist = converter.generate_ngspice_circuit(inductanceL1, inputVoltageIndex, operatingPointIndex);
        } else {
            OpenMagnetics::Cuk converter(inputsJson);
            auto designRequirements = converter.process_design_requirements();
            inductanceL1 = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
            if (!(inductanceL1 > 0)) throw std::runtime_error("Unable to calculate inductance");
            netlist = converter.generate_ngspice_circuit(inductanceL1, inputVoltageIndex, operatingPointIndex);
        }
        return netlist;
    } catch (const std::exception& exc) {
        return "Exception: " + std::string{exc.what()};
    }
}

std::string generate_zeta_ngspice_circuit(json inputsJson, size_t inputVoltageIndex, size_t operatingPointIndex) {
    try {
        bool isAdvanced = inputsJson.contains("desiredInductance");
        double inductanceL1;
        std::string netlist;
        if (isAdvanced) {
            OpenMagnetics::AdvancedZeta converter(inputsJson);
            inductanceL1 = converter.get_desired_inductance();
            netlist = converter.generate_ngspice_circuit(inductanceL1, inputVoltageIndex, operatingPointIndex);
        } else {
            OpenMagnetics::Zeta converter(inputsJson);
            auto designRequirements = converter.process_design_requirements();
            inductanceL1 = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
            if (!(inductanceL1 > 0)) throw std::runtime_error("Unable to calculate inductance");
            netlist = converter.generate_ngspice_circuit(inductanceL1, inputVoltageIndex, operatingPointIndex);
        }
        return netlist;
    } catch (const std::exception& exc) {
        return "Exception: " + std::string{exc.what()};
    }
}

std::string generate_four_switch_buck_boost_ngspice_circuit(json inputsJson, size_t inputVoltageIndex, size_t operatingPointIndex) {
    try {
        bool isAdvanced = inputsJson.contains("desiredInductance");
        double inductance;
        std::string netlist;
        if (isAdvanced) {
            OpenMagnetics::AdvancedFourSwitchBuckBoost converter(inputsJson);
            inductance = converter.get_desired_inductance();
            netlist = converter.generate_ngspice_circuit(inductance, inputVoltageIndex, operatingPointIndex);
        } else {
            OpenMagnetics::FourSwitchBuckBoost converter(inputsJson);
            auto designRequirements = converter.process_design_requirements();
            inductance = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
            if (!(inductance > 0)) throw std::runtime_error("Unable to calculate inductance");
            netlist = converter.generate_ngspice_circuit(inductance, inputVoltageIndex, operatingPointIndex);
        }
        return netlist;
    } catch (const std::exception& exc) {
        return "Exception: " + std::string{exc.what()};
    }
}

// --- Single-turns-ratio transformer topology (Weinberg) --------------------
std::string generate_weinberg_ngspice_circuit(json inputsJson, size_t inputVoltageIndex, size_t operatingPointIndex) {
    try {
        bool isAdvanced = inputsJson.contains("desiredInductance");
        double turnsRatio;
        double magnetizingInductance;
        std::string netlist;
        if (isAdvanced) {
            OpenMagnetics::AdvancedWeinberg converter(inputsJson);
            magnetizingInductance = converter.get_desired_inductance();
            turnsRatio = converter.get_desired_turns_ratio();
            netlist = converter.generate_ngspice_circuit(turnsRatio, magnetizingInductance, inputVoltageIndex, operatingPointIndex);
        } else {
            OpenMagnetics::Weinberg converter(inputsJson);
            auto designRequirements = converter.process_design_requirements();
            turnsRatio = designRequirements.get_turns_ratios()[0].get_nominal().value();
            magnetizingInductance = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
            if (!(magnetizingInductance > 0)) throw std::runtime_error("Unable to calculate magnetizing inductance");
            netlist = converter.generate_ngspice_circuit(turnsRatio, magnetizingInductance, inputVoltageIndex, operatingPointIndex);
        }
        return netlist;
    } catch (const std::exception& exc) {
        return "Exception: " + std::string{exc.what()};
    }
}

// --- Multi-turns-ratio transformer / resonant topologies -------------------
std::string generate_clllc_ngspice_circuit(json inputsJson, size_t inputVoltageIndex, size_t operatingPointIndex) {
    try {
        OpenMagnetics::Clllc converter(inputsJson);
        auto designRequirements = converter.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designRequirements.get_turns_ratios()) {
            if (tr.get_nominal()) turnsRatios.push_back(tr.get_nominal().value());
        }
        if (turnsRatios.empty()) throw std::runtime_error("Clllc: no turns ratios");
        double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
        if (!(magnetizingInductance > 0)) throw std::runtime_error("Clllc: no magnetizing inductance");
        return converter.generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, operatingPointIndex);
    } catch (const std::exception& exc) {
        return "Exception: " + std::string{exc.what()};
    }
}

std::string generate_vienna_ngspice_circuit(json inputsJson, size_t inputVoltageIndex, size_t operatingPointIndex) {
    try {
        bool isAdvanced = inputsJson.contains("desiredBoostInductance");
        std::unique_ptr<OpenMagnetics::Vienna> model;
        if (isAdvanced) model = std::make_unique<OpenMagnetics::AdvancedVienna>(inputsJson);
        else model = std::make_unique<OpenMagnetics::Vienna>(inputsJson);
        auto designRequirements = model->process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designRequirements.get_turns_ratios()) {
            if (tr.get_nominal()) turnsRatios.push_back(tr.get_nominal().value());
        }
        double boostInductance = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
        if (!(boostInductance > 0)) throw std::runtime_error("Vienna: no boost inductance");
        return model->generate_ngspice_circuit(turnsRatios, boostInductance, inputVoltageIndex, operatingPointIndex);
    } catch (const std::exception& exc) {
        return "Exception: " + std::string{exc.what()};
    }
}

std::string generate_pshb_ngspice_circuit(json inputsJson, size_t inputVoltageIndex, size_t operatingPointIndex) {
    try {
        OpenMagnetics::AdvancedPshb converter(inputsJson);
        auto inputs = converter.process();
        auto designRequirements = inputs.get_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designRequirements.get_turns_ratios()) {
            if (tr.get_nominal()) turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
        if (!(magnetizingInductance > 0)) throw std::runtime_error("PSHB: no magnetizing inductance");
        return converter.generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, operatingPointIndex);
    } catch (const std::exception& exc) {
        return "Exception: " + std::string{exc.what()};
    }
}

std::string generate_ahb_ngspice_circuit(json inputsJson, size_t inputVoltageIndex, size_t operatingPointIndex) {
    try {
        OpenMagnetics::AdvancedAsymmetricHalfBridge converter(inputsJson);
        auto inputs = converter.process();
        auto designRequirements = inputs.get_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designRequirements.get_turns_ratios()) {
            if (tr.get_nominal()) turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
        if (!(magnetizingInductance > 0)) throw std::runtime_error("AHB: no magnetizing inductance");
        return converter.generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, operatingPointIndex);
    } catch (const std::exception& exc) {
        return "Exception: " + std::string{exc.what()};
    }
}

// --- PFC (line-frequency model; distinct signature from the DC-DC family) ---
// PFC has no input-voltage / operating-point sweep indices — its SPICE model is
// parameterised by the boost inductance and circuit damping instead. Inductance
// is taken from an explicit "inductance" field or derived from the conduction
// mode, mirroring simulate_pfc_waveforms.
std::string generate_pfc_ngspice_circuit(json inputsJson, double dcResistance, double simulationTime, double timeStep) {
    try {
        OpenMagnetics::PowerFactorCorrection pfcInputs(inputsJson);
        double inductance;
        if (inputsJson.contains("inductance")) {
            inductance = inputsJson["inductance"].get<double>();
        } else {
            std::string mode = inputsJson.value("mode", "Continuous Conduction Mode");
            if (mode == "Continuous Conduction Mode") inductance = pfcInputs.calculate_inductance_ccm();
            else if (mode == "Critical Conduction Mode") inductance = pfcInputs.calculate_inductance_crcm();
            else inductance = pfcInputs.calculate_inductance_dcm();
        }
        if (!(inductance > 0)) throw std::runtime_error("PFC: unable to calculate inductance");
        double dcRes = inputsJson.value("dcResistance", dcResistance);
        return pfcInputs.generate_ngspice_circuit(inductance, dcRes, simulationTime, timeStep);
    } catch (const std::exception& exc) {
        return "Exception: " + std::string{exc.what()};
    }
}

void register_converter_bindings(py::module& m) {
    m.def("process_converter", &process_converter,
        "Process a converter topology specification to Inputs.",
        py::arg("topology_name"), py::arg("converter_json"), py::arg("use_ngspice") = true);
    
    m.def("design_magnetics_from_converter", &design_magnetics_from_converter,
        "Design magnetic components from a converter specification.",
        py::arg("topology_name"), py::arg("converter_json"), 
        py::arg("max_results") = 1, py::arg("core_mode_json") = "available cores",
        py::arg("use_ngspice") = true, py::arg("weights_json") = nullptr);
    
    m.def("process_flyback", &process_flyback, "Process Flyback converter.", py::arg("flyback"));
    m.def("process_buck", &process_buck, "Process Buck converter.", py::arg("buck"));
    m.def("process_boost", &process_boost, "Process Boost converter.", py::arg("boost"));
    m.def("process_single_switch_forward", &process_single_switch_forward, "Process Single-Switch Forward.", py::arg("forward"));
    m.def("process_two_switch_forward", &process_two_switch_forward, "Process Two-Switch Forward.", py::arg("forward"));
    m.def("process_active_clamp_forward", &process_active_clamp_forward, "Process Active Clamp Forward.", py::arg("forward"));
    m.def("process_push_pull", &process_push_pull, "Process Push-Pull converter.", py::arg("push_pull"));
    m.def("process_isolated_buck", &process_isolated_buck, "Process Isolated Buck.", py::arg("isolated_buck"));
    m.def("process_isolated_buck_boost", &process_isolated_buck_boost, "Process Isolated Buck-Boost.", py::arg("isolated_buck_boost"));
    m.def("process_current_transformer", &process_current_transformer, "Process Current Transformer.",
        py::arg("ct"), py::arg("turns_ratio"), py::arg("secondary_resistance") = 0.0);

    // 2026-05 additions: full coverage of MKF's 24 converter topologies.
    m.def("process_cuk", &process_cuk, "Process Cuk converter.", py::arg("cuk"));
    m.def("process_sepic", &process_sepic, "Process SEPIC converter.", py::arg("sepic"));
    m.def("process_zeta", &process_zeta, "Process Zeta converter.", py::arg("zeta"));
    m.def("process_four_switch_buck_boost", &process_four_switch_buck_boost,
        "Process Four-Switch Buck-Boost converter.", py::arg("converter"));
    m.def("process_asymmetric_half_bridge", &process_asymmetric_half_bridge,
        "Process Asymmetric Half-Bridge converter.", py::arg("converter"));
    m.def("process_weinberg", &process_weinberg, "Process Weinberg converter.", py::arg("converter"));
    m.def("process_vienna", &process_vienna, "Process Vienna Rectifier converter.", py::arg("converter"));
    m.def("process_clllc", &process_clllc, "Process CLLLC Resonant converter.", py::arg("converter"));
    m.def("process_src", &process_src, "Process Series Resonant converter (SRC).", py::arg("converter"));

    m.def("generate_ngspice_circuit", &generate_ngspice_circuit,
        "Return the canonical ngspice SPICE deck for the topology at a "
        "given operating point, sized from the magnetic design (turns "
        "ratios + magnetizing inductance). Wrapper around each Topology "
        "subclass's generate_ngspice_circuit method; processes the "
        "topology first so internal state (Lr/Cr for LLC, etc.) is "
        "populated. Returns {'netlist': '<spice>'} or {'error': '...'}.",
        py::arg("topology_name"), py::arg("converter_json"),
        py::arg("turns_ratios"), py::arg("magnetizing_inductance"),
        py::arg("vin_index") = 0, py::arg("op_index") = 0,
        py::arg("bridge_simulation_mode") = std::string(""),
        py::arg("spice_config") = nlohmann::json::object());

    m.def("get_extra_components_inputs", &get_extra_components_inputs,
        "Return the design requirements for extra components a topology brings "
        "alongside its main magnetic — resonant tank Lr/Cr for LLC, snubber "
        "components for PSFB, etc. Each entry is {kind: 'magnetic'|'capacitor', "
        "inputs: <MAS or CAS Inputs JSON>}. Mode 'IDEAL' returns archetypal "
        "requirements; 'REAL' refines using the supplied main magnetic's "
        "leakage / DCR / etc.",
        py::arg("topology_name"), py::arg("converter_json"),
        py::arg("mode") = "IDEAL", py::arg("magnetic_json") = nullptr);

    // Converter input parity aliases — same dispatch as process_converter but
    // exposed under the WASM-canonical name so frontend code that calls
    // `calculate_<topology>_inputs(...)` works against PyOM unchanged.
    #define BIND_CONVERTER_ALIAS(fn) m.def(#fn, &fn, "Topology inputs builder (WASM parity alias).", py::arg("inputs"))
    BIND_CONVERTER_ALIAS(calculate_flyback_inputs);
    BIND_CONVERTER_ALIAS(calculate_advanced_flyback_inputs);
    BIND_CONVERTER_ALIAS(calculate_buck_inputs);
    BIND_CONVERTER_ALIAS(calculate_advanced_buck_inputs);
    BIND_CONVERTER_ALIAS(calculate_boost_inputs);
    BIND_CONVERTER_ALIAS(calculate_advanced_boost_inputs);
    BIND_CONVERTER_ALIAS(calculate_single_switch_forward_inputs);
    BIND_CONVERTER_ALIAS(calculate_advanced_single_switch_forward_inputs);
    BIND_CONVERTER_ALIAS(calculate_two_switch_forward_inputs);
    BIND_CONVERTER_ALIAS(calculate_advanced_two_switch_forward_inputs);
    BIND_CONVERTER_ALIAS(calculate_active_clamp_forward_inputs);
    BIND_CONVERTER_ALIAS(calculate_advanced_active_clamp_forward_inputs);
    BIND_CONVERTER_ALIAS(calculate_push_pull_inputs);
    BIND_CONVERTER_ALIAS(calculate_advanced_push_pull_inputs);
    BIND_CONVERTER_ALIAS(calculate_isolated_buck_inputs);
    BIND_CONVERTER_ALIAS(calculate_advanced_isolated_buck_inputs);
    BIND_CONVERTER_ALIAS(calculate_isolated_buck_boost_inputs);
    BIND_CONVERTER_ALIAS(calculate_advanced_isolated_buck_boost_inputs);
    BIND_CONVERTER_ALIAS(calculate_cuk_inputs);
    BIND_CONVERTER_ALIAS(calculate_advanced_cuk_inputs);
    BIND_CONVERTER_ALIAS(calculate_sepic_inputs);
    BIND_CONVERTER_ALIAS(calculate_advanced_sepic_inputs);
    BIND_CONVERTER_ALIAS(calculate_zeta_inputs);
    BIND_CONVERTER_ALIAS(calculate_advanced_zeta_inputs);
    BIND_CONVERTER_ALIAS(calculate_four_switch_buck_boost_inputs);
    BIND_CONVERTER_ALIAS(calculate_advanced_four_switch_buck_boost_inputs);
    BIND_CONVERTER_ALIAS(calculate_weinberg_inputs);
    BIND_CONVERTER_ALIAS(calculate_advanced_weinberg_inputs);
    BIND_CONVERTER_ALIAS(calculate_clllc_inputs);
    BIND_CONVERTER_ALIAS(calculate_advanced_clllc_inputs);
    BIND_CONVERTER_ALIAS(calculate_vienna_inputs);
    BIND_CONVERTER_ALIAS(calculate_advanced_vienna_inputs);
    BIND_CONVERTER_ALIAS(calculate_src_inputs);
    BIND_CONVERTER_ALIAS(calculate_advanced_src_inputs);
    BIND_CONVERTER_ALIAS(calculate_llc_inputs);
    BIND_CONVERTER_ALIAS(calculate_advanced_llc_inputs);
    BIND_CONVERTER_ALIAS(calculate_cllc_inputs);
    BIND_CONVERTER_ALIAS(calculate_advanced_cllc_inputs);
    BIND_CONVERTER_ALIAS(calculate_dab_inputs);
    BIND_CONVERTER_ALIAS(calculate_advanced_dab_inputs);
    BIND_CONVERTER_ALIAS(calculate_pfc_inputs);
    BIND_CONVERTER_ALIAS(calculate_psfb_inputs);
    BIND_CONVERTER_ALIAS(calculate_advanced_psfb_inputs);
    BIND_CONVERTER_ALIAS(calculate_pshb_inputs);
    BIND_CONVERTER_ALIAS(calculate_advanced_pshb_inputs);
    BIND_CONVERTER_ALIAS(calculate_ahb_inputs);
    BIND_CONVERTER_ALIAS(calculate_advanced_ahb_inputs);
    #undef BIND_CONVERTER_ALIAS

    // CMC / DMC wrappers (parity with WebLibMKF). See definitions above.
    m.def("calculate_cmc_inputs", &calculate_cmc_inputs,
        "Build MAS Inputs (designRequirements + operatingPoints) for a CMC from "
        "wizard params (operatingVoltage, operatingCurrent, lineFrequency, "
        "numberOfWindings, and an impedance/insertion-loss/noise spec). Mirrors "
        "the WASM `calculate_cmc_inputs` used by el-choker's CMC wizard.",
        py::arg("cmc_inputs"));
    m.def("calculate_advanced_cmc_inputs", &calculate_advanced_cmc_inputs,
        "Advanced-mode CMC inputs builder (user supplies `desiredInductance` + "
        "`designFrequency` directly rather than a spec).",
        py::arg("cmc_inputs"));
    m.def("generate_cmc_ngspice_circuit", &generate_cmc_ngspice_circuit,
        "Generate the ngspice SPICE deck for a CMC at its design inductance.",
        py::arg("cmc_inputs"));
    m.def("simulate_cmc_lisn_waveforms", &simulate_cmc_lisn_waveforms,
        "Run the CISPR-standard LISN test simulation against the CMC at the "
        "requested impedance points. Returns the input MAS plus per-frequency "
        "waveforms (input voltage, winding currents, LISN voltage, CM attenuation).",
        py::arg("cmc_inputs"), py::arg("inductance"));
    m.def("simulate_cmc_ideal_waveforms", &simulate_cmc_ideal_waveforms,
        "Realistic CMC operating-point simulation: line voltage + switching "
        "noise (parametrised by parasitic capacitance / dV/dt). Mirrors the WASM "
        "function used by el-choker to produce operating points the MagneticAdviser "
        "needs to score loss-aware filters.",
        py::arg("cmc_inputs"), py::arg("inductance"),
        py::arg("parasitic_capacitance_pF") = 10.0, py::arg("dvdt_V_per_ns") = 50.0);
    m.def("calculate_dmc_inputs", &calculate_dmc_inputs,
        "Build MAS Inputs for a DMC from wizard params (differential-mode equivalent "
        "of calculate_cmc_inputs).",
        py::arg("dmc_inputs"));
    m.def("verify_dmc_attenuation", &verify_dmc_attenuation,
        "Verify a DMC + filter capacitor pair meets the user attenuation spec.",
        py::arg("dmc_inputs"), py::arg("inductance"), py::arg("capacitance") = 0.0);
    m.def("propose_dmc_design", &propose_dmc_design,
        "Propose a DMC inductance / capacitance pair satisfying the spec.",
        py::arg("dmc_inputs"));
    m.def("simulate_dmc_waveforms", &simulate_dmc_waveforms,
        "Simulate DMC time-domain waveforms (input/output voltage, inductor "
        "current, DM attenuation) at the EMI test frequencies.",
        py::arg("dmc_inputs"), py::arg("inductance"));
    m.def("generate_dmc_ngspice_circuit", &generate_dmc_ngspice_circuit,
        "Generate the ngspice SPICE deck for a DMC. Reads minimumInductance from "
        "the inputs JSON or falls back to propose_design().",
        py::arg("dmc_inputs"));

    // Phase B: simulate_*_ideal_waveforms (ngspice-based waveform + operating point extraction)
    m.def("simulate_flyback_ideal_waveforms", &simulate_flyback_ideal_waveforms,
        "Simulate Flyback ideal waveforms via ngspice.", py::arg("inputs"));
    m.def("simulate_flyback_with_magnetic", &simulate_flyback_with_magnetic,
        "Simulate Flyback with a pre-built magnetic JSON.", py::arg("inputs"), py::arg("magnetic"));
    m.def("simulate_buck_ideal_waveforms", &simulate_buck_ideal_waveforms,
        "Simulate Buck ideal waveforms via ngspice.", py::arg("inputs"));
    m.def("simulate_boost_ideal_waveforms", &simulate_boost_ideal_waveforms,
        "Simulate Boost ideal waveforms via ngspice.", py::arg("inputs"));
    m.def("simulate_sepic_ideal_waveforms", &simulate_sepic_ideal_waveforms,
        "Simulate SEPIC ideal waveforms via ngspice.", py::arg("inputs"));
    m.def("simulate_cuk_ideal_waveforms", &simulate_cuk_ideal_waveforms,
        "Simulate Cuk ideal waveforms via ngspice.", py::arg("inputs"));
    m.def("simulate_zeta_ideal_waveforms", &simulate_zeta_ideal_waveforms,
        "Simulate Zeta ideal waveforms via ngspice.", py::arg("inputs"));
    m.def("simulate_four_switch_buck_boost_ideal_waveforms", &simulate_four_switch_buck_boost_ideal_waveforms,
        "Simulate Four-Switch Buck-Boost ideal waveforms via ngspice.", py::arg("inputs"));
    m.def("simulate_forward_ideal_waveforms", &simulate_forward_ideal_waveforms,
        "Simulate Single-Switch Forward ideal waveforms via ngspice.", py::arg("inputs"));
    m.def("simulate_two_switch_forward_ideal_waveforms", &simulate_two_switch_forward_ideal_waveforms,
        "Simulate Two-Switch Forward ideal waveforms via ngspice.", py::arg("inputs"));
    m.def("simulate_active_clamp_forward_ideal_waveforms", &simulate_active_clamp_forward_ideal_waveforms,
        "Simulate Active-Clamp Forward ideal waveforms via ngspice.", py::arg("inputs"));
    m.def("simulate_push_pull_ideal_waveforms", &simulate_push_pull_ideal_waveforms,
        "Simulate Push-Pull ideal waveforms via ngspice.", py::arg("inputs"));
    m.def("simulate_isolated_buck_ideal_waveforms", &simulate_isolated_buck_ideal_waveforms,
        "Simulate Isolated Buck ideal waveforms via ngspice.", py::arg("inputs"));
    m.def("simulate_isolated_buck_boost_ideal_waveforms", &simulate_isolated_buck_boost_ideal_waveforms,
        "Simulate Isolated Buck-Boost ideal waveforms via ngspice.", py::arg("inputs"));
    m.def("simulate_weinberg_ideal_waveforms", &simulate_weinberg_ideal_waveforms,
        "Simulate Weinberg ideal waveforms via ngspice.", py::arg("inputs"));
    m.def("simulate_llc_ideal_waveforms", &simulate_llc_ideal_waveforms,
        "Simulate LLC ideal waveforms via ngspice.", py::arg("inputs"));
    m.def("simulate_cllc_ideal_waveforms", &simulate_cllc_ideal_waveforms,
        "Simulate CLLC ideal waveforms via ngspice.", py::arg("inputs"));
    m.def("simulate_clllc_ideal_waveforms", &simulate_clllc_ideal_waveforms,
        "Simulate CLLLC ideal waveforms via ngspice.", py::arg("inputs"));
    m.def("simulate_src_ideal_waveforms", &simulate_src_ideal_waveforms,
        "Simulate SRC ideal waveforms via ngspice.", py::arg("inputs"));
    m.def("simulate_dab_ideal_waveforms", &simulate_dab_ideal_waveforms,
        "Simulate DAB ideal waveforms via ngspice.", py::arg("inputs"));
    m.def("simulate_psfb_ideal_waveforms", &simulate_psfb_ideal_waveforms,
        "Simulate PSFB ideal waveforms via ngspice.", py::arg("inputs"));
    m.def("simulate_pshb_ideal_waveforms", &simulate_pshb_ideal_waveforms,
        "Simulate PSHB ideal waveforms via ngspice.", py::arg("inputs"));
    m.def("simulate_ahb_ideal_waveforms", &simulate_ahb_ideal_waveforms,
        "Simulate Asymmetric Half-Bridge ideal waveforms via ngspice.", py::arg("inputs"));
    m.def("simulate_vienna_ideal_waveforms", &simulate_vienna_ideal_waveforms,
        "Simulate Vienna Rectifier ideal waveforms via ngspice.", py::arg("inputs"));
    m.def("simulate_pfc_waveforms", &simulate_pfc_waveforms,
        "Simulate PFC waveforms via ngspice.", py::arg("inputs"));

    // Phase B: generate_*_ngspice_circuit (per-topology SPICE netlist generators)
    m.def("generate_flyback_ngspice_circuit", &generate_flyback_ngspice_circuit,
        "Generate Flyback ngspice circuit netlist.",
        py::arg("inputs"), py::arg("input_voltage_index") = 0, py::arg("operating_point_index") = 0);
    m.def("generate_buck_ngspice_circuit", &generate_buck_ngspice_circuit,
        "Generate Buck ngspice circuit netlist.",
        py::arg("inputs"), py::arg("input_voltage_index") = 0, py::arg("operating_point_index") = 0);
    m.def("generate_boost_ngspice_circuit", &generate_boost_ngspice_circuit,
        "Generate Boost ngspice circuit netlist.",
        py::arg("inputs"), py::arg("input_voltage_index") = 0, py::arg("operating_point_index") = 0);
    m.def("generate_sepic_ngspice_circuit", &generate_sepic_ngspice_circuit,
        "Generate SEPIC ngspice circuit netlist.",
        py::arg("inputs"), py::arg("input_voltage_index") = 0, py::arg("operating_point_index") = 0);
    m.def("generate_forward_ngspice_circuit", &generate_forward_ngspice_circuit,
        "Generate Single-Switch Forward ngspice circuit netlist.",
        py::arg("inputs"), py::arg("input_voltage_index") = 0, py::arg("operating_point_index") = 0);
    m.def("generate_two_switch_forward_ngspice_circuit", &generate_two_switch_forward_ngspice_circuit,
        "Generate Two-Switch Forward ngspice circuit netlist.",
        py::arg("inputs"), py::arg("input_voltage_index") = 0, py::arg("operating_point_index") = 0);
    m.def("generate_active_clamp_forward_ngspice_circuit", &generate_active_clamp_forward_ngspice_circuit,
        "Generate Active-Clamp Forward ngspice circuit netlist.",
        py::arg("inputs"), py::arg("input_voltage_index") = 0, py::arg("operating_point_index") = 0);
    m.def("generate_push_pull_ngspice_circuit", &generate_push_pull_ngspice_circuit,
        "Generate Push-Pull ngspice circuit netlist.",
        py::arg("inputs"), py::arg("input_voltage_index") = 0, py::arg("operating_point_index") = 0);
    m.def("generate_isolated_buck_ngspice_circuit", &generate_isolated_buck_ngspice_circuit,
        "Generate Isolated Buck ngspice circuit netlist.",
        py::arg("inputs"), py::arg("input_voltage_index") = 0, py::arg("operating_point_index") = 0);
    m.def("generate_isolated_buck_boost_ngspice_circuit", &generate_isolated_buck_boost_ngspice_circuit,
        "Generate Isolated Buck-Boost ngspice circuit netlist.",
        py::arg("inputs"), py::arg("input_voltage_index") = 0, py::arg("operating_point_index") = 0);
    m.def("generate_llc_ngspice_circuit", &generate_llc_ngspice_circuit,
        "Generate LLC ngspice circuit netlist.",
        py::arg("inputs"), py::arg("input_voltage_index") = 0, py::arg("operating_point_index") = 0);
    m.def("generate_cllc_ngspice_circuit", &generate_cllc_ngspice_circuit,
        "Generate CLLC ngspice circuit netlist.",
        py::arg("inputs"), py::arg("input_voltage_index") = 0, py::arg("operating_point_index") = 0);
    m.def("generate_src_ngspice_circuit", &generate_src_ngspice_circuit,
        "Generate SRC ngspice circuit netlist.",
        py::arg("inputs"), py::arg("input_voltage_index") = 0, py::arg("operating_point_index") = 0);
    m.def("generate_dab_ngspice_circuit", &generate_dab_ngspice_circuit,
        "Generate DAB ngspice circuit netlist.",
        py::arg("inputs"), py::arg("input_voltage_index") = 0, py::arg("operating_point_index") = 0);
    m.def("generate_psfb_ngspice_circuit", &generate_psfb_ngspice_circuit,
        "Generate PSFB ngspice circuit netlist.",
        py::arg("inputs"), py::arg("input_voltage_index") = 0, py::arg("operating_point_index") = 0);
    m.def("generate_cuk_ngspice_circuit", &generate_cuk_ngspice_circuit,
        "Generate Cuk ngspice circuit netlist.",
        py::arg("inputs"), py::arg("input_voltage_index") = 0, py::arg("operating_point_index") = 0);
    m.def("generate_zeta_ngspice_circuit", &generate_zeta_ngspice_circuit,
        "Generate Zeta ngspice circuit netlist.",
        py::arg("inputs"), py::arg("input_voltage_index") = 0, py::arg("operating_point_index") = 0);
    m.def("generate_four_switch_buck_boost_ngspice_circuit", &generate_four_switch_buck_boost_ngspice_circuit,
        "Generate Four-Switch Buck-Boost ngspice circuit netlist.",
        py::arg("inputs"), py::arg("input_voltage_index") = 0, py::arg("operating_point_index") = 0);
    m.def("generate_weinberg_ngspice_circuit", &generate_weinberg_ngspice_circuit,
        "Generate Weinberg ngspice circuit netlist.",
        py::arg("inputs"), py::arg("input_voltage_index") = 0, py::arg("operating_point_index") = 0);
    m.def("generate_clllc_ngspice_circuit", &generate_clllc_ngspice_circuit,
        "Generate CLLLC ngspice circuit netlist.",
        py::arg("inputs"), py::arg("input_voltage_index") = 0, py::arg("operating_point_index") = 0);
    m.def("generate_vienna_ngspice_circuit", &generate_vienna_ngspice_circuit,
        "Generate Vienna Rectifier ngspice circuit netlist.",
        py::arg("inputs"), py::arg("input_voltage_index") = 0, py::arg("operating_point_index") = 0);
    m.def("generate_pshb_ngspice_circuit", &generate_pshb_ngspice_circuit,
        "Generate PSHB ngspice circuit netlist.",
        py::arg("inputs"), py::arg("input_voltage_index") = 0, py::arg("operating_point_index") = 0);
    m.def("generate_ahb_ngspice_circuit", &generate_ahb_ngspice_circuit,
        "Generate Asymmetric Half-Bridge ngspice circuit netlist.",
        py::arg("inputs"), py::arg("input_voltage_index") = 0, py::arg("operating_point_index") = 0);
    m.def("generate_pfc_ngspice_circuit", &generate_pfc_ngspice_circuit,
        "Generate PFC ngspice circuit netlist (line-frequency model; no sweep indices).",
        py::arg("inputs"), py::arg("dc_resistance") = 0.1, py::arg("simulation_time") = 0.02,
        py::arg("time_step") = 1e-8);
}

} // namespace PyMKF
