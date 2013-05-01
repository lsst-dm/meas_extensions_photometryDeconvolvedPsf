// -*- lsst-c++ -*-
#ifndef LSST_MEAS_EXTENSIONS_PHOTOMETRY_DeconvolvedPsf_H
#define LSST_MEAS_EXTENSIONS_PHOTOMETRY_DeconvolvedPsf_H

#include "lsst/meas/algorithms/FluxControl.h"

namespace lsst { namespace meas { namespace extensions { namespace photometryDeconvolvedPsf {

/**
 *  @brief C++ control object for DeconvolvedPsf flux.
 *
 *  @sa DeconvolvedPsfFluxConfig.
 */
class DeconvolvedPsfFluxControl : public algorithms::FluxControl {
public:

    LSST_CONTROL_FIELD(fixed, bool,
                       "if true, use existing shape and centroid measurements instead of fitting");
    LSST_CONTROL_FIELD(nSigmaForRadius, double,
                       "Multiplier of rms size for aperture used to initially estimate the DeconvolvedPsf radius");
    LSST_CONTROL_FIELD(nIterForRadius, int, "Number of times to iterate when setting the DeconvolvedPsf radius");
    LSST_CONTROL_FIELD(nRadiusForFlux, double, "Number of DeconvolvedPsf radii for DeconvolvedPsf flux");
    LSST_CONTROL_FIELD(maxSincRadius, double,
                       "Largest aperture for which to use the slow, accurate, sinc aperture code");
    LSST_CONTROL_FIELD(minimumRadius, double,
                       "Minimum DeconvolvedPsf radius (if == 0.0 use PSF's DeconvolvedPsf radius). "
                       "Ignored if enforceMinimumRadius is false");
    LSST_CONTROL_FIELD(enforceMinimumRadius, bool, "If true check that the DeconvolvedPsf radius exceeds some minimum");
    LSST_CONTROL_FIELD(useFootprintRadius, bool,
                       "Use the Footprint size as part of initial estimate of DeconvolvedPsf radius");
    LSST_CONTROL_FIELD(smoothingSigma, double,
                       "Smooth image with N(0, smoothingSigma^2) Gaussian while estimating R_K");

    DeconvolvedPsfFluxControl() : 
        algorithms::FluxControl("flux.deconvolvedPsf"), fixed(false),
        nSigmaForRadius(6.0),
        nIterForRadius(1),
        nRadiusForFlux(2.5),
        maxSincRadius(10.0),
        minimumRadius(0.0),
        enforceMinimumRadius(true),
        useFootprintRadius(false),
        smoothingSigma(-1.0)
    {}

private:
    virtual PTR(algorithms::AlgorithmControl) _clone() const;
    virtual PTR(algorithms::Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

}}}} // namespace lsst::meas::extensions::photometryDeconvolvedPsf

#endif // !LSST_MEAS_EXTENSIONS_PHOTOMETRY_DeconvolvedPsf_H
