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

    LSST_CONTROL_FIELD(deconvolutionKernelSigma, double,
                       "Convolve with an N(0, deconvolutionKernelSigma^2) Gaussian to estimate kernel");
    LSST_CONTROL_FIELD(psfFlux, double, "Flux of the objects used to determine the PSF");
    LSST_CONTROL_FIELD(flux0, double, "Reference flux (see coeff)");
    LSST_CONTROL_FIELD(coeff, double, "Intensity smoothing is N(0, coeff*log10(apFlux/flux0))");
    LSST_CONTROL_FIELD(niter, int, "Maximum number of iterations");
    LSST_CONTROL_FIELD(rmsTol, double, "Tolerance in width when iterating");

    DeconvolvedPsfFluxControl() : 
        algorithms::FluxControl("flux.deconvolvedPsf"),
        deconvolutionKernelSigma(0.4),
        psfFlux(1e5),
        flux0(1e3),
        coeff(0.0),
        niter(15),
        rmsTol(1e-4)
    {}

private:
    virtual PTR(algorithms::AlgorithmControl) _clone() const;
    virtual PTR(algorithms::Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

}}}} // namespace lsst::meas::extensions::photometryDeconvolvedPsf

#endif // !LSST_MEAS_EXTENSIONS_PHOTOMETRY_DeconvolvedPsf_H
