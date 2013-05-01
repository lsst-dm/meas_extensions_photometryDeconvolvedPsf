// -*- LSST-C++ -*-
#include <numeric>
#include <cmath>
#include <functional>
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/image.h"
#include "lsst/afw/math/Integrate.h"
#include "lsst/meas/algorithms/Measure.h"

#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/FootprintFunctor.h"
#include "lsst/meas/algorithms/FluxControl.h"
#include "lsst/meas/algorithms/ScaledFlux.h"
#include "lsst/meas/extensions/photometryDeconvolvedPsf.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwDetection = lsst::afw::detection;
namespace afwGeom = lsst::afw::geom;
namespace afwImage = lsst::afw::image;
namespace afwMath = lsst::afw::math;

namespace lsst { namespace meas { namespace extensions { namespace photometryDeconvolvedPsf {

namespace {

/**
 * @brief A class that knows how to calculate fluxes using the PSF photometry algorithm
 * @ingroup meas/algorithms
 */
class DeconvolvedPsfFlux : public algorithms::FluxAlgorithm, public algorithms::ScaledFlux {
public:

    DeconvolvedPsfFlux(DeconvolvedPsfFluxControl const & ctrl, afw::table::Schema & schema) :
        FluxAlgorithm(ctrl, schema, "flux measured by a fit to the PSF model"),
        _fluxCorrectionKeys(ctrl.name, schema)
    {}

    virtual afw::table::KeyTuple<afw::table::Flux> getFluxKeys(int n=0) const {
        return FluxAlgorithm::getKeys();
    }

    virtual ScaledFlux::KeyTuple getFluxCorrectionKeys(int n=0) const {
        return _fluxCorrectionKeys;
    }

private:
    
    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(DeconvolvedPsfFlux);

    ScaledFlux::KeyTuple _fluxCorrectionKeys;
};

/**
 * Accumulate sum(x) and sum(x**2)
 */
template<typename T>
struct getSum2 {
    getSum2() : sum(0.0), sum2(0.0) {}
    
    getSum2& operator+(T x) {
        sum += x;
        sum2 += x*x;
        
        return *this;
    }
    
    double sum;                         // \sum_i(x_i)
    double sum2;                        // \sum_i(x_i^2)
};

template <typename TargetImageT, typename WeightImageT>
class FootprintWeightFlux : public afwDetection::FootprintFunctor<TargetImageT> {
public:
    FootprintWeightFlux(TargetImageT const& mimage, ///< The image the source lives in
                        PTR(WeightImageT) wimage    ///< The weight image
                       ) : afwDetection::FootprintFunctor<TargetImageT>(mimage),
                           _wimage(wimage),
                           _sum(0), _sumVar(0), _x0(0), _y0(0) {}
    
    /// @brief Reset everything for a new Footprint
    void reset() {}        
    void reset(afwDetection::Footprint const& foot) {
        _sumVar = _sum = 0.0;

        afwGeom::BoxI const& bbox(foot.getBBox());
        _x0 = bbox.getMinX();
        _y0 = bbox.getMinY();

        if (bbox.getDimensions() != _wimage->getDimensions()) {
            throw LSST_EXCEPT(pexExceptions::LengthErrorException,
                              (boost::format("Footprint at %d,%d -- %d,%d is wrong size "
                                             "for %d x %d weight image") %
                               bbox.getMinX() % bbox.getMinY() % bbox.getMaxX() % bbox.getMaxY() %
                               _wimage->getWidth() % _wimage->getHeight()).str());
        }
    }
    
    /// @brief method called for each pixel by apply()
    virtual void operator()(typename TargetImageT::xy_locator iloc, ///< locator pointing at the image pixel
                    int x,                                 ///< column-position of pixel
                    int y                                  ///< row-position of pixel
                   ) {
        _callImpl(iloc, x, y, typename TargetImageT::image_category());
    }

    /// Return the Footprint's flux
    double getSum() const { return _sum; }

    /// Return the variance of the Footprint's flux
    double getSumVar() const { return _sumVar; }

private:

    template <typename LocatorT>
    void _callImpl(LocatorT iloc, int x, int y, afw::image::detail::MaskedImage_tag) {
        double ival = iloc.image(0, 0);
        double vval = iloc.variance(0, 0);
        double wval = (*_wimage)(x - _x0, y - _y0);
        _sum += wval*ival;
        _sumVar += wval*wval*vval;        
    }

    template <typename LocatorT>
    void _callImpl(LocatorT iloc, int x, int y, afw::image::detail::Image_tag) {
        double ival = *iloc;
        double wval = (*_wimage)(x - _x0, y - _y0);
        _sum += wval * ival;
    }

    typename WeightImageT::Ptr const& _wimage;        // The weight image
    double _sum;                                      // our desired sum
    double _sumVar;
    int _x0, _y0;                                     // the origin of the current Footprint
};

template <typename TargetImageT>
std::pair<double,double> computeDeconvolvedPsfFlux(
    TargetImageT const image,
    PTR(afw::detection::Psf::Image) const & wimage,
    afw::geom::Point2D const & center
) {
    afwGeom::BoxI imageBBox(image.getBBox(afwImage::PARENT));
    FootprintWeightFlux<TargetImageT, afwDetection::Psf::Image> wfluxFunctor(image, wimage);
    // Build a rectangular Footprint corresponding to wimage
    afwDetection::Footprint foot(wimage->getBBox(afwImage::PARENT), imageBBox);
    wfluxFunctor.apply(foot);
    
    getSum2<afwDetection::Psf::Pixel> sum;
    sum = std::accumulate(wimage->begin(true), wimage->end(true), sum);
    
    double flux = wfluxFunctor.getSum()*sum.sum/sum.sum2;
    double fluxErr = ::sqrt(wfluxFunctor.getSumVar())*::fabs(sum.sum)/sum.sum2;
    return std::make_pair(flux, fluxErr);
}

/************************************************************************************************************/
/**
 * Calculate the desired psf flux
 */
template <typename PixelT>
void DeconvolvedPsfFlux::_apply(
    afw::table::SourceRecord & source, 
    afw::image::Exposure<PixelT> const& exposure,
    afw::geom::Point2D const & center
) const {
    source.set(getKeys().flag, true); // say we've failed so that's the result if we throw
    source.set(_fluxCorrectionKeys.psfFactorFlag, true);
    
    CONST_PTR(afwDetection::Psf) psf = exposure.getPsf();
    if (!psf) {
        throw LSST_EXCEPT(pexExceptions::RuntimeErrorException, "No PSF provided for PSF photometry");
    }

    PTR(afwDetection::Psf::Image) psfImage;
    try {
        psfImage = psf->computeImage(center);
    } catch (lsst::pex::exceptions::Exception & e) {
        LSST_EXCEPT_ADD(e, (boost::format("Computing PSF at (%.3f, %.3f)")
                            % center.getX() % center.getY()).str());
        throw e;
    }

    std::pair<double,double> result = computeDeconvolvedPsfFlux(exposure.getMaskedImage(), psfImage, center);
    source.set(getKeys().meas, result.first);
    source.set(getKeys().err, result.second);
    source.set(getKeys().flag, false);

    // The logic here could be a lot more efficient, because the weight image is the data image, but
    // for now it's more important to be absolutely certain we apply the exact same algorithm as
    // we did above.
    std::pair<double,double> psfResult = computeDeconvolvedPsfFlux(*psfImage, psfImage, center);
    source.set(_fluxCorrectionKeys.psfFactor, psfResult.first);
    source.set(_fluxCorrectionKeys.psfFactorFlag, false);

}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(DeconvolvedPsfFlux);

} // anonymous

PTR(algorithms::AlgorithmControl) DeconvolvedPsfFluxControl::_clone() const {
    return boost::make_shared<DeconvolvedPsfFluxControl>(*this);
}

PTR(algorithms::Algorithm) DeconvolvedPsfFluxControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const &
) const {
    return boost::make_shared<DeconvolvedPsfFlux>(*this, boost::ref(schema));
}

}}}}
