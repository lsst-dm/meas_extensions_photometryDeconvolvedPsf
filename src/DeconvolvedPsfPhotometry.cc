// -*- LSST-C++ -*-
#include <numeric>
#include <cmath>
#include <functional>
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/image.h"
#include "lsst/afw/math/ConvolveImage.h"
#include "lsst/afw/math/FunctionLibrary.h"
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
                              str(boost::format("Footprint at %d,%d -- %d,%d is wrong size "
                                                "for %d x %d weight image") %
                                  bbox.getMinX() % bbox.getMinY() % bbox.getMaxX() % bbox.getMaxY() %
                                  _wimage->getWidth() % _wimage->getHeight()));
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

template<typename ImageT>
double
calcRms(ImageT const& img)
{
    const int width = img.getWidth();
    const int height = img.getHeight();

    double const xcen = width/2;
    double const ycen = height/2;

    double sum = 0, sumrr = 0;
    for (int y = 0; y != height; ++y) {
        for (int x = 0; x != width; ++x) {
            double const val = img(x, y);
            sum += val;
            sumrr += val*(pow(x - xcen, 2) + pow(y - ycen, 2));
        }
    }
    sumrr /= sum;

    return ::sqrt(sumrr/2);
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
    DeconvolvedPsfFluxControl const & ctrl =
        static_cast<DeconvolvedPsfFluxControl const &>(this->getControl());

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
        LSST_EXCEPT_ADD(e, str(boost::format("Computing PSF at (%.3f, %.3f)")
                               % center.getX() % center.getY()));
        throw e;
    }
    //
    // Allow for the variation of the PSF with intensity
    //
    if (source.getApFluxFlag()) {
#if 1                                   // XXXXXXXXXXXXXXXXXXXXX
        throw LSST_EXCEPT(lsst::pex::exceptions::NotFoundException,
                          str(boost::format("Aperture flux is invalid for object at (%.3f, %.3f)")
                              % center.getX() % center.getY()));
#endif
    }

    double const objFlux = source.getPsfFlux(); // object's flux
    double const flux0 = ctrl.flux0;            // fiducial flux (see coeff)
    double const coeff = ctrl.coeff;            // intensity kernel is N(0, coeff*log10(flux/ctrl.flux0)^2)
    double const psfFlux = ctrl.psfFlux;        // characteristic flux of the PSF
    double const epsilon2 =
        ::pow(coeff*::log10(objFlux/flux0), 2) - 
        ::pow(coeff*::log10(psfFlux/flux0), 2); // (correction we need to apply)^2
 
    double const psfRms = calcRms(*psfImage);

    afwDetection::Psf::Image smoothed(psfImage->getDimensions());
    bool const normalize = true;
    bool const copyEdge = true;
    if (epsilon2 > 0) {                  // we need to smooth our PSF model to match the object
        double const psfRmsCorrected = ::sqrt(psfRms*psfRms + epsilon2);
        //
        // Iterate until we get the desired PSF width
        //
        double a = ::sqrt(epsilon2);    // width of needed smoothing kernel
        for (int i = 0; i < ctrl.niter; ++i) {
            a = ::hypot(a, 0.25);       // Undersampled kernels are essentially delta-functions
            afwMath::GaussianFunction1<double> gauss(a);
            int const kSize = 2*static_cast<int>(3*a + 1) + 1;
            afwMath::SeparableKernel const kernel(kSize, kSize, gauss, gauss);

            afwMath::convolve(smoothed, *psfImage, kernel, normalize, copyEdge);
            double const rms = calcRms(smoothed);
            if (rms >= psfRmsCorrected) {
                break;
            }

            *psfImage = smoothed;

            a = ::sqrt(::pow(psfRmsCorrected, 2) - ::pow(rms, 2));
            
            if (::fabs(a) < ctrl.rmsTol) { // good enough
                break;
            }
        }
    } else if (epsilon2 == 0.0) {
        ;
    } else {                            // we need to deconvolve the PSF model to match the object
        /*
         * We're going to deconvolve by using a linear expansion of the kernel 
         *
         * We start by convolving the PSF with a Gaussian, and expanding in a taylor series
         *           PSF \otimes N(delta^2) = PSF + delta*d(PSF)/d(width) + ...
         * so
         *           d(PSF)/d(width) ~ (PSF \otimes N(delta^2) - PSF)/delta
         *
         * we can then get a deconvolved PSF by evaluating
         *           PSF - |epsilon|*d(PSF)/d(width)
         *
         * Choosing delta's a black art;  we have a pixellated image so it can't be too small.
         *
         * In practice it works better to repeatedly subtract multiples of d(PSF)/d(width)
         * until our image has the desired RMS
         */
        double const psfRmsCorrected = ::sqrt(psfRms*psfRms + epsilon2); // n.b. epsilon2 < 0

        double delta = ctrl.deconvolutionKernelSigma; // smoothing length used to estimate deconvolution
        afwMath::GaussianFunction1<double> gauss(delta);
        int const kSize = 2*static_cast<int>(3*delta + 1) + 1;
        afwMath::SeparableKernel const kernel(kSize, kSize, gauss, gauss);
        {
            afwDetection::Psf::Image kim(kSize, kSize);
            kernel.computeImage(kim, true);
            delta = calcRms(kim);       // the width allowing for pixellation
        }

        afwMath::convolve(smoothed, *psfImage, kernel, normalize, copyEdge); // PSF \otimes N
        
        smoothed -= *psfImage;          // delta*d(psfImage)/d(width)
        smoothed /= delta;              // d(psfImage)/d(width)
        //
        // Iterate until we get the desired PSF width
        //
        {
            double a = ::sqrt(-epsilon2); // we know that we need to deconvolve, so take the absolute value
            for (int i = 0; i < ctrl.niter; ++i) {
                smoothed *= a;        //  ~ a*d(psfImage)/d(width)
                *psfImage -= smoothed;
                smoothed /= a;        //  ~ d(psfImage)/d(width)
                double const rms = calcRms(*psfImage);

                a = rms/psfRmsCorrected - 1;                
                
                if (::fabs(a) < ctrl.rmsTol) { // good enough
                    break;
                }
            }
        }
    }
#if 0
    psfImage = psf->computeImage(center); // debugging XXXXXXXXXXXXXXXXXXXXX
#endif
    //
    // OK, we have the proper PSF so we can proceed to make the measurement
    //
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
