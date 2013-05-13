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
#include "lsst/afw/math/Statistics.h"
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

PTR(afwMath::SeparableKernel)
makeGaussianKernel0(double alpha)
{
    afwMath::GaussianFunction1<double> gauss(alpha);
    int const kSize = 2*static_cast<int>(3*alpha + 1) + 1;
    return boost::make_shared<afwMath::SeparableKernel>(kSize, kSize, gauss, gauss);
}

PTR(afwMath::SeparableKernel)
makeGaussianKernel(double alpha, int const niter=10, double const tol=1e-4)
{
    if (niter <= 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                          str(boost::format("niter must be >= 1; saw %d") % niter));
    }

    PTR(afwMath::SeparableKernel) kernel; // the desired kernel

    double rms0 = alpha;                                        // desired rms
    double epsilon_min = alpha, epsilon_max = ::hypot(alpha, 1); // range of possible values
    double rms_min, rms_max;            // rms of kernel with epsilon = epsilon_{min,max}
    rms_min = rms_max = std::numeric_limits<double>::quiet_NaN(); // unknown
    double rms;                                                   // measured rms

    for (int i = 0; i != niter; ++i) {
        kernel = makeGaussianKernel0(alpha);
        {
            afwDetection::Psf::Image kim(kernel->getDimensions());
            kernel->computeImage(kim, true);

            rms = calcRms(kim);
        }
        //printf("%2d %.4f <= %.4f <= %.4f rms: %.4f\n", i, epsilon_min, alpha, epsilon_max, rms);

        if (::fabs(rms - rms0) < tol) {
            break;
        } else if (rms > rms0) {
            if (alpha <= epsilon_min) {
                break;
            }
            if (!utils::isfinite(rms_max)) {
                double tmp = alpha;
                alpha = 0.5*(alpha + epsilon_min);
                epsilon_max = tmp;
            } else {
                double nalpha = epsilon_max - (rms_max - rms0)*(epsilon_max - alpha)/(rms_max - rms);
                if (nalpha < epsilon_min) {
                    nalpha = 0.5*(alpha + epsilon_min);
                }
                epsilon_max = alpha;
                alpha = nalpha;
            }
            rms_max = rms;
        } else {
            if (alpha == epsilon_max) {
                break;
            }
            if (!utils::isfinite(rms_min)) {
                double const tmp = alpha;
                alpha = 0.5*(alpha + epsilon_max);
                epsilon_min = tmp;
            } else {
                double nalpha = epsilon_min + (rms0 - rms_min)*(alpha - epsilon_min)/(rms - rms_min);
                if (nalpha > epsilon_max) {
                    nalpha = 0.5*(alpha + epsilon_max);
                }
                epsilon_min = alpha;
                alpha = nalpha;
            }
            rms_min = rms;
        }
    }
    //std::cout << "alpha = " << alpha << std::endl;

    return kernel;
}

/************************************************************************************************************/
/*
 * Lucy-Richardson deconvolution of img with an N(0, epsilon^2) kernel
 */
template<typename ImageT>
PTR(ImageT)
deconvolveLR(ImageT const& phi_,        // image to deconvolve
             double epsilon,            // with an N(0, epsilon^2) kernel
             int niter=10)              // number of L-R iterations
{
    ImageT phi = phi_;                  // we need a writeable copy as we need to add a bias to make it > 0
    afw::math::Statistics const stats =
        afw::math::makeStatistics(phi, afw::math::MIN | afw::math::MAX | afw::math::MEAN);
    double const min = ::fabs(stats.getValue(afw::math::MIN));
    double const max = ::fabs(stats.getValue(afw::math::MAX));
    double const bkgd = min + 1e-3*(max - min);
    phi += bkgd;

    PTR(afw::math::SeparableKernel const) kernel = makeGaussianKernel(epsilon); // N(0, epsilon^2)
    /*
     * The input data's phi; the deconvolved data's psi
     */
    PTR(ImageT) psi_r = boost::make_shared<ImageT>(phi.getDimensions()); // previous estimate of psi
    PTR(ImageT) psi_rp1 = boost::make_shared<ImageT>(phi.getDimensions()); // current estimate of psi

    ImageT phi_r(phi.getDimensions());  // current estimate of observed data (== psi otimes kernel)
    ImageT tmp(phi.getDimensions());                                       // temp space
    
    *psi_r = stats.getValue(afw::math::MEAN); // initial estimate is mean of input data
    for (int i = 0; i != niter; ++i) {
        bool const normalize = true;
        bool const copyEdge = true;
        afwMath::convolve(phi_r, *psi_r, *kernel, normalize, copyEdge);

        tmp = phi;
        tmp /= phi_r;
        afwMath::convolve(*psi_rp1, tmp, *kernel, normalize, copyEdge);
        *psi_rp1 *= *psi_r;

#if 0
        psi_r.swap(psi_rp1);            // update
#else
        *psi_r = *psi_rp1;              // update
#endif
    }
    *psi_r -= bkgd;

    return psi_r;
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
    if (source.getPsfFluxFlag()) {
        throw LSST_EXCEPT(lsst::pex::exceptions::NotFoundException,
                          str(boost::format("Aperture flux is invalid for object at (%.3f, %.3f)")
                              % center.getX() % center.getY()));
    }

    double const objFlux = source.getPsfFlux(); // object's flux
    double const flux0 = ctrl.flux0;            // fiducial flux (see coeff)
    double const coeff = ctrl.coeff;            // intensity kernel is N(0, coeff*log10(flux/ctrl.flux0)^2)
    double const psfFlux = ctrl.psfFlux;        // characteristic flux of the PSF
    double const epsilon2 =
        ::pow(coeff*::log10(objFlux/flux0), 2) - 
        ::pow(coeff*::log10(psfFlux/flux0), 2); // (correction we need to apply)^2
 
    afwDetection::Psf::Image smoothed(psfImage->getDimensions());
    if (epsilon2 > 0) {                  // we need to smooth our PSF model to match the object
        PTR(afwMath::SeparableKernel const) kernel = makeGaussianKernel(::sqrt(epsilon2));
        bool const normalize = true;
        bool const copyEdge = true;
        afwMath::convolve(smoothed, *psfImage, *kernel, normalize, copyEdge);
        *psfImage = smoothed;
    } else if (epsilon2 == 0.0) {
        ;
    } else {                            // we need to deconvolve the PSF model to match the object
        psfImage = deconvolveLR(*psfImage, ::sqrt(-epsilon2), ctrl.niter);
    }

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
