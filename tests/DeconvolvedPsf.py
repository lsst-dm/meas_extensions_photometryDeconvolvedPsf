#!/usr/bin/env python
"""
Tests for measuring things

Run with:
   python DeconvolvedPsf.py
or
   python
   >>> import DeconvolvedPsf; DeconvolvedPsf.run()
"""

import math, os, sys, unittest
import numpy as np
import lsst.utils.tests as tests
import lsst.pex.exceptions as pexExceptions
import lsst.pex.logging as pexLogging
import lsst.pex.policy as pexPolicy
import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as afwEllipses
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg
import lsst.meas.extensions.photometryDeconvolvedPsf as DeconvolvedPsf

try:
    type(verbose)
except NameError:
    verbose = 1
    display = False
    ds9Frame = 0
pexLogging.Trace_setVerbosity("meas.photometry.deconvolvedPsf", verbose)

import lsst.afw.display.ds9 as ds9
import lsst.afw.display.utils as displayUtils

def calcRms(img):
    xcen = img.getWidth()//2
    ycen = img.getHeight()//2

    sumrr = 0
    for y in range(img.getHeight()):
        for x in range(img.getWidth()):
            sumrr += float(img[x, y])*((x - xcen)**2 + (y - ycen)**2)
    sumrr /= img.getArray().sum()

    return math.sqrt(sumrr/2)


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class DeconvolvedPsfPhotometryTestCase(unittest.TestCase):
    """A test case for measuring DeconvolvedPsf quantities"""

    def setUp(self):
        self.width, self.height = 70, 70
        self.objImg = None

        self.psfFlux = 1e5              # Flux corresponding to the PSF
        self.flux0 = 1e3                # Characteristic flux (see coeff)
        self.coeff = 0.1                # Smooth with N(0, self.coeff*log10(flux/flux0)^2)
        
    def tearDown(self):
        if self.objImg:
            del self.objImg

    def makePsf(self, alpha, b, objFlux):
        ksize = 25                      # size of desired kernel

        epsilon = self.coeff*math.log10(objFlux/self.flux0) if self.flux0 != 0.0 else 0.0
        psf = measAlg.DoubleGaussianPsf(ksize, ksize,
                                        math.hypot(alpha,   epsilon),
                                        math.hypot(2*alpha, epsilon), b)

        if False:
            kim = afwImage.ImageD(ksize, ksize)
            kim = psf.computeImage()
            print "alpha =", alpha, "eps =", epsilon, "rms =", calcRms(kim)

        return psf

    def makeAndMeasure(self, objFlux, alpha, b, dx=0.0, dy=0.0):
        """Make and measure a PSF"""

        xcen, ycen = 0.5*self.width + 11 + dx, 0.5*self.height + 12 + dy
        #
        # Create the PSF
        #
        psf = self.makePsf(alpha, b, self.psfFlux)
        #
        # Make the object
        #
        self.objImg = None
        if not self.objImg:
            gal = afwImage.ImageF(self.width, self.height)
            gal.setXY0(10, 10)

            obj = self.makePsf(alpha, b, objFlux).computeImage(afwGeom.PointD(xcen, ycen))
            obj *= objFlux/obj.getArray().sum()

            if False:               # requires support for gal[obj.getBBox(), afwImage.PARENT]
                gal[obj.getBBox(afwImage.PARENT), afwImage.PARENT] = obj.convertF()
            else:
                gal.Factory(gal, obj.getBBox(afwImage.PARENT), afwImage.PARENT)[:] <<= obj.convertF()

            self.objImg = afwImage.makeExposure(afwImage.makeMaskedImage(gal))
            self.objImg.setPsf(psf)

            self.objImg.getMaskedImage().getVariance()[:] = 1.0

            if display:
                ds9.mtv(self.objImg, frame=ds9Frame, title="%g %g" % (alpha, b))

                ds9.dot("+", xcen - self.objImg.getX0(), ycen - self.objImg.getY0(),
                        size=1, ctype=ds9.RED, frame=ds9Frame)
                ds9.pan(xcen - self.objImg.getX0(), ycen - self.objImg.getY0(), frame=ds9Frame)
        #
        # Do the measuring
        #
        msConfig = measAlg.SourceMeasurementConfig()
        msConfig.algorithms.names.add("flux.sinc")
        msConfig.algorithms.names.add("flux.psf")
        msConfig.algorithms.names.add("flux.deconvolvedPsf")
        msConfig.algorithms.names.remove("correctfluxes")
        msConfig.slots.apFlux = "flux.sinc"

        msConfig.algorithms["flux.deconvolvedPsf"].priority = 4 # i.e. run after other flux algorithms
        #msConfig.algorithms["flux.deconvolvedPsf"].deconvolutionKernelSigma = 0.4
        msConfig.algorithms["flux.deconvolvedPsf"].coeff = self.coeff
        msConfig.algorithms["flux.deconvolvedPsf"].psfFlux = self.psfFlux
        msConfig.algorithms["flux.deconvolvedPsf"].flux0 = self.flux0
        #msConfig.algorithms["flux.deconvolvedPsf"].niter = 15
        #msConfig.algorithms["flux.deconvolvedPsf"].rmsTol = 1e-4
        
        schema = afwTable.SourceTable.makeMinimalSchema()
        ms = msConfig.makeMeasureSources(schema) # add our fields
        
        table = afwTable.SourceTable.make(schema)
        msConfig.slots.setupTable(table)
        source = table.makeRecord()

        ss = afwDetection.FootprintSet(self.objImg.getMaskedImage(), afwDetection.Threshold(0.1))
        feet = ss.getFootprints()
        assert(len(feet) > 0)
        fp = ss.getFootprints()[0]
        source.setFootprint(fp)

        center = afwGeom.Point2D(xcen, ycen)
        ms.apply(source, self.objImg, center)

        flux = source.get("flux.deconvolvedPsf")
        fluxErr = source.get("flux.deconvolvedPsf.err")
        flags = source.get("flux.deconvolvedPsf.flags")

        if display:
            xc, yc = xcen - self.objImg.getX0(), ycen - self.objImg.getY0()
            ds9.dot("x", xc, yc, ctype=ds9.MAGENTA, size=1, frame=ds9Frame)
            displayUtils.drawFootprint(fp, XY0=self.objImg.getXY0())

            shape = source.getShape()

        return flux, fluxErr, flags, source.get("flux.psf")

    def testEllipticalGaussian(self):
        """Test measuring the DeconvolvedPsf quantities of an elliptical Gaussian"""
        #
        # Make and the objects
        #
        ab_vals = [(1.0, 0.0), (2.0, 0.0), (2.0, 0.1), ]
        nIter = 2
        for dx in (0.0, 0.5,):
            for dy in (0.0, 0.5,):
                for alpha, b in ab_vals:
                    for objFlux in (1e2, 1e3, 1e4, 1e5, 1e6, ):
                        flux, fluxErr, flags, psfFlux = self.makeAndMeasure(objFlux, alpha, b, dx=dx, dy=dy)
                        
                        failFlux =  math.isnan(flux) or flags or abs(flux/objFlux - 1) > 2.5e-3
                        
                        ID = "alpha,b %4.1f, %5.2f  dx,dy = %.1f,%.1f " % (alpha, b, dx, dy)
                        msg = "%s  flux_DeconvolvedPsf: %9.4g v. exact value %9.4g (error %5.2f%%) (psfFlux error %5.2f%%)" % \
                            (ID, flux, objFlux, 100*(flux/objFlux - 1), 100*(psfFlux/objFlux - 1))

                        if False:
                            print msg
                            #continue

                        self.assertFalse(failFlux, msg)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(DeconvolvedPsfPhotometryTestCase)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
