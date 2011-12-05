/**
***************************************************************************
* @file brick/computerVision/utilities.hh
*
* Header file declaring utility functions for brick::computerVision
* library.
*
* Copyright (C) 2006, 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_UTILITIES_HH
#define BRICK_COMPUTERVISION_UTILITIES_HH

#include <brick/computerVision/colorspaceConverter.hh>
#include <brick/computerVision/image.hh>
#include <brick/numeric/transform2D.hh>

namespace brick {

  namespace computerVision {
    
    /**
     * This function returns by reference an image which either shares
     * or copies the data from the input array.
     *
     * If possible, this function returns an Image which references the
     * same memory as the input array, but in which each pixel is the
     * aggregate of the appropriate number of elements from the array.
     * For example, if this function is called with template argument
     * RGB8, and an Nx(3*M) array of UnsignedInt8 is passed as the
     * argument, then the return value will be an NxM Image<RGB8> which
     * references the same memory as the argument.  Imagine that the
     * first three elements of the first row of the argument are 12, 14,
     * and 72.  In this case, the upper-left RGB value of the returned
     * image image is {12, 14, 72}.  This memory sharing only works if
     * the compiler does not "pad" the pixel struct to byte align its
     * members.
     *
     * If memory sharing is not possible, then this function tests the
     * size of argument outputImage, reinitializes it only if the size
     * does not match the size of the input array, and then copies the
     * data from inputArray to outputImage.
     *
     * WARNING: If data sharing is possible, the returned image does no
     * memory management.  It is only valid until the original input
     * data is destroyed.
     *
     * If you use this function in conjunction with
     * dissociateColorComponents(), you can simply ignore the question
     * of whether the data sharing works or not.  For example, you might
     * use associateColorComponents to get an Image to operate on, do
     * all of your image processing, and then use
     * dissocateColorComponents to get back to a flat array.  If data
     * sharing is possible, all of your operations will have taken place
     * in the original array.  If data sharing isn't possible, then the
     * two functions will take care of copying back and forth between
     * the image and the flat array.
     * 
     * @param inputArray This argument is the array from which to
     * construct the return image.
     * 
     * @param outputImage This argument is the image which will be
     * modified or copied to.
     *
     * @return The return value is true if the returned image shares
     * data with the input array, false otherwise.
     */
    template<ImageFormat FORMAT>
    bool
    associateColorComponents(
      brick::numeric::Array2D<typename ImageFormatTraits<FORMAT>::ComponentType>& inputArray,
      Image<FORMAT>& outputImage);


    /**
     * @deprecated {Please use the two-argument version of
     * associateColorComponents() instead.}
     *
     * This function tries to return an Image which references the
     * same memory as the input array, but in which each pixel is the
     * aggregate of the appropriate number of elements from the array.
     * If it is not possible to do so, then this function will throw a
     * LogicException.  This function is deprecated.  Instead, ues the
     * two-argument form of associateColorComponents().
     * 
     * WARNING: The returned image does no memory management.  It is only
     * valid until the original input image is destroyed.
     * 
     * @param inputArray This argument is the array from which to
     * construct the return image.
     * 
     * @return The return value is an image pointing to the input array
     * data, but in which the color components are lumped into pixels.
     */
    template<ImageFormat FORMAT>
    Image<FORMAT>
    associateColorComponents(
      brick::numeric::Array2D<typename ImageFormatTraits<FORMAT>::ComponentType>& inputArray);


    /** 
     * This function takes an image in one colorspace and generates a
     * corresponding image in a second colorspace.  Use this function as
     * follows:
     *
     *   Image<GRAY8> = convertColorspace<GRAY8, RGB16>(myRGB16Image);
     *
     * or equivalently, the second template argument can be left implicit,
     *
     *   Image<GRAY8> = convertColorspace<GRAY8>(myRGB16Image);
     *
     * This function only works for image formats in which there's a
     * one-to-one mapping between input and output pixel types.  When
     * converting between formats which don't match this requirement,
     * such as when converting from RGB8 to YUV420, please use a
     * different routine.
     *
     * @param inputImage This argument is the image to be converted.
     * 
     * @return The return value is an image in the converted colorspace.
     */
    template<ImageFormat OUTPUT_FORMAT, ImageFormat INPUT_FORMAT>
    Image<OUTPUT_FORMAT>
    convertColorspace(const Image<INPUT_FORMAT>& inputImage);


    /**
     * This function returns by reference an array which either shares
     * or copies the data from the input image.
     *
     * If possible, this function returns an Array2D which references
     * the same memory as the input image, but in which each pixel has
     * been "flattened" so that the returned array has a separate
     * element for each color component of each pixel.  For example, if
     * this function is called with an NxM Image<RGB8> as its argument,
     * the returned value will be an Nx(3*M) array of 8-bit unsigned
     * ints which references the same memorty.  Imagine that the
     * upper-left RGB value of the image is {12, 14, 72}.  In this case,
     * the first three elements of first row of the returned image will
     * be 12, 14, and 72.  This memory sharing only works if the
     * compiler does not "pad" the pixel struct to byte align its
     * members.
     *
     * If memory sharing is not possible, then this function tests the
     * size of argument outputArray, reinitializes it only if the size
     * does not match the size of the input image, and then copies the
     * data from inputImage to outputArray.
     *
     * WARNING: If data sharing is possible, the returned array does no
     * memory management.  It is only valid until the original input
     * data is destroyed.
     *
     * If you use this function in conjunction with
     * associateColorComponents(), you can simply ignore the question
     * of whether the data sharing works or not.  For example, you might
     * use associateColorComponents to get an Image to operate on, do
     * all of your image processing, and then use
     * dissocateColorComponents to get back to a flat array.  If data
     * sharing is possible, all of your operations will have taken place
     * in the original array.  If data sharing isn't possible, then the
     * two functions will take care of copying back and forth between
     * the image and the flat array.
     * 
     * @param inputImage This argument is the image from which to
     * construct the return array.
     * 
     * @outputArray This argument is the array which will be modified or
     * copied to.
     *
     * @return The return value is true if the returned array shares
     * data with the input image, false otherwise.
     */
    template<ImageFormat FORMAT>
    bool
    dissociateColorComponents(
      Image<FORMAT>& inputImage,
      brick::numeric::Array2D<typename ImageFormatTraits<FORMAT>::ComponentType>& outputArray);


    /** 
     * @deprecated {Please use the two-argument version of
     * dissociateColorComponents() instead.}
     *
     * This function tries to return an Array2D which references the
     * same memory as the input image, but in which each pixel has
     * been "flattened" so that the returned array has a separate
     * element for each color component of each pixel.  If it is not
     * possible to do so, then this function will throw a
     * LogicException.  This function is deprecated.  Instead, ues the
     * two-argument form of dissociateColorComponents().
     * 
     * WARNING: The returned array does no memory management.  It is only
     * valid until the original input image is destroyed.
     * 
     * @param inputImage This argument is the image from which to
     * construct the return array.
     * 
     * @return The return value is a 2D array pointing to the input
     * image data, but in which the color components are not lumped into
     * pixels.
     */
    template<ImageFormat FORMAT>
    brick::numeric::Array2D<typename ImageFormatTraits<FORMAT>::ComponentType>
    dissociateColorComponents(Image<FORMAT>& inputImage);


    /** 
     * Return the transform that, in the least-squares sense, most
     * nearly transforms each point in fromPoints to the corresponding
     * point in toPoints.  All of the iterator arguments must come
     * from sequences of numeric::Vector2D<Type>.  This function
     * currently only works with Type set to Float64 (double).
     * 
     * @param toPointsBegin This argument and the next define a
     * sequence of Vector2D<Type> to which fromPoints correspond.
     * 
     * @param toPointsEnd This argument and the previous define a
     * sequence of Vector2D<Type> to which fromPoints correspond.
     * 
     * @param fromPointsBegin This argument begins a sequence of
     * Vector2D<Type> that we imagine can be mapped onto toPoints by
     * an affine transform.  This sequence must have at least the same
     * number of elements as toPoints.
     * 
     * @return The return value is the affine transform that most
     * nearly transforms each element of fromPoints to the
     * corresponding element of toPoints.
     */
    template <class Type, class Iter0, class Iter1>
    brick::numeric::Transform2D<Type>
    estimateAffineTransform(Iter0 toPointsBegin, Iter0 toPointsEnd,
                            Iter1 fromPointsBegin);

    
    /**
     * This function subsamples its input to create a new, smaller
     * image.  The number of rows in the resulting image is the
     * largest integer, r, such that (rowStep * (r - 1) + 1) <=
     * inputImage.rows().  The number of columns in the resulting
     * image is the largest integer, c, such that (columnStep * (c -
     * 1) + 1) <= inputImage.columns().  The pixel values in the
     * resulting image are the values at the intersection of every
     * rowStep-th row and columnStep-th column of inputImage, starting
     * with the value at (row, column) = (0, 0), which is always
     * copied directly in to pixel (0, 0) of the result image.  In
     * other words, element (0, 0) of the resulting image is equal to
     * element (0, 0) of inputImage, element (0, 1) of the resulting
     * image is equal to element (0, columnStep) of inputImage,
     * element (0, 2) of the resulting image is equal to element (0, 2
     * * columnStep) of inputImage, element (2, 4) of the resulting
     * image is equal to element (2 * rowStep, 4 * columnStep) in
     * inputImage.
     *
     * @param inputImage This argument is the image to be subsampled.
     *
     * @param rowStep This argument controls which rows of inputImage
     * contribute to the result.
     *
     * @param columnStep This argument controls which columns of inputImage
     * contribute to the result.
     *
     * @return The return value is the subsampled image.
     */
    template<ImageFormat Format>
    Image<Format>
    subsample(const Image<Format>& inputImage,
              size_t rowStep = 2,
              size_t columnStep = 2);


    /**
     * This function interpolates its input to create a new, larger
     * image.  The number of rows in the resulting image is one fewer
     * than twice the number of rows in the input image.  The number
     * of columns in the resulting image is one fewer than twice the
     * number of columns in the input image.  Where output pixels
     * overlap input pixels, they are copied directly from the input
     * image.  Where output pixels fall in between input pixels, they
     * are constructed by bilinear interpolation.  The pixel at (0, 0)
     * in the output image falls directly on top of pixel (0, 0) in
     * the input image.  (2, 2) in the output image falls on top of
     * (1, 1) in the input image, and so forth.
     *
     * Template argument IntermediateFormat specifies the type of
     * pixel that should be used during calculation of the
     * interpolated values.  It should have sufficient precision to
     * contain the sum of four input pixels.
     *
     * @param inputImage This argument is the image to be interpolated.
     *
     * @return The return value is the supersampled image.
     */
    template<ImageFormat InputFormat,
             ImageFormat OutputFormat,
             ImageFormat IntermediateFormat>
    Image<OutputFormat>
    supersample(const Image<InputFormat>& inputImage);
  
    
    /** 
     * This function creates a new array and copies into it the pixel
     * data from the input image.  If the image format is one which has
     * multiple interleaved color components, such as RGB8, then the
     * returned array will have individual elements for each component
     * of each pixel.  For example, if the first row of an HSV8 image is
     * [{h0, s0, v0}, {h1, s1, v1}, {h2, s2, v2}, ... ], then the first
     * row of the array returned by toArray will be [h0, s0, v0, h1, s1,
     * v1, h2, s2, v2, ...].
     * 
     * @param inputImage This argument is the image to be copied.
     * 
     * @return The return value is an array containing the copied data.
     */
    template <class Type, ImageFormat FORMAT>
    brick::numeric::Array2D<Type>
    toArray(const Image<FORMAT>& inputImage);

  } // namespace computerVision    

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/utilities_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_UTILITIES_HH */
