/**
***************************************************************************
* @file brick/iso12233/demos/analyzePatch.cc
*
* Demo program to compute the MTF of an image patch.
*
* Copyright (C) 2020 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/exception.hh>

#include <brick/computerVision/image.hh>
#include <brick/computerVision/imageIO.hh>

#include <brick/iso12233/oecf.hh>
#include <brick/iso12233/iso12233.hh>

#include <brick/numeric/array1D.hh>
#include <brick/numeric/subArray2D.hh>

#include <brick/utilities/optionParser.hh>

namespace bcm = brick::common;
namespace bcv = brick::computerVision;
namespace bis = brick::iso12233;
namespace bnm = brick::numeric;

namespace {

  struct AppConfig
  {
    std::string inputFileName;
    std::string outputImageFileName;
    std::size_t topLeftRow = 0;
    std::size_t topLeftColumn = 0;
    std::size_t patchWidth = 64;
    std::size_t patchHeight = 64;
    std::size_t windowSize = 0;
  };

  double estimateMtf50(bnm::Array1D<double> const& mtf);
  void parseArguments(AppConfig& appConfig, int argc, char **argv);
    
} // namespace


int main(int argc, char* argv[])
{
  // Record any user input and configuration.
  AppConfig appConfig;
  parseArguments(appConfig, argc, argv);

  // Get file input.
  std::string comment;

#if HAVE_LIBPNG
  bcv::Image<bcv::GRAY8> inputImage = bcv::readPNG<bcv::GRAY8>(
    appConfig.inputFileName, comment);
#else
  bcv::Image<bcv::GRAY8> inputImage = bcv::readPGM8(
    appConfig.inputFileName, comment);
#endif

  // Condition input arguments based on image size.
  appConfig.topLeftRow = std::min(appConfig.topLeftRow,
                                  inputImage.rows() - 2);
  appConfig.topLeftColumn = std::min(appConfig.topLeftColumn,
                                  inputImage.columns() - 2);
  if(appConfig.topLeftRow + appConfig.patchHeight > inputImage.rows() ||
     appConfig.topLeftColumn + appConfig.patchHeight > inputImage.columns()) {
    BRICK_THROW(bcm::IndexException, "main()",
                "Input patch extends past image edge.");
  }

  // Snip out image patch.
  std::size_t bottomLeftColumn = appConfig.topLeftColumn + appConfig.patchWidth;
  std::size_t bottomLeftRow = appConfig.topLeftRow + appConfig.patchHeight;
  bnm::Array2D<bcv::Image<bcv::GRAY8>::PixelType> patchArray = 
    bnm::subArray(inputImage, bnm::Slice(appConfig.topLeftRow, bottomLeftRow),
                  bnm::Slice(appConfig.topLeftColumn, bottomLeftColumn));
  bcv::Image<bcv::GRAY8> inputPatch(patchArray);
  if(appConfig.outputImageFileName != "") {
#if HAVE_LIBPNG
    bcv::writePNG(appConfig.outputImageFileName, inputPatch, "");
#else
    bcv::writePGM8(appConfig.outputImageFileName, inputPatch, "");
#endif
  }
  
  // Try to process the image and output the result.
  bnm::Array1D<double> mtf = bis::iso12233<double>(
    inputPatch, appConfig.windowSize, bis::ShiftOECF<double,
    bcv::GRAY8>(inputPatch));

  // Compute some statistics.
  auto mtf50 = estimateMtf50(mtf);

  std::cout << "MTF: " << mtf << std::endl;
  std::cout << "MTF50: " << mtf50 << " lp/pix" << std::endl;

  return 0;
}


namespace {

  double estimateMtf50(bnm::Array1D<double> const& mtf)
  {
    // Find the first element less than 0.5.
    std::size_t ii = 0;
    while(ii < mtf.size()) {
      if(mtf[ii] < 0.5) {
        break;
      }
      ++ii;
    }

    // Sanity check.
    if(ii == 0 || ii >= mtf.size()) {
      BRICK_THROW(bcm::ValueException, "estimateMtf50()",
                  "Input MTF doesn't drop below 0.5 in the expected range.");
    }

    // Interpolate to find the index at which we imagine the MTF drops
    // below 0.5.  Using linear interpolation, we have:
    //   alpha * mtf[ii] + (1 - alpha) * mtf[ii - 1] = 0.5
    //   alpha * (mtf[ii] - mtf[ii - 1]) = 0.5 - mtf[ii - 1]
    //   alpha = (0.5 - mtf[ii - 1]) / (mtf[ii] - mtf[ii - 1])
    auto alpha = (0.5 - mtf[ii - 1]) / (mtf[ii] - mtf[ii - 1]);

    // Sanity check.
    if(alpha < 0.0 || alpha > 1.0) {
      BRICK_THROW(bcm::LogicException, "estimateMtf50()",
                  "Interpolation error.");
    }

    // Now convert to line-pairs per pixel.
    double mtf50Index = (ii - 1) + alpha;
    return mtf50Index / mtf.size();
  }

  void parseArguments(AppConfig& appConfig, int argc, char **argv)
  {
    // Create an OptionParser that automatically generates help
    // messages, and that prints usage and calls exit(65) if it
    // encounters an inappropriate command line.
    brick::utilities::OptionParser optionParser(65);

    // Specify required positional arguments.
    optionParser.addPositionalArgument(
      "PNG_FILE",
      "Image containing input patch.  This must be an 8-bit grayscale PNG "
      "image from which we can snip out a patch that contains a single "
      "\"slanted edge.\"  The patch should be significantly wider (say, more "
      "than 2x) than argument WINDOW_SIZE.  The patch should have plenty "
      "(say, more than 10, fewer than 1000) rows.", true);

    // Specify positional arguments that are not required.
    {
      std::ostringstream defaultValueStream;
      defaultValueStream << appConfig.topLeftColumn;
      optionParser.addPositionalArgument(
        "TOP_LEFT_COLUMN",
        "Column coordinate of the upper left corner of the image patch to be "
        "analyzed.", false, defaultValueStream.str());
    }

    {
      std::ostringstream defaultValueStream;
      defaultValueStream << appConfig.topLeftRow;
      optionParser.addPositionalArgument(
        "TOP_LEFT_ROW",
        "Row coordinate of the upper left corner of the image patch to be "
        "analyzed.", false, defaultValueStream.str());
    }
    
    // Specify options which do not take arguments.
    // optionParser.addOption(
    // "TEXT_TARGET", "-t", "--text-target",
    // "Setting this option generates old-style datasets in which the target "
    // "values are written to a text file, not an image.");

    // Specify options that require values.
    optionParser.addOptionWithValue(
      "OUTPUT_PNG", "-o", "--output-image", appConfig.outputImageFileName,
      "If this option is not set to the empty string, the selected image "
      "patch will be saved to the specified PNG output file.");

    optionParser.addOptionWithValue(
      "PATCH_HEIGHT", "-h", "--patch-height", appConfig.patchHeight,
      "This argument should generally be between 10 and 1000.  It, along with "
      "PATCH_WIDTH, controlsmust be signficantly larger (like 2x or so) than "
      "WINDOW_SIZE.  It, along with PATCH_HEIGHT, controls the size of the "
      "image patch that will be extracted and passed to the e-SFR code.  Try "
      "to make this patch include the relevant slanted edge only.");

    optionParser.addOptionWithValue(
      "PATCH_WIDTH", "-w", "--patch-width", appConfig.patchWidth,
      "This argument must be signficantly larger (like 2x or so) than "
      "WINDOW_SIZE.  It, along with PATCH_HEIGHT, controls the size of the "
      "image patch that will be extracted and passed to the e-SFR code.  Try "
      "to make this patch include the relevant slanted edge only.");

    optionParser.addOptionWithValue(
      "WINDOW_SIZE", "-s", "--window-size", appConfig.windowSize,
      "This argument must (for now) be a power of two, and is the size of the "
      "horizontal window surrounding the slanted edge that should be "
      "evaluated.  Choose a number high enough to be sure that any lingering "
      "effects of being near an edge have died out by the time you are "
      "windowSize / 2 pixels from the edge itself.  The larger this number, "
      "the less your result will be affected by FFT artifacts, so shoot for "
      "10 times the width of your blurred edge.  The smaller this number, "
      "the less likely the window will overlap the edge of the selected "
      "patch, so don't make it unnecessarily large.  If set to zero here, "
      "WINDOW_SIZE will be automatically adjusted to the largest power of "
      "two that is smaller than (PATCH_WIDTH / 2).  If set to a number that "
      "is not a power of two, WINDOW_SIZE will be rounded down to the nearest "
      "power of two.");

    // Parse program arguments.
    optionParser.parseCommandLine(argc, argv);

    // Recover results of parse.
    appConfig.inputFileName = optionParser.getValue("PNG_FILE");
    appConfig.outputImageFileName = optionParser.getValue("OUTPUT_PNG");
    appConfig.topLeftColumn = optionParser.convertValue<std::size_t>(
      "TOP_LEFT_COLUMN");
    appConfig.topLeftRow = optionParser.convertValue<std::size_t>(
      "TOP_LEFT_ROW");
    appConfig.patchHeight = optionParser.convertValue<std::size_t>(
      "PATCH_HEIGHT");
    appConfig.patchWidth = optionParser.convertValue<std::size_t>(
      "PATCH_WIDTH");
    appConfig.windowSize = optionParser.convertValue<std::size_t>(
      "WINDOW_SIZE");


    // If windowSize wasn't set, pick a reasonable value here.
    if(appConfig.windowSize == 0) {
      appConfig.windowSize = appConfig.patchWidth / 2;
    }

    // Round windowSize down to the nearest power of two.
    std::size_t powerOfTwo = 2;
    while(powerOfTwo <= appConfig.windowSize) {
      powerOfTwo <<= 1;
    }
    appConfig.windowSize = (powerOfTwo >> 1);
  }

} // namespace
