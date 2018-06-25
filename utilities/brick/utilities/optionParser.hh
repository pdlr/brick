/**
***************************************************************************
* @file brick/utilities/optionParser.hh
*
* Header file declaring the OptionParser class.
*
* Copyright (C) 2006-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_UTILITIES_OPTIONPARSER_HH
#define BRICK_UTILITIES_OPTIONPARSER_HH

#include <map>
#include <string>
#include <vector>
#include <brick/utilities/optionDescription.hh>
#include <brick/utilities/stringManipulation.hh>

namespace brick {

  namespace utilities {

    /**
     ** The OptionParser class parses program options.  It supports
     ** the following features.
     **
     **   - Short style options: "du -s -k"
     **   - Long gnu-style options: "gnumeric --version --no-splash"
     **   - Options with values: "xclock -g 200x200 -rgb5"
     **   - Named positional arguments: "ls foo.txt"
     **   - Repeated options: "gcc -I/usr/X11/include -I. -I../include"
     **   - Default values for unspecified options.
     **   - Automatic option documentation for usage messages.
     **   - Short option stacking (if desired): "ls -ltr"
     **   - Automatic program exit and printing of usage messages (if desired).
     **
     ** Here is example code:
     **
     ** @code
     **   // Create an OptionParser which automatically generates
     **   // help messages, and which prints usage and calls exit(65)
     **   // if it encounters an inappropriate command line.
     **   OptionParser optionParser(65);
     **
     **   // Specify options which do not take arguments.
     **   optionParser.addOption(
     **     "NOLINK", "-c", "--no_link",
     **     "Don't run the linker.  Instead, generate .o files.");
     **   optionParser.addOption(
     **     "DEBUG", "-g", "--debug",
     **     "Include debugging symbols in generated object code.");
     **
     **   // Specify options that require values.
     **   optionParser.addOptionWithValue(
     **     "OPTIMIZATION_LEVEL", "-O", "--optimization_level", "0",
     **     "Set optimization level.");
     **   optionParser.addOptionWithValue(
     **     "INCLUDE_DIR", "-I", "--include_dir", "",
     **     "Add INCLUDE_DIR to include search path.");
     **   // ...Note that default values can be non-string types, as long
     **   // ...as they can be converted to strings via the output stream
     **   // ...operator.
     **   optionParser.addOptionWithValue(
     **     "WARNING_LEVEL", "-W", "--warning_level", "3",
     **     "Add INCLUDE_DIR to include search path.");
     **
     **   // Specify required positional arguments.
     **   optionParser.addPositionalArgument(
     **     "INPUT_FILE", "File to be compiled.", true);
     **
     **   // Specify positional arguments which are not required.
     **   optionParser.addPositionalArgument(
     **     "SECOND_FILE", "Additional file to be compiled.", false);
     **
     **   // Parse program arguments.
     **   optionParser.parseCommandLine(argc, argv);
     **
     **   // Recover results of parse.
     **   bool debug = optionParser.getCount("DEBUG") != 0;
     **   int  optimizationLevel =
     **     optionParser.convertValue<int>("OPTIMIZATION_LEVEL");
     **   std::string inputFile0 = optionParser.getValue("INPUT_FILE");
     **   std::string inputFile1 = optionParser.getValue("SECOND_FILE");
     **   std::vector<std::string> includeDirs;
     **   for(size_t includeIndex = 0;
     **       includeIndex < optionParser.getCount("INCLUDE_DIR");
     **       ++includeIndex) {
     **     includeDirs.push_back(
     **       optionParser.getValue("INCLUDE_DIR", includeIndex));
     **   }
     ** @endcode
     **
     ** A command line argument is considered to match the short
     ** version of an option if it begins with the entire text of the
     ** short version.  For example, each of "-p", "-pdf", and "-pp"
     ** would match an option with a short version of "-p".  Arguments
     ** which match the short version of an option, but have
     ** additional text (such as "-pdf" in the example above) are
     ** interpreted in one of two ways.  If constructor argument
     ** allowStackedShortOptions was set to true, and the short
     ** version of the matched option starts with the character '-',
     ** and the matched option does not take an value, then the
     ** remaining text will be prepended with a '-' character and
     ** re-parsed.  In this case, the "-pdf" from the previous example
     ** would be parsed as if it were typed "-p -df".  In all other
     ** cases, the remaining text will be interpreted as an argument
     ** to the matched option, that is "-pdf" would be parsed as "-p"
     ** with the argument "df".
     **
     ** A command line argument is considered to match the long
     ** version of an option if one of three things is true: first,
     ** the text of the argument matches the text of the long version
     ** of the option exactly; second the option was specified with
     ** "allowPartialMatch" set to true, and the entire text of the
     ** argument matches an initial substring of the long version of
     ** the option; or third, the argument begins with a string that
     ** matches under one of the other two cases, and then continues
     ** with a string that starts with an equals sign.  For an option
     ** with long version "--optimization_level", these three cases
     ** are illustrated by the three argument strings
     ** "--optimization_level", "--opt", and "--opt=3".  In the final
     ** case, the text following the equals sign is interpreted as a
     ** value being specified for the option.
     **/
    class OptionParser {
    public:

      /**
       * Default constructor.  OptionParser instances constructed in
       * this way will throw an IOException from parseCommandLine() if
       * the command line doesn't match the specified options.  You
       * can catch the exception and the print usage, etc.
       *
       * @param allowExtraArguments This argument specifies whether
       * extra positional argument should be permitted.  If
       * allowExtraArguments is false, then member function
       * parseCommandLine() will throw an IOException whenever the
       * number of positional arguments on the command line exceeds
       * the number specified using member function
       * addPositionalArgument().  If allowExtraArguments is true,
       * then any extra positional arguments will be recorded without
       * throwing an exception.  These arguments can then be recovered
       * using member function getExtraPositionalArguments().
       *
       * @param allowStackedShortOptions This argument specifies
       * whether it is permissible to decompose arguments such as
       * "-plr" into "-p -l -r".  If allowStackedShortOptions is true,
       * then this type of decomposition will be carried out for short
       * options which do not accept values.  If option "-p" does
       * accept a value, then "-plr" will (as always) be interpreted
       * as "-p" with the argument "lr".  Note that this type of
       * decomposition may lead to ambiguity.  For example, if you
       * have an option "-p" and an option "-pl", then it's hard to
       * know how "-plr" should parse.  If any of your short options
       * are longer than a dash plus a single character, then you may
       * wish disable stacked short options by setting
       * allowStackedShortOptions to false.
       *
       * @param allowOptionishArguments This argument specifies how to
       * handle command line arguments which start with "-", but which
       * don't correspond to any of the specified options.  If
       * allowOptionishArguments is true, then such arguments will be
       * simply recorded as positional arguments.  If
       * allowOptionishArguments is false, then such arguments will
       * cause an IOException to be thrown for "unrecognized option."
       */
      OptionParser(bool allowExtraArguments = true,
                   bool allowStackedShortOptions = true,
                   bool allowOptionishArguments = false);


      /**
       * This constructor specifies that if a malformed commandline is
       * parsed, the OptionParser should exit instead of throwing an
       * exception, optionally printing a usage message first.
       *
       * @param exitCode This argument specifies the exit code which
       * should be used to report a bad command line.
       *
       * @param handleMinusMinusHelp This argument specifies whether
       * the parser should automatically watch of "--help" arguments.
       * If handleMinusMinusHelp is true, the parser will print a
       * usage message and then exit(0) when a "--help" option is
       * encountered.
       *
       * @param printUsage This argument specifies whether, when a
       * malformed commandline is parsed, a usage message should be
       * printed before exiting with the specified exit code.
       *
       * @param allowExtraArguments Please refer to the documentation
       * for the default constructor to find a description of this
       * argument.
       *
       * @param allowStackedShortOptions Please refer to the
       * documentation for the default constructor to find a
       * description of this argument.
       *
       * @param allowOptionishArguments Please refer to the
       * documentation for the default constructor to find a
       * description of this argument.
       */
      explicit
      OptionParser(int exitCode,
                   bool handleMinusMinusHelp = true,
                   bool printUsage = true,
                   bool allowExtraArguments = true,
                   bool allowStackedShortOptions = true,
                   bool allowOptionishArguments = false);


      /**
       * Destructor.
       */
      ~OptionParser();


      /**
       * Specifies summary text that will be printed as part of the
       * usage message after the program synopsis, but before any
       * option or argument descriptions.
       *
       * @param description This argument should briefly describe the
       * function of the program.
       */
      void
      addDescription(std::string const& description) {
        m_description = description;
      }


      /**
       * This member function adds an option to be recognized during
       * command line parsing.
       *
       * @param name This argument specifies a unique name for the
       * option.  It will be used later to refer to this option when
       * checking parsing results.
       *
       * @param shortVersion This argument specifies a short version
       * of the option, such as "-p" or "-rgb".  Setting this argument
       * to "" indicates that this option does not have a short
       * version.
       *
       * @param longVersion This argument specifies a long version of
       * the option, such as "--do_calibration".  Setting this
       * argument to "" indicates that this option does not have a
       * long version.
       *
       * @param docString This argument specifies a breif description
       * for the option.  If the automatic documentation features of
       * the optionParser are not being used, then this argument can
       * be safely set to "".
       *
       * @param allowPartialMatch This argument specifies whether to
       * allow partial matches for the long version of this option.
       * For more information, please see the class comment for
       * OptionParser.
       */
      void
      addOption(std::string const& name,
                std::string const& shortVersion,
                std::string const& longVersion,
                std::string const& docString,
                bool allowPartialMatch = true);


      /**
       * This member function adds an option to be recognized during
       * command line parsing.  It is like the 5-argument version of
       * addOption(), but allows greater control over getUsage()
       * output.
       *
       * @param name This argument specifies a unique name for the
       * option.  It will be used later to refer to this option when
       * checking parsing results.
       *
       * @param shortVersion This argument specifies a short version
       * of the option, such as "-p" or "-rgb".  Setting this argument
       * to "" indicates that this option does not have a short
       * version.
       *
       * @param longVersion This argument specifies a long version of
       * the option, such as "--do_calibration".  Setting this
       * argument to "" indicates that this option does not have a
       * long version.
       *
       * @param docString This argument specifies a breif description
       * for the option.  If the automatic documentation features of
       * the optionParser are not being used, then this argument can
       * be safely set to "".
       *
       * @param sectionName This specifies in which section of the
       * getUsage() output this option should be described.  Set this
       * to the empty string if you are not using sections.
       *
       * @param prefix This string will be prepended to the value of
       * argument name before the option names are alphabetized for in
       * getUsage() output.  Use it if you want to control the order
       * in which options are shown to the user.
       *
       * @param allowPartialMatch This argument specifies whether to
       * allow partial matches for the long version of this option.
       * For more information, please see the class comment for
       * OptionParser.
       */
      void
      addOption(std::string const& name,
                std::string const& shortVersion,
                std::string const& longVersion,
                std::string const& docString,
                std::string const& sectionName,
                std::string const& prefix = "",
                bool allowPartialMatch = true);


      /**
       * This member function adds an option which requires a value to
       * the list which will be recognized during command line
       * parsing.  It is like the 8-argument version of
       * addOption(), but allows greater control over getUsage()
       * output.
       *
       * @param name This argument specifies a unique name for the
       * option.  It will be used later to refer to this option when
       * checking parsing results.
       *
       * @param shortVersion This argument specifies a short version
       * of the option, such as "-p" or "-rgb".  Setting this argument
       * to "" indicates that this option does not have a short
       * version.
       *
       * @param longVersion This argument specifies a long version of
       * the option, such as "--do_calibration".  Setting this
       * argument to "" indicates that this option does not have a
       * long version.
       *
       * @param defaultValue This argument specifies the default value
       * for the option.  If the option is not found on the command
       * line, then the default value will be used.  If the default
       * value is not a std:string, it will be converted to a string
       * using a std::ostringstream and the stream output operator.
       *
       * @param docString This argument specifies a breif description
       * for the option.  If the automatic documentation features of
       * the optionParser are not being used, then this argument can
       * be safely set to "".
       *
       * @param sectionName This specifies in which section of the
       * getUsage() output this option should be described.  Set this
       * to the empty string if you are not using sections.
       *
       * @param prefix This string will be prepended to the value of
       * argument name before the option names are alphabetized for in
       * getUsage() output.  Use it if you want to control the order
       * in which options are shown to the user.
       *
       * @param requireArgument This argument is currently not supported.
       *
       * @param allowPartialMatch This argument specifies whether to
       * allow partial matches for the long version of this option.
       * For more information, please see the class comment for
       * OptionParser.
       *
       * @param allowOptionishValue This argument specifies whether or
       * not the supplied value (on the command line) is permitted to
       * begin with the character '-'.
       */
      template<class Type>
      void
      addOptionWithValue(std::string const& name,
                         std::string const& shortVersion,
                         std::string const& longVersion,
                         Type const& defaultValue,
                         std::string const& docString,
                         std::string const& sectionName = "",
                         std::string const& prefix = "",
                         bool requireArgument = true,
                         bool allowPartialMatch = true,
                         bool allowOptionishValue = false);


      /**
       * This member function adds a named positional argument.  The
       * order in which positional arguments are specified is
       * important.  This is most clearly illustrated with an example:
       *
       * @code
       *   int argc = 4;
       *   char* argv = {"progname", "arg0", "arg1", "arg2"};
       *   OptionParser optionParser;
       *   optionParser.addPositionalArgument("FOO", "");
       *   optionParser.addPositionalArgument("BAR", "");
       *   optionParser.addPositionalArgument("BAZ", "");
       *   optionParser.parseCommandLine(argc, argv);
       *   // The next line returns "arg0".
       *   std::string firstArgument = optionParser.getValue("FOO");
       *   // The next line returns "arg1".
       *   std::string secondArgument = optionParser.getValue("BAR");
       *   // The next line returns "arg2".
       *   std::string thirdArgument = optionParser.getValue("BAZ");
       * @endcode
       *
       * @param name This argument specifies a unique name for the
       * positional argument.  It will be used later to refer to this
       * argument when checking parsing results.
       *
       * @param docString This argument specifies a breif description
       * for the argument.  If the automatic documentation features of
       * the optionParser are not being used, then this argument can
       * be safely set to "".
       *
       * @param isRequired This argument specifies whether the absence
       * of this positional argument should cause the parsing of a
       * command line to fail.  If isRequired is set to false, and
       * this positional argument is missing from the command line,
       * then the default value specified by argument defaultValue
       * will be used.
       */
      void
      addPositionalArgument(std::string const& name,
                            std::string const& docString,
                            bool isRequired = false,
                            std::string const& defaultValue = "");


      /**
       * Use this function to group options into sections.  This can
       * make it easier for the user to understand the text returned
       * by getUsage().
       *
       * @param sectionName This argument identifies the section, but
       * will only be printed in the usage message if argument
       * sectionDescription is an empty string.  Sections will be
       * described in alphabetical order by sectionName.
       *
       * @param sectionDescription This argument describes the
       * section, and will be printed at the beginning of the section.
       * You should make this string include the heading name -- as
       * you'd like it printed -- and any other important information.
       * If sectionDescription is specified as the empty string,
       * argument sectionName will be used instead.
       */
      void
      addSection(std::string const& sectionName,
                 std::string const& sectionDescription = "");


      /**
       * This member function is like getValue(std::string const&),
       * but attempts to convert the returned string to the specified
       * type.  It works for all built-in types, and for user-defined
       * types which define a stream input operator.  If the
       * conversion fails, and you used the OptionParser constructor
       * which specifies an exit code, then the OptionParser will
       * behave as if the command line parsing failed, printing usage
       * and calling exit() as appropriate.  If the conversion fails
       * and the non-exit code constructor was used, then getValue()
       * will throw a ConversionException.
       *
       * @param name This argument specifies the option or positional
       * argument whose value is to be queried.  If name does not
       * match the name argument of a previous call to addOption(),
       * addOptionWithValue(), or addPositionalArgument(), then the
       * conversion will fail.
       *
       * @return The return value is the value for the requested
       * option or positional argument.
       */
      template<class Type>
      Type
      convertValue(std::string const& name);
      template<class Type>


      /**
       * This member function is just like convertValue(std::string
       * const&), except that it allows the user to specify min and
       * max values for the converted value.  If argument clampResult
       * is true, then the converted value will be clamped to the
       * range [lowerBound, upperBound].  If argument clampResult is
       * false, then out-of-range values will result variously exit,
       * print usage and exit, or throw a ConversionException, as if
       * command-line parsing had failed, depending on which
       * constructer arguments were specified.  Note that a converted
       * value equal to upperBound is acceptable.
       *
       * @param name This argument specifies the option or positional
       * argument whose value is to be queried.  If name does not
       * match the name argument of a previous call to addOption(),
       * addOptionWithValue(), or addPositionalArgument(), then the
       * conversion will fail.
       *
       * @param lowerBound This argument is the lowest acceptable
       * value for the converted value.
       *
       * @param upperBound This argument is the highest acceptable
       * value for the converted value.
       *
       * @param clampResult This argument specifies whether
       * out-of-range values should be clamped to the acceptable
       * range, or whether a ConversionException should be thrown.
       *
       * @return The return value is the converted value of the
       * requested command-line argument.
       */
      Type
      convertValue(std::string const& name,
                   Type const& lowerBound,
                   Type const& upperBound,
                   bool clampResult = true);


      /**
       * This member function is like getValue(std::string const&, int),
       * but attempts to convert the returned string to the specified
       * type.  It works for all built-in types, and for user-defined
       * types which define a stream input operator.  If the
       * conversion fails, and you used the OptionParser constructor
       * which specifies an exit code, then the OptionParser will
       * behave as if the command line parsing failed, printing usage
       * and calling exit() as appropriate.  If the conversion fails
       * and the non-exit code constructor was used, then getValue()
       * will throw a ConversionException.
       *
       * @param name This argument is as described in the
       * documentation for getValue(std::string const&, int).
       *
       * @param valueIndex This argument is as described in the
       * documentation for getValue(std::string const&, int).
       *
       * @return The return value is the value for the requested
       * option or positional argument.
       */
      template<class Type>
      Type
      convertValue(std::string const& name, int valueIndex);


      /**
       * This member function is just like convertValue(std::string
       * const&, int), except that it allows the user to specify min
       * and max values for the converted value.  If argument
       * clampResult is true, then the converted value will be clamped
       * to the range [lowerBound, upperBound].  If argument clampResult
       * is false, then out-of-range values will result variously exit,
       * print usage and exit, or throw a ConversionException, as if
       * command-line parsing had failed, depending on which
       * constructer arguments were specified.  Note that a converted
       * value equal to upperBound is acceptable.
       *
       * @param name This argument specifies the option or positional
       * argument whose value is to be queried.  If name does not
       * match the name argument of a previous call to addOption(),
       * addOptionWithValue(), or addPositionalArgument(), then the
       * conversion will fail.
       *
       * @param lowerBound This argument is the lowest acceptable
       * value for the converted value.
       *
       * @param upperBound This argument is the highest acceptable
       * value for the converted value.
       *
       * @param clampResult This argument specifies whether
       * out-of-range values should be clamped to the acceptable
       * range, or whether a ConversionException should be thrown.
       *
       * @return The return value is the converted value of the
       * requested command-line argument.
       */
      template<class Type>
      Type
      convertValue(std::string const& name, int valueIndex,
                   Type const& lowerBound,
                   Type const& upperBound,
                   bool clampResult = true);


      /**
       * This member function indicates how many times the specified
       * option or positional argument was specified in the most
       * recently parsed command line.  For positional arguments, the
       * return value is always 0 or 1.  For options, the return value
       * may be higher, since options may be specified more than once.
       *
       * @param name This argument is as described in the
       * documentation for getValue(std::string const&, int).
       *
       * @return The return value is the value specified on the
       * commandline for the indicated option or positional argument.
       * If the command line did not provide that option or argument,
       * a default value will be returned.
       */
      size_t
      getCount(std::string const& name) const;


      /**
       * If constructor argument allowExtraArguments was set to true,
       * and if number of positional arguments on the most recently
       * parsed command line exceeds the number of named positional
       * arguments specified using addPositionalArgument(), then this
       * member function will return the extra positional arguments,
       * in the order they were found on the command line, as a vector
       * of strings.
       *
       * @return The return value is a vector of extra positional
       * arguments.
       */
      std::vector<std::string>
      getExtraPositionalArguments();


      /**
       * This member function returns a formatted string describing
       * the available command line options, as specified using
       * addOption() and addOptionWithValue().  It is useful for
       * building your own usage messages if you used the default
       * constructor.  If you used the constructor which specifies an
       * exit code, you probably have no use for
       * getOptionsDescription().  The string returned by
       * getOptionsDescription() will be much more informative if the
       * docString arguments of the addOption() and
       * addOptionWithValue() calls were not empty.
       *
       * @return The return value is a formatted string describing the
       * configured options.
       */
      std::string
      getOptionsDescription(std::string const& sectionName = "");


      /**
       * This member function returns a formatted string suitable for
       * a usage message.  It includes the string returned by
       * getOptionsDescription(), so if you call getUsage(), you
       * probably don't need to call getOptionsDescription().  If you
       * used the constructor which specifies an exit code, you
       * probably have no use for getUsage().  The string returned by
       * getOptionsDescription() will be much more informative if the
       * docString arguments of addOption(), addOptionWithValue(), and
       * addPositionalArgument() calls were not empty.
       *
       * @return The return value is a formatted usage message.
       */
      std::string
      getUsage();


      /**
       * Following a successful parse of a command line, this member
       * function returns the requested option value or positional
       * argument value as a string.  If an option was specified more
       * than once, this member function returns the value specified
       * latest on the command line.  This member function is
       * equivalent to calling getValue<Type>(const std::string, int)
       * with the int argument set to -1.  If name does not match any
       * specified option or positional argument, then a
       * ValueException will be thrown.
       *
       * @param name This argument specifies the option or positional
       * argument whose value is to be queried.  If name does not
       * match the name argument of a previous call to addOption(),
       * addOptionWithValue(), or addPositionalArgument(), then the
       * returned value will be the empty string.
       *
       * @return The return value is the value for the requested
       * option or positional argument, or the empty string if no such
       * value was specified on the most recently parsed command line.
       */
      std::string
      getValue(std::string const& name);


      /**
       * Following a successful parse of a command line, this member
       * function returns the (valueIndex)th value of the requested
       * option or positional argument, where 0 means first, 1 means
       * second, etc.
       *
       * @param name This argument specifies the option or positional
       * argument whose value is to be queried.  If name does not
       * match the name argument of a previous call to addOption(),
       * addOptionWithValue(), or addPositionalArgument(), then the
       * returned value will be the empty string.
       *
       * @param valueIndex This argument specifies which occurrence of
       * the option/argument is being queried.  If valueIndex is set
       * to -1, than the last occurrence of the option will be
       * returned.  ValueIndex must always be less than the number
       * returned by getCount(name), and therefore must always be 0 or
       * -1 for positional arguments.
       *
       * @return The return value is the value for the requested
       * option or positional argument, or the empty string if no such
       * value was specified on the most recently parsed command line.
       */
      std::string
      getValue(std::string const& name, int valueIndex);


      /**
       * This member function should be called after all calls to
       * addOption(), addOptionWithValue(), and
       * addPositionalArgument() are completed in order to parse a
       * command line and store the result for subsequent access using
       * the getCount() and getValue() methods.  If the command line
       * doesn't match the expected options & positional arguments,
       * one of several things will happen, depending on which
       * constructor & constructor arguments were used to create the
       * OptionParser.  The default behavior is to throw an
       * IOException, but the OptionParser can also be configured to
       * exit with a specific error code, or to print a usage message
       * before exiting with a specific error code.  If the
       * allowMinusMinusHelp constructor argument was specified, and a
       * "--help" option is encountered, then the OptionParser
       * instance will print a usage message and then call exit(0).
       *
       * @param argc This argument can be taken directly from the argc
       * argument to main().
       *
       * @param argv This argument can be taken directly from the argv
       * argument to main().
       */
      void
      parseCommandLine(int argc, const char* argv[]);


      /**
       * This member function is identical to parseCommandLine(int
       * argc, const char* argv[]), except that it does not require a
       * const qualifier on its second argument.
       *
       * @param argc This argument can be taken directly from the argc
       * argument to main().
       *
       * @param argv This argument can be taken directly from the argv
       * argument to main().
       */
      void
      parseCommandLine(int argc, char* argv[]);


    private:

      // Private typedefs.
      typedef std::map<std::string, OptionDescription> OptionMap;


      // Private enum for specifying how to handle illegal commandlines.
      enum ErrorBehavior {ExitOnError, ThrowOnError, UsageAndExitOnError};


      // Private member functions
      bool
      findOptionDescription(std::string const& argument,
                            OptionDescription& optionDescription,
                            size_t& typedLength,
                            bool& isShortMatch);


      void
      recordPositionalArgument(std::string const& argument);


      // Depending on constructor arguments, either exit, or print usage
      // and exit, or do nothing.
      void
      usageAndExitIfAppropriate(std::string const& what);

      void
      usageAndExitIfAppropriate(brick::common::Exception const& ee) {
        this->usageAndExitIfAppropriate(std::string(ee.what()));
      }



      // Private data members.
      bool m_allowExtraArguments;
      bool m_allowOptionishArguments;
      bool m_allowStackedShortOptions;
      ErrorBehavior m_errorBehavior;
      int m_exitCode;
      std::vector<std::string> m_extraArgumentValues;
      bool m_handleHelp;
      size_t m_numberOfPosArgRequired;
      size_t m_numberOfPosArgParsed;
      std::string m_description;
      std::map<std::string, size_t> m_optionCounts;
      std::map<std::string, OptionDescription> m_optionDescriptions;
      std::map<std::string, OptionMap> m_optionDescriptionsBySection;
      std::map< std::string, std::vector<std::string> > m_optionValues;
      std::map< std::string, std::string > m_sectionDescriptions;
      std::vector<std::string> m_positionalArgumentNames;
      std::vector<std::string> m_positionalArgumentDefaultValues;
      std::vector<std::string> m_positionalArgumentDocs;
      std::string m_programName;

    }; // class OptionParser

  } // namespace utilities

} // namespace brick

/* ========= Definitions of inline and template functions below. ========== */

#include <sstream>

namespace brick {

  namespace utilities {


    // This member function adds an option which requires a value to
    // the list which will be recognized during command line
    // parsing.  This template is specialized for std::string.
    template<>
    void
    OptionParser::
    addOptionWithValue<std::string>(std::string const& name,
                                    std::string const& shortVersion,
                                    std::string const& longVersion,
                                    std::string const& defaultValue,
                                    std::string const& docString,
                                    std::string const& sectionName,
                                    std::string const& prefix,
                                    bool requireArgument,
                                    bool allowPartialMatch,
                                    bool allowOptionishValue);


    // This member function adds an option which requires a value to
    // the list which will be recognized during command line
    // parsing.
    template<class Type>
    void
    OptionParser::
    addOptionWithValue(std::string const& name,
                       std::string const& shortVersion,
                       std::string const& longVersion,
                       Type const& defaultValue,
                       std::string const& docString,
                       std::string const& sectionName,
                       std::string const& prefix,
                       bool requireArgument,
                       bool allowPartialMatch,
                       bool allowOptionishValue)
    {
      std::ostringstream buffer;
      buffer << defaultValue;
      this->addOptionWithValue(
        name, shortVersion, longVersion, buffer.str(), docString,
        sectionName, prefix, requireArgument, allowPartialMatch,
        allowOptionishValue);
    }


    template<class Type>
    Type
    OptionParser::
    convertValue(std::string const& name)
    {
      std::string resultString = this->getValue(name);

      try {
        return convertString<Type>(resultString);
      } catch(ConversionException const&) {
        std::ostringstream message;
        message << "Error when parsing option " << name << ": "
                << "can't convert \"" << resultString
                << "\" to appropriate type.";
        this->usageAndExitIfAppropriate(message.str());
        BRICK_THROW(ConversionException, "OptionParser::convertValue()",
                    message.str().c_str());
      }
    }


    template<class Type>
    Type
    OptionParser::
    convertValue(std::string const& name,
                 Type const& lowerBound,
                 Type const& upperBound,
                 bool clampResult)
    {
      Type result = this->convertValue<Type>(name);

      if(clampResult) {
        result = std::max(result, lowerBound);
        result = std::min(result, upperBound);
      } else if(result < lowerBound || result > upperBound) {
        std::ostringstream message;
        message << "Error when parsing option " << name << ": "
                << "result (" << result << ") is outside of the "
                << "specified range ([" << lowerBound << ", " << upperBound
                << "]).";
        this->usageAndExitIfAppropriate(message.str());
        BRICK_THROW(ConversionException, "OptionParser::convertValue()",
                    message.str().c_str());
      }
      return result;
    }


    template<class Type>
    Type
    OptionParser::
    convertValue(std::string const& name, int valueIndex)
    {
      std::string resultString = this->getValue(name, valueIndex);

      try {
        return convertString<Type>(resultString);
      } catch(ConversionException&) {
        std::ostringstream message;
        message << "Error when parsing option " << name << ": "
                << "can't convert \"" << resultString
                << "\" to appropriate type.";
        this->usageAndExitIfAppropriate(message.str());
        BRICK_THROW(ConversionException, "OptionParser::convertValue()",
                  message.str().c_str());
      }
    }


    template<class Type>
    Type
    OptionParser::
    convertValue(std::string const& name, int valueIndex,
                 Type const& lowerBound,
                 Type const& upperBound,
                 bool clampResult)
    {
      Type result = this->convertValue<Type>(name, valueIndex);

      if(clampResult) {
        result = std::max(result, lowerBound);
        result = std::min(result, upperBound);
      } else if(result < lowerBound || result > upperBound) {
        std::ostringstream message;
        message << "Error when parsing option " << name << ": "
                << "result (" << result << ") is outside of the "
                << "specified range ([" << lowerBound << ", " << upperBound
                << "]).";
        this->usageAndExitIfAppropriate(message.str());
        BRICK_THROW(ConversionException, "OptionParser::convertValue()",
                    message.str().c_str());
      }
      return result;
    }

  } // namespace utilities

} // namespace brick

#endif /* #ifndef BRICK_UTILITES_OPTIONPARSER_HH */
