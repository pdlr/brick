/**
***************************************************************************
* @file brick/utilities/optionParser.cc
*
* Source file defining the OptionParser class.
*
* Copyright (C) 2006-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <algorithm>
#include <sstream>
#include <brick/common/exception.hh>
#include <brick/utilities/optionParser.hh>
#include <brick/utilities/stringManipulation.hh>

namespace brick {

  namespace utilities {

    // Default constructor.
    OptionParser::
    OptionParser(bool allowExtraArguments,
                 bool allowStackedShortOptions,
                 bool allowOptionishArguments)
      : m_allowExtraArguments(allowExtraArguments),
        m_allowOptionishArguments(allowOptionishArguments),
        m_allowStackedShortOptions(allowStackedShortOptions),
        m_errorBehavior(ThrowOnError),
        m_exitCode(0),
        m_extraArgumentValues(),
        m_handleHelp(false),
        m_numberOfPosArgRequired(0),
        m_numberOfPosArgParsed(0),
        m_description(),
        m_optionCounts(),
        m_optionDescriptions(),
        m_optionDescriptionsBySection(),
        m_optionValues(),
        m_positionalArgumentNames(),
        m_positionalArgumentDefaultValues(),
        m_positionalArgumentDocs(),
        m_programName()
    {
      // Make sure there's a section for global options.
      this->addSection("", "Options");
    }


    // Constructor.
    OptionParser::
    OptionParser(int exitCode,
                 bool handleMinusMinusHelp,
                 bool printUsage,
                 bool allowExtraArguments,
                 bool allowStackedShortOptions,
                 bool allowOptionishArguments)
      : m_allowExtraArguments(allowExtraArguments),
        m_allowOptionishArguments(allowOptionishArguments),
        m_allowStackedShortOptions(allowStackedShortOptions),
        m_errorBehavior(printUsage ? UsageAndExitOnError : ExitOnError),
        m_exitCode(exitCode),
        m_extraArgumentValues(),
        m_handleHelp(handleMinusMinusHelp),
        m_numberOfPosArgRequired(0),
        m_numberOfPosArgParsed(0),
        m_description(),
        m_optionCounts(),
        m_optionDescriptions(),
        m_optionDescriptionsBySection(),
        m_optionValues(),
        m_positionalArgumentNames(),
        m_positionalArgumentDefaultValues(),
        m_positionalArgumentDocs(),
        m_programName()
    {
      // Make sure there's a section for global options.
      this->addSection("", "Options:");
    }


    // Destructor.
    OptionParser::
    ~OptionParser()
    {
      // Empty.
    }


    // This member function adds an option to be recognized during
    // command line parsing.
    void
    OptionParser::
    addOption(std::string const& name,
              std::string const& shortVersion,
              std::string const& longVersion,
              std::string const& docString,
              bool allowPartialMatch)
    {
      this->addOption(name, shortVersion, longVersion, docString,
                      "", "", allowPartialMatch);
    }


    // This member function adds an option to be recognized during
    // command line parsing.
    void
    OptionParser::
    addOption(std::string const& name,
              std::string const& shortVersion,
              std::string const& longVersion,
              std::string const& docString,
              std::string const& sectionName,
              std::string const& prefix,
              bool allowPartialMatch)
    {
      // Does this option already exist?
      if(m_optionDescriptions.find(name) != m_optionDescriptions.end()) {
        std::ostringstream message;
        message << "Option \"" << name << "\" has already been specified.";
        BRICK_THROW(brick::common::StateException, "OptionParser::addOption()",
                    message.str().c_str());
      }

      // Figure out to which section this option belongs.
      auto optionMapIter = m_optionDescriptionsBySection.find(sectionName);
      if(optionMapIter == m_optionDescriptionsBySection.end()) {
        this->addSection(sectionName);
        optionMapIter = m_optionDescriptionsBySection.find(sectionName);
      }

      // On with the show.
      std::string alphabetizeAs = prefix + name;
      OptionDescription newOption(
        name, shortVersion, longVersion, false, allowPartialMatch, docString);

      // Add the new optionDescription to both
      // m_optionDescriptionsBySection (via optionMapIter) and
      // m_optionDescriptions.  This lets us access it conveniently
      // both by section and by name.
      optionMapIter->second.insert(std::make_pair(alphabetizeAs, newOption));
      m_optionDescriptions.insert(std::make_pair(name, newOption));
      m_optionCounts.insert(std::make_pair(name, 0));
    }


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
                                    bool allowOptionishValue)
    {
      if(requireArgument == false) {
        BRICK_THROW(brick::common::NotImplementedException,
                    "OptionParser::addOptionWithValue()",
                    "Options with optional arguments are not yet supported.");
      }
      if(allowOptionishValue == true) {
        BRICK_THROW(brick::common::NotImplementedException,
                    "OptionParser::addOptionWithValue()",
                    "Optionish option values are not yet supported.");
      }

      // Does this option already exist?
      if(m_optionDescriptions.find(name) != m_optionDescriptions.end()) {
        std::ostringstream message;
        message << "Option \"" << name << "\" has already been specified.";
        BRICK_THROW(brick::common::StateException,
                    "OptionParser::addOptionWithValue()",
                    message.str().c_str());
      }

      // Figure out to which section this option belongs.
      auto optionMapIter = m_optionDescriptionsBySection.find(sectionName);
      if(optionMapIter == m_optionDescriptionsBySection.end()) {
        this->addSection(sectionName);
        optionMapIter = m_optionDescriptionsBySection.find(sectionName);
      }

      // On with the show.
      std::string alphabetizeAs = prefix + name;
      OptionDescription newOption(
        name, shortVersion, longVersion, true, allowPartialMatch, docString,
        defaultValue);

      // Add the new optionDescription to both
      // m_optionDescriptionsBySection (via optionMapIter) and
      // m_optionDescriptions.  This lets us access it conveniently
      // both by section and by name.
      optionMapIter->second.insert(std::make_pair(alphabetizeAs, newOption));
      m_optionDescriptions.insert(std::make_pair(name, newOption));
      m_optionCounts.insert(std::make_pair(name, 0));
    }


    void
    OptionParser::
    addPositionalArgument(std::string const& name,
                          std::string const& docString,
                          bool isRequired,
                          std::string const& defaultValue)
    {
      m_positionalArgumentNames.push_back(name);
      m_positionalArgumentDocs.push_back(docString);
      m_positionalArgumentDefaultValues.push_back(defaultValue);
      if(isRequired) {
        m_numberOfPosArgRequired = m_positionalArgumentNames.size();
      }
    }


    // Use this function to group options into sections.
    void
    OptionParser::
    addSection(std::string const& sectionName,
               std::string const& sectionDescription)
    {
      if(m_sectionDescriptions.find(sectionName)
         != m_sectionDescriptions.end()) {
        std::ostringstream message;
        message << "Section \"" << sectionName
                << "\" has already been specified.";
        BRICK_THROW(brick::common::StateException,
                    "OptionParser::addSection()",
                    message.str().c_str());
      }

      if(sectionDescription == "") {
        // If no sectionDescription is available, substitute the
        // section name.
        m_sectionDescriptions.insert(
          std::make_pair(sectionName, sectionName));
      } else {
        m_sectionDescriptions.insert(
          std::make_pair(sectionName, sectionDescription));
      }

      m_optionDescriptionsBySection.insert(
        std::make_pair(sectionName, OptionMap()));
    }


    size_t
    OptionParser::
    getCount(std::string const& name) const
    {
      std::map<std::string, size_t>::const_iterator iter =
        m_optionCounts.find(name);
      if(iter == m_optionCounts.end()) {
        std::ostringstream message;
        message << "Invalid option name: " << name << ".";
        BRICK_THROW(brick::common::ValueException, "OptionParser::getCount()",
                    message.str().c_str());
      }
      return iter->second;
    }


    std::vector<std::string>
    OptionParser::
    getExtraPositionalArguments()
    {
      return m_extraArgumentValues;
    }


    std::string
    OptionParser::
    getOptionsDescription(std::string const& sectionName)
    {
      std::ostringstream optionsStream;
      auto optionMapIter = m_optionDescriptionsBySection.find(sectionName);
      if(optionMapIter == m_optionDescriptionsBySection.end()) {
        std::ostringstream message;
        message << "Unrecognized section name: " << sectionName;
        BRICK_THROW(brick::common::ValueException,
                    "OptionParser::getOptionsDescription()",
                    message.str().c_str());
      }

      if(!optionMapIter->second.empty()) {
        auto optionIter = optionMapIter->second.begin();
        while(optionIter != optionMapIter->second.end()) {
          optionsStream << optionIter->second << "\n";
          ++optionIter;
        }
      }
      return optionsStream.str();
    }


    std::string
    OptionParser::
    getUsage()
    {
      std::string indentString = "        ";
      std::ostringstream usageStream;
      usageStream << "Usage:\n"
                  << "  " << m_programName << " ";
      if(!m_optionDescriptions.empty()) {
        usageStream << "[options] ";
      }
      size_t argIndex = 0;
      while(argIndex < m_numberOfPosArgRequired) {
        usageStream << m_positionalArgumentNames[argIndex] << " ";
        ++argIndex;
      }
      while(argIndex < m_positionalArgumentNames.size()) {
        usageStream << "[" << m_positionalArgumentNames[argIndex] << "] ";
        ++argIndex;
      }

      if(m_description != "") {
        usageStream << "\n\n" << m_description << "\n";
      }

      if(!m_positionalArgumentNames.empty()) {
        usageStream << "\n\nPositional arguments:\n";
        for(argIndex = 0; argIndex < m_numberOfPosArgRequired; ++argIndex) {
          usageStream
            << "  " << m_positionalArgumentNames[argIndex] << "\n"
            << wrapString(indentString + m_positionalArgumentDocs[argIndex],
                          indentString)
            << "\n";
        }
        while(argIndex < m_positionalArgumentNames.size()) {
          usageStream
            << "  " << m_positionalArgumentNames[argIndex] << "\n"
            << wrapString(indentString + m_positionalArgumentDocs[argIndex],
                          indentString) << "\n"
            << wrapString((indentString + "Default value: \""
                           + m_positionalArgumentDefaultValues[argIndex]),
                          indentString) << "\""
            << "\n";
          ++argIndex;
        }
      }

      if(!m_optionDescriptions.empty()) {
        for(auto sectionIter = m_sectionDescriptions.begin();
            sectionIter != m_sectionDescriptions.end(); ++sectionIter) {
          usageStream << "\n" << sectionIter->second << "\n";
          usageStream << this->getOptionsDescription(sectionIter->first);
        }
      }

      return usageStream.str();
    }


    std::string
    OptionParser::
    getValue(std::string const& name)
    {
      return this->getValue(name, -1);
    }


    std::string
    OptionParser::
    getValue(std::string const& name, int valueIndex)
    {
      // Look for the specified name in the recorded option/argument values.
      std::map< std::string, std::vector<std::string> >::const_iterator vIter =
        m_optionValues.find(name);
      if(vIter != m_optionValues.end()) {
        // Found an appropriate entry in the recorded values.
        if(valueIndex >= static_cast<int>((vIter->second).size())) {
          std::ostringstream message;
          message << "Option/argument " << name << " was specified "
                  << (vIter->second).size() << " times, "
                  << "yet getValue() was called with valueIndex set to "
                  << valueIndex << ".";
          BRICK_THROW(brick::common::IndexException, "OptionParser::getValue()",
                      message.str().c_str());
        }
        if(valueIndex < 0
           || valueIndex >= static_cast<int>((vIter->second).size())) {
          valueIndex = (vIter->second).size() - 1;
        }
        return (vIter->second)[valueIndex];
      }

      // Nothing appropriate found in the recorded values.  Check to
      // see if the specified name corresponds to an option, and if
      // so, return the default value for that option.
      std::map<std::string, OptionDescription>::const_iterator odIter =
        m_optionDescriptions.find(name);
      if(odIter != m_optionDescriptions.end()) {
        return (odIter->second).getDefaultValue();
      }

      // No appropriate option found.  Check to see if the specified
      // name corresponds to a positional argument, and if so, return
      // the default value for that positional argument.
      std::vector<std::string>::const_iterator panIter =
        std::find(m_positionalArgumentNames.begin(),
                  m_positionalArgumentNames.end(),
                  name);
      if(panIter != m_positionalArgumentNames.end()) {
        size_t argIndex = panIter - m_positionalArgumentNames.begin();
        return m_positionalArgumentDefaultValues[argIndex];
      }

      // // No appropriate positional argument found.  Return the empty
      // // string.
      // return "";

      // No appropriate positional argument found.
      std::ostringstream message;
      message << "The value of Option/argument " << name << " was requested, "
              << "but no option or argument with this name exists.";
      BRICK_THROW(brick::common::ValueException, "OptionParser::getValue()",
                  message.str().c_str());
    }


    void
    OptionParser::
    parseCommandLine(int argc, const char* argv[])
    {
      // An enum representing the possible states of the parser below.
      enum ParseState {
        BRICK_OP_GETTING_NEW_ARG,
        BRICK_OP_IDENTIFYING_ARG,
        BRICK_OP_RECORDING_POS_ARG,
        BRICK_OP_RECORDING_OPTION,
        BRICK_OP_RECORDING_VALUE,
        BRICK_OP_FINISHING_SHORT_OPTION,
        BRICK_OP_FINISHING_LONG_OPTION,
        BRICK_OP_FINISHED,
        BRICK_OP_ERROR
      };

      if(argc < 1) {
        BRICK_THROW(brick::common::ValueException,
                    "OptionParser::parseCommandLine()",
                    "Argument argc must be greater than or equal to 1.");
      }

      // Reset state for new parse.
      m_optionValues.clear();
      std::map<std::string, size_t>::iterator iter = m_optionCounts.begin();
      while(iter != m_optionCounts.end()) {
        iter->second = 0;
        ++iter;
      }
      m_extraArgumentValues.clear();

      // Record "special" arguments.
      m_programName = argv[0];

      // Parse the remaining command line, catching any errors.
      try {
        int argIndex = 1;
        std::string currentArgument = "";
        OptionDescription optionDescription;
        size_t typedLength;
        bool isShortMatch;
        std::ostringstream errorMessage;

        ParseState currentState = BRICK_OP_GETTING_NEW_ARG;
        while(currentState != BRICK_OP_FINISHED) {
          switch(currentState) {
          case BRICK_OP_GETTING_NEW_ARG:
            if(argIndex >= argc) {
              currentState = BRICK_OP_FINISHED;
            } else {
              currentArgument = argv[argIndex];
              ++argIndex;
              currentState = BRICK_OP_IDENTIFYING_ARG;
            }
            break;
          case BRICK_OP_IDENTIFYING_ARG:
            // Special case for --help arguments.
            if(m_handleHelp && (currentArgument == "--help")) {
              std::cout << this->getUsage() << std::endl;
              exit(0);
            }

            // See if any of the user-specified options match.
            {
              bool optionFound = this->findOptionDescription(
                currentArgument, optionDescription, typedLength, isShortMatch);
              if(!optionFound) {
                currentState = BRICK_OP_RECORDING_POS_ARG;
              } else {
                currentState = BRICK_OP_RECORDING_OPTION;
              }
            }
            break;
          case BRICK_OP_RECORDING_POS_ARG:
            // This argument doesn't match any of our options.  It must
            // be a positional argument.
            if(!m_allowOptionishArguments && currentArgument[0] == '-') {
              errorMessage << "Unrecognized option \"" << currentArgument
                           << "\".";
              currentState = BRICK_OP_ERROR;
            } else {
              this->recordPositionalArgument(currentArgument);
              currentState = BRICK_OP_GETTING_NEW_ARG;
            }
            break;
          case BRICK_OP_RECORDING_OPTION:
            // Increment the count for this option.
            ++((this->m_optionCounts.find(optionDescription.getName()))
               ->second);

            // And make sure we do the appropriate followup.
            if(optionDescription.requiresValue()) {
              currentState = BRICK_OP_RECORDING_VALUE;
            } else if(isShortMatch) {
              currentState = BRICK_OP_FINISHING_SHORT_OPTION;
            } else {
              currentState = BRICK_OP_FINISHING_LONG_OPTION;
            }
            break;
          case BRICK_OP_RECORDING_VALUE:
            // This argument does match one of our options.
            if(typedLength < currentArgument.size()) {
              this->m_optionValues[optionDescription.getName()].push_back(
                currentArgument.substr(typedLength, std::string::npos));
              currentState = BRICK_OP_GETTING_NEW_ARG;
            } else if(argIndex < argc) {
              this->m_optionValues[optionDescription.getName()].push_back(
                argv[argIndex]);
              ++argIndex;
              currentState = BRICK_OP_GETTING_NEW_ARG;
            } else {
              errorMessage << "Missing value for option \"" << currentArgument
                           << "\".";
              currentState = BRICK_OP_ERROR;
            }
            break;
          case BRICK_OP_FINISHING_SHORT_OPTION:
            if(typedLength >= currentArgument.size()) {
              currentState = BRICK_OP_GETTING_NEW_ARG;
            } else if(m_allowStackedShortOptions
                      && typedLength >= 2
                      && currentArgument[0] == '-') {
              currentArgument =
                "-" + currentArgument.substr(typedLength, std::string::npos);
              currentState = BRICK_OP_IDENTIFYING_ARG;
            } else {
              errorMessage
                << "Option \""
                << currentArgument.substr(0, typedLength)
                << "\" does not take a value, yet was specified as \""
                << currentArgument << ".";
              currentState = BRICK_OP_ERROR;
            }
            break;
          case BRICK_OP_FINISHING_LONG_OPTION:
            if(typedLength >= currentArgument.size()) {
              currentState = BRICK_OP_GETTING_NEW_ARG;
            } else {
              errorMessage
                << "Option \""
                << currentArgument.substr(0, typedLength)
                << "\" does not take a value, yet was specified as \""
                << currentArgument << ".";
              currentState = BRICK_OP_ERROR;
            }
            break;
          case BRICK_OP_ERROR:
            BRICK_THROW(brick::common::IOException,
                        "OptionParser::parseCommandLine()",
                        errorMessage.str().c_str());
            break;
          default:
            // Should never get here.
            break;
          }
        }

        // Make sure we got the right number of positional arguments.
        if(m_numberOfPosArgParsed < m_numberOfPosArgRequired) {
          std::ostringstream message;
          message << "Expected " << m_numberOfPosArgRequired
                  << " positional argument(s), but only found "
                  << m_numberOfPosArgParsed << ".";
          BRICK_THROW(brick::common::IOException,
                      "OptionParser::parseCommandLine()",
                      message.str().c_str());
        }
        if(!m_allowExtraArguments && (m_extraArgumentValues.size() != 0)) {
          std::ostringstream message;
          message << "Expected no more than "
                  << m_positionalArgumentNames.size()
                  << " positional arguments, but found "
                  << m_numberOfPosArgParsed << ".";
          BRICK_THROW(brick::common::IOException,
                      "OptionParser::parseCommandLine()",
                      message.str().c_str());
        }

      } catch(const brick::common::IOException& caughtException) {
        // If we don't exit, simply re-throw the exception.
        this->usageAndExitIfAppropriate(caughtException);
        throw;
      }
    }


    // This member function is identical to parseCommandLine(int
    // argc, const char* argv[]), except that it does not require a
    // const qualifier on its second argument.
    void
    OptionParser::
    parseCommandLine(int argc, char* argv[])
    {
      const char** newArgv = new const char*[argc];
      for(int index0 = 0; index0 < argc; ++index0) {
        newArgv[index0] = argv[index0];
      }
      try {
        this->parseCommandLine(argc, newArgv);
      } catch(...) {
        delete[] newArgv;
        throw;
      }
      delete[] newArgv;
    }


    bool
    OptionParser::
    findOptionDescription(std::string const& argument,
                          OptionDescription& optionDescription,
                          size_t& typedLength,
                          bool& isShortMatch)
    {
      // Ugly algorithm here makes arg parsing be O(N^2) where N the
      // number of args.  For now, the 80/20 rule dictates that we not
      // sweat about it.
      typedLength = 0;
      bool foundOption = false;
      std::map<std::string, OptionDescription>::const_iterator iter =
        m_optionDescriptions.begin();
      while(iter != m_optionDescriptions.end()) {
        bool localIsShortMatch;
        size_t localTypedLength =
          (iter->second).getMatchLength(argument, localIsShortMatch);
        if(localTypedLength != 0) {
          if(foundOption) {
            std::ostringstream message;
            message << "Argument \"" << argument << "\" is ambiguous.";
            BRICK_THROW(brick::common::IOException,
                        "OptionParser::findOptionDescription()",
                        message.str().c_str());
          } else {
            typedLength = localTypedLength;
            foundOption = true;
            optionDescription = iter->second;
            isShortMatch = localIsShortMatch;
          }
        }
        ++iter;
      }
      return foundOption;
    }


    void
    OptionParser::
    recordPositionalArgument(std::string const& argument)
    {
      if(m_numberOfPosArgParsed < m_positionalArgumentNames.size()) {
        std::string argumentName =
          m_positionalArgumentNames[m_numberOfPosArgParsed];
        m_optionValues[argumentName].push_back(argument);
      } else {
        m_extraArgumentValues.push_back(argument);
      }
      ++m_numberOfPosArgParsed;
    }


    // Depending on constructor arguments, either exit, or print usage
    // and exit, or do nothing.
    void
    OptionParser::
    usageAndExitIfAppropriate(std::string const& what)
    {
      if(m_errorBehavior == UsageAndExitOnError) {
        std::cout << what << "\n\n";
        std::cout << this->getUsage() << std::endl;
        exit(m_exitCode);
      } else if(m_errorBehavior == ExitOnError) {
        exit(m_exitCode);
      }
    }

  } // namespace utilities

} // namespace brick
