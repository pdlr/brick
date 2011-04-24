/**
***************************************************************************
* @file brick/utilities/optionDescription.cc
* 
* Source file defining the OptionDescription class.
*
* Copyright (C) 2006-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <sstream>
#include <brick/utilities/stringManipulation.hh>
#include <brick/utilities/optionDescription.hh>

namespace brick {

  namespace utilities {
  
    // Constructor.
    OptionDescription::
    OptionDescription()
      : m_name(),
        m_shortAppearance(),
        m_longAppearance(),
        m_requiresValue(),
        m_allowPartialMatch(),
        m_docString(),
        m_defaultValue()
    {
      // Empty.
    }
    
      
    // Constructor.
    OptionDescription::
    OptionDescription(std::string const& name,
                      std::string const& shortAppearance,
                      std::string const& longAppearance,
                      bool requiresValue,
                      bool allowPartialMatch,
                      std::string const& docString,
                      std::string const& defaultValue)
      : m_name(name),
        m_shortAppearance(shortAppearance),
        m_longAppearance(longAppearance),
        m_requiresValue(requiresValue),
        m_allowPartialMatch(allowPartialMatch),
        m_docString(docString),
        m_defaultValue(defaultValue)
    {
      // Empty.
    }


    // Destructor.
    OptionDescription::
    ~OptionDescription()
    {
      // Empty.
    }


    size_t
    OptionDescription::
    getMatchLength(std::string const& argument,
                   bool& isShortMatch) const
    {
      if(m_shortAppearance.size() != 0) {
        if((argument.size() == m_shortAppearance.size())
           && (argument == m_shortAppearance)) {
          isShortMatch = true;
          return argument.size();
        }
        if((argument.size() > m_shortAppearance.size())
           && (argument.compare(0, m_shortAppearance.size(), m_shortAppearance)
               == 0)) {
          isShortMatch = true;
          return m_shortAppearance.size();
        }
      }
      if(m_longAppearance.size() != 0) {
        // Split arguments of the form "--foobar=baz".
        std::vector<std::string> argumentParts =
          splitString(argument, "=", true, 1);
        if(argumentParts[0].size() == m_longAppearance.size()
           && argumentParts[0] == m_longAppearance) {
          if(argumentParts.size() == 1) {
            isShortMatch = false;
            return m_longAppearance.size();
          } else {
            isShortMatch = false;
            return m_longAppearance.size() + 1;
          }
        }
        if(m_allowPartialMatch
           && argumentParts[0].size() > 0
           && argumentParts[0].size() < m_longAppearance.size()
           && (m_longAppearance.compare(0, argumentParts[0].size(), argument)
               == 0)) {
          if(argumentParts.size() == 1) {
            isShortMatch = false;
            return argumentParts[0].size();
          } else {
            isShortMatch = false;
            return argumentParts[0].size() + 1;
          }
        }
      }
      isShortMatch = false;
      return 0;
    }


    std::string
    OptionDescription::
    getOptionDoc() const
    {
      std::string indentString = "        ";
      std::ostringstream docStream;
      if(m_shortAppearance != "") {
        if(m_requiresValue) {
          docStream << m_shortAppearance << " " << m_name;
        } else {
          docStream << m_shortAppearance;
        }
      }
      if(m_shortAppearance != "" && m_longAppearance != "") {
        docStream << ", ";
      }
      if(m_longAppearance != "") {
        if(m_requiresValue) {
          docStream << m_longAppearance << "=" << m_name;
        } else {
          docStream << m_longAppearance;
        }
      }
      docStream << "\n"
                << wrapString(indentString + m_docString, indentString);
      if(m_defaultValue != "") {
        docStream << "\n"
                  << wrapString(indentString + "Default value: \""
                                + m_defaultValue + "\"", indentString);
      }
      return docStream.str();
    }

    
    bool
    OptionDescription::
    isMatch(std::string const& argument) const
    {
      bool dummy;
      return (this->getMatchLength(argument, dummy) != 0);
    }

    
    std::ostream&
    operator<<(std::ostream& stream,
               const OptionDescription& optionDescription)
    {
      stream << "  " << optionDescription.getOptionDoc();
      return stream;
    }
    
      
  } // namespace utilities

} // namespace brick
