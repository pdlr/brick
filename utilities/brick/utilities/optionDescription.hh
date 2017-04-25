/**
***************************************************************************
* @file brick/utilities/optionDescription.hh
*
* Header file declaring the OptionDescription class.
*
* Copyright (C) 2006-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_UTILITIES_OPTIONDESCRIPTION_HH
#define BRICK_UTILITIES_OPTIONDESCRIPTION_HH

#include <ostream>
#include <string>

namespace brick {

  namespace utilities {
    
    /**
     ** The OptionDescription class
     **/
    class OptionDescription {
    public:

      /** 
       * Default constructor.
       */
      OptionDescription();

      
      /** 
       * Constructor.
       */
      OptionDescription(std::string const& name,
                        std::string const& shortAppearance,
                        std::string const& longAppearance,
                        bool requireValue,
                        bool allowPartialMatch,
                        std::string const& docString,
                        std::string const& defaultValue = "");
      
      /**
       * Destructor.
       */
      ~OptionDescription();

      std::string const&
      getAlphabetizationKey() const {return m_key;}
      

      std::string const&
      getDefaultValue() const {return m_defaultValue;}

      
      std::string const&
      getName() const {return m_name;}


      size_t
      getMatchLength(std::string const& argument,
                     bool& isShortMatch) const;

      
      std::string
      getOptionDoc() const;
      

      bool
      isMatch(std::string const& argument) const;
      
      
      bool
      requiresValue() const {return m_requiresValue;}

    private:

      std::string m_key;
      std::string m_name;
      std::string m_shortAppearance;
      std::string m_longAppearance;
      bool m_requiresValue;
      bool m_allowPartialMatch;
      std::string m_docString;
      std::string m_defaultValue;
    };


    std::ostream&
    operator<<(std::ostream& stream,
               const OptionDescription& optionDescription);

  } // namespace utilities

} // namespace brick

#endif /* #ifndef BRICK_UTILITES_OPTIONDESCRIPTION_HH */
