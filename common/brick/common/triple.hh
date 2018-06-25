/**
***************************************************************************
* @file brick/common/triple.hh
*
* Header file declaring Triple class.
*
* Copyright (C) 2003-2010 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_COMMON_TRIPLE_HH
#define BRICK_COMMON_TRIPLE_HH

namespace brick {

  namespace common {

    /**
     ** The Triple class provides a convenient way to pass groups of
     ** three things around.  It's intended to be just like std::pair,
     ** but with three elements.  It is based on Stroustrup's implementation
     ** of class pair, taken from "The C++ Programming Language," 3rd
     ** Edition, Addison Wesley, 1997-2000.
     **/
    template <class Type0, class Type1, class Type2>
    class Triple
    {
    public:
      /* ============ Public typedefs ============ */

      /// This describes the type of the first element of the triple.
      typedef Type0 first_type;

      /// This describes the type of the second element of the triple.
      typedef Type1 second_type;

      /// This describes the type of the Third element of the triple.
      typedef Type2 third_type;


      /* ============ Public member functions ============ */
      /**
       * Default constructor initializes all members to default values.
       */
      Triple() : first(Type0()), second(Type1()), third(Type2()) {}


      /**
       * This constructor copies each of its arguments into the triple.
       *
       * @param element0 This argument will be copied to the first element
       * of the new Triple.
       * @param element1 This argument will be copied to the first element
       * of the new Triple.
       * @param element2 This argument will be copied to the first element
       * of the new Triple.
       */
      Triple(const Type0& element0, const Type1& element1, const Type2& element2)
        : first(element0), second(element1), third(element2) {}


      /**
       * This constructor copies each element of another triple.  Note
       * that the types of the elements of the copied triple will be
       * explicitly cast to types of the elements of the the constructed
       * triple.
       *
       * @param source This argument is the triple to be copied.
       */
      template <class OtherType0, class OtherType1, class OtherType2>
      Triple(const Triple<OtherType0, OtherType1, OtherType2>& other)
        : first(static_cast<Type0>(other.first)),
          second(static_cast<Type1>(other.second)),
          third(static_cast<Type2>(other.third)) {}


      /**
       * Empty destructor.
       */
      ~Triple() {}


      /* ============ Public data members ============ */

      /**
       ** This public data member holds the first element of the triple.
       **/
      Type0 first;

      /**
       ** This public data member holds the second element of the triple.
       **/
      Type1 second;

      /**
       ** This public data member holds the third element of the triple.
       **/
      Type2 third;
    };


    /**
     * makeTriple is a convenience function for creating Triples.  It is
     * intended to be just like std::make_pair().
     *
     * @param element0 This argument specifies the first element of the triple.
     * @param element1 This argument specifies the second element of the triple.
     * @param element2 This argument specifies the third element of the triple.
     * @return A Triple instance containing the three arguments.
     */
    template <class Type0, class Type1, class Type2>
    inline Triple<Type0, Type1, Type2>
    makeTriple(const Type0& element0, const Type1& element1,
               const Type2& element2) {
      return Triple<Type0, Type1, Type2>(element0, element1, element2);
    }

  } // namespace common

} // namespace brick

#endif /* #ifndef BRICK_COMMON_TRIPLE_HH */
