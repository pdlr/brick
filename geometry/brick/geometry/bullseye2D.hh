/**
***************************************************************************
* @file brick/geometry/bullseye2D.hh
*
* Header file declaring the Bullseye2D class template.
*
* Copyright (C) 2008-2013 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_GEOMETRY_BULLSEYE2D_HH
#define BRICK_GEOMETRY_BULLSEYE2D_HH

#include <iostream>
#include <brick/geometry/ellipse2D.hh>
#include <brick/numeric/array1D.hh>
#include <brick/numeric/vector2D.hh>

namespace brick {

  namespace geometry {
    
    /**
     ** The Bullseye2D class represents an bullseye in 2D space.  The
     ** bullseye is characterized by a center point, and then a number
     ** of concentric rings surrounding the center point.  These rings
     ** may be elliptical, as if viewing a circular bullseye through
     ** oblique orthogonal projection.
     **/
    template <class Type>
    class Bullseye2D {
    public:
      
      /** 
       * The default constructor initializes to a 3-ring circular
       * bullseye, with rings at radii of 0.5, 1.0, and 1.5.
       */
      inline
      Bullseye2D();

      
      /** 
       * This constructor initializes to a circular bullseye, with a
       * user-specified number of rings at radii of 0.5, 1.0, 1.5,
       * etc.
       * 
       * @param numberOfRings This argument must be 1 or greater, and
       * specifies the number of rings in the bullseye.
       */
      explicit inline
      Bullseye2D(brick::common::UInt32 numberOfRings);

      
      /** 
       * Construct a bullseye, explicitly setting parameters.  TBD:
       * provide example documentation here.
       * 
       * @param ellipse This argument specifies an ellipse that
       * describes each of the rings of the bullseye.  Each ring will
       * simply be a scaled version of this ellipse (centered on the
       * same origin).
       * 
       * @param scalesBegin This argument is the beginning of of a
       * sequence of Type scale parameters, each of which
       * specifies the size of one of the bulleye rings as a
       * proportion of the ellipse specified by the other constructor
       * argument.  The bullseye will have as many rings as there are
       * elements in this sequence.
       * 
       * @param scalesEnd This argument is the end iterator for the
       * sequence that starts with argument scalesBegin.
       */
      template <class Iter>
      Bullseye2D(Ellipse2D<Type> const& ellipse,
                 Iter scalesBegin, Iter scalesEnd);
      
      
      /** 
       * The copy constructor deep copies its argument.
       * 
       * @param source This argument is the class instance to be copied.
       */
      inline
      Bullseye2D(Bullseye2D<Type> const& source);


      /** 
       * Destructor.
       */
      ~Bullseye2D() {}


      /** 
       * The assignment operator deep copies its argument.
       * 
       * @param source This argument is the class instance to be copied.
       * 
       * @return The return value is a reference to *this.
       */
      inline Bullseye2D<Type>&
      operator=(Bullseye2D<Type> const& source);


      /** 
       * Estimate bullseye parameters from a series of points on the
       * bullseye.  After calling this member function, *this will
       * match the input points as closely as possible.  The
       * estimation performed in this member is an extension of the
       * closed form ellipse estimation algorithm of Halir and Flusser
       * [1], which is based on Fitzgibbon's work [2].
       *
       * [1] Halir, R., and Flusser, J., "Numerically Stable Direct
       * Least Squares Fitting Of Bullseyes." 1998.
       *
       * [2] Fitzgibbon, A. W., Pilu, M, and Fischer, R. B., "Direct
       * Least Squares Fitting of Bullseyes." Proc. of the 13th
       * International Conference on Pattern Recognition, pp 253–257,
       * 1996.
       *
       * @param pointsBeginIter This argument is an iterator pointing
       * to the beginning of a sequence of Vector2D<Type> instances.
       * Points belonging to the first ring of the bullseye should
       * preceed points belonging to the second ring.  Points
       * belonging to the second ring should proceed points belonging
       * to the third ring, and so forth.
       * 
       * @param pointsEndIter This argument is an iterator pointing to the
       * end of the sequence of Vector2D<Type> instances.
       * 
       * @param countsBeginIter This argument is an iterator pointing
       * to the beginning of a sequence of counts indicating how many
       * of the input points belong to the first ring of the bullseye,
       * how many to the second, and so on.  The number of rings in
       * the bullseye will be updated to match the length of this
       * sequence.
       * 
       * @param countsEndIter This argument is an iterator pointing
       * to the end of the sequence of counts
       * 
       * @param computeResidual This argument indicates whether or not
       * an overal residual value should be computed.
       * 
       * @return If computeResidual is true, then the return value is
       * the computed algebraic residual, else the return value is 0.
       */
      template<class PointsIterType, class CountsIterType>
      Type
      estimate(PointsIterType pointsBeginIter, PointsIterType pointsEndIter,
               CountsIterType countsBeginIter, CountsIterType countsEndIter,
               bool computeResidual = true);


      /** 
       * This member function duplicates the functionality of the
       * other estimate() member function, but returns algebraic
       * residuals for each input point.
       *
       * @param pointsBeginIter This argument is an iterator pointing
       * to the beginning of a sequence of Vector2D<Type> instances.
       * Points belonging to the first ring of the bullseye should
       * preceed points belonging to the second ring.  Points
       * belonging to the second ring should proceed points belonging
       * to the third ring, and so forth.
       * 
       * @param pointsEndIter This argument is an iterator pointing to the
       * end of the sequence of Vector2D<Type> instances.
       * 
       * @param countsBeginIter This argument is an iterator pointing
       * to the beginning of a sequence of counts indicating how many
       * of the input points belong to the first ring of the bullseye,
       * how many to the second, and so on.  The number of rings in
       * the bullseye will be updated to match the length of this
       * sequence.
       * 
       * @param countsEndIter This argument is an iterator pointing
       * to the end of the sequence of counts
       * 
       * @param residualIter This argument points to the beginning of
       * a sequence through which the residual associated with each
       * input will be returned.
       *
       * @param computeResidual Setting this argument to false will
       * prevent the residual from being computed.  Only worthwhile if
       * you're _really_ in a hurry.
       */
      template<class PointsIterType, class CountsIterType, class ResidualIter>
      void
      estimate(PointsIterType pointsBeginIter, PointsIterType pointsEndIter,
               CountsIterType countsBeginIter, CountsIterType countsEndIter,
               ResidualIter residualIter, bool computeResidual = true);



      /** 
       * Returns an ellipse that describes each of the rings of the
       * bullseye.  Each ring is simply be a scaled version of this
       * ellipse (centered on the same origin).  The scales of the
       * rings can be found by calling member function getScales().
       * 
       * @return The return value is an ellipse describing the rings
       * of the bullseye.
       */
      Ellipse2D<Type>
      getEllipse() const;
      

      /** 
       * Return the number of concentric rings this bullseye has.
       * 
       * @return The return value is the number of rings.
       */
      unsigned int
      getNumberOfRings() const {return m_scales.size();}

      
      /** 
       * This member function returns the geometric center of the bullseye.
       * 
       * @return The return value is the point at the centroid of the
       * bullseye.
       */
      brick::numeric::Vector2D<Type> const&
      getOrigin() const {return m_origin;}
      

      /** 
       * This member function returns a sequence of scale parameters,
       * each of which specifies the size of one of the bulleye rings
       * as a proportion of the ellipse returned by member function
       * getEllipse().  The bullseye has as many rings as there are
       * elements in this sequence. TBD: provide example documentation
       * here.
       * 
       * @return The return value is a sequence of scale parameters.
       */
      std::vector<Type> const&
      getScales() const {return m_scales;}

      
      /** 
       * This member function returns a vector pointing from the
       * center of the bullseye to the point on a ring of the bullseye
       * that is farthest from the center.  Note that there are two
       * such farthest points on opposite sides of each ring.  The
       * vector returned by this member function will remain
       * consistent for the life of the bullseye, and the semimajor
       * axis of each ring will point in the same direction.
       * 
       * @param ringNumber This argument specifies which ring of the
       * bullseye will be considered.  Setting this to zero indicates
       * the first ring, setting it to one indicates the second, and
       * so on.  This number must be less than the number of rings in
       * the bullseye.
       * 
       * @return The return value is a vector pointing along the
       * semimajor axis of the specified ring of the bullseye.
       */
      brick::numeric::Vector2D<Type>
      getSemimajorAxis(unsigned int ringNumber) const {
        return m_scales[ringNumber] * m_semimajorAxis;
      }
      

      /** 
       * This member function returns a vector pointing from the
       * center of the bullseye to the point on a ring of the bullseye
       * that is closest to the center.  Note that there are two
       * such farthest points on opposite sides of each ring.  The
       * vector returned by this member function will remain
       * consistent for the life of the bullseye, and the semiminor
       * axis of each ring will point in the same direction.
       * 
       * @param ringNumber This argument specifies which ring of the
       * bullseye will be considered.  Setting this to zero indicates
       * the first ring, setting it to one indicates the second, and
       * so on.  This number must be less than the number of rings in
       * the bullseye.
       * 
       * @return The return value is a vector pointing along the
       * semiminor axis of the specified ring of the bullseye.
       */
      brick::numeric::Vector2D<Type>
      getSemiminorAxis(unsigned int ringNumber) const {
        return m_scales[ringNumber] * m_semiminorAxis;
      }


      /** 
       * Sets the position of the bullseye center.
       * 
       * @param origin This argument is a Vector2D instance describing
       * the desired origin.
       */
      void
      setOrigin(brick::numeric::Vector2D<Type> const& origin) {
        m_origin = origin;
      }

    private:
      // Private member functions.

      // Convert from implicit bullseye representation to trigonometric
      void
      convertAlgebraicToTrigonometric(
        brick::numeric::Array1D<Type> const& algebraicParameters,
        brick::numeric::Vector2D<Type>& origin,
        brick::numeric::Vector2D<Type>& semimajorAxis,
        brick::numeric::Vector2D<Type>& semiminorAxis,
        std::vector<Type>& scales);
      

      // Private data members.
      brick::numeric::Vector2D<Type> m_origin;
      brick::numeric::Vector2D<Type> m_semimajorAxis;
      brick::numeric::Vector2D<Type> m_semiminorAxis;
      std::vector<Type> m_scales;

    }; // class Bullseye2D



    /* ======= Non-member functions. ======= */

    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream, Bullseye2D<Type> const& bullseye);
    
  } // namespace geometry
    
} // namespace brick


// Include definitions of inline and template functions.
#include <brick/geometry/bullseye2D_impl.hh>

#endif /* #ifndef BRICK_GEOMETRY_BULLSEYE2D_HH */
