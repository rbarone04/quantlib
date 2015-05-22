/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Ferdinando Ametrano

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file compositediscountcurve.hpp
    \brief Composite discount curve
*/

#ifndef quantlib_composite_discount_curve_hpp
#define quantlib_composite_discount_curve_hpp

#include <ql/termstructures/yieldtermstructure.hpp>

namespace QuantLib {

    //! Composite discount curve joining two different curves at a given date
    /*! The CompositeDiscountCurve inherits reference date, DayCounter,
        Calendar, and settlement days from the first curve,
        max date from the second one.

        \note This term structure will remain linked to the original
              structures, i.e., any changes in the latters will be
              reflected in this structure as well.

        \ingroup yieldtermstructures

        \test
        - todo
    */
    class CompositeDiscountCurve : public YieldTermStructure {
      public:
        CompositeDiscountCurve(const Handle<YieldTermStructure>& first,
                               const Handle<YieldTermStructure>& second,
                               const Date& joinDate,
                               bool allowExtrapolatedJunction,
                               bool allowExtrapolation);
        //! the date at which discount = 1.0
        const Date& referenceDate() const;
        //! \name CompositeDiscountCurve inspectors
        //@{
        const Handle<YieldTermStructure>& first() const;
        const Handle<YieldTermStructure>& second() const;
        Date joinDate() const;
        bool allowExtrapolatedJunction() const;
        bool allowExtrapolation() const;
        //@}
        //! \name YieldTermStructure interface
        //@{
        DayCounter dayCounter() const;
        Calendar calendar() const;
        Natural settlementDays() const;
        Date maxDate() const;
      protected:
        DiscountFactor discountImpl(Time) const;
        //@}
      private:
        Handle<YieldTermStructure> first_, second_;
        Date joinDate_;
        bool allowExtrapolatedJunction_, allowExtrapolation_;
    };


    // inline definitions

    inline
    const Handle<YieldTermStructure>& CompositeDiscountCurve::first() const {
        return first_;
    }

    inline
    const Handle<YieldTermStructure>& CompositeDiscountCurve::second() const {
        return second_;
    }

    inline Date CompositeDiscountCurve::joinDate() const {
        return joinDate_;
    }

    inline bool CompositeDiscountCurve::allowExtrapolatedJunction() const {
        return allowExtrapolatedJunction_;
    }

    inline bool CompositeDiscountCurve::allowExtrapolation() const {
        return allowExtrapolation_;
    }

    inline const Date& CompositeDiscountCurve::referenceDate() const {
        return first_->referenceDate();
    }

    inline DayCounter CompositeDiscountCurve::dayCounter() const {
        return first_->dayCounter();
    }

    inline Calendar CompositeDiscountCurve::calendar() const {
        return first_->calendar();
    }

    inline Natural CompositeDiscountCurve::settlementDays() const {
        return first_->settlementDays();
    }

    inline Date CompositeDiscountCurve::maxDate() const {
        return second_->maxDate();
    }

}

#endif
