/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2008 Roland Lichters

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

#include <numeric>

#include <ql/experimental/credit/integralntdengine.hpp>
#include <ql/cashflows/fixedratecoupon.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/experimental/credit/basket.hpp>

namespace QuantLib {

    void IntegralNtdEngine::calculate() const {
        Date today = Settings::instance().evaluationDate();

        results_.errorEstimate = Null<Real>();
        results_.value = 0.0;
        results_.premiumValue = 0.0;
        results_.upfrontPremiumValue = 0.;
        Real accrualValue = 0.0;
        Real claimValue = 0.0;
        Date d, d0;
        for (Size i = 0; i < arguments_.premiumLeg.size(); i++) {
            boost::shared_ptr<FixedRateCoupon> coupon =
                boost::dynamic_pointer_cast<FixedRateCoupon>(
                    arguments_.premiumLeg[i]);
            Date d = arguments_.premiumLeg[i]->date();
            if (d > discountCurve_->referenceDate()) {
                std::vector<Probability> probsTriggering =
                    arguments_.basket->probsBeingNthEvent(arguments_.ntdOrder, 
                        d);
                Probability defaultProb = 
                    std::accumulate(probsTriggering.begin(), 
                    probsTriggering.end(), 0.);

                results_.premiumValue += arguments_.premiumLeg[i]->amount()
                    * discountCurve_->discount(d)
                    * (1.0 - defaultProb);

                if (coupon->accrualStartDate() >= 
                    discountCurve_->referenceDate())
                    d = coupon->accrualStartDate();
                else
                    d = discountCurve_->referenceDate();

                // do steps of specified size
                d0 = d;
                Period stepSize = integrationStepSize_;
                probsTriggering =
                    arguments_.basket->probsBeingNthEvent(arguments_.ntdOrder, 
                        d0);
                Probability defProb0 = std::accumulate(probsTriggering.begin(), 
                    probsTriggering.end(), 0.);
                do {
                    DiscountFactor disc = discountCurve_->discount(d);
                    std::vector<Probability> probsTriggering1 =
                        arguments_.basket->probsBeingNthEvent(
                            arguments_.ntdOrder, d);
                    Probability defProb1 = 
                        std::accumulate(probsTriggering1.begin(), 
                            probsTriggering1.end(), 0.);
                    Probability dcfdd = defProb1 - defProb0;
                    defProb0 = defProb1;

                    if (arguments_.settlePremiumAccrual)
                        accrualValue += coupon->accruedAmount(d)*disc*dcfdd;
                    //Recoveries might differ along names, depending on which 
                    //name is triggering the contract the loss will be different
                    for(Size iName=0; iName<arguments_.basket->remainingSize(); 
                        iName++) {
                        claimValue -= (probsTriggering1[iName]-
                            probsTriggering[iName])
                            * arguments_.basket->claim()->amount(d, 
                                arguments_.notional, 
                                arguments_.basket->recoveryRate(d, iName))
                            * disc;
                    }

                    probsTriggering = probsTriggering1;

                    d0 = d;
                    d = d0 + stepSize;
                    // reduce step size ?
                    if (stepSize != 1*Days && d > coupon->accrualEndDate()) {
                        stepSize = 1*Days;
                        d = d0 + stepSize;
                    }
                }
                while (d <= coupon->accrualEndDate());
            }
        }

		// The upfront might be due before the curve ref date...
        if (!arguments_.premiumLeg[0]->hasOccurred(today))
            results_.upfrontPremiumValue =
                arguments_.basket->remainingNotional() 
                    * arguments_.upfrontRate
                    * discountCurve_->discount(
                        boost::dynamic_pointer_cast<FixedRateCoupon>(
                            arguments_.premiumLeg[0])->accrualStartDate());
        if (arguments_.side == Protection::Buyer) {
            results_.premiumValue *= -1;
            accrualValue *= -1;
            claimValue *= -1;
            results_.upfrontPremiumValue *= -1;
        }

        results_.value = results_.premiumValue + accrualValue + claimValue + 
            results_.upfrontPremiumValue;

        results_.fairPremium = -arguments_.premiumRate * claimValue 
            / (results_.premiumValue + accrualValue);
        // alternatively use results buffers and omit locals.
        results_.protectionValue = claimValue;

        results_.additionalResults["fairPremium"] = results_.fairPremium;
        results_.additionalResults["premiumLegNPV"] = 
            results_.premiumValue + results_.upfrontPremiumValue;
        results_.additionalResults["protectionLegNPV"] = 
            results_.protectionValue;
    }

}
