/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*!
Copyright (C) 2002, 2003 Sadruddin Rejeb
Copyright (C) 2004 Ferdinando Ametrano
Copyright (C) 2005, 2006, 2007 StatPro Italia srl
Copyright (C) 2015 Riccardo Barone

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

#include <ql/quantlib.hpp>

// This makes it easier to use array literals (alas, no std::vector literals)
#define LENGTH(a) (sizeof(a)/sizeof(a[0]))

#ifdef BOOST_MSVC
/* Uncomment the following lines to unmask floating-point
exceptions. Warning: unpredictable results can arise...

See http://www.wilmott.com/messageview.cfm?catid=10&threadid=9481
Is there anyone with a definitive word about this?
*/
// #include <float.h>
// namespace { unsigned int u = _controlfp(_EM_INEXACT, _MCW_EM); }
#endif

#include <boost/timer.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ql/utilities/dataformatters.hpp>

using namespace QuantLib;

#if defined(QL_ENABLE_SESSIONS)
namespace QuantLib {

    Integer sessionId() { return 0; }

}
#endif
//Define more comfortable data structures 
//for different interest rate instruments
struct Datum {
    Integer settlementDays;
    Period maturity;
    Rate rate;
};

struct FraDatum {
    Integer settlementDays;
    Integer nStart;
    Integer nMaturity;
    Rate rate;
};

struct SwapDatum {
    Integer settlementDays;
    Period maturity;
    Rate rate;
};
//Define three different selection strategies for calibration
enum SelectionMethod { total, coterminal, swapTenor };
//Method to select the proper instruments according to the selection strategy
std::vector<Size> selection(const std::vector<Integer> optExpiry,
    const std::vector<Integer> swLength, const SelectionMethod select,
    const Integer c){
    std::vector<Size> ct;
    switch (select)
    {
    case total:
        for (Size j = 0; j < swLength.size(); j++)
        {
            for (Size i = 0; i < optExpiry.size(); i++)
            {
                ct.push_back(j*optExpiry.size() + i);
            }
        }
        break;
    case coterminal:
        QL_REQUIRE(c <= optExpiry.back() + swLength.back(), 
                   "coterminal " <<c<< "Y greater than max quoted coterminal " 
                    << optExpiry.back() + swLength.back() << "Y");
        for (Size j = 0; j < swLength.size(); j++)
        {
            for (Size i = 0; i < optExpiry.size(); i++)
            {
                if (optExpiry[i] + swLength[j] == c)
                {
                    ct.push_back(j*optExpiry.size() + i);
                }
            }
        }
        break;
    case swapTenor:
        QL_REQUIRE(c <= swLength.back(), 
                   "Underlying swap length " << c 
                   << "Y greater than max quoted underlying swap length " 
                   <<swLength.back() << "Y");
        for (Size j = 0; j < swLength.size(); j++)
        {
            if (swLength[j] == c)
            {
                for (Size i = 0; i < optExpiry.size(); i++)
                {
                    ct.push_back(j*optExpiry.size() + i);
                }
            }
        }
        break;
    default:
        QL_FAIL("Unknown selection type");
        break;
    }
    
    
    return ct;
}

/*******************************
** SWAPTION VOLATILITY MATRIX **
*******************************/
//Defining the relevant periond in the Swaption Volatility matrix
std::vector<Integer> optionExpiry = {
    1, 2, 3, 4, 5, 7, 10, 15, 20, 25, 30 };

std::vector<Integer> swapLenghts = {
     2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30 };

//Defining the Swaption volatility matrix (swaptionLength X optionExpiry)
std::vector<Volatility> swaptionVols = {
   1.606, 1.144, 0.869, 0.692, 0.568, 0.424, 0.364, 0.356, 0.426, 0.470, 0.438,
   1.107, 0.904, 0.735, 0.605, 0.517, 0.411, 0.361, 0.360, 0.433, 0.461, 0.429,
   0.917, 0.763, 0.639, 0.548, 0.481, 0.399, 0.358, 0.365, 0.437, 0.452, 0.420,
   0.798, 0.673, 0.577, 0.509, 0.459, 0.391, 0.357, 0.368, 0.437, 0.442, 0.411,
   0.717, 0.614, 0.536, 0.483, 0.441, 0.384, 0.357, 0.373, 0.441, 0.437, 0.409,
   0.664, 0.578, 0.509, 0.464, 0.428, 0.379, 0.357, 0.377, 0.441, 0.432, 0.407,
   0.621, 0.548, 0.489, 0.450, 0.418, 0.376, 0.357, 0.383, 0.441, 0.427, 0.404,
   0.593, 0.526, 0.474, 0.439, 0.411, 0.374, 0.359, 0.387, 0.438, 0.421, 0.402,
   0.573, 0.511, 0.464, 0.432, 0.406, 0.372, 0.360, 0.390, 0.437, 0.417, 0.403,
   0.487, 0.444, 0.414, 0.389, 0.371, 0.350, 0.349, 0.372, 0.394, 0.379, 0.381,
   0.459, 0.425, 0.400, 0.383, 0.369, 0.350, 0.346, 0.355, 0.368, 0.363, 0.378,
   0.459, 0.429, 0.405, 0.387, 0.373, 0.353, 0.346, 0.353, 0.369, 0.374, 0.381,
   0.460, 0.431, 0.407, 0.388, 0.372, 0.353, 0.344, 0.355, 0.373, 0.372, 0.377
};

/*
Method to calibrate a given model with a customizable instrument
selection. It displays the relative errors of the whole matrix and stores other
useful data in .txt files
*/
void calibrateModel(
    const boost::shared_ptr<ShortRateModel>& model,
    const std::vector<boost::shared_ptr<CalibrationHelper> >& selectedHelpers,
    const std::vector<boost::shared_ptr<CalibrationHelper> >& allHelpers) {
    std::cout << "Starting calibration on "
              << selectedHelpers.size() << " elements" << std::endl;
    std::ofstream modelPrice;
    modelPrice.open("modelPrice.txt",std::fstream::app);
    std::ofstream absErr;
    absErr.open("absErr.txt");
    std::ofstream relErr;
    relErr.open("relErr.txt");
    boost::timer t;

    //Calibration with Levenberg-Marquardt method

    LevenbergMarquardt om;
    model->calibrate(selectedHelpers, om,
        EndCriteria(10000, 500, 1.0e-8, 1.0e-8, 1.0e-8));

    //Output the calibration time
    double seconds = t.elapsed();
    Integer hours = int(seconds / 3600);
    seconds -= hours * 3600;
    Integer minutes = int(seconds / 60);
    seconds -= minutes * 60;
    std::cout << "Calibration completed in ";
    if (hours > 0)
        std::cout << hours << " h ";
    if (hours > 0 || minutes > 0)
        std::cout << minutes << " m ";
    std::cout << std::fixed << std::setprecision(0)
        << seconds << " s\n" << std::endl;
    
    // Output the repricing relative errors
    std::cout << "Calibration (relative) errors:\n" << std::endl;
    modelPrice << "\n--------------------------------------------\n";
    for (Size j = 0; j < swapLenghts.size(); j++)
    {
        //std::cout << "|";
        for (Size i = 0; i < optionExpiry.size(); i++)
        {
            Size k = j*optionExpiry.size() + i;
            Real g2Price = allHelpers[k]->modelValue();
            //Volatility implied = allHelpers[k]->impliedVolatility(npv, 1e-12,
            //    10000, swaptionVols[k] - 0.3, swaptionVols[k] + 0.3);
            
            const Real lowerPrice = allHelpers[k]->blackPrice(0.001);
            const Real upperPrice = allHelpers[k]->blackPrice(10);
            Volatility implied;
            if (g2Price <= lowerPrice)
                implied = 0.001;
            else
                if (g2Price >= upperPrice)
                    implied = 10.0;
                else
                    implied = allHelpers[k]->impliedVolatility(
                    g2Price, 1e-12, 5000, 0.001, 10);

            Volatility err = (implied - swaptionVols[k]);
            Real marketPrice = allHelpers[k]->blackPrice(swaptionVols[k]);
            //Real err = g2Price - marketPrice;
            std::cout << " " << std::fixed << std::setprecision(5)
                << err / swaptionVols[k] << " ";
            absErr << std::fixed << std::setprecision(4)
                   << g2Price - marketPrice << ", ";
            relErr << std::fixed << std::setprecision(4)
                   << (g2Price - marketPrice) / marketPrice << ", ";
            modelPrice << std::fixed << std::setprecision(4)
                      << g2Price << ", ";
        }
        //std::cout << " |\n";
        absErr << "\n";
        relErr << "\n";
        modelPrice << "\n";
        std::cout << std::endl;
    }
    modelPrice.close();
    absErr.close();
    relErr.close();
}


int main(int, char*[]) {

    try {

        boost::timer timer;
        std::cout << std::endl;

        /******************
        *** MARKET DATA ***
        ******************/

        Calendar calendar = TARGET();
        Date today(16, September, 2015);
        Integer fixingDays = 2;
        Date evaluationDate = calendar.adjust(today);

        Settings::instance().evaluationDate() = evaluationDate;

        Date settlementDate = calendar.advance(evaluationDate, 
                                               fixingDays, 
                                               Days);

        std::cout << "Today: " << evaluationDate.weekday()
            << ", " << evaluationDate << std::endl;
        //
        std::cout << "Evaluation date: " << evaluationDate.weekday()
            << ", " << evaluationDate << std::endl;
        //
        std::cout << "Settlement date: " << settlementDate.weekday()
            << ", " << settlementDate << std::endl;

        Integer tenor = 6;

        // DEPOs
        Datum depoDataON[] = {
            //{settlementDays,Maturity, Rate}
            { 0, 1 * Days, -0.001300 },
            { 1, 1 * Days, -0.001300 },
            { 2, 1 * Days, -0.001301 }
        };
        Datum depoDataSixMonth[] = {
            //{settlementDays,Maturity, Rate}
            { 1, 1 * Days, 0.000353 }, { 2, 1 * Weeks, 0.000353 },
            { 2, 2 * Weeks, 0.000353 }, { 2, 3 * Weeks, 0.000353 },
            { 2, 1 * Months, 0.000353 }, { 2, 2 * Months, 0.000367 },
            { 2, 3 * Months, 0.000372 }, { 2, 4 * Months, 0.000366 },
            { 2, 5 * Months, 0.000376 }, { 2, 6 * Months, 0.000376 }
        };

        // FRAs
        FraDatum fraDataSixMonth[] = {
            //{settlementDays,Start,End,Rate}
            { 2, 1, 7, 0.00037 }, { 2, 2, 8, 0.00039 },
            { 2, 3, 9, 0.00042 }, { 2, 4, 10, 0.00043 },
            { 2, 5, 11, 0.00045 }, { 2, 6, 12, 0.00046 },
            { 2, 7, 13, 0.0005 }, { 2, 8, 14, 0.00053 },
            { 2, 9, 15, 0.00056 }, { 2, 10, 16, 0.00063 },
            { 2, 11, 17, 0.00068 }, { 2, 12, 18, 0.00074 },
            { 2, 13, 19, 0.00084 }, { 2, 14, 20, 0.00093 },
            { 2, 15, 21, 0.00103 }, { 2, 16, 22, 0.00115 },
            { 2, 17, 23, 0.00127 }, { 2, 18, 24, 0.00139 }

        };

        // SWAPs
        SwapDatum swapDataSixMonth[] = {
            //{settlementDays, Maturity, Rate}
            { 2, 3 * Years, 0.00147 }, { 2, 4 * Years, 0.00247 },
            { 2, 5 * Years, 0.00365 }, { 2, 6 * Years, 0.00495 },
            { 2, 7 * Years, 0.00629 }, { 2, 8 * Years, 0.00761 },
            { 2, 9 * Years, 0.00881 }, { 2, 10 * Years, 0.00986 },
            { 2, 11 * Years, 0.01081 }, { 2, 12 * Years, 0.01165 },
            { 2, 13 * Years, 0.01238 }, { 2, 14 * Years, 0.01300 },
            { 2, 15 * Years, 0.01352 }, { 2, 16 * Years, 0.01397 },
            { 2, 17 * Years, 0.01435 }, { 2, 18 * Years, 0.01467 },
            { 2, 19 * Years, 0.01492 }, { 2, 20 * Years, 0.01512 },
            { 2, 21 * Years, 0.01527 }, { 2, 22 * Years, 0.01538 },
            { 2, 23 * Years, 0.01545 }, { 2, 24 * Years, 0.01550 },
            { 2, 25 * Years, 0.01554 }, { 2, 26 * Years, 0.01557 },
            { 2, 27 * Years, 0.01560 }, { 2, 28 * Years, 0.01563 },
            { 2, 29 * Years, 0.01566 }, { 2, 30 * Years, 0.01568 },
            { 2, 31 * Years, 0.01570 }, { 2, 32 * Years, 0.01573 },
            { 2, 33 * Years, 0.01575 }, { 2, 34 * Years, 0.01577 },
            { 2, 35 * Years, 0.01578 }, { 2, 36 * Years, 0.01579 },
            { 2, 37 * Years, 0.01579 }, { 2, 38 * Years, 0.015796 },
            { 2, 39 * Years, 0.01580 }, { 2, 40 * Years, 0.01580 },
            { 2, 50 * Years, 0.01539 }, { 2, 60 * Years, 0.01524 }
        };

        SwapDatum swapDataON[] = {
            { 2, 1 * Weeks, -0.00131 },
            { 2, 2 * Weeks, -0.00132 }, { 2, 3 * Weeks, -0.00131 },
            { 2, 1 * Months, -0.00132 }, { 2, 2 * Months, -0.00133 },
            { 2, 3 * Months, -0.00135 }, { 2, 4 * Months, -0.00137 },
            { 2, 5 * Months, -0.00139 }, { 2, 6 * Months, -0.0014 },
            { 2, 7 * Months, -0.00142 }, { 2, 8 * Months, -0.00143 },
            { 2, 9 * Months, -0.00144 }, { 2, 10 * Months, -0.00146 },
            { 2, 11 * Months, -0.00147 }, { 2, 1 * Years, -0.00148 },
            { 2, 13 * Months, -0.00149 }, { 2, 14 * Months, -0.0015 },
            { 2, 15 * Months, -0.00151 }, { 2, 16 * Months, -0.00151 },
            { 2, 17 * Months, -0.00151 }, { 2, 18 * Months, -0.00151 },
            { 2, 19 * Months, -0.0015 }, { 2, 20 * Months, -0.00149 },
            { 2, 21 * Months, -0.00147 }, { 2, 22 * Months, -0.00145 },
            { 2, 23 * Months, -0.00143 }, { 2, 2 * Years, -0.00141 },
            { 2, 3 * Years, -0.00092 }, { 2, 4 * Years, -0.00008 },
            { 2, 5 * Years, 0.00099 }, { 2, 6 * Years, 0.00228 },
            { 2, 7 * Years, 0.00362 }, { 2, 8 * Years, 0.00496 },
            { 2, 9 * Years, 0.00619 }, { 2, 10 * Years, 0.0073 },
            { 2, 11 * Years, 0.00831 }, { 2, 12 * Years, 0.00918 },
            { 2, 13 * Years, 0.00995 }, { 2, 14 * Years, 0.01062 },
            { 2, 15 * Years, 0.01121 }, { 2, 16 * Years, 0.01172 },
            { 2, 17 * Years, 0.01216 }, { 2, 18 * Years, 0.01252 },
            { 2, 19 * Years, 0.01283 }, { 2, 20 * Years, 0.01307 },
            { 2, 21 * Years, 0.01326 }, { 2, 22 * Years, 0.01341 },
            { 2, 23 * Years, 0.01353 }, { 2, 24 * Years, 0.01362 },
            { 2, 25 * Years, 0.0137 }, { 2, 26 * Years, 0.01377 },
            { 2, 27 * Years, 0.01383 }, { 2, 28 * Years, 0.01389 },
            { 2, 29 * Years, 0.01395 }, { 2, 30 * Years, 0.014 },
            { 2, 35 * Years, 0.01423 }, { 2, 40 * Years, 0.01433 },
            { 2, 50 * Years, 0.01407 }, { 2, 60 * Years, 0.01403 }
        };

        /*********************
        ***  BUILD EUR ON  ***
        *********************/
        std::cout << std::endl;
        std::cout << " ---------- " << std::endl;
        std::cout << "|  EUR ON  |" << std::endl;
        std::cout << " ---------- " << std::endl;

        /*** QUOTES ***/

        // DEPOs
        std::vector<boost::shared_ptr<Quote>> depoONQuote;
        for (Size i = 0; i < LENGTH(depoDataON); i++)
        {
            boost::shared_ptr<Quote>quote(new SimpleQuote(depoDataON[i].rate));
            depoONQuote.push_back(quote);
        }
        // SWAPs
        std::vector<boost::shared_ptr<Quote>> swapONQuote;
        for (Size i = 0; i < LENGTH(swapDataON); i++)
        {
            boost::shared_ptr<Quote>quote(new SimpleQuote(swapDataON[i].rate));
            swapONQuote.push_back(quote);
        }

        /*** RATE HELPERS ***/

        DayCounter depoOnDayCounter = Actual360();
        BusinessDayConvention depoOnDayConvention = Following;
        //DEPOs
        std::vector<boost::shared_ptr<RateHelper>> rateHelperON;

        for (Size i = 0; i < LENGTH(depoDataON); i++)
        {
            boost::shared_ptr<RateHelper> instrON(new
                DepositRateHelper(Handle<Quote>(depoONQuote.at(i)),
                depoDataON[i].maturity,
                depoDataON[i].settlementDays,
                calendar,
                depoOnDayConvention,
                true,
                depoOnDayCounter));
            rateHelperON.push_back(instrON);
        }
        //SWAPs
        boost::shared_ptr<Eonia> eonia(new Eonia);

        for (Size i = 0; i < LENGTH(swapDataON); i++)
        {
            boost::shared_ptr<RateHelper> instrON(new
                OISRateHelper(swapDataON[i].settlementDays,
                swapDataON[i].maturity,
                Handle<Quote>(swapONQuote.at(i)),
                eonia));
            rateHelperON.push_back(instrON);
        }

        /*** CURVE BUILDING ***/

        DayCounter termStructureDayCounter = Actual365Fixed();

        double tolerance = 1.0e-15;

        Handle<YieldTermStructure> EURONTS(
            boost::shared_ptr<YieldTermStructure>(
            new PiecewiseYieldCurve<Discount, LogCubic>(
            evaluationDate, rateHelperON,
            termStructureDayCounter,
            tolerance,
            LogCubic(
            CubicInterpolation::Spline, true,
            CubicInterpolation::SecondDerivative, 0.0,
            CubicInterpolation::SecondDerivative, 0.0))));

        EURONTS->enableExtrapolation();

        //Discount at start date
        std::cout << "Discount at " << evaluationDate;
        std::cout << " : " << std::fixed
            << EURONTS->discount(evaluationDate, false) << std::endl;
        //Discount at max date
        std::cout << "Discount at " << EURONTS->maxDate();
        std::cout << " : " << std::fixed
            << EURONTS->discount(EURONTS->maxDate(), false) << std::endl;

        /*********************
        ***  BUILD EUR 6M  ***
        *********************/
        std::cout << std::endl;
        std::cout << " ---------- " << std::endl;
        std::cout << "|  EUR 6M  |" << std::endl;
        std::cout << " ---------- " << std::endl;
        boost::shared_ptr<IborIndex> indexON(new
            Euribor6M(EURONTS));

        /*** QUOTES ***/

        // DEPOs
        std::vector<boost::shared_ptr<Quote>> depoSixMonthQuote;
        for (Size i = 0; i < LENGTH(depoDataSixMonth); i++)
        {
            boost::shared_ptr<Quote> quote(
                new SimpleQuote(depoDataSixMonth[i].rate));
            depoSixMonthQuote.push_back(quote);
        }
        // FRAs
        std::vector<boost::shared_ptr<Quote>> fraSixMonthQuote;
        for (Size i = 0; i < LENGTH(fraDataSixMonth); i++)
        {
            boost::shared_ptr<Quote> quote(
                new SimpleQuote(fraDataSixMonth[i].rate));
            fraSixMonthQuote.push_back(quote);
        }
        // SWAPs
        std::vector<boost::shared_ptr<Quote>> swapSixMonthQuote;
        for (Size i = 0; i < LENGTH(swapDataSixMonth); i++)
        {
            boost::shared_ptr<Quote> quote(
                new SimpleQuote(swapDataSixMonth[i].rate));
            swapSixMonthQuote.push_back(quote);
        }


        /*** RATE HELPERS ***/

        DayCounter depo6mDayCounter = Actual360();
        BusinessDayConvention depo6mDayConvention = Following;
        //DEPOs
        std::vector<boost::shared_ptr<RateHelper>> rateHelper6M;

        for (Size i = 0; i < LENGTH(depoDataSixMonth); i++)
        {
            boost::shared_ptr<RateHelper> instr6M(new
                DepositRateHelper(Handle<Quote>(depoSixMonthQuote.at(i)),
                depoDataSixMonth[i].maturity,
                fixingDays,
                calendar,
                depo6mDayConvention,
                true,
                depo6mDayCounter));
            rateHelper6M.push_back(instr6M);
        }
        //FRAs
        BusinessDayConvention fra6mDayConvention = ModifiedFollowing;

        for (Size i = 0; i < LENGTH(fraDataSixMonth); i++)
        {
            boost::shared_ptr<RateHelper> instr6M(new
                FraRateHelper(Handle<Quote>(fraSixMonthQuote.at(i)),
                fraDataSixMonth[i].nStart,
                fraDataSixMonth[i].nMaturity,
                fixingDays,
                calendar,
                fra6mDayConvention,
                true,
                depo6mDayCounter));
            rateHelper6M.push_back(instr6M);
        }
        //SWAPs
        Frequency sw6mFixedLegFreq = Annual;
        BusinessDayConvention sw6mFixedLegConv = ModifiedFollowing;
        DayCounter sw6mFixedLegDayCounter = Thirty360(Thirty360::European);
        boost::shared_ptr<IborIndex> euribor6M(new Euribor6M);


        for (Size i = 0; i < LENGTH(swapDataSixMonth); i++)
        {
            boost::shared_ptr<RateHelper> instr6M(new
                SwapRateHelper(Handle<Quote>(swapSixMonthQuote.at(i)),
                swapDataSixMonth[i].maturity,
                calendar,
                sw6mFixedLegFreq,
                sw6mFixedLegConv,
                sw6mFixedLegDayCounter,
                euribor6M));
            rateHelper6M.push_back(instr6M);
        }


        /*** CURVE BUILDING ***/

        Handle<YieldTermStructure> EURSixMonthTS(
            boost::shared_ptr<YieldTermStructure>(
            new PiecewiseYieldCurve<Discount, LogCubic>(
            settlementDate, rateHelper6M,
            termStructureDayCounter,
            tolerance,
            LogCubic(
            CubicInterpolation::Spline, true,
            CubicInterpolation::SecondDerivative, 0.0,
            CubicInterpolation::SecondDerivative, 0.0))));

        EURSixMonthTS->enableExtrapolation();

        //Discount at start date
        std::cout << "Discount at " << evaluationDate;
        std::cout << " : " << std::fixed
            << EURSixMonthTS->discount(settlementDate, false)
            << std::endl;
        //Discount at max date
        std::cout << "Discount at " << EURONTS->maxDate();
        std::cout << " : " << std::fixed
            << EURSixMonthTS->discount(EURONTS->maxDate(), false)
            << std::endl;

        /*******************
        *** ATM SWAPTION ***
        *******************/
        std::cout << std::endl;
        std::cout << " -------------- " << std::endl;
        std::cout << "| ATM SWAPTION |" << std::endl;
        std::cout << " -------------- " << std::endl;

        //Forward Swap pricing
        std::cout << "Forward Swap pricing" << std::endl;

        Frequency fixedLegFreq = Annual;
        Frequency floatingLegFreq = Semiannual;

        BusinessDayConvention fixedLegConv = Unadjusted;
        BusinessDayConvention floatingLegConv = ModifiedFollowing;

        DayCounter fixedLegDayCounter = Thirty360(Thirty360::European);

        VanillaSwap::Type type = VanillaSwap::Payer;
        Rate dummyFixedRate = 0.03;
        boost::shared_ptr<IborIndex> indexSixMonths(new
            Euribor6M(EURSixMonthTS));

        Date startDate = calendar.advance(settlementDate, 1, Years,
            floatingLegConv);
        Date maturity = calendar.advance(startDate, 5, Years, floatingLegConv);

        Schedule fixedSchedule(startDate,
            maturity,
            Period(fixedLegFreq),
            calendar,
            fixedLegConv,
            fixedLegConv,
            DateGeneration::Forward,
            false);

        Schedule floatSchedule(startDate,
            maturity,
            Period(floatingLegFreq),
            calendar,
            floatingLegConv,
            floatingLegConv,
            DateGeneration::Forward,
            false);


        boost::shared_ptr<VanillaSwap> swap(new VanillaSwap(
            type,
            1000.0,
            fixedSchedule,
            dummyFixedRate,
            fixedLegDayCounter,
            floatSchedule,
            indexSixMonths,
            0.0,
            indexSixMonths->dayCounter()));

        boost::shared_ptr<PricingEngine> swapEngine(
            new DiscountingSwapEngine(EURONTS));
        swap->setPricingEngine(swapEngine);

        Rate fixedATMRate = swap->fairRate();

         std::cout << "1Y Fwd Swap5Y Fair rate: "
             << io::rate(fixedATMRate) << std::endl;

        boost::shared_ptr<VanillaSwap> atmSwap(new VanillaSwap(
            type,
            1000.0,
            fixedSchedule,
            fixedATMRate,
            fixedLegDayCounter,
            floatSchedule,
            indexSixMonths,
            0.0,
            indexSixMonths->dayCounter()));

        std::vector<boost::shared_ptr<CalibrationHelper> > swaptions;
        //Define a time grid (used for pricing with trees)
        //std::list<Time> times;
        /*** CALIBRATION HELPERS ***/
        //Create calibration helpers
        for (Size i = 0; i < swaptionVols.size(); i++)
        {
            boost::shared_ptr<Quote> vol(new SimpleQuote(swaptionVols.at(i)));
            swaptions.push_back(boost::shared_ptr<CalibrationHelper>(new
                SwaptionHelper(
                Period(optionExpiry.at(i%optionExpiry.size()), Years),
                Period(swapLenghts.at(i / optionExpiry.size()), Years),
                Handle<Quote>(vol),
                indexSixMonths,
                indexSixMonths->tenor(),
                indexSixMonths->dayCounter(),
                indexSixMonths->dayCounter(),
                EURSixMonthTS,
                CalibrationHelper::CalibrationErrorType::ImpliedVolError)));
           // swaptions.back()->addTimesTo(times);
        }

        // Building time-grid
        //TimeGrid grid(times.begin(), times.end(), 50);
        
        //Print default Model parameters
        boost::shared_ptr<G2> modelG2(new G2(EURSixMonthTS,
                                              0.62981346,
                                              0.02399075,
                                              0.03582301,
                                              0.01115666,
                                              -0.97563256));
        std::cout << "Default G2++ parameters:\n"
            << "a     = " << modelG2->params()[0] << ", "
            << "sigma = " << modelG2->params()[1] << "\n"
            << "b     = " << modelG2->params()[2] << ", "
            << "eta   = " << modelG2->params()[3] << "\n"
            << "rho   = " << modelG2->params()[4]
            << std::endl << std::endl;

        /*** Model calibrations ***/
        std::cout << "Model Calibration\n" << std::endl;
        
        for (Size i = 0; i < swaptions.size(); i++)
        {
            swaptions.at(i)->setPricingEngine(boost::shared_ptr<PricingEngine>(
                new G2SwaptionEngine(modelG2, 6.0, 60)));
        }
        for (SelectionMethod method = total ; method <= swapTenor; 
             method=static_cast<SelectionMethod>(method+1))
        {
            std::vector<Size> c = selection(optionExpiry, 
                                            swapLenghts, 
                                            method, 
                                            10);
            std::vector<boost::shared_ptr<CalibrationHelper> >selectedSwaption;

            for (Size i = 0; i < c.size(); i++)
            {
                selectedSwaption.push_back(swaptions[c[i]]);
            }

            if (method == 1)
            {
                std::vector<Size> d = selection(optionExpiry,
                    swapLenghts,
                    method,
                    40);
                for (Size i = 0; i < c.size(); i++)
                {
                    selectedSwaption.push_back(swaptions[d[i]]);
                }
            }
            std::cout << "Calibration method = " << method + 1 << std::endl;
            calibrateModel(modelG2, selectedSwaption, swaptions);

            std::cout << "\nCalibrated G2++ paramenters:\n"
                << "a     = "<<std::setprecision(8)<<modelG2->params()[0]<<", "
                << "sigma = "<<std::setprecision(8)<<modelG2->params()[1]<<"\n"
                << "b     = "<<std::setprecision(8)<<modelG2->params()[2]<<", "
                << "eta   = "<<std::setprecision(8)<<modelG2->params()[3]<<"\n"
                << "rho   = "<<std::setprecision(8)<<modelG2->params()[4]
                << std::endl <<std::endl;
        }
        

        double seconds = timer.elapsed();
        Integer hours = int(seconds / 3600);
        seconds -= hours * 3600;
        Integer minutes = int(seconds / 60);
        seconds -= minutes * 60;
        std::cout << " \nRun completed in ";
        if (hours > 0)
            std::cout << hours << " h ";
        if (hours > 0 || minutes > 0)
            std::cout << minutes << " m ";
        std::cout << std::fixed << std::setprecision(0)
            << seconds << " s\n" << std::endl;

        return 0;
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "unknown error" << std::endl;
        return 1;
    }
}