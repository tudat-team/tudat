/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      WARNING: There seems to be a bug in the computation of the L3 location!
 *
 */

#ifndef TUDAT_LIBRATION_POINT_FUNCTIONS_H
#define TUDAT_LIBRATION_POINT_FUNCTIONS_H

#include <cmath>

#include <memory>

#include <Eigen/Core>

#include "tudat/math/basic/basicFunction.h"

namespace tudat
{

namespace circular_restricted_three_body_problem
{

//! Simple definition of a test function, so that it can be used by all root-finder unit tests.
struct LibrationPointLocationFunction
{
    //! Default destructor.
    virtual ~LibrationPointLocationFunction( ) { }

    //! Expected true location of the root.
    virtual double getTrueRootLocation( ) = 0;

    //! Accuracy of the true value of the root.
    virtual double getTrueRootAccuracy( ) = 0;

    //! Get a reasonable initial guess of the root location.
    virtual double getInitialGuess( ) = 0;

    //! Get a reasonable lower boundary for the root location.
    virtual double getLowerBound( ) = 0;

    //! Get a reasonable upper boundary for the root location.
    virtual double getUpperBound( ) = 0;
};


struct LibrationPointLocationFunction1 : public LibrationPointLocationFunction,
        public tudat::basic_mathematics::BasicFunction< double, double >
{
    //! Maximum order of the derivative before throwing an exception.
    unsigned int maximumDerivativeOrder;

    //! Create a function, where aMaximumDerivativeOrder is the maximum order of the derivative.
    LibrationPointLocationFunction1( unsigned int aMaximumDerivativeOrder, const double massParameter )
        : maximumDerivativeOrder( aMaximumDerivativeOrder ), massParameter_( massParameter )
    { }

    //! Mathematical test function.
    double evaluate( const double inputValue )
    {
        // Define Mathematical function: f(x) =
        return pow(inputValue, 5.0) - (3.0 - massParameter_) * pow(inputValue, 4.0)
               + (3.0 - 2.0 * massParameter_) * pow(inputValue, 3.0) - massParameter_ * pow(inputValue, 2.0)
               + 2.0 * massParameter_ * inputValue - massParameter_;
    }

    //! Derivatives of mathematical test function.
    double computeDerivative( const unsigned int order, const double inputValue )
    {
        // Sanity check.
        if ( order > maximumDerivativeOrder )
        {
            throw std::runtime_error( "The root-finder should not evaluate higher derivatives!" );
        }

        // Return the analytical expression for the derivatives.
        if ( order == 0 )
        {
            // Return the function value: y =
            return evaluate( inputValue );
        }

        else if ( order == 1 )
        {
            // Return the first derivative function value: y =
            return 5.0*pow(inputValue, 4.0) - 4.0*(3.0 - massParameter_)*pow(inputValue, 3.0)
                   + 3.0*(3.0 - 2.0 * massParameter_) * pow(inputValue, 2.0) - 2.0*massParameter_ * inputValue
                   + 2.0 * massParameter_;
        }

        else if ( order == 2 )
        {
            // Return the second derivative function value: y = .
            return 20.0*pow(inputValue, 3.0) - 12.0*(3.0 - massParameter_)*pow(inputValue, 2.0)
                   + 6.0*(3.0 - 2.0 * massParameter_) * inputValue - 2.0*massParameter_;
        }

        else
        {
            throw std::runtime_error(
                        "An error occured when evaluating the order of the derivative." );
        }
    }

    //! Crash on integration as root_finders should not execute these.
    double computeDefiniteIntegral( const unsigned int order, const double lowerBound,
                                    const double upperbound )
    {
        throw std::runtime_error( "The root-finder should not evaluate integrals!" );
    }

    //! Get the expected true location of the root.
    /*!
     * Returns the expected true location of the function root, here \f$1\f$.
     *
     * \return True location of the root.
     */
    double getTrueRootLocation( )
    {
        return 1.0;
    }

    //! Get the accuracy of the true location of the root.
    /*!
     * Returns the accuracy of the true location of the function root, here 1e-15.
     *
     * \return Accuracy of the true location of the root.
     */
    double getTrueRootAccuracy( ) { return 1.0e-15; }

    //! Get a reasonable initial guess of the root location.
    /*!
     * Returns a reasonable initial guess for the true location of the function root, here 4.
     *
     * \return Initial guess for the true location of the function root.
     */
    double getInitialGuess( ) { return 1.0; }

    //! Get a reasonable lower boundary for the root location.
    /*!
     * Returns a reasonable lower bound for the true location of the function root, here -1.
     *
     * \return Lower bound for the true location of the function root.
     */
    double getLowerBound( ) { return 0.0; }

    //! Get a reasonable upper boundary for the root location.
    /*!
     * Returns a reasonable upper bound for the true location of the function root, here 4.
     *
     * \return Upper bound for the true location of the function root.
     */
    double getUpperBound( ) { return 1.0; }

protected:

    double massParameter_;
private:
};

struct LibrationPointLocationFunction2 : public LibrationPointLocationFunction,
        public tudat::basic_mathematics::BasicFunction< double, double >
{
    //! Maximum order of the derivative before throwing an exception.
    unsigned int maximumDerivativeOrder;

    //! Create a function, where aMaximumDerivativeOrder is the maximum order of the derivative.
    LibrationPointLocationFunction2( unsigned int aMaximumDerivativeOrder, const double massParameter )
        : maximumDerivativeOrder( aMaximumDerivativeOrder ), massParameter_( massParameter )
    { }

    //! Mathematical test function.
    double evaluate( const double inputValue )
    {
        // Define Mathematical function: f(x) =
        return pow(inputValue, 5.0) + (3.0 - massParameter_) * pow(inputValue, 4.0)
               + (3.0 - 2.0 * massParameter_) * pow(inputValue, 3.0) - massParameter_ * pow(inputValue, 2.0)
               - 2.0 * massParameter_ * inputValue - massParameter_;
    }

    //! Derivatives of mathematical test function.
    double computeDerivative( const unsigned int order, const double inputValue )
    {
        // Sanity check.
        if ( order > maximumDerivativeOrder )
        {
            throw std::runtime_error( "The root-finder should not evaluate higher derivatives!" );
        }

        // Return the analytical expression for the derivatives.
        if ( order == 0 )
        {
            // Return the function value: y =
            return evaluate( inputValue );
        }

        else if ( order == 1 )
        {
            // Return the first derivative function value: y =
            return 5.0*pow(inputValue, 4.0) + 4.0*(3.0 - massParameter_)*pow(inputValue, 3.0)
                   + 3.0*(3.0 - 2.0 * massParameter_) * pow(inputValue, 2.0) - 2.0*massParameter_ * inputValue
                   - 2.0 * massParameter_;
        }

        else if ( order == 2 )
        {
            // Return the second derivative function value: y = .
            return 20.0*pow(inputValue, 3.0) + 12.0*(3.0 - massParameter_)*pow(inputValue, 2.0)
                   + 6.0*(3.0 - 2.0 * massParameter_) * inputValue - 2.0*massParameter_;
        }

        else
        {
            throw std::runtime_error(
                        "An error occured when evaluating the order of the derivative." );
        }
    }

    //! Crash on integration as root_finders should not execute these.
    double computeDefiniteIntegral( const unsigned int order, const double lowerBound,
                                    const double upperbound )
    {
        throw std::runtime_error( "The root-finder should not evaluate integrals!" );
    }

    //! Get the expected true location of the root.
    /*!
     * Returns the expected true location of the function root, here \f$1\f$.
     *
     * \return True location of the root.
     */
    double getTrueRootLocation( )
    {
        return 1.0;
    }

    //! Get the accuracy of the true location of the root.
    /*!
     * Returns the accuracy of the true location of the function root, here 1e-15.
     *
     * \return Accuracy of the true location of the root.
     */
    double getTrueRootAccuracy( ) { return 1.0e-15; }

    //! Get a reasonable initial guess of the root location.
    /*!
     * Returns a reasonable initial guess for the true location of the function root, here 4.
     *
     * \return Initial guess for the true location of the function root.
     */
    double getInitialGuess( ) { return 1.0; }

    //! Get a reasonable lower boundary for the root location.
    /*!
     * Returns a reasonable lower bound for the true location of the function root, here -1.
     *
     * \return Lower bound for the true location of the function root.
     */
    double getLowerBound( ) { return 0.0; }

    //! Get a reasonable upper boundary for the root location.
    /*!
     * Returns a reasonable upper bound for the true location of the function root, here 4.
     *
     * \return Upper bound for the true location of the function root.
     */
    double getUpperBound( ) { return 1.0; }

protected:

    double massParameter_;

private:
};


} // namespace circular_restricted_three_body_problem
} // namespace tudat

#endif //
