/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/aerodynamics/rarefiedFlowInteractionModel.h"

#include <Eigen/Core>

#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/math/basic/mathematicalConstants.h"

#include <array>

namespace tudat
{
namespace aerodynamics
{
/* Computes aerodynamic force coefficients for a single panel
* 
*/
Eigen::Vector3d RarefiedFlowInteractionModel::computePanelForceCoefficientVector( 
    double CosineOfNormalDragAngle, //gammai in Doornbos
    double CosineOfNormalLiftAngle, //li in Doornbos
    double panelSurfaceArea,
    Eigen::Vector3d u_lift,
    Eigen::Vector3d u_drag,
    double Vinf,
    double T_atm,
    std::array number_densities,
    double total_number_density,
    double Aref)
{
    
    // To be moved to correct location (mathematicalConstants.h maybe?)
    std::array<double, 8> nrl_msise00_species_m = {4.002602, 15.999, 28.0134, 31.9988, 39.948, 1.008, 14.0067, 15.999}; // [g/mol] atomic mass of species in NRLMSISE-00 model

    // Initialize force coefficient vector
    Eigen::Vector3d forceCoefficientVector = Eigen::Vector3d::Zero();

    rhoO = number_densities.get(7); // atomic oxygen density

    if (currentCosineOfNormalDragAngle > 0){ // check if panel is pointing into the flow
        for (int j_species = 0; j_species < nrl_msise00_species_m.size(); j_species++) {
            Cdij = get_Cd_ij(Vinf, T_atm, species_m[j_species], CosineOfNormalDragAngle, panelSurfaceArea, Aref, rhoO) * number_densities[j_species] / total_number_density;
            Clij = get_Cl_ij(Vinf, T_atm, species_m[j_species], CosineOfNormalDragAngle, CosineOfNormalLiftAngle, panelSurfaceArea, Aref, rhoO) * number_densities[j_species] / total_number_density;

            panelForceCoefficientVector += Cdij * u_drag + Clij * u_lift;
        }
    }
    
    return panelForceCoefficientVector;
}

Eigen::Vector3d RarefiedFlowInteractionModel::computePanelMomentCoefficientVector( 
    Eigen::Vector3d panelForceCoefficientVector,
    Eigen::Vector3d panelPositionVector,
    double lref,
    )
{
    return panelMomentCoefficientVector = panelPositionVector.cross(panelForceCoefficientVector) / lref;
}


double RarefiedFlowInteractionModel::get_Cd_ij(double Vinf, double Tinf, double mj, double gammai, double Ai, double Aref, double rhoO)
    {
    // Function to calculate the drag coefficient for a single species j on a single panel i
    // Vinf: freestream velocity
    // Tinf: freestream temperature
    // mj: molecular mass of species j
    // gammai: cosine of the angle between the freestream velocity and the panel normal
    // Ai: area of panel i
    // Aref: reference area
    // rhoO: atomic oxygen density

    double Sj = get_Sj(Vinf, Tinf, mj);
    double Pij = get_Pij(Sj, gammai);
    double Gj = get_Gj(Sj);
    double Qj = 1.0 + Gj;
    double Zij = get_Zij(Sj, gammai);
    double alpha = get_alpha(rhoO, Tinf);
    double Vre_Vinf = get_Vre_Vinf(alpha, Tinf, Vinf);

    double Cd_ij = Pij / pow(pi, 0.5);
    Cd_ij += gammai * Qj * Zij;
    Cd_ij += gammai / 2.0 * Vre_Vinf * (gammai * pow(pi, 0.5) * Zij + Pij);
    Cd_ij *= Ai / Aref;

    return Cd_ij;
}

double RarefiedFlowInteractionModel::get_Cl_ij(double Vinf, double Tinf, double mj, double gammai, double li, double Ai, double Aref, double rhoO){
    // Function to calculate the drag coefficient for a single species j on a single panel i
    // Vinf: freestream velocity
    // Tinf: freestream temperature
    // mj: molecular mass of species j
    // gammai: cosine of the angle between the freestream velocity and the panel normal
    // li: sine of the angle between the freestream velocity and the panel normal
    // Ai: area of panel i
    // Aref: reference area
    // rhoO: atomic oxygen density

    double Sj = get_Sj(Vinf, Tinf, mj);
    double Pij = get_Pij(Sj, gammai);
    double Gj = get_Gj(Sj);
    double Zij = get_Zij(Sj, gammai);
    double alpha = get_alpha(rhoO, Tinf);
    double Vre_Vinf = get_Vre_Vinf(alpha, Tinf, Vinf);

    double Cl_ij = li * Gj * Zij;
    Cl_ij += li / 2.0 * Vre_Vinf * (gammai * pow(pi, 0.5) * Zij + Pij);
    Cl_ij *= Ai / Aref;

    return Cl_ij;
}

double RarefiedFlowInteractionModel::get_Sj(double Vinf, double Tinf, double mj){
    return Vinf / pow(2.0 * k * Tinf/mj, 0.5);
}

double RarefiedFlowInteractionModel::get_Gj(double Sj){
    return 1.0 / (2.0*pow(Sj, 2.0));;
}

double RarefiedFlowInteractionModel::get_Pij(double Sj, double gammai){
    return (1.0/Sj) * exp(-pow(gammai, 2.0) * pow(Sj, 2.0));
}

double RarefiedFlowInteractionModel::get_Zij(double Sj, double gammai){
    return 1.0 + erf(gammai * Sj);
}

double RarefiedFlowInteractionModel::get_alpha(double rhoO, double Tinf){
    // Function to calculate the energy accomodation coefficient, Miyata et. al. 2018
    // rhoO: atomic oxygen density
    // Tinf: freestream temperature
    // return (7.5e-17 * rhoO * Tinf) / (1 + 7.5e-17 * rhoO * Tinf);
    return 1.0;
}

double RarefiedFlowInteractionModel::get_Vre_Vinf(double alpha, double Ti, double Vinf){
    // Function to calculate the ratio of rebound velocity to freestream velocity
    // alpha: energy accomodation coefficient
    // Ti: temperature of panel i
    // Vinf: freestream velocity
    return pow(0.5 * ( 1.0 + alpha * ( (4.0*R*Ti)/pow(Vinf, 2.0) - 1.0 ) ), 0.5);
}


} // tudat
} // electromagnetism
