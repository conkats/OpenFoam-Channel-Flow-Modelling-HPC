/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "RefEBRSM.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"
#include "fixedInternalValueFvPatchField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(RefEBRSM, 0);
addToRunTimeSelectionTable(RASModel, RefEBRSM, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

RefEBRSM::RefEBRSM
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel,
    const word& turbulenceModelName,
    const word& modelName 
)
:
    RASModel(typeName, U, phi, lamTransportModel),
    GenElliptic(U, phi, lamTransportModel),

    alpha_
    (
        IOobject
        (
            "alpha",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    crossTurbDiffusion_(coeffDict_.lookupOrAddDefault<Switch>("crossTurbDiffusion", false)),
    CMu
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CMu",
            coeffDict_,
            0.21
        )
    ),
    CSEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CSEps",
            coeffDict_,
            0.1826
        )
    ),
    CSR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CSR",
            coeffDict_,
            0.21
        )
    ),
    A1p_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A1p",
            coeffDict_,
            0.065
        )
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    // The reference model does not use SGDH 
    SGDHEps_(coeffDict_.lookupOrAddDefault("SGDHEps", false)),
    SGDHR_(coeffDict_.lookupOrAddDefault("SGDHR", false))

{

    volVectorField n = fvc::grad(alpha_);
    n /= mag(n) + dimensionedScalar("nsmall", n.dimensions(), VSMALL);
    volSymmTensorField nn = symm(n*n);
    volScalarField Ts("T", T());

    // nut_ from Sylvains paper . Not needed for the reference model. However, it is needed if SGDH is used 
    nut_ = ( (1-pow(alpha_, 3)) * (R_&&nn) + pow(alpha_, 3) * k_) *Cmu_ *Ts;
    nut_.correctBoundaryConditions();
    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool RefEBRSM::read()
{
    if (RASModel::read())
    {
        CSEps_.readIfPresent(coeffDict());
        CSR_.readIfPresent(coeffDict());
        A1p_.readIfPresent(coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}


void RefEBRSM::correct()
{
    GenElliptic::correct();

    if (!turbulence_)
    {
        return;
    }

    volSymmTensorField P = -twoSymm(R_ & fvc::grad(U_));
    volScalarField G("RASModel::G", 0.5*mag(tr(P)));

    volScalarField Ts("T", T());
 
    dimensionedScalar epsilonSmall_
	(
			"epsilonSmall_",
			dimensionSet(0,2,-3,0,0,0,0),
			scalar(ROOTVSMALL)
	);

      volScalarField tauwlim=  k_/(epsilon_ + epsilonSmall_);
      volScalarField tauhlim=  k_/(epsilon_ + epsilonSmall_);
      volScalarField taunolim=  k_/(epsilon_ + epsilonSmall_);  

    // Direction of alpha
    volVectorField n = fvc::grad(alpha_);
    n /= mag(n) + dimensionedScalar("nsmall", n.dimensions(), VSMALL);
    volSymmTensorField nn = symm(n*n);

   #include "../include/epsilonWallI2.H" // set patch internal eps values

    // split R_ into normal diffusion and cross diffusion terms
    volSymmTensorField Rdiag = R_;
    dimensionedScalar kzero = kMin_ * 0.0;
    Rdiag.replace(symmTensor::XY, kzero);
    Rdiag.replace(symmTensor::YZ, kzero);
    Rdiag.replace(symmTensor::XZ, kzero);
    volSymmTensorField Rupper = R_ - Rdiag;

    symmTensor  minDiagR = gMin(Rdiag);
    if(
        minDiagR.xx() < 0.0 ||
        minDiagR.yy() < 0.0 ||
        minDiagR.zz() < 0.0 )
    {
        Info << "muDurbin::correct():: Warning! " << nl
             << "negative diagonal for R. I will probably fail soon! Rdiag.min = "
             << minDiagR << endl;
    }

    surfaceScalarField Tsf = fvc::interpolate(Ts, "interpolate(T)");
    surfaceSymmTensorField Rdiagf  = fvc::interpolate(Rdiag, "interpolate(R)");
    surfaceSymmTensorField Rupperf = fvc::interpolate(Rupper, "interpolate(R)");

    // weight of the homogeneous contribution of the pressure strain
    // indicating wall distance

    volScalarField fa = pow(alpha_, 3);

    volScalarField C1Modify = 1.0 + A1p_*(1.0 - fa)
        * (G/(epsilon_ + epsilonSmall_)); 

    // Dissipation equation 
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
        + fvm::div(phi_, epsilon_)
        - fvm::Sp(fvc::div(phi_), epsilon_)
        - fvm::laplacian(nu(), epsilon_)
        - ( // Turbulent diffusion..................................................... 
            SGDHEps_ == true
          ? fvm::laplacian(nut_/sigmaEps_, epsilon_)
          : fvm::laplacian(CSEps_ * Tsf * Rdiagf, epsilon_)
          )
      ==
        C1_ * G/Ts * C1Modify
        - fvm::Sp(C2_/Ts, epsilon_)
    );
  
    if(crossTurbDiffusion_)
    {
        epsEqn() -= fvc::laplacian(CSEps_ * Tsf * Rupperf, epsilon_);
    }

    epsEqn().relax();
    epsEqn().boundaryManipulate(epsilon_.boundaryField());

    solve(epsEqn);
    bound(epsilon_, epsilonMin_);

    // wall contribution of pressure strain (weight included). Note:
    // repeated typo in Eqs. (9) and (A.10) in (Manceau & Hanjalic 2002)
    volSymmTensorField psw =  (-5.0) * (1.0-fa) /tauwlim *
        ( 
            twoSymm(R_ & nn) 
            - 0.5 * (R_ && nn) * (nn + I)
        );

    // used LRR-IP as default; SSG is also option.    
    volScalarField imSrc = Clrr1_*fa*epsilon_/k_;
    volSymmTensorField exSrc = (Clrr1_ * twoThirdsI) * fa * epsilon_ - fa * Clrr2_*dev(P);
    
    if(SSG_) // SSG model. Disabled by default
    {
        volSymmTensorField bij("bij", 0.5 * dev(R_/k_));
        volTensorField fbij(symm2full(bij));
        volTensorField     gradU = fvc::grad(U_);
        volSymmTensorField Sij("Sij", dev(symm(gradU)));
        volTensorField fSij(symm2full(Sij));
        volTensorField     Wij =  skew(gradU);

    // IMPORTANT used                      not Cg4_ * k_ * dev( twoSymm(fbij & fSij) ) ..........................................................................................
        imSrc = 0.5 * fa * (Cg1_ / taunolim + Cg1s_ * G / k_);
        exSrc = fa *
            (  (Cg1_ * k_ / taunolim + Cg1s_ * G) * oneThirdI
            +    Cg2_ * k_/taunolim * dev(symm(bij & bij))
            +   (Cg3_ - Cg3s_ * sqrt(bij && bij)) * k_ * Sij
            +    Cg4_ * k_ * ( twoSymm(fbij & fSij) - 2.0/3.0 * I * (fbij && fSij))
            +    Cg5_ * k_ *   twoSymm(bij & Wij) );
    }
   
    // Reynolds stress equation
    tmp<fvSymmTensorMatrix> REqn
        (
            fvm::ddt(R_)
            + fvm::div(phi_, R_)
            - fvm::Sp(fvc::div(phi_), R_)
            - ( // Turbulent diffusion..................................................... 
                 SGDHR_ == true
                 ? fvm::laplacian(nut_/sigmaK_, R_)
                 : fvm::laplacian(CSR_ * Tsf * Rdiagf, R_) // Daly-Harlow diffusion model
              ) 
            - fvm::laplacian(nu(), R_)
            + fvm::Sp(imSrc, R_)                      // pressure-rate-of-strain, implicit (RHS) source
            + fvm::Sp((1-fa)/tauwlim, R_)         // wall contribution of epsilon
            ==                                        // ==
            P                                         // production tensor
            - twoThirdsI * k_/tauhlim * fa              // homogeneous contribution of epsilon
            + exSrc                                   // pressure-rate-of-strain, explicit (LHS) source
            + psw                                     // pressure-rate-of-strain, wall contribution
        );

    if(crossTurbDiffusion_)
    {
        REqn() -= fvc::laplacian(CSR_*Ts*Rupper, R_, "laplacian(R)");
    }
   
    REqn().relax();
    solve(REqn);
    
    R_.max
    (
    dimensionedSymmTensor
            (
                "zero",
                R_.dimensions(),
                symmTensor
                (
                    kMin_.value(), -GREAT, -GREAT,
                    kMin_.value(), -GREAT,
                    kMin_.value()
                )
           )
    );

    k_ = 0.5*tr(R_);

    bound(k_, kMin_);
    
    volScalarField L_ = L();
    Ts = T(); // re-compute time scale


    tmp<fvScalarMatrix> alphaEqn
    (
        fvm::laplacian(alpha_)
     ==
        fvm::Sp(1.0/sqr(L_), alpha_)
        - scalar(1.0)/( sqr(L_) )
    );
    
    alphaEqn().relax();
    solve(alphaEqn);

    // Re-calculate viscosity
    nut_ = ( (1-pow(alpha_, 3)) * (R_&&nn) + pow(alpha_, 3) * k_) *Cmu_ *Ts;
    nut_.correctBoundaryConditions();


    // Finally, re-calculate turbulent viscosity according to V2F model
    // nut_ = 2.0/3.0 * Cmu_ * max((R_ && nn), k0_) * Ts;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
