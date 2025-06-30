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

#include "RefEBRSMEterm1.H"
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

defineTypeNameAndDebug(RefEBRSMEterm1, 0);
addToRunTimeSelectionTable(RASModel, RefEBRSMEterm1, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

RefEBRSMEterm1::RefEBRSMEterm1
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
            0.22
        )
    ),
    CSEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CSEps",
            coeffDict_,
            0.183
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
            0.1
        )
    ),
    A1Eterm_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A1Eterm",
            coeffDict_,
            0.085
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
    SGDHEps_(coeffDict_.lookupOrAddDefault("SGDHEps", true)),
    SGDHR_(coeffDict_.lookupOrAddDefault("SGDHR", true)),
    nutCMuk2_(coeffDict_.lookupOrAddDefault("nutCMuk2", false))

{

    volVectorField n = fvc::grad(alpha_);
    n /= mag(n) + dimensionedScalar("nsmall", n.dimensions(), VSMALL);
    volSymmTensorField nn = symm(n*n);
    volScalarField Ts("T", T());

    dimensionedScalar epsilonSmall_
	(
			"epsilonSmall_",
			dimensionSet(0,2,-3,0,0,0,0),
			scalar(ROOTVSMALL)
	);

    if(nutCMuk2_)
    {
    nut_ = CMu*sqr(k_)/(epsilon_+ epsilonSmall_);
    }
    else
    {
    nut_ = ( (1-pow(alpha_, 3)) * (R_&&nn) + pow(alpha_, 3) * k_) *Cmu_ *Ts;
    }
    nut_.correctBoundaryConditions();
    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool RefEBRSMEterm1::read()
{
    if (RASModel::read())
    {
        CSEps_.readIfPresent(coeffDict());
        CSR_.readIfPresent(coeffDict());
        A1p_.readIfPresent(coeffDict());
        A1Eterm_.readIfPresent(coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}


void RefEBRSMEterm1::correct()
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
    //time-scale limiters ............................................................................................
//    volScalarField tauwlim=  max( k_/(epsilon_ + epsilonSmall_),  10.0 * pow(alpha_, 3) * sqrt(nu()/(epsilon_ + epsilonSmall_)));
//    volScalarField tauhlim=  max( k_/(epsilon_ + epsilonSmall_),  10.0 * sqrt(nu()/(epsilon_ + epsilonSmall_)));
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

    // Dissipation equation 
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
        + fvm::div(phi_, epsilon_)
        - fvm::Sp(fvc::div(phi_), epsilon_)
        - fvm::laplacian(nu(), epsilon_, "laplacian(epsilon)")
        - ( // Turbulent diffusion..................................................... 
            SGDHEps_ == true
          ? fvm::laplacian(nut_/sigmaEps_, epsilon_, "laplacian(epsilon)")
          : fvm::laplacian(CSEps_ * Tsf * Rdiagf, epsilon_, "laplacian(epsilon)")
          )
      ==
        C1_ * G/Ts 
        - fvm::Sp(C2_/Ts, epsilon_)
        + A1Eterm_* nu() * (R_&&nn) * k_/(epsilon_ + epsilonSmall_) * (1-fa) * fvc::magSqrGradGrad(U_)       //E term ...........................
    );

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
                 ? fvm::laplacian(nut_/sigmaK_, R_, "laplacian(R)")
                 : fvm::laplacian(CSR_ * Tsf * Rdiagf, R_, "laplacian(R)") // Daly-Harlow diffusion model
              ) 
            - fvm::laplacian(nu(), R_, "laplacian(R)")
            + fvm::Sp(imSrc, R_)                      // pressure-rate-of-strain, implicit (RHS) source
            + fvm::Sp((1-fa)/tauwlim, R_)         // wall contribution of epsilon
            ==                                        // ==
            P                                         // production tensor
            - twoThirdsI * k_/tauhlim * fa              // homogeneous contribution of epsilon
            + exSrc                                   // pressure-rate-of-strain, explicit (LHS) source
            + psw                                     // pressure-rate-of-strain, wall contribution
        );

    
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
    if(nutCMuk2_)
    {
    nut_ = CMu*sqr(k_)/(epsilon_+ epsilonSmall_);
    }
    else
    {
    nut_ = ( (1-pow(alpha_, 3)) * (R_&&nn) + pow(alpha_, 3) * k_) *Cmu_ *Ts;
    }
    nut_.correctBoundaryConditions();
   
/**
       // Correct wall shear stresses

    const fvPatchList& patches = mesh_.boundary();
 
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            symmTensorField& Rw = R_.boundaryField()[patchi];

            const scalarField& nutw = nut_.boundaryField()[patchi];

            const vectorField snGradU(U_.boundaryField()[patchi].snGrad());

            const vectorField& faceAreas
                = mesh_.Sf().boundaryField()[patchi];

            const scalarField& magFaceAreas
                = mesh_.magSf().boundaryField()[patchi];

            forAll(curPatch, facei)
            {
                // Calculate near-wall velocity gradient
                tensor gradUw
                    = (faceAreas[facei]/magFaceAreas[facei])*snGradU[facei];

                // Calculate near-wall shear-stress tensor
                tensor tauw = -nutw[facei]*2*symm(gradUw);

                // Reset the shear components of the stress tensor
                Rw[facei].xy() = tauw.xy();
                Rw[facei].xz() = tauw.xz();
                Rw[facei].yz() = tauw.yz();
            }
        }
    }
**/

    // Finally, re-calculate turbulent viscosity according to V2F model
    // nut_ = 2.0/3.0 * Cmu_ * max((R_ && nn), k0_) * Ts;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
