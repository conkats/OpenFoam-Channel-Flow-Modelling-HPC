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

Class
    Foam::incompressible::RASModels::RefEBRSMEterm1

Description
    S Lardeau and R. Manceau paper provided in https://hal.archives-ouvertes.fr/hal-01051799/document (Computations of canonical and complex flow
configurations using a robust formulation of the elliptic-blending Reynolds-Stress model)

SourceFiles
    RefEBRSMEterm1.C

\*---------------------------------------------------------------------------*/

#ifndef RefEBRSMEterm1_H
#define RefEBRSMEterm1_H

#include "GenElliptic.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class RefEBRSMEterm1 Declaration
\*---------------------------------------------------------------------------*/

class RefEBRSMEterm1
:
    public GenElliptic
{
    // Private data
      // Fields
            volScalarField alpha_;
            Switch crossTurbDiffusion_; // turn on cross diffusion due to Reynolds stresses
            dimensionedScalar CMu;
            dimensionedScalar CSEps_; 
            dimensionedScalar CSR_; 
            dimensionedScalar A1p_; 
            dimensionedScalar A1Eterm_;
            volScalarField nut_; 

    //...................................................................................................
    //      bool LimitTS; // limit the timescales tauwlim and tauhlim.......................................................
            bool SGDHEps_;
            bool SGDHR_;
            bool nutCMuk2_;
           
public:

    //- Runtime type information
    TypeName("RefEBRSMEterm1");

    // Constructors

        //- Construct from components
        RefEBRSMEterm1
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~RefEBRSMEterm1()
    {}


    // Member Functions
        
        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
