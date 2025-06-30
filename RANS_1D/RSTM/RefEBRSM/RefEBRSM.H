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
    Foam::incompressible::RASModels::RefEBRSM

Description
    The reference model provided in Appendix C of http://ac.els-cdn.com/S0142727X14001179/1-s2.0-S0142727X14001179-main.pdf?_tid=1e73b44c-7c30-11e7-be3f-00000aacb360&acdnat=1502193357_2298de9f18e5288eaed1c49d6db57303

SourceFiles
    RefEBRSM.C

\*---------------------------------------------------------------------------*/

#ifndef RefEBRSM_H
#define RefEBRSM_H

#include "GenElliptic.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class RefEBRSM Declaration
\*---------------------------------------------------------------------------*/

class RefEBRSM
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
            volScalarField nut_; 

    //...................................................................................................
    //      bool LimitTS; // limit the timescales tauwlim and tauhlim.......................................................
            bool SGDHEps_;
            bool SGDHR_;
           
public:

    //- Runtime type information
    TypeName("RefEBRSM");

    // Constructors

        //- Construct from components
        RefEBRSM
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~RefEBRSM()
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
