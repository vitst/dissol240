/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.


Description
    contact info: vitaliy.starchenko@gmail.com
    
    Write a description

SourceFiles
    ellipse.C

\*---------------------------------------------------------------------------*/

#ifndef ellipse_H
#define ellipse_H

#include "calcType.H"
#include <map>
#include <list>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace calcTypes
{

/*---------------------------------------------------------------------------*\
                           Class ellipse Declaration
\*---------------------------------------------------------------------------*/

class ellipse
:
    public calcType
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        ellipse(const ellipse&);

        //- Disallow default bitwise assignment
        void operator=(const ellipse&);
        
        word processingType_;
        label N_, M_, K_;
        label N1_, M1_, K1_;
        
        scalar dx, dy;

        scalar maxPosX, minPosX; // maximum and minimum in X direction
        scalar maxPosZ, minPosZ; // maximum and minimum in Z direction
        
        scalarField pointsXY_2D_X;
        scalarField pointsXY_2D_Z;
        
        pointField pointsXYonsurface;
        std::map<label, label> oppositeWallUniform;
        
protected:

    // Member Functions

        // Calculation routines

            //- Initialise - typically setting static variables,
            //  e.g. command line arguments
            virtual void init();

            //- Pre-time loop calculations
            virtual void preCalc
            (
                const argList& args,
                const Time& runTime,
                const fvMesh& mesh
            );

            //- Time loop calculations
            virtual void calc
            (
                const argList& args,
                const Time& runTime,
                const fvMesh& mesh
            );

            void build_surface_points
            (
              const fvMesh& mesh
            );
            
            scalar primitive_simpson_integration
            (
              scalarField& x,
              scalarField& y  
            );

        // I-O
            //- Write component fields
            void write_temp
            (
              const fvMesh& mesh,
              const Time& runTime
            );
            void write_aspect
            (
              const fvMesh& mesh,
              const Time& runTime
            );
            void write_q
            (
              const fvMesh& mesh,
              const Time& runTime
            );
            void write_p
            (
              const fvMesh& mesh,
              const Time& runTime
            );
            void write_Ccup
            (
              const fvMesh& mesh,
              const Time& runTime
            );

public:

    //- Runtime type information
    TypeName("ellipse");

    // Constructor
        ellipse();

    //- Destructor
        virtual ~ellipse();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace calcTypes
} // End namespace Foam

#endif

// ************************************************************************* //
