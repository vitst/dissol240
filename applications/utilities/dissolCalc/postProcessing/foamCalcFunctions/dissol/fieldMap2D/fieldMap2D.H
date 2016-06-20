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
    fieldMap2D.C

\*---------------------------------------------------------------------------*/

#ifndef fieldMap2D_H
#define fieldMap2D_H

#include "calcType.H"
#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace calcTypes
{

/*---------------------------------------------------------------------------*\
                           Class fieldMap2D Declaration
\*---------------------------------------------------------------------------*/

class fieldMap2D
:
    public calcType
{
    // Private Member Functions
  
        const static int maxNumProcPoints = 1000000; // max size of the pattern
        
        int thisTimeSize;
        int totNumLoop;
        
        int sizeAA;
        
        int curNum;
        
        int curBlock;
        
        int curEnd;
        
        //- Disallow default bitwise copy construct
        fieldMap2D(const fieldMap2D&);

        //- Disallow default bitwise assignment
        void operator=(const fieldMap2D&);
        
        word processingType_;
        label N_, M_, K_;
        label N1_, M1_, K1_;
        
        int N1M1;
        
        scalar dx, dy;

        scalar maxPosX, minPosX; // maximum and minimum in X direction
        scalar maxPosY, minPosY; // maximum and minimum in Y direction
        scalar maxPosZ, minPosZ; // maximum and minimum in Z direction
        
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
            void write_h
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

            void write_all
            (
              const fvMesh& mesh,
              const Time& runTime
            );
            
public:

    //- Runtime type information
    TypeName("fieldMap2D");

    // Constructor
        fieldMap2D();

    //- Destructor
        virtual ~fieldMap2D();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace calcTypes
} // End namespace Foam

#endif

// ************************************************************************* //
