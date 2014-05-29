/* DIFFUSION ASSIGNATION HEADER */

// Guard
#ifndef _H_DIFFUSION_
#define _H_DIFFUSION_

// Standard Libraries
#include <utility>
#include <cstdio>
#include <iostream>
#include <cmath>
//inSilico
#include <base/io/vtk/LegacyWriter.hpp>
#include <base/io/Format.hpp>
#include <base/mesh/Size.hpp>
#include <base/linearAlgebra.hpp>
#include <base/post/evaluateField.hpp>

namespace diffusion
{
	const double coorTol = 10e-6;
	const unsigned maxIter = 10;

    // Initial condition
    template<typename VECDIM, typename DOF>
    static void initialState( const VECDIM& x, DOF* doFPtr )
    {
        const double value = 0.;
        doFPtr -> setValue( 0, value );
    }


	// Set Boundaries conditions
    template<typename VECDIM, typename DOF, typename POOL>
   	static void boundary ( const VECDIM& x, DOF* doFPtr, POOL& pool, int comp )
    {
		double conc = pool.medium.getConc ( comp );
		doFPtr -> constrainValue( 0, conc );
    }
	



	/* LINK TO MESH */
	// ID Finder
	std::size_t idLocator ( double * position, int * size, int numElements )
	{
		std::size_t id;
		double idCalc[3];
		double elemSize [3];

		for ( int i=0; i < 3; ++i )
		{
			elemSize[i] = (double) size[i] / (double) numElements;
			idCalc[i] = std::floor ( position[i] / elemSize[i] ); 
		}
		id = (int)( idCalc[2] * numElements * numElements ) + 
			 (int)( idCalc[1] * numElements ) + (int) idCalc[0];
		
		return id;
	}

	// Diffusion Constant setting
	template<typename ELEMENT, typename VECTOR>
	double diffusionConstant ( const ELEMENT* geomEp,
						  	   const typename ELEMENT::GeomFun::VecDim& xi,
							   VECTOR & elemDiff )
	{
		std::size_t id = geomEp -> getID();	
		return ( elemDiff[id] );
	}

	

	// Consumption Rate setting 
	template<typename ELEMENT, typename FIELD ,typename VECTOR>
	base::Vector<1>::Type comsRate ( const ELEMENT* geomEp,
						   				const typename ELEMENT::GeomFun::VecDim& xi,
                                      	const FIELD& field,
						   				VECTOR & elemRate )
	{
		std::size_t id = geomEp -> getID();	
		    
   		const typename FIELD::Element* fieldEp = field.elementPtr( geomEp -> getID() );
    
    	const base::Vector<1>::Type A = base::post::evaluateField( geomEp, fieldEp, xi );

		return ( elemRate[id] * A );
	}
	

	
	// Send discrete ifno to continous
	template <typename POOL, typename MESH, typename VECTOR>
	void diffusionLink ( POOL& pool, MESH& mesh, VECTOR& elemDiff, VECTOR& elemRate, int numElements )
	{
		int poolSize = pool.cells.size();
		int* scaffoldSize = pool.scaffold.getSize();	
		double elemSize = base::mesh::Size<typename MESH::Element>::apply( 
													mesh.elementPtr( 0 ) ); 
		for ( int i=0; i < poolSize; ++i )
		{
			if ( pool.cells[i].getStatus() != 0 )
			{
				std::size_t idPos;
				double * position = pool.cells[i].getPosition(); 
				
				idPos = idLocator ( position, scaffoldSize, numElements );	
			
				double volumeProp = pow ( pool.cells[i].getRadius(), 3 ) / 
									pow ( elemSize, 3 ); 

				if ( volumeProp > 1.0 ) { volumeProp = 1.0; };

				// Calculate Diffusion Constant and Set
				double newDiff = elemDiff.at ( idPos ) +
								 volumeProp * ( pool.scaffold.getDiffConst ( 2 ) -
												elemDiff.at ( idPos ) ); 
					
				elemDiff.at ( idPos ) = newDiff;

				// Calculate Comsumption Rate and Set
				double newRate = volumeProp * ( pool.cells[i].getRate() ) + 
												elemRate.at (idPos);
				elemRate.at ( idPos ) = -newRate;
			}
		}
	};

	// Get concentration and process Cell Health
	template <typename POOL, typename MESH, typename FIELD>
	void healthLink ( POOL& pool, MESH& mesh, FIELD& concentration, int numElements )
	{
		int poolSize = pool.cells.size();
		int* scaffoldSize = pool.scaffold.getSize();	
		double elemSize = base::mesh::Size<typename MESH::Element>::apply( 
													mesh.elementPtr( 0 ) ); 
		for ( int i=0; i < poolSize; ++i )
		{
			if ( pool.cells[i].getStatus() != 0 )
			{
				std::size_t idPos;
				double * position = pool.cells[i].getPosition(); 
				idPos = idLocator ( position, scaffoldSize, numElements );	
				
				//Get concentration
				base::Vector<1,double>::Type concElem;
				concElem = base::post::evaluateField ( 	
									mesh.elementPtr ( idPos ),
									concentration.elementPtr ( idPos ),
									base::ShapeCentroid<MESH::Element::shape>::apply() 
										  			 );
											
				//Check Health
				pool.procHealth ( i, concElem[0] );
				
			}
		}
	}







	/* OUTPUT */
	// Output to a VTK file
	template<typename MESH, typename CONC>
	void writeVTKFile( const std::string& baseName,
					   const unsigned     step,
					   const MESH&        mesh,
					   const CONC&        conc )
	{
		// VTK Legacy
		const std::string vtkFile =
			"./results/diffusion/" + baseName + "_" + base::io::leadingZeros( step ) 
					 + ".vtk";
		std::ofstream vtk( vtkFile.c_str() );
		base::io::vtk::LegacyWriter vtkWriter( vtk );
		vtkWriter.writeUnstructuredGrid( mesh );
		base::io::vtk::writePointData( vtkWriter, mesh, conc, "concentration" );
		vtk.close();
	}




}

#endif
