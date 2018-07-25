/* 
 * File:   rhsvector.h
 * Author: Author: Nadish Saini
 *
 * Created on March 1, 2016, 7:00 PM
 */

#ifndef RHSVECTOR_H
#define	RHSVECTOR_H
#include <fstream>
#include <iostream>

using namespace std;

/****Function to calculate Flux on Right Face of each cell****/
void rightfaceflux(vector< vector<double> > uelem, int degree, double a, vector<double> &faceflux, int nelem)
{
    
    if (degree == 1)
    {
        for (int i=0; i<nelem;i++)
        {
            double fluxrface;
            if (a > 0)
            {
                fluxrface= a*(uelem[i][0]);
            }
            if (a < 0)
            {
                fluxrface= a*(uelem[i+1][0]);
            }
            
            faceflux[i]=fluxrface;
        }
    }
    if (degree == 2)
    {
        for (int i=0; i<nelem;i++)
        {
            double fluxrface;
            if (a > 0)
            {
                fluxrface= a*(uelem[i][0]+uelem[i][1]);
            }
            if (a < 0)
            {
                fluxrface= a*(uelem[i+1][0] - uelem[i+1][1]);
            }
            
            faceflux[i]=fluxrface;
        }
    }
    
    if (degree == 3)
    {
        for (int i=0;i<nelem;i++)
        {
            double fluxrface;
            if(a > 0)
            {
                fluxrface= a*(uelem[i][0]+uelem[i][1]+uelem[i][2]);
            }
            if(a < 0)
            {
                fluxrface= a*(uelem[i+1][0]-uelem[i+1][1]+uelem[i+1][2]);
            }
            faceflux[i]=fluxrface;
        }
    }
}


/***Function to calculate Boundary Integral for each Element***/
void boundintegral(vector<double> faceflux, vector< vector<double> > &rhsbound, double a, int degree, int nelem)
{
    if (degree == 1)
    {
        for (int i=1;i < nelem-1;i++)
        {
            rhsbound[i][0]= faceflux[i] - faceflux[i-1];
        }
    }
    if (degree == 2)
    {
        for (int i=1;i < nelem-1;i++)
        {
            rhsbound[i][0]= faceflux[i] - faceflux[i-1];
            rhsbound[i][1]=faceflux[i] + faceflux[i-1];
        }
    }
    
    if (degree == 3)
    {
        for (int i=1;i<nelem-1;i++)
        {
            rhsbound[i][0]= faceflux[i] - faceflux[i-1];
            rhsbound[i][1]=faceflux[i] + faceflux[i-1];
            rhsbound[i][2]=faceflux[i] - faceflux[i-1];
        }
    }
}

/***Function to calculate Domain Integral for each Element***/
void domainintegral(vector< vector<double> > &rhsdomain, vector< vector<double> > uelem, int degree, int nelem, double a)
{
    if (degree == 1)
    {
        for (int i=1; i< nelem-1; i++)
        {
            rhsdomain[i][0]=0;
        }
    }
    if (degree == 2)
    {
        for (int i=1; i< nelem-1; i++)
        {
            rhsdomain[i][0]=0;
            rhsdomain[i][1]=2*a*uelem[i][0];
        }
    }
    
    if (degree == 3)
    {
        for (int i=1; i< nelem-1; i++)
        {
            rhsdomain[i][0]=0;
            rhsdomain[i][1]=2*a*uelem[i][0];
            rhsdomain[i][2]=2*a*uelem[i][1];
        }
    }
}

void rhsvector(vector< vector<double> > &R, int nelem, vector< vector<double> > uelem, int degree, double a, int fluxscheme, double deltax)
{
    
    /*****Calculate Fluxes at right face of each cell****/
    /**Note this will calculate fluxes for all elements on the right face*/
    vector<double> faceflux (nelem,0);
    rightfaceflux(uelem, degree, a, faceflux, nelem);   
    /**************************************************/
    
    /****Calculate Boundary Integral****/
    vector< vector<double> > rhsbound(nelem, vector<double> (degree, 0));    
    boundintegral(faceflux, rhsbound, a, degree, nelem);   
    /******************************************/
    
    
    
    /************Calculate Domain Integral**************/
    vector< vector<double> > rhsdomain(nelem, vector<double> (degree, 0));
    domainintegral(rhsdomain, uelem,  degree,  nelem,  a);
    /**************************************************/
    
    /**********Finally R.H.S. Vector******/
    for (int i=1; i< nelem-1;i++)
    {
        for (int j=0; j<degree; j++)
        {
            R[i][j] = rhsbound[i][j] - rhsdomain[i][j];
            R[i][j] = -R[i][j];
        }
    }
      
    
}

#endif	/* RHSVECTOR_H */
