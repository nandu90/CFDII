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

double minmod(double a, double b)
{
    double result;
    if(a*b > 0)
    {
        if (fabs(a) < fabs(b))
        {
            result= a;
        }
        else if (fabs (a) >= fabs(b))
        {
            result= b;
        }
    }
    else if (a*b <= 0)
    {
        result= 0;
    }
    return result;
}

double mc_limit(double a, double b, double c)
{
    /*double result;
    if (a*b > 0 && a*c > 0)
    {
        double sign=b/fabs(b);
        
        double num;
        if(fabs(a) < fabs(c))
        {
            num=2*fabs(a);
        }
        else if(fabs(a) >= fabs(c))
        {
            num=2*fabs(c);
        }
        
        if(fabs(b) < num)
        {
            result = fabs(b)*sign;
        }
        else if(fabs(b) >= num)
        {
            result = num*sign;
        }
    }*/
    
    double num=minmod(2*a,2*c);
    double result =minmod(b,num);
    return result;
}

/****Function to calculate Flux on Right Face of each cell****/
void rightfaceflux(vector< vector<double> > uelem, int degree, double a, vector<double> &faceflux, int nelem, int reconstruct, double deltax)
{
    
    if (degree == 1)
    {
        for (int i=0; i<nelem-1;i++)
        {
            double uL,uR;
            int k=reconstruct;
            if(i==0)
            {
                uL=uelem[i][0]+0.25*(1-k)*(uelem[i][0]-uelem[nelem-2][0])+0.25*(1+k)*(uelem[i+1][0]-uelem[i][0]);
            }
            else
            {
            uL=uelem[i][0]+0.25*(1-k)*(uelem[i][0]-uelem[i-1][0])+0.25*(1+k)*(uelem[i+1][0]-uelem[i][0]);
            }
            
            if(i==nelem-2)
            {
                uR=uelem[i+1][0]-0.25*(1+k)*(uelem[i+1][0]-uelem[i][0])-0.25*(1-k)*(uelem[2][0]-uelem[i+1][0]);
            }
            else
            {
                uR=uelem[i+1][0]-0.25*(1+k)*(uelem[i+1][0]-uelem[i][0])-0.25*(1-k)*(uelem[i+2][0]-uelem[i+1][0]);
            }
            
            double fluxrface;
            
            if(k==2)
            {
            if (a > 0)
            {
                fluxrface= a*(uelem[i][0]);
            }
            if (a < 0)
            {
                fluxrface= a*(uelem[i+1][0]);
            }
            }
            
            else if(k == 3) //minmod limiter
            {
                double uface=0;
                double num1=0;
                double num2=0;
                
                if(i == 0)  //consideration for periodic BC
                {
                    num1=(uelem[i][0]-uelem[nelem-2][0])/deltax;
                }
                else
                {
                    num1=(uelem[i][0]-uelem[i-1][0])/deltax;
                }
                num2=(uelem[i+1][0]-uelem[i][0])/deltax;
                uface= uelem[i][0] + minmod(num1,num2)*deltax/2;
                
                fluxrface=a*uface;
            }
            
            else if(k == 4) //MC limiter
            {
                double uface=0;
                double num1=0;
                double num2=0;
                double num3=0;
                
                if(i == 0)  //consideration for periodic BC
                {
                    num1=(uelem[i][0]-uelem[nelem-2][0])/deltax;
                    num2=(uelem[i+1][0]-uelem[nelem-2][0])/(2*deltax);
                }
                else
                {
                    num1=(uelem[i][0]-uelem[i-1][0])/deltax;
                    num2=(uelem[i+1][0]-uelem[i-1][0])/(2*deltax);
                }
                num3=(uelem[i+1][0]-uelem[i][0])/deltax;
                uface= uelem[i][0] + mc_limit(num1,num2,num3)*deltax/2;
                
                fluxrface=a*uface;
            }
            
            else// if (k == 0 || k==-1 || k==1)
            {
            if (a > 0)
            {
                fluxrface= a*(uL);
            }
            if (a < 0)
            {
                fluxrface= a*(uR);
            }
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

void rhsvector(vector< vector<double> > &R, int nelem, vector< vector<double> > uelem, int degree, double a, int fluxscheme, double deltax, int reconstruct)
{
    
    /*****Calculate Fluxes at right face of each cell****/
    /**Note this will calculate fluxes for all elements on the right face*/
    vector<double> faceflux (nelem,0);
    rightfaceflux(uelem, degree, a, faceflux, nelem, reconstruct,deltax);   
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
