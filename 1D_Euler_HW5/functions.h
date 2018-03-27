/* 
 * File:   functions.h
 * Author: nash
 *
 * Created on March 11, 2016, 4:08 PM
 */

#ifndef FUNCTIONS_H
#define	FUNCTIONS_H

#include <vector>
#include <math.h>
#include <omp.h>
using namespace std;



/*****Function to construct velocity vector for element from Inital velocity profile*/
void elemVector (vector< vector< vector<double> > > &Uelem, vector<vector <double> > U, int degree, int nelem, vector<vector <double> > Ucenter, int neqns)
{
    if (degree == 1)
    {
        for (int j=0; j<nelem; j++)
        {
            for (int i=0; i< neqns; i++)
            {
                Uelem[j][i][0]=Ucenter[j][i];
            }
        }
    }
    if (degree == 2)
    {
        for (int j=0;j<nelem;j++)
        {
            for (int i=0; i< neqns; i++)
            {
               Uelem[j][i][0]=(U[j][i]+U[j+1][i])*0.5;
               Uelem[j][i][1]=(U[j+1][i]-U[j][i])*0.5;
            }
        }
    }
    
    if (degree == 3)
    {
        for (int j=0; j<nelem; j++)
        {
            for (int i=0;i<neqns;i++)
            {
               Uelem[j][i][0]=(U[j][i]+U[j+1][i]+4*Ucenter[j][i])/6;
               Uelem[j][i][1]=(U[j][i]-U[j+1][i])/2;
               Uelem[j][i][2]=(U[j][i]+U[j+1][i]-2*Ucenter[j][i])/3;
            }
        }
    }
}

/***Function to update Ghost Velocities***/
void ghostupdate(vector< vector< vector<double> > > &Uelem, int nelem, int degree, int neqns, int BC)
{
    for (int j=0; j<neqns; j++)
    {
        for (int k=0; k<degree; k++)
        {
            if (BC ==1)
            {
            Uelem[0][j][k]=Uelem[nelem-2][j][k];
        
            Uelem[nelem-1][j][k]=Uelem[1][j][k];
            }
            
            else if (BC == 2)
            {
                Uelem[0][j][k]=Uelem[1][j][k];
                Uelem[nelem-1][j][k]=Uelem[nelem-2][j][k];
            }
        }        
    }
}

/****Function to Generate Ghost Cells***/
void ghost(vector <vector <double> > &W, int &nelem, int degree, int neqns, vector<vector <double> > &Wcenter, int BC)
{
    vector <vector <double> > tempW(nelem+1, vector<double> (neqns,0));
    
    for (int j=0; j< neqns; j++)
    {
       for (int i=0; i< nelem+1;i++)
       {
           tempW[i][j]=W[i][j];
       }
    }
    nelem=nelem+2;
    W.resize(nelem+1);
    W[nelem].resize(neqns);
    W[nelem-1].resize(neqns);
    
               
    for (int j=0;j<neqns;j++)
    {
        for(int i=1;i<nelem;i++)
        {
            W[i][j]=tempW[i-1][j];
        }
    }
    
    //////
    if(BC ==1)
    {
    for (int j=0;j<neqns;j++)
    {        
        W[0][j]=tempW[nelem-3][j];
        W[nelem][j] = tempW[1][j];
    }
    for(int j=0;j<neqns;j++)
    {
        Wcenter[0][j]=Wcenter[nelem-2][j];
        Wcenter[nelem-1][j]=Wcenter[1][j];
    }
    }
    
    else if (BC==2)
    {
        for (int j=0;j<neqns;j++)
        {
            W[0][j] = tempW[0][j];
            W[nelem][j]=tempW[nelem-2][j];
            
            Wcenter[0][j]=Wcenter[1][j];
            Wcenter[nelem-1][j]=Wcenter[nelem-2][j];
        }
    }
}

/***Function to find element properties (length)***/
void elemprop(vector< vector<double> > &geoel, int nelem, double deltax, double x0)
{        
        for (int i=0;i<nelem;i++)
        {
            double x1=x0+(i-1)*deltax;
            double x2=x1+deltax;
            
            geoel[i].resize(1);
            geoel[i][0]=(x2-x1);
            //cout<<geoel[i][0]<<endl;
        }
}

/***Function to construct Mass matrix for elements***/
void massmatrix(vector< vector<double> > &M, vector< vector<double> > geoel, int nelem, int degree)
{
    if (degree == 1)
    {
        for (int i=0; i<nelem;i++)
        {
            M[i][0]=geoel[i][0];
        }
    }
    if (degree == 2)
    {
        for (int i=0; i<nelem;i++)
        {
            M[i][0]=geoel[i][0];
            M[i][1]=geoel[i][0]/3;
        }
    }    
    else if (degree == 3)
    {
        for (int i=0; i<nelem;i++)
        {
            M[i][0]=geoel[i][0];
            M[i][1]=geoel[i][0]/3;
            M[i][2]=geoel[i][0]/5;
        }
    }
}

/***Function for Time integration by forward Euler Method***/
void euler (vector< vector< vector<double> > > &Uelem, vector< vector<double> > M, vector< vector< vector<double> > > R, double deltat, int nelem, int degree, int neqns)
{
    #pragma omp parallel for schedule(dynamic)
    for ( int i=1;i<nelem-1;i++)
    {
        for (int j=0;j<neqns;j++)
        {
            for (int k=0;k<degree;k++)
            {
                Uelem[i][j][k] = Uelem[i][j][k] + deltat*R[i][j][k]/M[i][k];
            }
        }
    }
}


#endif	/* FUNCTIONS_H */

