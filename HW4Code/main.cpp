/* 
 * File:   main.cpp
 * Author: Author: Nadish Saini
 *
 * Created on March 1, 2016, 12:38 PM
 */
#include <iostream>
# include <stdlib.h>
# include <math.h>
# include <fstream>
#include <cmath>
#include <time.h>
#include <algorithm>
#include <sstream>
#include <string>
#include <vector>
#include <fenv.h>
#include <stdio.h>
#include <limits.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <omp.h>
#include <time.h>
#include <cctype>
#include <cstring>
#include "controlinput.h"
#include "rhsvector.h"

using namespace std;


/*****Function to construct velocity vector for element from Inital velocity profile*/
void elemVel (vector< vector<double> > &uelem, vector<double> u, int degree, int nelem, vector<double> ucenter)
{
    if (degree == 1)
    {
        for (int i=0; i< nelem; i++)
        {
            uelem[i][0]=ucenter[i];
        }
    }
    if (degree == 2)
    {
        for (int i=0; i< nelem; i++)
        {
            uelem[i][0]=(u[i]+u[i+1])*0.5;
            uelem[i][1]=(u[i+1]-u[i])*0.5;
        }
    }
    
    if (degree == 3)
    {
        for (int i=0;i<nelem;i++)
        {
            uelem[i][0]=(u[i]+u[i+1]+4*ucenter[i])/6;
            uelem[i][1]=(u[i+1]-u[i])/2;
            uelem[i][2]=(u[i]+u[i+1]-2*ucenter[i])/3;
        }
    }
}

/***Function to update Ghost Velocities***/
void ghostupdate(vector <vector <double> > &uelem, int nelem, int degree)
{
        uelem[0][0]=uelem[nelem-2][0];
        uelem[0][1]=uelem[nelem-2][1];
        
        uelem[nelem-1][0]=uelem[1][0];
        uelem[nelem-1][1]=uelem[1][1];
}

/****Function to Generate Ghost Cells***/
void ghost(vector <double> &u, int &nelem, int degree)
{
    vector <double> tempu(nelem+1);
    
    for (int i=0; i< nelem+1;i++)
    {
        tempu[i]=u[i];
    }
    
        nelem=nelem+2;
        u.resize(nelem+1);
        
        u[0]=tempu[nelem-2];
        u[nelem] = tempu[1];
        for(int i=1;i<nelem;i++)
        {
            u[i]=tempu[i-1];
        }
    
}

/***Function to find element properties (length)***/
void elemprop(vector< vector<double> > &geoel, int nelem, int degree, double deltax, double x0)
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
void euler (vector< vector<double> > &uelem,vector< vector<double> > M, vector< vector<double> > R, double deltat, int nelem, int degree)
{
    for (int i=1;i<nelem-1;i++)
    {
        for (int j=0;j<degree;j++)
        {
            uelem[i][j] = uelem[i][j] + deltat*R[i][j]/M[i][j];
        }
    }
}



int main()
{
    /****Initializing Variables*****/
    double x0,xn;
    int nelem;
    int degree;
    double cfl;
    double a;
    vector<double> u;
    int fluxscheme=0;
    double totaltime;
    int stages;
    vector<double> ucenter;
    /******************************/
    
    /*****Read Control File*****/    
    control(x0,xn,nelem,u,degree,cfl,a,fluxscheme,totaltime,stages,ucenter);
    /****************************/
       
    
    /*****Generate Ghost Cells*****/
    ghost(u,nelem,degree);
    vector<double> ui (nelem+1,0);
    /**Store Initial Soln in a separate matrix*/
    for (int i=1;i<nelem;i++)
    {
        ui[i]=u[i];
    }
    /****************************/
    
    
    
    
    /****Get element properties****/
    double deltax=(xn-x0)/(nelem-2);  //Length of cell
    vector< vector<double> > geoel;
    geoel.resize(nelem);
    elemprop(geoel, nelem, degree, deltax, x0);
    
    /********************************/
    
    
    /***Compute Mass Matrix****/
    vector< vector<double> > M (nelem, vector<double> (degree,0));
    massmatrix(M,geoel, nelem, degree);  
    
    /*************************/
    
    
    /****Construct Velocity vector for each element****/
    vector< vector<double> > uelem(nelem, vector<double> (degree,0));
    elemVel (uelem, u, degree, nelem,ucenter);
    
    /**************************************************/
    
    
            
    double deltat=cfl*deltax/fabs(a);
    double t=0;
    vector< vector<double> > R (nelem, vector<double> (degree,0));
    
    
    /****Time loop*****/
    while (t<=totaltime)
    {
        if (stages == 1)
        {
        rhsvector(R, nelem, uelem, degree, a, fluxscheme, deltax);
        euler(uelem, M, R, deltat, nelem, degree);
        ghostupdate(uelem,nelem,degree);
        }
        
        
        
        if (stages ==2)
        {
            vector< vector<double> > utemp(nelem, vector<double> (degree,0));
            for (int i=1;i<nelem-1;i++)
            {
                for (int j=0;j<degree;j++)
                {
                    utemp[i][j]=uelem[i][j];
                }
            }
            rhsvector(R, nelem, uelem, degree, a, fluxscheme, deltax);
            euler(uelem, M, R, deltat, nelem, degree);
            ghostupdate(uelem,nelem,degree);
            
            rhsvector(R, nelem, uelem, degree, a, fluxscheme, deltax);
            euler(uelem, M, R, deltat, nelem, degree);
            for (int i=1;i<nelem-1;i++)
            {
                for (int j=0;j<degree;j++)
                {
                    uelem[i][j]= 0.5*utemp[i][j]+0.5*uelem[i][j];
                }
            }
            ghostupdate(uelem,nelem,degree);
            
        }
        
        else if (stages == 3)
        {
            vector< vector<double> > utemp(nelem, vector<double> (degree,0));
            for (int i=1;i<nelem-1;i++)
            {
                for (int j=0;j<degree;j++)
                {
                    utemp[i][j]=uelem[i][j];
                }
            }
            rhsvector(R, nelem, uelem, degree, a, fluxscheme, deltax);
            euler(uelem, M, R, deltat, nelem, degree);
            ghostupdate(uelem,nelem,degree);
            
            rhsvector(R, nelem, uelem, degree, a, fluxscheme, deltax);
            euler(uelem, M, R, deltat, nelem, degree);
            ghostupdate(uelem,nelem,degree);
            for (int i=1;i<nelem-1;i++)
            {
                for (int j=0;j<degree;j++)
                {
                    uelem[i][j]= (3.0/4.0)*utemp[i][j]+(1.0/4.0)*uelem[i][j];
                }
            }
            ghostupdate(uelem,nelem,degree);
            
            rhsvector(R, nelem, uelem, degree, a, fluxscheme, deltax);
            euler(uelem, M, R, deltat, nelem, degree);
            for (int i=1;i<nelem-1;i++)
            {
                for (int j=0;j<degree;j++)
                {
                    uelem[i][j]= (1.0/3.0)*utemp[i][j]+(2.0/3.0)*uelem[i][j];
                }
            }
            ghostupdate(uelem,nelem,degree);
        }
        
        t=t+deltat;
    }
    
    /*****Get the final Nodal Velocities***********/
    /***Reconstruct using original function****/
    for(int i=1; i<nelem;i++)
    {
        if (degree == 1)
        {
            u[i]= uelem[i][0];
        }
        if (degree == 2)
        {
        u[i]= uelem[i][0]-uelem[i][1];
        }
        
        if (degree == 3)
        {
            u[i]= uelem[i][0]-uelem[i][1]+uelem[i][2];
        }
    }
    
    
    ofstream output;
    output.open("out.txt");
    for (int i=1; i<nelem;i++)        
    {
        output<<x0+(i-1)*deltax<<" "<<ui[i]<<" "<<u[i]<<endl;
    }
    output.close();
    
    
    
    
    
}
         

