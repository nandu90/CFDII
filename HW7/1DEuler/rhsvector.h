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
#include <sstream>
#include <omp.h>
#include "fluxes.h"
#include "higherorderfluxes.h"

/****Function to calculate Flux on Right Face of each cell****/
void rightfaceflux(vector< vector<vector<double> > > Uelem, int degree, vector<vector <double> > &faceflux, int nelem, int neqns, double g, int fluxscheme, double t, double deltax)
{
    
    if (degree == 1)
    {
        for (int i=0; i<nelem-1;i++)
        {
            double uL=Uelem[i][1][0]/Uelem[i][0][0];
            double pL=(g-1)*(Uelem[i][2][0]-0.5*Uelem[i][1][0]*uL);
            double cL=sqrt(fabs(g*pL/Uelem[i][0][0]));
            double ML=uL/cL;
            
            double uR=Uelem[i+1][1][0]/Uelem[i+1][0][0];
            double pR=(g-1)*(Uelem[i+1][2][0]-0.5*Uelem[i+1][1][0]*uR);
            double cR=sqrt(fabs(g*pR/Uelem[i+1][0][0]));
            double MR=uR/cR;
            //cout<<ML<<" "<<MR<<endl;
            
            if (fluxscheme == 1)
            {
                vanLeer(ML, cL, Uelem[i][0][0], MR, cR, Uelem[i+1][0][0], faceflux[i], g, Uelem[i][2][0],Uelem[i+1][2][0],pL,pR);
            }
            else if (fluxscheme == 2)
            {
                stegerWarming(ML, cL, Uelem[i][0][0], MR, cR, Uelem[i+1][0][0], faceflux[i], g, Uelem[i][2][0],Uelem[i+1][2][0],pL,pR);
            }
            else if (fluxscheme == 3)
            {
                godunov(ML, cL, Uelem[i][0][0], MR, cR, Uelem[i+1][0][0], faceflux[i], g, Uelem[i][2][0],Uelem[i+1][2][0],pL,pR,t);
            }
            else// if (fluxscheme == 4 || fluxscheme == 5 ||fluxscheme)
            {
                if (i == 0)  //extrapolated BC
                {
                    higherordervanLeer(Uelem[i], Uelem[i], Uelem[i+1], Uelem[i+2], faceflux[i], g,fluxscheme, deltax, neqns);
                }
                else if (i == nelem-2)
                {
                    higherordervanLeer(Uelem[i], Uelem[i-1], Uelem[i+1], Uelem[i+1], faceflux[i], g,fluxscheme, deltax, neqns);
                }
                else
                {
                    higherordervanLeer(Uelem[i], Uelem[i-1], Uelem[i+1], Uelem[i+2], faceflux[i], g,fluxscheme, deltax, neqns);
                }
            }
            //cout<<faceflux[i][0]<<" "<<faceflux[i][1]<<" "<<faceflux[i][2]<<endl;
        }
    }
    
    else if (degree == 2)
    {
        
        for (int i=0; i<nelem-1;i++)
        {
            double UhL[neqns];
            double UhR[neqns];
            for (int j=0;j<neqns;j++)
            {
                UhL[j]=Uelem[i][j][0]+Uelem[i][j][1];
                UhR[j]=Uelem[i+1][j][0]-Uelem[i+1][j][1];
            }
            double uL=UhL[1]/UhL[0];
            double pL=(g-1)*(UhL[2]-0.5*UhL[1]*uL);
            double cL=sqrt(fabs(g*pL/UhL[0]));
            double ML=uL/cL;
            
            double uR=UhR[1]/UhR[0];
            double pR=(g-1)*(UhR[2]-0.5*UhR[1]*uR);
            double cR=sqrt(fabs(g*pR/UhR[0]));
            double MR=uR/cR;
            //cout<<ML<<" "<<MR<<endl;
            
            if (fluxscheme == 1)
            {
                vanLeer(ML, cL, UhL[0], MR, cR, UhR[0], faceflux[i], g,UhL[2],UhR[2],pL,pR);
            }
            else if (fluxscheme == 2)
            {
                stegerWarming(ML, cL, UhL[0], MR, cR, UhR[0], faceflux[i], g,UhL[2],UhR[2],pL,pR);
            }
            else if (fluxscheme == 3)
            {
                godunov(ML, cL, UhL[0], MR, cR, UhR[0], faceflux[i], g,UhL[2],UhR[2],pL,pR,t);
            }
        }
    }
    
    else if (degree == 3)
    {
        for(int i=0;i<nelem-1;i++)
        {
            double UhL[neqns];
            double UhR[neqns];
            for (int j=0;j<neqns;j++)
            {
                UhL[j]=Uelem[i][j][0]+Uelem[i][j][1]+Uelem[i][j][2];
                UhR[j]=Uelem[i+1][j][0]-Uelem[i+1][j][1]+Uelem[i+1][j][2];
            }
            double uL=UhL[1]/UhL[0];
            double pL=(g-1)*(UhL[2]-0.5*UhL[1]*uL);
            double cL=sqrt(fabs(g*pL/UhL[0]));
            double ML=uL/cL;
            
            double uR=UhR[1]/UhR[0];
            double pR=(g-1)*(UhR[2]-0.5*UhR[1]*uR);
            double cR=sqrt(fabs(g*pR/UhR[0]));
            double MR=uR/cR;
            
            if (fluxscheme == 1)
            {
                vanLeer(ML, cL, UhL[0], MR, cR, UhR[0], faceflux[i], g,UhL[2],UhR[2],pL,pR);
            }
            
            else if(fluxscheme == 2)
            {
                stegerWarming(ML, cL, UhL[0], MR, cR, UhR[0], faceflux[i], g,UhL[2],UhR[2],pL,pR);
            }
            else if(fluxscheme == 3)
            {
                godunov(ML, cL, UhL[0], MR, cR, UhR[0], faceflux[i], g,UhL[2],UhR[2],pL,pR,t);
            }
        }
    }
    
}


/***Function to calculate Boundary Integral for each Element***/
void boundintegral(vector<vector <double> > faceflux, vector< vector<vector <double> > > &rhsbound, int degree, int nelem, int neqns)
{
    if (degree == 1)
    {
        for (int i=1;i < nelem-1;i++)
        {
            for (int j=0; j< neqns;j++)
            {
                rhsbound[i][j][0]= faceflux[i][j] - faceflux[i-1][j];
            }
        }
    }
    if (degree == 2)
    {
        for (int i=1;i < nelem-1;i++)
        {
            for (int j=0; j<neqns; j++)
            {
                rhsbound[i][j][0]= faceflux[i][j] - faceflux[i-1][j];
                rhsbound[i][j][1]=faceflux[i][j] + faceflux[i-1][j];
            }
        }
    }
    
    if (degree == 3)
    {
        for (int i=1;i<nelem-1;i++)
        {
            for(int j=0;j<neqns;j++)
            {
                rhsbound[i][j][0]= faceflux[i][j] - faceflux[i-1][j];
                rhsbound[i][j][1]=faceflux[i][j] + faceflux[i-1][j];
                rhsbound[i][j][2]=faceflux[i][j] - faceflux[i-1][j];
            }
        }
    }
}

/***Function to calculate Domain Integral for each Element***/
void domainintegral(vector< vector<vector <double> > > &rhsdomain, vector< vector<vector<double> > > Uelem, int degree, int nelem, int neqns, double g,double deltax)
{
    if (degree == 1)
    {
        for (int i=1; i< nelem-1; i++)
        {
            for (int j=0;j <neqns;j++)
            {
                rhsdomain[i][j][0]=0;
            }
        }
    }
    if (degree == 2)
    {
        for (int i=1; i< nelem-1; i++)
        {
            double ua1,ua2,ub1,ub2,uc1,uc2;
            ua1=Uelem[i][0][0];
            ua2=Uelem[i][0][1];
            
            ub1=Uelem[i][1][0];
            ub2=Uelem[i][1][1];
            
            uc1=Uelem[i][2][0];
            uc2=Uelem[i][2][1];
            
            for (int j=0;j <neqns;j++)
            {
                rhsdomain[i][j][0]=0;
            }
            
            double gpoint=0;
            
            
            rhsdomain[i][0][1]=ub1*2.0;
            rhsdomain[i][1][1]=(g-1.0)*(uc1-((ub1*ub1)*(1.0/2.0))/ua1)*2.0+((ub1*ub1)*2.0)/ua1;
            rhsdomain[i][2][1]=(ub1*(uc1+(g-1.0)*(uc1-((ub1*ub1)*(1.0/2.0))/ua1))*2.0)/ua1;
        }
    }
    
    if (degree == 3)
    {
        for (int i=1; i< nelem-1; i++)
        {
            double ua1,ua2,ub1,ub2,uc1,uc2,ua3,ub3,uc3;
            ua1=Uelem[i][0][0];
            ua2=Uelem[i][0][1];
            ua3=Uelem[i][0][2];
            
            ub1=Uelem[i][1][0];
            ub2=Uelem[i][1][1];
            ub3=Uelem[i][1][2];
            
            uc1=Uelem[i][2][0];
            uc2=Uelem[i][2][1];
            uc3=Uelem[i][2][2];
            
            for (int j=0;j <neqns;j++)
            {
                rhsdomain[i][j][0]=0;
            }
            
            rhsdomain[i][0][1]=ub1*2.0;
            rhsdomain[i][1][1]=(ua1*(ub1*ub1)*9.0+ua1*(ub2*ub2)*3.0-(ua1*ua1)*uc1*6.0+(ua2*ua2)*uc1*2.0-g*ua1*(ub1*ub1)*3.0-g*ua1*(ub2*ub2)+g*(ua1*ua1)*uc1*6.0-g*(ua2*ua2)*uc1*2.0-ua2*ub1*ub2*6.0+g*ua2*ub1*ub2*2.0)/((ua1*ua1)*3.0-ua2*ua2);
            rhsdomain[i][2][1]=-((ub1-sqrt(3.0)*ub2*(1.0/3.0))*(-uc1+(g-1.0)*(-uc1+sqrt(3.0)*uc2*(1.0/3.0)+pow(ub1-sqrt(3.0)*ub2*(1.0/3.0),2.0)/(ua1*2.0-sqrt(3.0)*ua2*(2.0/3.0)))+sqrt(3.0)*uc2*(1.0/3.0)))/(ua1-sqrt(3.0)*ua2*(1.0/3.0))+((ub1+sqrt(3.0)*ub2*(1.0/3.0))*(uc1+sqrt(3.0)*uc2*(1.0/3.0)+(g-1.0)*(uc1+sqrt(3.0)*uc2*(1.0/3.0)-pow(ub1+sqrt(3.0)*ub2*(1.0/3.0),2.0)/(ua1*2.0+sqrt(3.0)*ua2*(2.0/3.0)))))/(ua1+sqrt(3.0)*ua2*(1.0/3.0));
             
                       
            rhsdomain[i][0][2]=ub2*2.0;
            rhsdomain[i][1][2]=sqrt(3.0)*sqrt(5.0)*((g-1.0)*(uc1+uc3*(2.0/5.0)-pow(ub1+ub3*(2.0/5.0)-sqrt(3.0)*sqrt(5.0)*ub2*(1.0/5.0),2.0)/(ua1*2.0+ua3*(4.0/5.0)-sqrt(3.0)*sqrt(5.0)*ua2*(2.0/5.0))-sqrt(3.0)*sqrt(5.0)*uc2*(1.0/5.0))*3.0+(pow(ub1+ub3*(2.0/5.0)-sqrt(3.0)*sqrt(5.0)*ub2*(1.0/5.0),2.0)*3.0)/(ua1+ua3*(2.0/5.0)-sqrt(3.0)*sqrt(5.0)*ua2*(1.0/5.0)))*(-1.0/9.0)+sqrt(3.0)*sqrt(5.0)*((g-1.0)*(uc1+uc3*(2.0/5.0)-pow(ub1+ub3*(2.0/5.0)+sqrt(3.0)*sqrt(5.0)*ub2*(1.0/5.0),2.0)/(ua1*2.0+ua3*(4.0/5.0)+sqrt(3.0)*sqrt(5.0)*ua2*(2.0/5.0))+sqrt(3.0)*sqrt(5.0)*uc2*(1.0/5.0))*3.0+(pow(ub1+ub3*(2.0/5.0)+sqrt(3.0)*sqrt(5.0)*ub2*(1.0/5.0),2.0)*3.0)/(ua1+ua3*(2.0/5.0)+sqrt(3.0)*sqrt(5.0)*ua2*(1.0/5.0)))*(1.0/9.0);
            rhsdomain[i][2][2]=(sqrt(1.5E1)*(ub1+ub3*(2.0/5.0)+sqrt(1.5E1)*ub2*(1.0/5.0))*(uc1+uc3*(2.0/5.0)+sqrt(1.5E1)*uc2*(1.0/5.0)+(g-1.0)*(uc1+uc3*(2.0/5.0)+sqrt(1.5E1)*uc2*(1.0/5.0)-(pow(ub1*5.0+ub3*2.0+sqrt(1.5E1)*ub2,2.0)*(1.0/1.0E1))/(ua1*5.0+ua3*2.0+sqrt(1.5E1)*ua2)))*(5.0/3.0))/(ua1*5.0+ua3*2.0+sqrt(1.5E1)*ua2)-(sqrt(1.5E1)*(ub1+ub3*(2.0/5.0)-sqrt(1.5E1)*ub2*(1.0/5.0))*(uc1+uc3*(2.0/5.0)-sqrt(1.5E1)*uc2*(1.0/5.0)+(g-1.0)*(uc1+uc3*(2.0/5.0)-sqrt(1.5E1)*uc2*(1.0/5.0)-(pow(ub1*5.0+ub3*2.0-sqrt(1.5E1)*ub2,2.0)*(1.0/5.0))/(ua1*1.0E1+ua3*4.0-sqrt(1.5E1)*ua2*2.0)))*(5.0/3.0))/(ua1*5.0+ua3*2.0-sqrt(1.5E1)*ua2);
            
        }
    }
}

void rhsvector(vector< vector<vector<double> > > &R, int nelem, vector< vector<vector<double> > > Uelem, int degree, int fluxscheme, double deltax, int neqns, double g, double deltat)
{
    
    /*****Calculate Fluxes at right face of each cell****/
    /**Note this will calculate fluxes for all elements on the right face*/
    vector<vector <double> > faceflux (nelem,vector<double> (neqns,0));
    rightfaceflux(Uelem, degree, faceflux, nelem, neqns,g, fluxscheme,deltat,deltax);   
    /**************************************************/
    
    /****Calculate Boundary Integral****/
    vector< vector<vector <double> > > rhsbound(nelem, vector<vector <double> >(neqns, vector<double>(degree,0)));    
    boundintegral(faceflux, rhsbound , degree, nelem,neqns);   
    /******************************************/
    
    
    
    /************Calculate Domain Integral**************/
    vector< vector<vector <double> > > rhsdomain(nelem, vector<vector <double> >(neqns, vector<double>(degree,0)));
    domainintegral(rhsdomain,Uelem, degree, nelem,neqns,g,deltax);
    /**************************************************/
    
    
    /**********Finally R.H.S. Vector******/
    for (int i=1; i< nelem-1;i++)
    {
        for (int k=0; k<neqns; k++)
        {
            for (int j=0; j<degree; j++)
            {
                R[i][k][j] = rhsbound[i][k][j] - rhsdomain[i][k][j];
                R[i][k][j] = -R[i][k][j];
            }
        }
    }
    /*ofstream check;
    check.open("check.txt",ios::trunc);
    for(int i=1;i<nelem-1;i++)
    {
        check<<i<<" ";
        for(int j=0;j<neqns;j++)
        {
            check<<rhsdomain[i][j][0]<<" "<<rhsdomain[i][j][1]<<" "<<rhsdomain[i][j][2]<<" ";
        }
        check<<endl;
    }
    check.close();  */
    
}

#endif	/* RHSVECTOR_H */
