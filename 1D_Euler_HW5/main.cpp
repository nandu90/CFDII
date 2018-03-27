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
#include "functions.h"





int main()
{
    
    feenableexcept(FE_INVALID | FE_OVERFLOW |FE_DIVBYZERO);
    
    time_t t1,t2;
    t1=time(0);
    
    
    /****Initializing Variables*****/
    double x0,xn;
    int nelem;
    int degree;
    int neqns;
    double cfl;
    vector<vector <double> > W;  ///Array containing rho,v,p
    int fluxscheme=0;
    double totaltime;
    int stages;
    double deltat;
    double g;
    vector<vector <double> > Wcenter;
    int BC;
    int tcontrol;
    /******************************/
    
    /*****Read Control File*****/    
    control(x0,xn,nelem,W,degree,cfl,fluxscheme,totaltime,stages,Wcenter,neqns,deltat,g,BC,tcontrol);
    /****************************/
    cout<<degree<<" "<<deltat<<" "<<tcontrol<<" "<<neqns<<" "<<BC<<" "<<fluxscheme<<" "<<g<<endl;
    double deltax=(xn-x0)/(nelem);  //Length of cell
    
    
    
    /*****Generate Ghost Cells*****/
    ghost(W,nelem,degree,neqns,Wcenter,BC);
    
    /****************************/
    
    
    
    
    /****Get element properties****/
    
    vector< vector<double> > geoel;
    geoel.resize(nelem);
    elemprop(geoel, nelem, deltax, x0);
    
    /********************************/
    
    
    /***Compute Mass Matrix****/
    vector< vector<double> > M (nelem, vector<double> (degree,0));
    massmatrix(M,geoel, nelem, degree);  
    
    /*************************/
    
    /****Construct Solution Vector at Nodes*****/
    vector<vector <double> > U (nelem+1, vector<double> (neqns,0));
    //vector<vector <double> > Leigen (neqns, vector<double> (neqns,0));
    for (int i=0;i<nelem+1;i++)
    {
        U[i][0]=W[i][0];
        U[i][1]=W[i][0]*W[i][1];
        U[i][2]=(W[i][2]/(g-1)) + 0.5*W[i][0]*pow(W[i][1],2);
    }
    /***Also Construct Solution Vector at cell centers****/
    vector<vector <double> > Ucenter (nelem, vector<double> (neqns,0));
    for (int i=0;i<nelem;i++)
    {
        Ucenter[i][0]=Wcenter[i][0];
        Ucenter[i][1]=Wcenter[i][0]*Wcenter[i][1];
        Ucenter[i][2]=(Wcenter[i][2]/(g-1)) + 0.5*Wcenter[i][0]*pow(Wcenter[i][1],2);
    }
    
    /****Now Construct Solution vector for each element****/
    vector< vector< vector<double> > > Uelem(nelem, vector< vector <double> >(neqns, vector<double> (degree,0)));
    elemVector (Uelem, U, degree, nelem,Ucenter, neqns);
    
    
    /**************************************************/
    /*ofstream check;
    check.open("check.txt",ios::trunc);
    for(int i=0;i<nelem;i++)
    {
        check<<i<<" ";
        for(int j=0;j<neqns;j++)
        {
            check<<Uelem[i][j][0]<<" "<<Uelem[i][j][1]<<" "<<Uelem[i][j][2]<<" ";
        }
        check<<endl;
    }
    check.close();*/  
    
            
    
    double t=0;
    vector< vector< vector<double> > > R(nelem, vector< vector <double> >(neqns, vector<double> (degree,0)));
    
    
    /****Time loop*****/
    while (t<=totaltime)
    {
        if(tcontrol == 2)
        {
            double tempdeltat;
            double Uh[neqns];
            for(int i=1;i<nelem-1;i++)
            {
                for(int j=0;j<neqns;j++)
                {
                   if(degree  == 1)
                   {
                       Uh[j]=Uelem[i][j][0];
                   }
                   else if(degree == 2)
                   {
                       Uh[j]=Uelem[i][j][0];
                   }
                   else if(degree == 3)
                   {
                       Uh[j]=Uelem[i][j][0]-Uelem[i][j][2]/2;
                   }
                }
                double u=Uh[1]/Uh[0];
                double p=(g-1)*(Uh[2] - 0.5*Uh[1]*u);
                double c=sqrt(fabs(g*p/Uh[0]));
                
                if(i==1)
                {
                   deltat=cfl*deltax/((2*(degree-1)+1)*(fabs(u)+c));
                }
                else
                {
                    tempdeltat=cfl*deltax/((2*(degree-1)+1)*(fabs(u)+c));
                    if(tempdeltat <= deltat)
                    {
                        deltat=tempdeltat;
                    }
                }
            }
        }
        
        if (stages == 1)
        {
        rhsvector(R, nelem, Uelem, degree, fluxscheme, deltax, neqns, g, deltat);
        euler(Uelem, M, R, deltat, nelem, degree,neqns);
        ghostupdate(Uelem,nelem,degree,neqns,BC);
        }
        
        
        
        if (stages ==2)
        {
            vector< vector< vector<double> > > Utemp(nelem, vector< vector <double> >(neqns, vector<double> (degree,0)));
            for (int i=1;i<nelem-1;i++)
            {
                for (int k=0;k<neqns;k++)
                {
                   for (int j=0;j<degree;j++)
                   {
                       Utemp[i][k][j]=Uelem[i][k][j];
                   }
                }
            }
            rhsvector(R, nelem, Uelem, degree, fluxscheme, deltax, neqns, g, deltat/2);
            euler(Uelem, M, R, deltat, nelem, degree,neqns);
            ghostupdate(Uelem,nelem,degree,neqns,BC);
            
            rhsvector(R, nelem, Uelem, degree, fluxscheme, deltax, neqns, g, deltat/2);
            euler(Uelem, M, R, deltat, nelem, degree,neqns);
            for (int i=1;i<nelem-1;i++)
            {
                for (int k=0;k<neqns;k++)
                {
                    for (int j=0;j<degree;j++)
                    {
                        Uelem[i][k][j]= 0.5*Utemp[i][k][j]+0.5*Uelem[i][k][j];
                    }
                }
            }
            ghostupdate(Uelem,nelem,degree,neqns,BC);
            
        }
        
        else if (stages == 3)
        {
            vector< vector< vector<double> > > Utemp(nelem, vector< vector <double> >(neqns, vector<double> (degree,0)));
            
            for (int i=1;i<nelem-1;i++)
            {
                for (int k=0;k<neqns;k++)
                {
                   for (int j=0;j<degree;j++)
                   {
                       Utemp[i][k][j]=Uelem[i][k][j];
                   }
                }
            }
            rhsvector(R, nelem, Uelem, degree, fluxscheme, deltax, neqns, g, deltat/3);
            euler(Uelem, M, R, deltat, nelem, degree,neqns);
            ghostupdate(Uelem,nelem,degree,neqns,BC);
            
            rhsvector(R, nelem, Uelem, degree, fluxscheme, deltax, neqns, g, deltat/3);
            euler(Uelem, M, R, deltat, nelem, degree,neqns);
            ghostupdate(Uelem,nelem,degree,neqns,BC);
            
            for (int i=1;i<nelem-1;i++)
            {
                for (int k=0;k<neqns;k++)
                {
                    for (int j=0;j<degree;j++)
                    {
                        Uelem[i][k][j]= (3.0/4.0)*Utemp[i][k][j]+(1.0/4.0)*Uelem[i][k][j];
                    }
                }
            }
            ghostupdate(Uelem,nelem,degree,neqns,BC);
            
            rhsvector(R, nelem, Uelem, degree, fluxscheme, deltax, neqns, g, deltat/3);
            euler(Uelem, M, R, deltat, nelem, degree,neqns);
            
            for (int i=1;i<nelem-1;i++)
            {
                for (int k=0;k<neqns;k++)
                {
                    for (int j=0;j<degree;j++)
                    {
                        Uelem[i][k][j]= (1.0/3.0)*Utemp[i][k][j]+(2.0/3.0)*Uelem[i][k][j];
                    }
                }
            }
            ghostupdate(Uelem,nelem,degree,neqns,BC);
        }
        cout<<t<<endl;
        t=t+deltat;
    }
    
    
    
    /*****Get the final Nodal Velocities***********/
    /***Reconstruct using original function****/
    vector< vector<vector <double> > > Wnew(nelem, vector<vector <double> > (degree, vector<double> (neqns,0)));
    
    vector< vector <double> >  Machnew(nelem, vector <double>  (degree, 0));
    vector< vector <double> >  Machold(nelem, vector <double>  (degree, 0));
    double x;
    double deltaxprint;
    
    if(degree == 1)
    {
        x=x0+deltax/2;
        deltaxprint=deltax;
        for(int i=1;i<nelem-1;i++)
        {
            vector<double> Uh(neqns,0);
            for(int j=0;j<neqns;j++)
            {
                Uh[j]=Uelem[i][j][0];
            }
            Wnew[i][0][0]=Uh[0];
                Wnew[i][0][1]=Uh[1]/Uh[0];
                Wnew[i][0][2]=(g-1)*(Uh[2]-0.5*Uh[1]*Wnew[i][0][1]);
            double c=sqrt(g*Wnew[i][0][2]/Wnew[i][0][0]);
            Machnew[i][0]=Wnew[i][0][1]/c;
            Machold[i][0]=Wcenter[i][1]/sqrt(g*Wcenter[i][2]/Wcenter[i][0]);
        }
    }
    else if(degree == 2)
    {
        x=x0;
        deltaxprint=deltax;
        for(int i=1;i<nelem-1;i++)
        {
            vector< vector<double> > Uh(degree, vector<double>(neqns,0));
            for (int j=0;j<neqns;j++)
            {
                Uh[0][j]=Uelem[i][j][0]-Uelem[i][j][1];
                Uh[1][j]=Uelem[i][j][0]+Uelem[i][j][1];
            }
            
            for(int j=0;j<degree;j++)
            {
                Wnew[i][j][0]=Uh[j][0];
                Wnew[i][j][1]=Uh[j][1]/Uh[j][0];
                Wnew[i][j][2]=(g-1)*(Uh[j][2]-0.5*Uh[j][1]*Wnew[i][j][1]);
                double c=sqrt(fabs(g*Wnew[i][j][2]/Wnew[i][j][0]));
                Machnew[i][j]=Wnew[i][j][1]/c;
            }
            Machold[i][0]=W[i][1]/sqrt(g*W[i][2]/W[i][0]);
            Machold[i][1]=W[i+1][1]/sqrt(g*W[i+1][2]/W[i+1][0]);

        }
    }
    else if(degree == 3)
    {
        x=x0;
        deltaxprint=deltax/2;
        
        for(int i=1;i<nelem-1;i++)
        {
            vector< vector<double> > Uh(degree, vector<double>(neqns,0));
            for (int j=0;j<neqns;j++)
            {
                Uh[0][j]=Uelem[i][j][0]-Uelem[i][j][1]+Uelem[i][j][2];
                Uh[1][j]=Uelem[i][j][0]-Uelem[i][j][2]/2;
                Uh[2][j]=Uelem[i][j][0]+Uelem[i][j][1]+Uelem[i][j][2];
            }
            for(int j=0;j<degree;j++)
            {
                Wnew[i][j][0]=Uh[j][0];
                Wnew[i][j][1]=Uh[j][1]/Uh[j][0];
                Wnew[i][j][2]=(g-1)*(Uh[j][2]-0.5*Uh[j][1]*Wnew[i][j][1]);
                double c=sqrt(fabs(g*Wnew[i][j][2]/Wnew[i][j][0]));
                Machnew[i][j]=Wnew[i][j][1]/c;
            }
            Machold[i][0]=W[i][1]/sqrt(g*W[i][2]/W[i][0]);
            Machold[i][1]=Wcenter[i][1]/sqrt(g*Wcenter[i][2]/Wcenter[i][0]);
            Machold[i][2]=W[i+1][1]/sqrt(g*W[i+1][2]/W[i+1][0]);
        }
    }
    
    ofstream output;
    output.open("out.txt");
    if(degree == 1)
    {
        for(int i=1;i<nelem-1;i++)
        {
            output<<x;
            for(int j=0;j<neqns;j++)
            {
                output<<" "<<Wnew[i][0][j]<<" "<<Wcenter[i][j];
            }
            output<<" "<<Machnew[i][0]<<" "<<Machold[i][0];
            output<<endl;
            x=x+deltax;
        }
    }
    
    else if(degree == 2)
    {
        for(int i=1;i<nelem-1;i++)
        {
            for(int k=0;k<degree;k++)
            {
                if(k == 0){output<<x;}
                else {output<<x+deltax;}
            
            for(int j=0;j<neqns;j++)
            {
                output<<" "<<Wnew[i][k][j]<<" "<<W[i+k][j];
            }
            output<<" "<<Machnew[i][k]<<" "<<Machold[i][k];
            output<<endl;
            
            }
            x=x+deltax;
        }
        
    }
    
    else if(degree == 3)
    {
        for(int i=1;i<nelem-1;i++)
        {
            for (int k=0;k<degree;k++)
            {
            if(k == 0){output<<x;}
                else if (k==1){output<<x+deltaxprint;}
                else if (k==2){output<<x+2*deltaxprint;}
                
            for(int j=0;j<neqns;j++)
            {
                if (k == 0)
                {
                    output<<" "<<Wnew[i][k][j]<<" "<<W[i][j];
                }
                else if (k==1)
                {
                    output<<" "<<Wnew[i][k][j]<<" "<<Wcenter[i][j];
                }
                else if (k==2)
                {
                    output<<" "<<Wnew[i][k][j]<<" "<<W[i+1][j];
                }
            }
            output<<" "<<Machnew[i][k]<<" "<<Machold[i][k];
            output<<endl;
            }
            x=x+deltax;
        }
    }
    
    output.close();
    
    
    t2=time(0);
    
    double seconds= difftime(t2,t1);
    
    cout<<"Total run time: "<<seconds<<" secs"<<endl;
    
    
}
         

