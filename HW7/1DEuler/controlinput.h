/* 
 * File:   controlinput.h
 * Author: Nadish Saini
 *
 * Created on March 1, 2016, 12:40 PM
 */

#ifndef CONTROLINPUT_H
#define	CONTROLINPUT_H
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include "evalexpression.h"

using namespace std;



void control(double &x0, double &xn, int &ncells, vector<vector <double> > &W, int &degree, double &cfl, int &fluxscheme, double &totaltime, int &stages, vector<vector <double> > &Wcenter, int &neqns, double &deltat, double &g, int &BC, int &tcontrol)
{
    
    
    
    vector<vector <string> > expression;
    vector<vector <double> > xbreak;
    vector<vector <bool> > xbreakleft;
    vector<int> divisions;
    
    ifstream controlfile;
    controlfile.open("control.txt");
    string line;
    while(getline(controlfile,line))
    {
        string junk;
        string junk2;
        istringstream ss(line);
        ss>>junk; 
        if(junk.compare("Domain") == 0)
        {
            getline(controlfile,line);
            istringstream ss2(line);
            ss2>>junk2>>x0;   
            
            getline(controlfile,line);
            istringstream ss3(line);
            ss3>>junk2>>xn; 
            
            getline(controlfile,line);
            istringstream ss4(line);
            ss4>>junk2>>ncells; 
        }
        else if(junk.compare("InitialConditions") == 0)
        {   
            ss>>neqns;
            divisions.resize(neqns);
            expression.resize(neqns);
            xbreak.resize(neqns);
            xbreakleft.resize(neqns);
            while(getline(controlfile,line))
            {
                //cout<<"here"<<endl;
                istringstream ss2(line);
                ss2>>junk2;
                if (junk2.compare("End_Initial_Conditions") == 0)
                {
                    break;
                }
                else
                {
                    bool read=false;
                    int div=0;
                    if(junk2.compare("rho") == 0)
                    {
                         div=0;
                         read=true;
                    }
                    else if(junk2.compare("p") == 0)
                    {
                         div=2;
                         read=true;
                    }
                    else if(junk2.compare("ux") ==0)
                    {
                        div=1;
                        read=true;
                    }
                    
                    if(read==true)
                    {
                    getline(controlfile,line);
                    istringstream ss3(line);
                    ss3>>divisions[div];
                    expression[div].resize(divisions[div]);
                    xbreak[div].resize(divisions[div]-1);
                    xbreakleft[div].resize(divisions[div]-1);
                    for(int i=0;i<xbreakleft.size();i++)
                    {
                        xbreakleft[div][i]=false;
                    }
                    
                    for (int i=0; i<divisions[div];i++)
                    {
                        getline(controlfile,line);
                        istringstream ss4(line);
                        ss4>>expression[div][i];
                        if(i < divisions[div]-1)
                        {
                            ss4>>junk2>>junk2>>junk2>>xbreak[div][i]>>junk2;
                            if(junk2.compare("]") == 0)
                            {
                                xbreakleft[div][i]=true;
                            }
                        }
                    }
                    }
                }
                
            }
            
        }
        else if(junk.compare("PolyDegree") == 0)
        {
            ss >> degree;
        }
        else if(junk.compare("CFL") == 0)
        {
            ss >> cfl;
            tcontrol = 2;
        }
        else if(junk.compare("Van-Leer") == 0)
        {
            fluxscheme=1;
        }
        else if(junk.compare("Steger-Warming") == 0)
        {
            fluxscheme=2;
        }
        else if(junk.compare("Godunov") == 0)
        {
            fluxscheme=3;
        }
        else if(junk.compare("Minmod-Primitive") == 0)
        {
            fluxscheme=4;
        }
        else if(junk.compare("Minmod-Conservative") == 0)
        {
            fluxscheme=5;
        }
        else if(junk.compare("MC-Primitive") == 0)
        {
            fluxscheme=6;
        }
        else if(junk.compare("MC-Conservative") == 0)
        {
            fluxscheme=7;
        }
        else if(junk.compare("totaltime") == 0)
        {
            ss>>totaltime;
        }
        else if(junk.compare("Stages") == 0)
        {
            ss>>stages;
        }
        else if(junk.compare("timestep") == 0)
        {
            ss>>deltat;
            tcontrol=1;
        }
        else if(junk.compare("gamma") == 0)
        {
            ss>>g;
        }
        else if(junk.compare("periodic") == 0)
        {
            BC=1;
        }
        else if(junk.compare("extrapolated") ==0)
        {
            BC=2;
        }
            
        
    }
    controlfile.close();
    
    W.resize(ncells+1);
    Wcenter.resize(ncells+2);
    for(int i=0;i<ncells+1;i++)
    {
        W[i].resize(neqns);
    }
    for (int i=0;i<ncells+2;i++)
    {
        Wcenter[i].resize(neqns);
    }
    
    //cout<<divisions[0]<<" "<<divisions[1]<<" "<<divisions[2]<<endl;
    int k=0;
    
    double deltax=(xn-x0)/ncells;
    double x=x0;
    
    /***Find Nodal Values at grid point using function specified in Control File***/
    for (int j=0;j<neqns;j++)
    {
        x=x0;
        k=0;
        for (int i=0;i<ncells+1;i++)
        {
        if (k < divisions[j]-1)
        {
        if(x < xbreak[j][k])
        {
            evaluate<double>(expression[j][k],x,W[i][j]);
        }
        else if(fabs(x-xbreak[j][k]) < 1e-5)
        {
            if(xbreakleft[j][k] == true)
            {
                evaluate<double>(expression[j][k],x,W[i][j]);
                k++;
            }
            else if (xbreakleft[j][k] == false)
            {
                k++;
                evaluate<double>(expression[j][k],x,W[i][j]);
            }
        }
        else
        {
            k++;
            evaluate<double>(expression[j][k],x,W[i][j]);
        }
        }
        else
        {
            evaluate<double>(expression[j][k],x,W[i][j]);
        }
        x=x+deltax;        
        }
    }
    
    
    /**Find Cell Center Values using function specified in Control File*/
    
    for (int j=0; j< neqns;j++)
    {
        x=x0+deltax/2;
        k=0;
        for (int i=1; i< ncells+1;i++)
        {
        if (k < divisions[j]-1)
        {
        if(x < xbreak[j][k])
        {
            evaluate<double>(expression[j][k],x,Wcenter[i][j]);
        }
        else if(fabs(x-xbreak[j][k]) < 1e-5)
        {
            if(xbreakleft[j][k] == true)
            {
                evaluate<double>(expression[j][k],x,Wcenter[i][j]);
                k++;
            }
            else if (xbreakleft[j][k] == false)
            {
                k++;
                evaluate<double>(expression[j][k],x,Wcenter[i][j]);
            }
        }
        else
        {
            k++;
            evaluate<double>(expression[j][k],x,Wcenter[i][j]);
        }
        }
        else
        {
            evaluate<double>(expression[j][k],x,Wcenter[i][j]);
        }
        x=x+deltax;
        }
    }
    
}


#endif	/* CONTROLINPUT_H */

