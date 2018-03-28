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
#include "exprtk.hpp"
#include <vector>

using namespace std;

template <typename T>
void evaluate(string expression_string, double xvalue, double &u)
{
   typedef exprtk::symbol_table<T> symbol_table_t;
   typedef exprtk::expression<T>     expression_t;
   typedef exprtk::parser<T>             parser_t;
   
   T x = T(xvalue);
   
   symbol_table_t symbol_table;
   symbol_table.add_variable("x",x);
   symbol_table.add_constants();
   
   expression_t expression;
   expression.register_symbol_table(symbol_table);

   parser_t parser;
   parser.compile(expression_string,expression);
   
   T result = expression.value();
   
   
   
   u=result;
}

void control(double &x0, double &xn, int &ncells, vector<double> &u, int &degree, double &cfl, double &a, int &fluxscheme, double &totaltime, int &stages, vector<double> &ucenter, int &reconstruct)
{
    
    
    
    vector<string> expression;
    vector<double> xbreak;
    int divisions;
    
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
            
            getline(controlfile,line);
            istringstream ss2(line);
            ss2>>divisions;
            expression.resize(divisions);
            xbreak.resize(divisions-1);
            
            for(int i=0;i<divisions;i++)
            {
            getline(controlfile,line);
            istringstream ss3(line);
            ss3>>expression[i]>>junk2>>junk2;
            if(i < divisions-1)
            {
                ss3>>xbreak[i];
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
        }
        else if(junk.compare("a") == 0)
        {
            ss>>a;
        }
        else if(junk.compare("Upwind") == 0)
        {
            fluxscheme=1;
        }
        else if(junk.compare("Lax-Friedrichs") == 0)
        {
            fluxscheme=2;
        }
        else if(junk.compare("CentralDifference") == 0)
        {
            fluxscheme=3;
        }
        else if(junk.compare("totaltime") == 0)
        {
            ss>>totaltime;
        }
        else if(junk.compare("Stages") == 0)
        {
            ss>>stages;
        }
        
        else if(junk.compare("Fromm") == 0)
        {
            reconstruct=0;
        }
        else if(junk.compare("Beam-Warming") == 0)
        {
            reconstruct=-1;
        }
        else if(junk.compare("Lax-Wendroff") == 0)
        {
            reconstruct=1;
        }
        else if(junk.compare("Minmod") == 0)
        {
            reconstruct=3;
        }
        else if(junk.compare("MC_limiter") == 0)
        {
            reconstruct=4;
        }
            
        
    }
    controlfile.close();
    
    u.resize(ncells+1);
    int k=0;
    
    double deltax=(xn-x0)/ncells;
    double x=x0;
    
    /***Find Nodal Values at grid point using function specified in Control File***/
    for (int i=0;i<ncells+1;i++)
    {
        if (k < divisions-1)
        {
        if(x <= xbreak[k])
        {
            evaluate<double>(expression[k],x,u[i]);
        }
        else
        {
            k++;
        }
        }
        else
        {
            evaluate<double>(expression[k],x,u[i]);
        }
        x=x+deltax;        
    }
    
    ucenter.resize(ncells+2);
    x=x0+deltax/2;
    k=0;
    
    /**Find Cell Center Values using function specified in Control File*/
    for (int i=1; i< ncells+1;i++)
    {
        if (k < divisions-1)
        {
        if(x <= xbreak[k])
        {
            evaluate<double>(expression[k],x,ucenter[i]);
        }
        else
        {
            k++;
        }
        }
        else
        {
            evaluate<double>(expression[k],x,ucenter[i]);
        }
        x=x+deltax;
    }
    
}


#endif	/* CONTROLINPUT_H */

