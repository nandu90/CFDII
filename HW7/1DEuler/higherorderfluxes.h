/* 
 * File:   fluxes.h
 * Author: nash
 *
 * Created on March 14, 2016, 9:19 PM
 */

#ifndef HIGHERORDERFLUXES_H
#define	HIGHERORDERFLUXES_H

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

void cons_to_prim(vector<double> &W, vector< vector<double> > U, double g)
{
    W[0]=U[0][0];
    W[1]=U[1][0]/U[0][0];
    W[2]=(g-1)*(U[2][0] - 0.5*pow(U[1][0],2)/U[0][0]);
}



void higherordervanLeer(vector< vector<double> >Ui, vector< vector<double> >Uiminus, vector< vector<double> >Uiplus, vector< vector<double> >Uiplus2, vector <double> &flux, double g,int fluxscheme, double deltax, int neqns)
{
    double rhoL,uL,cL,ML,uEL,pL;
    double rhoR,uR,cR,MR,uER,pR;
    
    if(fluxscheme == 4 || fluxscheme == 6) //based on primitive variables
    {
        
        /*Convert vectors to primitive for all 3 cells*/
        vector<double> Wi(neqns,0);
        vector<double> Wiminus(neqns,0);
        vector<double> Wiplus(neqns,0);
        
        cons_to_prim(Wi, Ui,g);
        cons_to_prim(Wiminus, Uiminus,g);
        cons_to_prim(Wiplus, Uiplus,g);
        /****************/
        vector<double> Wface(neqns,0);
        
        if (fluxscheme == 4)
        {
        for(int j=0; j< neqns; j++)
        {
            double num1=0;
            double num2=0;  
            
            num1=(Wi[j]-Wiminus[j])/deltax;
            num2=(Wiplus[j]-Wi[j])/deltax;
            Wface[j]= Wi[j] + minmod(num1,num2)*(deltax*0.5);
        }
        }
        
        else if(fluxscheme == 6)
        {
            for(int j=0; j< neqns; j++)
            {
            double num1=0;
            double num2=0;  
            double num3=0;
            
            num1=(Wi[j]-Wiminus[j])/deltax;
            num2=(Wiplus[j]-Wiminus[j])/(2*deltax);
            num3=(Wiplus[j]-Wi[j])/deltax;
            
            Wface[j]= Wi[j] + mc_limit(num1,num2,num3)*deltax/2.0;
            }
        }
        
        rhoL=Wface[0];
        uL=Wface[1];
        pL=Wface[2];
        cL=sqrt(fabs(g*pL/rhoL));
        ML=uL/cL;
        uEL=pL/(g-1) + 0.5*rhoL*uL*uL;
        
        vector<double> Wiplus2(neqns,0);
        cons_to_prim(Wiplus2, Uiplus2,g);
        
        if (fluxscheme == 4)
        {
        for(int j=0; j< neqns; j++)
        {
            double num1=0;
            double num2=0;  
            
            num1=(Wiplus[j]-Wi[j])/deltax;
            num2=(Wiplus2[j]-Wiplus[j])/deltax;
            Wface[j]= Wiplus[j] - minmod(num1,num2)*(deltax*0.5);
        }
        }
        
        else if (fluxscheme == 6)
        {
            for(int j=0; j< neqns; j++)
            {
            double num1=0;
            double num2=0;  
            double num3=0;
            
            num1=(Wiplus[j]-Wi[j])/deltax;
            num2=(Wiplus2[j]-Wi[j])/(2*deltax);
            num3=(Wiplus2[j]-Wiplus[j])/deltax;
            
            Wface[j]= Wiplus[j] - mc_limit(num1,num2,num3)*deltax/2.0;
            }
        }
        
        rhoR=Wface[0];
        uR=Wface[1];
        pR=Wface[2];
        cR=sqrt(fabs(g*pL/rhoL));
        MR=uR/cR;
        uER=pR/(g-1) + 0.5*rhoR*uR*uR;
        
    }
    
    else if (fluxscheme == 5 || fluxscheme == 7) //based on Conservative Variables
    {
        vector<double> Uface(neqns,0);
        
        if (fluxscheme == 5)
        {
        for(int j=0; j< neqns; j++)
        {
            double num1=0;
            double num2=0;  
            
            num1=(Ui[j][0]-Uiminus[j][0])/deltax;
            num2=(Uiplus[j][0]-Ui[j][0])/deltax;
            Uface[j]= Ui[j][0] + minmod(num1,num2)*(deltax*0.5);
        }
        }
        
        else if (fluxscheme == 7)
        {
            for(int j=0; j< neqns; j++)
            {
            double num1=0;
            double num2=0;  
            double num3=0;
            
            num1=(Ui[j][0]-Uiminus[j][0])/deltax;
            num2=(Uiplus[j][0]-Uiminus[j][0])/(2*deltax);
            num3=(Uiplus[j][0]-Ui[j][0])/deltax;
            
            Uface[j]= Ui[j][0] + mc_limit(num1,num2,num3)*deltax/2.0;
            }
        }
        
        rhoL=Uface[0];
        uL=Uface[1]/Uface[0];
        pL=(g-1)*(Uface[2] - 0.5*pow(Uface[1],2)/Uface[0]);
        cL=sqrt(fabs(g*pL/rhoL));
        ML=uL/cL;
        uEL=Uface[2];
        
        if(fluxscheme == 5)
        {
        for(int j=0; j< neqns; j++)
        {
            double num1=0;
            double num2=0;  
            
            num1=(Uiplus[j][0]-Ui[j][0])/deltax;
            num2=(Uiplus2[j][0]-Uiplus[j][0])/deltax;
            Uface[j]= Uiplus[j][0] - minmod(num1,num2)*(deltax*0.5);
        }
        }
        
        else if (fluxscheme == 7)
        {
            for(int j=0; j< neqns; j++)
            {
            double num1=0;
            double num2=0;  
            double num3=0;
            
            num1=(Uiplus[j][0]-Ui[j][0])/deltax;
            num2=(Uiplus2[j][0]-Ui[j][0])/(2*deltax);
            num3=(Uiplus2[j][0]-Uiplus[j][0])/deltax;
            
            Uface[j]= Uiplus[j][0] - mc_limit(num1,num2,num3)*deltax/2.0;
            }
        }
        
        rhoR=Uface[0];
        uR=Uface[1]/Uface[0];
        pR=(g-1)*(Uface[2] - 0.5*pow(Uface[1],2)/Uface[0]);
        cR=sqrt(fabs(g*pR/rhoR));
        MR=uR/cR;
        uER=Uface[2];
    }
    
    double massL,momL,energyL;
    double massR,momR,energyR;
    if (ML >= 1)
    {
        massL= rhoL*cL*ML;
        momL=rhoL*pow(cL*ML,2)+pL;
        energyL=cL*ML*(uEL+ pL);
    }
    else if (ML <= -1)
    {
        massL=0;
        momL=0;
        energyL=0;
    }
    else if (ML > -1 && ML < 1)
    {
        massL=rhoL*cL*pow(0.5*(ML+1),2);
        momL=massL*((g-1)*ML*cL+2*cL)/g;
        energyL=massL*pow((g-1)*ML*cL+2*cL, 2)/(2*(g+1)*(g-1));
        
    }
    
    
    if (MR >= 1)
    {
        massR=0;
        momR=0;
        energyR=0;
    }
    else if (MR <= -1)
    {
        massR= rhoR*cR*MR;
        momR=rhoR*pow(cR*MR,2)+pR;
        energyR=cR*MR*(uER + pR);
    }
    else if (MR > -1 && MR < 1)
    {
        massR=-rhoR*cR*pow(0.5*(MR-1),2);
        momR=massR*((g-1)*MR*cR-2*cR)/g;
        energyR=massR*pow((g-1)*MR*cR-2*cR, 2)/(2*(g+1)*(g-1));
    }
    
    
    flux[0]=massL+massR;
    flux[1]=momL+momR;
    flux[2]=energyL+energyR;
}


#endif	/*HIGHERORDERFLUXES_H*/

