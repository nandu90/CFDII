/* 
 * File:   fluxes.h
 * Author: nash
 *
 * Created on March 14, 2016, 9:19 PM
 */

#ifndef FLUXES_H
#define	FLUXES_H

using namespace std;
void stegerWarming(double ML, double cL, double rhoL, double MR, double cR, double rhoR, vector <double> &flux, double g, double uEL, double uER, double pL, double pR)
{
    double massL,momL,energyL;
    double massR,momR,energyR;
    double uL=ML*cL;
    double uR=MR*cR;
    
    if (ML >= 1)
    {
        massL= rhoL*cL*ML;
        momL=rhoL*pow(cL*ML,2)+pL;
        energyL=cL*ML*(uEL+ pL);
    }
    else if (ML <1 && ML >= 0)
    {
        massL=((g-1)/g)*rhoL*uL + (rhoL/(2*g))*(uL+cL);
        momL=((g-1)/g)*rhoL*uL*uL+ (rhoL/(2*g))*(uL+cL)*(uL+cL);
        energyL= ((g-1)/g)*rhoL*uL*pow(uL,2)/2 + (rhoL/(2*g))*(uL+cL)*(uL*uL/2 + cL*cL/(g-1) + cL*uL);
    }
    else if (ML < 0 && ML > -1)
    {
        massL=(rhoL/(2*g))*(uL+cL);
        momL=massL*(uL+cL);
        energyL=massL*(uL*uL/2 + cL*cL/(g-1) + cL*uL);
    }
    else if(ML <= -1)
    {
        massL=0;
        momL=0;
        energyL=0;
    }
    
    if (MR >= 1)
    {
        massR= 0;
        momR=0;
        energyR=0;
    }
    else if (MR <1 && MR >= 0)
    {
        massR=(rhoR/(2*g))*(uR-cR);
        momR=massR*(uR-cR);
        energyR=massR*(uR*uR/2 + cR*cR/(g-1) - cR*uR);
    }
    else if (MR < 0 && MR > -1)
    {
        massR=((g-1)/g)*rhoR*uR + (rhoR/(2*g))*(uR-cR);
        momR=((g-1)/g)*rhoR*uR*uR+ (rhoR/(2*g))*(uR-cR)*(uR-cR);
        energyR= ((g-1)/g)*rhoR*uR*pow(uR,2)/2 + (rhoR/(2*g))*(uR-cR)*(uR*uR/2 + cR*cR/(g-1) - cR*uR);
    }
    else if(MR <= -1)
    {
        massR= rhoR*cR*MR;
        momR=rhoR*pow(cR*MR,2)+pR;
        energyR=cR*MR*(uER+ pR);
    }
    
    
    flux[0]=massL+massR;
    flux[1]=momL+momR;
    flux[2]=energyL+energyR;
}

void vanLeer(double ML, double cL, double rhoL, double MR, double cR, double rhoR, vector <double> &flux, double g, double uEL, double uER, double pL, double pR)
{
    double massL,momL,energyL;
    double massR,momR,energyR;
    double uL=ML*cL;
    double uR=MR*cR;
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



void func_and_diff(double &func , double &diff, double rhol, double rhor, double vl, double vr, double pl, double pr, double g, double P, double cl, double cr)
{    
    
    func = 1-(vr-vl + (cr*(P-1)/g)/(sqrt(fabs(1+(g+1)*(P-1)/(2*g)))))*((g-1)/(2*cl))-pow(fabs(P),(g-1)/(2*g))*pow(fabs(pr/pl),(g-1)/(2*g));
    diff = - ((cr/pow(fabs(g*(((P - 1)*(g + 1))/(2*g) + 1)), 0.5) - (cr*(P - 1)*(g + 1))/(4*g*g*pow(fabs((((P - 1)*(g + 1))/(2*g) + 1)),3/2)))*(g - 1))/(2*cl) - (pow(fabs(P),((g - 1)/(2*g) - 1))*pow(fabs(pr/pl),((g - 1)/(2*g)))*(g - 1))/(2*g);
}

void region4(double &v4, double &c4, double &p4, double &rho4, double x, double cl, double g, double t, double vl,double pl,double x0)
{
    
    v4=(1/(g+1))*(vl*(g-1)+2*cl+2*(x-x0)/t);
    c4=v4-(x-x0)/t;
    p4=pl*pow(c4/cl,2*g/(g-1));
    rho4=g*p4/pow(c4,2);
}

void godunov(double ML, double cL, double rhoL, double MR, double cR, double rhoR, vector <double> &flux, double g, double uEL, double uER, double pL, double pR, double t)
{
    double uL=ML*cL;
    double uR=MR*cR;
    
    
    double Pi=(pL+pR)/(2*pR);
    
    if (pR < 0)
    {
        pR=0.01;
    }
    if (pR > 1)
    {
        pR=1;
    }
    if(pL > 1)
    {
        pL=1;
    }
    if(pL < 0)
    {
        pL=0.01;
    }
    double P;
    double err;
    do
    {
        double func,diff;
        func_and_diff(func,diff,rhoL,rhoR,uL,uR,pL,pR,g,Pi,cL,cR);
        P = Pi - (func/diff);
        err=fabs(P-Pi);
        Pi = P;
    }while (err/Pi > 0.00001);
    
    double x0=0;
    
    /********Region 2 solution********/
    double p2=P*pR;
    double v2=uR+(cR/g)*(P-1)/sqrt(fabs(1+(g+1)*(P-1)/(2*g)));
    double alpha=(g+1)/(g-1);
    double rho2=rhoR*(1+alpha*P)/(alpha+P);
    double c2=sqrt(fabs(g*p2/rho2));
    
    /*******Region 3 solution*********/
    double p3=p2;
    double v3=v2;
    double rho3=rhoL*pow(fabs(p3/pL),(1/g));
    double c3=sqrt(fabs(g*p3/rho3));
    
    /*****Shock speed******/
    double vs;
    if(v2 == uR)
    {
        vs=v2;
    }
    else
    {
        vs=uR+(P-1)*pow(cR,2)/(g*(v2-uR));
    }
    
    /*****Shock location****/
    double xs=vs*t+x0;
    /*****Contact discontinuity location ******/
    double xc=v2*t+x0;
    /*****Expansion wave end location******/
    double x4=(v3-c3)*t+x0;
    /*****Left region location****(Expansion wave beginning)*/
    double xl=(uL-cL)*t+x0;
    
    double rho,u,c,p;
    
    if(xl>=x0)
    {
        rho=rhoL;
        u=uL;
        c=cL;
        p=pL;
    }
    else if(xl < x0 && x4 >= x0)
    {
        region4(u, c, p, rho, x0, cL, g, t, uL,pL,x0);
    }
    else if(x4 < x0 && xc >= x0)
    {
        u=v3;
        rho=rho3;
        p=p3;
        c=c3;
    }
    else if(xc < x0 && xs >= x0)
    {
        u=v2;
        rho=rho2;
        p=p2;
        c=c2;
    }
    else if(xs < x0)
    {
        u=uR;
        rho=rhoR;
        p=pR;
        c=cR;
    }
    
    double rhoE= p/(g-1) + 0.5*rho*u*u;
    flux[0]=rho*u;
    flux[1]=rho*u*u + p;
    flux[2]=u*(rhoE + p);
}


#endif	/* FLUXES_H */

