/* 
 * File:   rungeKutta.h
 * Author: nash
 *
 * Created on March 2, 2016, 5:57 PM
 */

#ifndef RUNGEKUTTA_H
#define	RUNGEKUTTA_H
#include <vector>

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

#endif	/* RUNGEKUTTA_H */

