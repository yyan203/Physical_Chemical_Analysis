//
// Created by Yongjian Yang on 3/22/20.
// Neutron Scattering Analysis based on book/paper by
// Hansen & McDonald, Theory of simple liquids, p110 Eq. 4.1.29
// Binder & Kob book, <glassy materials and disordered solids: an introduction to their statistical mechanics> p50, Eq.2.39
// the above are used for calculation of the neutron structure factor F(k)
// David A. Keen, Journal of Applied Crystallography, 2001, 34, 172-177
//
#include <map>
#include <vector>
#include <utility>
#include <math.h>
#include "neutron.h"
#include "stdio.h"
#include "stdlib.h"
#include <cstdio>
//#include <xraylib.h>
#include "x-ray.h"


class scatter_analysis {
private:
    std::vector<double> Fq;   // F(Q) total-scattering structure factor
    std::vector<double> Fxq;  // F_X-Ray(Q) total-scattering structure factor of X-ray
    std::vector<double> qi;   // i.e. Q, scattering vector
    std::vector<double> Gr;   // G(r) total-radial distribution function
    std::vector<double> Gxr;  // G(r) total-radial distribution function
    std::vector<double> ri;   // r distance
    std::map< std::pair<int, int>, double* > Aij_Q; // Aij(Q), Fabe-Ziman partial structure factors
public:
    void getStructFactor(std::map< std::pair<int, int>, double* > &g, int NumberOfAtomtype,
                         int NumberOfrBins, int NumberOfqBins, double qmax, double Lmin,
                         double rho, double c[], Neutron &neu,
                         std::vector<std::tuple<double,double,double,double> > &FQ,
                         std::vector<std::tuple<double,double,double,double,double> > &GR){
        double dq = (qmax-0.0)/NumberOfqBins;
        double dr = 0.5*Lmin/NumberOfrBins;
        double R = 0.5 * Lmin; // half length of the smallest dimension
        double fij_Q;
        Xray xray;
        //for (auto i = 0;i<250;i++) printf("%lf %lf\n",0.1*i,xray.getfq(1,0.1*i));
        printf("Lmin=%lf qmax=%lf dq=%lf dr=%lf\n",Lmin,qmax,dq,dr);
        for(int i = 0; i < NumberOfqBins; i++) {
            qi.push_back((i+0.5)*dq);
            Fq.push_back(0.0);
            Fxq.push_back(0.0);
        }
        for(int i = 0; i < NumberOfrBins; i++) {
            Gr.push_back(0.0);
            Gxr.push_back(0.0);
            ri.push_back((i + 0.5) * dr);
        }

        // hold total without consider element difference
        Aij_Q[std::make_pair(-1, -1)] = new double[NumberOfqBins];
        for (int k = 0; k < NumberOfqBins; k++) {
            Aij_Q[std::make_pair(-1, -1)][k] = 0.0; // debug
            for (int l = 0; l < NumberOfrBins; l++) {
                Aij_Q[std::make_pair(-1, -1)][k] +=
                        ri[l] * (g[std::make_pair(-1, -1)][l] - 1.0) * sin(qi[k] * ri[l]) / qi[k] * sin(M_PI * ri[l] / R) / (M_PI * ri[l] / R);
                //if(k==0 && l < 201) printf("---->%d %lf %lf %lf\n",l,Aij_Q[std::make_pair(-1, -1)][k],g[std::make_pair(-1, -1)][l],qi[k]);
            }
            Aij_Q[std::make_pair(-1,-1)][k] *= dr * 4.0 * M_PI * rho;
            //printf("AijQ= %lf ", Aij_Q[std::make_pair(-1,-1)][k]);
        }
        // for computation of Gxr // see David A. Keen Jounral of Applied Crystallography, 2001, Eq.63 and 64
        double cK_total = 0.0; // SUM(ci*Ki)
        for (int i = 1;i < NumberOfAtomtype; i++) cK_total += c[i] * xray.getK(i);
        cK_total *= cK_total;

        // See the description on top of this code for the reference I used for calculation of Neutron structure factor Fq[k] using Kob's definition
        // all other calculation are based on  David A. Keen, Journal of Applied Crystallography, 2001, 34, 172-177
        double delta_fun = 0.0;
        double prefactor = 0.0; // Binder & Kob book, p50, Eq.2.39 the prefactor (actually, its inverse) in front of the summation
        for (int i = 1; i < NumberOfAtomtype + 1; i++)
         {prefactor += c[i] * neu.getb(i) * neu.getb(i);}
        prefactor = 1.0 / prefactor;

        for (int i = 1; i < NumberOfAtomtype + 1; i++)
            for (int j = 1; j < NumberOfAtomtype + 1; j++) {
                Aij_Q[std::make_pair(i, j)] = new double[NumberOfqBins];
                for (int k = 0; k < NumberOfqBins; k++) {
                    Aij_Q[std::make_pair(i, j)][k] = 0.0;

                    for (int l = 0; l < NumberOfrBins; l++) {
                        Aij_Q[std::make_pair(i, j)][k] +=
                                // NOT good! has delta_function at k->0 ri[l] * (g[std::make_pair(i,j)][l]) * sin(qi[k] * ri[l]) / qi[k] * sin(M_PI * ri[l] / R) / (M_PI * ri[l] / R) ;
                                //ri[l] * (g[std::make_pair(i,j)][l] - 1.0) * sin(qi[k] * ri[l]) / qi[k];
                                 ri[l] * (g[std::make_pair(i,j)][l] - 1.0) * sin(qi[k] * ri[l]) / qi[k] * sin(M_PI * ri[l] / R) / (M_PI * ri[l] / R) ;
                        if(isnan(Aij_Q[std::make_pair(i, j)][k]) != 0){
                            printf("%lf %lf %lf %d %d %d\n",ri[l],g[std::make_pair(i,j)][l],qi[k], i, j, l);
                            exit(0);
                        }
                        // G(r) = SUM(ci*cj*bi*bj*[gij(r) - 1])

                        // Gr[l] += c[i] * c[j] * neu.getb(i) * neu.getb(j) * (g[std::make_pair(i, j)][l] - 1);
                        // Gxr[l] += c[i] * c[j] * xray.getK(i) * xray.getK(j) / cK_total * (g[std::make_pair(i, j)][l] - 1);
                    }
                    Aij_Q[std::make_pair(i, j)][k] *= dr * 4.0 * M_PI * rho * c[i] * c[j]; // Hansen book, p110, Eq. 4.1.29
                    if (i == j) {delta_fun = 1.0;}
                    else {delta_fun = 0;}
                    Aij_Q[std::make_pair(i, j)][k] += delta_fun * c[i];  // Hansen book, p110, Eq. 4.1.29
                    Fq[k] += prefactor * neu.getb(i) * neu.getb(j) * Aij_Q[std::make_pair(i, j)][k]; // Binder & Kob book, p50, Eq.2.39

                    // for X-ray scattering
                    fij_Q = 0.0;
                    if(qi[k] <= 25){for (auto m =1;m<=NumberOfAtomtype;m++) {
                            fij_Q += c[m] * xray.getfq(m, qi[k]);
                        }
                        fij_Q *= fij_Q;
                        fij_Q = xray.getfq(i, qi[k]) * xray.getfq(j, qi[k]) / fij_Q;
                    }
                    if(qi[k] <= 25){Fxq[k] += c[i] * c[j] * fij_Q * (Aij_Q[std::make_pair(i, j)][k] - 1);}
                    if(k==5) printf("(k=5) c%d:%lf c%d:%lf f%d%d_Q:%lf f%d:%lf A%d%d_Q:%lf\n",
                                    i, c[i],j, c[j],i,j,fij_Q,j, xray.getfq(j,qi[k]),i,j,Aij_Q[std::make_pair(i, j)][k] - 1);
                }
        }
        for (auto i = 0; i < Fq.size(); i++) {
            //printf("%lf %lf %lf\n", qi[i], Fq[i], Aij_Q[std::make_pair(-1, -1)][i]);
            FQ.push_back(std::make_tuple(qi[i], Fq[i], Fxq[i], Aij_Q[std::make_pair(-1, -1)][i]));
        }
        for (auto i = 0; i < Gr.size(); i++) {
            //printf("%lf %lf %lf\n", qi[i], Fq[i], Aij_Q[std::make_pair(-1, -1)][i]);
            GR.push_back(std::make_tuple(ri[i], Gr[i], 4 * M_PI * ri[i] * rho * Gr[i], Gxr[i], 4 * M_PI * ri[i] * rho * Gxr[i])); // last two terms are D(r) for neutron and D(r) for X-ray
        }
    }
};


