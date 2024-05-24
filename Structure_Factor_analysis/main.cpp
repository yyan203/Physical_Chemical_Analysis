/****************************************************
 * calculate radial distribution function g(r)
 * of ensemble configurations in "Lata4olivia file"
 * Yongjian Yang, Penn State University, 2018
 *****************************************/
#include "stdio.h"
#include "stdlib.h"
#include <math.h>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <utility>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>
#include <map>
#include <tuple>


#include "scatter_analysis.h"


#define MAXConfig 10000
#define MAXEqui 10000
#define MAXNumberOfParticles 1114200
#define MAXNumberOfBins 10000

#define SQR(x) ((x)*(x))
#define CUB(x) ((x)*(x)*(x))

int main(int argc, char **argv)
{
    char inputfilename[100], outputfilename[100];
    int NumberOfConfig, NumberOfEqui, NumberOfEnd, NumberOfParticles, NumberOfAtomtype;
    int SampleNumber;
    int snapshot,line;
    //double L;
    double Lx,Ly,Lz,Lmin;
    double xlo,xhi,ylo,yhi, zlo, zhi;
    int id;  // atom id
    int nmols;  // atom id
    double q; // charge
    //float ix,iy,iz;
    //int itype;
    double rho;
    double qmax;
    double dr;
    //double ** g[MAXNumberOfBins];
    int NumberOfBins, NumberOfqBins;
    int i,j,box_length_flag;
    double dx,dy,dz;
    double rij;
    char AtomID;
    std::ifstream fpxyz;
    std::string str;
    std::stringstream ss;
    FILE *fpoutput;
    FILE *ffq;
    FILE *ggr; // G(r) and D(r)

    if(argc != 12) {
        // Provide the correct syntax of this analysis code
        printf("Correct syntax: exe inputdumpfile outputfile NumConfig NumEqui(start) NumEnd(end) NumAtom NumAtomType NumOfrBins NumOfqBins qmax box_length_flag<1 or 0>\n");
        printf("NumEqui<start_frame> NumEnd<end> NumBin<r Bin> box_length_flag<1, use Lx,Ly,Lz at NumEqui, else update Lx,y,z every frame>\n");
        printf("Assuming Particles of each species are not decreasing or increasing, otherwise please change c[i]\n");
        printf("Check lata4olivia file, 6-8 lines which contain xlo, xly, xlz should only contain two columns, sometimes, these file may have 3 column\n");
        exit(0);
    }
    else {
        //Note that argv[0] is "md.analysis.cxx"
        sscanf(argv[1], "%s", inputfilename);
        sscanf(argv[2], "%s", outputfilename);
        sscanf(argv[3], "%d", &NumberOfConfig);
        sscanf(argv[4], "%d", &NumberOfEqui);
        sscanf(argv[5], "%d", &NumberOfEnd);
        sscanf(argv[6], "%d", &NumberOfParticles);
        sscanf(argv[7], "%d", &NumberOfAtomtype);
        sscanf(argv[8], "%d", &NumberOfBins);
        sscanf(argv[9], "%d", &NumberOfqBins);
        sscanf(argv[10], "%lf", &qmax); // maximum Q value
        sscanf(argv[11], "%d", &box_length_flag); // if 1, only use length at NumberOfEqui, else, update Lx,Ly,Lz every frame
        //scanf(argv[7], "%f", &Lx);
        //scanf(argv[8], "%f", &Ly);
        //scanf(argv[9], "%f", &Lz);
    }
    int PairNum = NumberOfAtomtype * NumberOfAtomtype - NumberOfAtomtype;
    //double* x[NumberOfAtomtype + 1],*y[NumberOfAtomtype + 1],*z[NumberOfAtomtype + 1];
    double* x,* y,* z;
    x = new double[MAXNumberOfParticles];
    y = new double[MAXNumberOfParticles];
    z = new double[MAXNumberOfParticles];
    int * type;
    type = new int[MAXNumberOfParticles];
    // double x[MAXNumberOfParticles],y[MAXNumberOfParticles],z[MAXNumberOfParticles];
    // int type[MAXNumberOfParticles];

    double c[NumberOfAtomtype+1]; // composition
    double rho_all[NumberOfAtomtype+1]; //number density of each species
    for (auto i=0;i<=NumberOfAtomtype;i++){c[i]=0;rho_all[i]=0;}

    std::map< std::pair<int, int>, double* > g;
    g[std::make_pair(-1, -1)] = new double[MAXNumberOfBins];
    for (int i = 1; i < NumberOfAtomtype + 1; i++) {
        for (int j = 1; j < NumberOfAtomtype + 1; j++) {
            g[std::make_pair(i, j)] = new double[MAXNumberOfBins];
            for (int k = 0; k < MAXNumberOfBins; k++) {
                g[std::make_pair(i, j)][k] = 0.0;
                g[std::make_pair(-1,-1)][k] = 0.0; // for debug, used to store total G(r)
                //if (i==1 && j==3 && k==0) printf("g=%lf-----------\n",g[std::make_pair(1, 3)][0] = 0.0);
            }
        }
    }
    rho = 0.0;
    //rho = NumberOfParticles/CUB(L);
     Lmin = 10000;
    SampleNumber = 0;
    fpxyz.open(inputfilename);
    if (fpxyz.is_open()) {printf("file is opened successfully\n");}
    else {printf("cannot open file\n"); exit(0);}
    for(snapshot=0;snapshot<NumberOfConfig;snapshot++)
    {
        //for(line=0;line<=NumberOfParticles;line++) // N+1 lines
        for(line=0;line<=NumberOfParticles+8;line++) // N+2 lines
        {
            getline(fpxyz, str);
            std::stringstream ss(str);
            if(line == 3) {
                ss >> nmols; // printf("%d atoms in snapshot %d\n",nmols, snapshot);
            } // getline(fpxyz, str);
            else if(line == 5)
                    {ss >> xlo >> xhi;} // fscanf(fpxyz,"%lf %lf\n",&xlo,&xhi);
            else if(line == 6)
                    {ss >> ylo >> yhi;} // fscanf(fpxyz,"%lf %lf\n",&ylo,&yhi);
            else if(line == 7)
                    {ss >> zlo >> zhi; } // fscanf(fpxyz,"%lf %lf\n",&zlo,&zhi);
            else if(line > 8){
                ss>>id>>type[line-9]>>q>>x[line-9]>>y[line-9]>>z[line-9]; // fscanf(fpxyz,"%*s %d %*lf %lf %lf %lf\n",&type[line-9],&x[line-9],&y[line-9],&z[line-9]);
                if(snapshot == NumberOfEqui) c[type[line-9]]++;
                if(type[line-9]<0 || type[line-9]>5){printf("snapshot:%d line:%d atom type %d wrong (not between 0 and 5)\n",snapshot,line,type[line-9]); exit(0);}
            }
            else
                ; // fscanf(fpxyz,"%*[^\n]\n",NULL);
        }
        if(box_length_flag==0){Lx = xhi - xlo; Ly = yhi - ylo; Lz = zhi - zlo;}
        else {
            if(snapshot == NumberOfEqui) {Lx = xhi - xlo; Ly = yhi - ylo; Lz = zhi - zlo;}
        }
        //Lx = 31.21; Ly = 31.54; Lz = 37.365;
        printf("Number of particles type %d in the system\n",NumberOfAtomtype);
        if(snapshot == NumberOfEqui) {
            printf("Lx %lf Ly %lf Lz %lf\n",Lx, Ly, Lz);
            for (auto i=1;i<=NumberOfAtomtype;i++){
                rho_all[i]=c[i]/Lx/Ly/Lz;
                c[i]/=NumberOfParticles;
            printf("type %d count %lf rho %lf\n",i,c[i] * NumberOfParticles/8,rho_all[i]);
        }
    }
    if(snapshot >= NumberOfEqui && snapshot <= NumberOfEnd) {
        rho += NumberOfParticles / (Lx * Ly * Lz); /////////////////////
        //rho += NumberOfParticles/(31.21 * 31.54 * 37.365);
        Lmin = std::min(std::min(Lx, Ly), Lz);
        dr = 0.5 * Lmin / NumberOfBins;
        SampleNumber++;
        for (i = 0; i < NumberOfParticles; i++){
            //printf("Frame %d: %lf percent finished!\n",snapshot,double(i)/ double(NumberOfParticles) * 100.00);
            for (j = i + 1; j < NumberOfParticles; j++) {
                dx = x[i] - x[j];
                dx = dx - Lx * round(dx / Lx);
                dy = y[i] - y[j];
                dy = dy - Ly * round(dy / Ly);
                dz = z[i] - z[j];
                dz = dz - Lz * round(dz / Lz);
                rij = sqrt(SQR(dx) + SQR(dy) + SQR(dz));
                if (rij < 0.5 * Lmin) {
                    if(type[i]>5 || type[j]>5 || type[i]<0 || type[j]<0) printf("%d %d %d %d %f %f\n",i,j,type[i],type[j],rij,dr);
                    g[std::make_pair(type[i], type[j])][(int) (rij / dr)] += 1.0;
                    g[std::make_pair(type[j], type[i])][(int) (rij / dr)] += 1.0;
                    g[std::make_pair(-1, -1)][(int) (rij / dr)] += 2.0;
                }
            }
        }
    }
}
fpxyz.close();

printf("density %lf SampleNumber %d\n",rho,SampleNumber);
rho/=double(SampleNumber); // average num density

fpoutput = fopen(outputfilename,"w");
// fprintf(fpoutput, "r = %lf\tg(r) = %lf\n", (i + 0.5) * dr, g[i] / NumberOfParticles);
fprintf(fpoutput, "#r g(r):");
int index=2;
for (int i = 1; i < NumberOfAtomtype + 1; i++)
    for (int j = 1; j < NumberOfAtomtype + 1; j++) {
        fprintf(fpoutput, "\t%d-%d(%d)", i,j,index);
        index++;
    }
fprintf(fpoutput, "\tTotal_G(r)(%d)", index);
fprintf(fpoutput, "\n");

for(auto k=0;k<NumberOfBins;k++){
    fprintf(fpoutput, "%lf", (k + 0.5) * dr);
    g[std::make_pair(-1, -1)][k] /= (SampleNumber * 4.0 * M_PI / 3.0 * (CUB(k + 1) - CUB(k)) * CUB(dr) * rho); // Total G(r), for debug
    for (int i = 1; i < NumberOfAtomtype + 1; i++)
        for (int j = 1; j < NumberOfAtomtype + 1; j++) {
            g[std::make_pair(i, j)][k] /= SampleNumber;
            if(rho_all[j] !=0 ) g[std::make_pair(i, j)][k] /=
                                        4.0 * M_PI / 3.0 * (CUB(k + 1) - CUB(k)) * CUB(dr) * rho_all[j];
            if (c[i] > 0 && rho_all[j]){
                    //fprintf(fpoutput, "\t%lf", g[std::make_pair(i, j)][k] / c[i]);
                    g[std::make_pair(i, j)][k] /= (NumberOfParticles * c[i]);
                    fprintf(fpoutput, "\t%lf", g[std::make_pair(i, j)][k]);
                }
                else
                   fprintf(fpoutput, "\t%lf",g[std::make_pair(i, j)][k]);
            }
        g[std::make_pair(-1, -1)][k] /= NumberOfParticles;
        fprintf(fpoutput, "\t%lf", g[std::make_pair(-1, -1)][k]);
        fprintf(fpoutput, "\n");
    }
    fclose(fpoutput);
    printf("density %lf\n",rho);
    printf("sampling number %d\n",SampleNumber);
    if(box_length_flag==1)
        printf("use fixed box length\n");
    else
        printf("use updated box length every frame\n");

    printf("*******************************************************\n");

    if(isnan(g[std::make_pair(1, 2)][0]) != 0){
        printf("g wrong 2%lf\n",g[std::make_pair(1, 3)][0]);
        exit(0);
    }

    scatter_analysis analysis;
    Neutron neu;
    std::vector<std::tuple<double,double,double,double> > Fq;
    std::vector<std::tuple<double,double,double,double,double> > Gr;
    analysis.getStructFactor(g,NumberOfAtomtype,NumberOfBins,NumberOfqBins,qmax,Lmin,rho,c,neu,Fq,Gr);
    ffq = fopen("Fq.dat","w");
    fprintf(ffq, "#q\tF(Q)-neutron\tF(Q)-Xray\tA_ij(Q)_total(without considering different species)\n");
    for (auto i=0; i < Fq.size();i++) fprintf(ffq, "%lf\t%lf\t%lf\t%lf\n",std::get<0>(Fq[i]),std::get<1>(Fq[i]),std::get<2>(Fq[i]),std::get<3>(Fq[i]));
    fclose(ffq);
    // print G(r) and D(r) differential correlation function
    ggr = fopen("Gr.dat","w");
    fprintf(ggr, "#r\tG(r)\tD(r)(neutron)\tGx(r)\tDx(r)(X-ray)\n");
    for (auto i=0; i < Gr.size();i++) fprintf(ggr, "%lf\t%lf\t%lf\t%lf\t%lf\n",std::get<0>(Gr[i]),std::get<1>(Gr[i]),std::get<2>(Gr[i]),std::get<3>(Gr[i]),std::get<4>(Gr[i]));
    fclose(ggr);
    return 0;
}
