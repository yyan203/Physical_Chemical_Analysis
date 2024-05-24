//
// Created by Yongjian Yang on 3/23/18.
//

#ifndef STRUCTURE_FACTOR_3_X_RAY_H
#define STRUCTURE_FACTOR_3_X_RAY_H
#include <map>
#include <math.h>

class Xray {
private:
    std::map< int,double* > a; // Xray scattering factor table, from table in http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php
    std::map< int,double* > b; // Xray scattering factor table
    std::map< int,double > c; // Xray scattering factor table
    std::map< int,double > K; // K effective number of species, approximatly equal to atomic number Z

public:
    Xray (){
        /* debug */
        // for Si-O system
        K[1] = 8.0 ; // O
        K[2] = 14.0 ; // Si

        // K[1] = 6.0 ; // C
        // K[2] = 1.0 ; // H
        // K[3] = 8.0 ; // O
        // K[4] = 7.0 ; // N
        // K[5] = 30.0 ; // Zn

        for (int i=0;i<=5;i++) {
            a[i] = new double[6];
            b[i] = new double[6];
        }
        // C H O N Zn2+
        // O Si
        a[1][1]=3.0485	 ;b[1][1]=13.2771	;a[1][2]= 2.2868  ;b[1][2]=5.7011	;a[1][1]=1.5463	;b[1][1]=0.3239 ;a[1][4]=0.867 ;b[1][4]=32.9089;c[1]=0.2508;
        a[2][1]=6.2915   ;b[2][1]=2.4386 	;a[2][2]= 3.0353  ;b[2][2]=32.3337  ;a[2][3]=1.9891 ;b[2][3]=0.6785 ;a[2][4]=1.541 ;b[2][4]=81.6937;c[2]=1.1407;
        // a[1][1]=2.31 ;b[1][1]=20.8439	;a[1][2]= 1.02	  ;b[1][2]=10.2075 ;a[1][3]=1.5886	;b[1][3]=0.5687 ;a[1][4]=0.865   ;b[1][4]=51.6512;c[1]=0.2156;
        // a[2][1]=0.489918 ;b[2][1]=20.6593	;a[2][2]= 0.262003;b[2][2]=7.74039;a[2][3]=0.196767 ;b[2][3]=49.5519 ;a[2][4]=0.049879 ;b[2][4]=2.20159;c[2]=0.001305;
        // a[4][1]=12.2126 ;b[4][1]=0.0057   ;a[4][2]= 3.1322  ;b[4][2]=9.8933	;a[4][3]=2.0125	;b[4][3]=28.9975 ;a[4][4]=1.1663   ;b[4][4]=0.5826;c[4]=-11.529;
        // a[5][1]=11.9719 ;b[5][1]=2.9946   ;a[5][2]= 7.3862  ;b[5][2]=0.2031	;a[5][3]=6.4668	;b[5][3]=7.0826 ;a[5][4]=1.394   ;b[5][4]=18.0995;c[5]=0.7807;

    }
    // e is element and q is scattering vector
    double getfq(int e,double q) {
        if (e < 1 || e > 5) {printf("index error!"); exit(0);}
        double fq=0.0;
        for (int i=1;i<5;i++) fq += a[e][i] * exp(-1*b[e][i] * (q / M_PI / 4) * (q / M_PI / 4));
        fq += c[e];
        return fq;
    }
    double getK(int e) {
       return K[e];
    }
};

#endif //STRUCTURE_FACTOR_3_X_RAY_H
