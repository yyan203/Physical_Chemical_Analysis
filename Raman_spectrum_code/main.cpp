// This is an example c++ (c-style code for analyzing MD results
// The basic components are essentially the same as the perl script lmps.analysis.pl. However, note the differences in implementing and speed of these two languages.

//Tasks
//(1) Calculate RAMAN using PRB 2018, 97, 054106, Shcheblanov et al.

//Please Remember to cite Yongjian's work on Raman spectrum if you use this code.

//Memory check in CCNI
// valgrind --leak-check=full ./md.cell dump.melt haha 1 3

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"
#include "mkl_lapacke.h"
#include "mkl_solvers_ee.h"
#include "math.h"
#include "string.h"
#include <vector>
#include <queue>
#include <map>

//#include <complex>
#include <iostream>
#define ONELINEMAX 2000
#define MAXDUMP 10
#define MaxCellAtoms 4000
#define MaxCellNbr 100
#define PI 3.1415926
#define VERBOSE 1

/////////////////////////////////
// Global variables defination //
/////////////////////////////////

class Atom {  // this is about the information of every single atom
public:
    int id;
    int type;
    float x[3];
    float v[3];
    float f[3];
    void PrintInfo() {
        printf("Atom: id[%d] type[%d] x[%f,%f,%f]\n", id,type,x[0],x[1],x[2]);
    };
};



void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
    MKL_INT i, j;
    printf("\n %s\n", desc);
    FILE * fp;
    //fp = fopen(desc, "w");
    //std::vector<double> sum_eigenvector = std::vector<double>(n, 0.0);
    for (i = 0; i < m; i++) {
        for (j = n-1; j >= 0; j--) {
            printf(" %.5f", a[i * lda + j]);
        //if (m >= 1) sum_eigenvector[j]+= a[i* lda + j];
        }
        printf("\n");
    }
    //for (int i = 0; i < sum_eigenvector.size(); i++) fprintf(fp, "%lg\n", sum_eigenvector[i]);
    //fclose(fp);
}
  // get eigenvalue
  void LAPACK_dsyev(int N, double *a, double *w){
      MKL_INT  n = N, lda = N, info;
          /* Executable statements */
          printf( "\n#########################################\nDiagonalizing Dynamical Matrix using LAPACKE_dsyev, please wait for a few seconds ... ...\n\n\n" );
          /* Solve eigenproblem */
          info = LAPACKE_dsyev( LAPACK_ROW_MAJOR, 'V', 'U', n, a, lda, w );
          if( info > 0 ) {
              printf( "The algorithm failed to compute eigenvalues.\n" );
              exit( 1 );
          }
          /* Print eigenvalues */
          //print_matrix( "Eigenvalues", 1, n, w, 1 );
          /* Print eigenvectors */
          //print_matrix( "Eigenvectors (stored columnwise)", n, n, a, lda );
      printf("Print eigenvectors a[][1]\n");
      for (int i=0; i < 10; i++){
          printf("%lg ", a[n*i+1]);
      }
      printf("\n");
  }

void calculate_RAMAN(double *E, double *A, double ***D_alpha,
                     int N, double V_sqrt, Atom **myatom,
                     std::vector<int> &id2index, char * output_prefix, int timestep,
                     double unit2_Rad_per_Second,
                     double TEMPER // Temperature in Kelvin;
                     ){
    printf("Print eigenvectors a[][1]\n");
    for (int i=0; i < 10; i++){
        printf("%lg ", A[N*i+1]);
    }
    printf("\n");
    // N is number of eigenvalues
    //#################### do Raman calculation here:
    std::cout<<"Calculating Raman ...........\n"
            "Number of mode is: "<< N <<"\n"
            "Temperature is: "<< TEMPER << "\n"
            "unit conversion factor to Rad/s is: " << unit2_Rad_per_Second << std::endl;
    //for(int i=1;i <= num_of_atom; i++) printf("%lg ", D_alpha[i]->D[0][1][2]); // i from 1 to num_of_atom

    std::vector<double> I_n = std::vector<double>(N, 0.0);
    std::vector<double> I_P = std::vector<double>(1500, 0.0);
    std::vector<double> I_P_factor = std::vector<double>(1500, 0.0);
    std::vector<int> I_P_HH_count = std::vector<int>(1500, 0);
    std::vector<double> I_P_HH = std::vector<double>(1500, 0.0);
    std::vector<double> I_P_a_n_square = std::vector<double>(1500, 0.0);
    std::vector<double> I_P_b_n_square = std::vector<double>(1500, 0.0);
    std::vector<double> I_P_HV = std::vector<double>(1500, 0.0);
    std::vector<double> w;
    std::vector<double> I;
    double a_n, b_n_square, w_n_THz, w_n_cm, w_L, w_n_sqr, R_ii, mode, M_I; // M_I is mass of Ith atom
    int round_w, myatom_id, atomtype;
    double n_w;
    double h = 4.135668e-15; // eV*s; Planck constant
    double k_B = 8.617e-5; // boltzmann constant  eV/K;
    double Pi= 3.1415926;
    double C = 2.998e8; // light speed, 299800000 m/s
    double V = V_sqrt * V_sqrt;
    double R_11, R_22, R_33, R_12, R_13, R_23; // for b_n_square
    double I_n_HH, I_n_HV;
    std::vector<double> I_n_HH_sum=std::vector<double>(1500,0.0);
    MKL_INT NN=N;
    MKL_INT IND=1;

  double *eigenvector = new double[N];
  for (int eigenvalue = 0; eigenvalue < N; eigenvalue++) {  // eigenvalue is colume number
      //std::cout<<"calculate "<< eigenvalue << " th eigenvalue\n" << std::endl;
      a_n = 0.0; b_n_square = 0.0;
      w_n_sqr = E[eigenvalue];
      if (w_n_sqr <= 0) continue;

      for (int i=0; i < N; i++){  // i is row number
          //if(eigenvalue==1 && i < 10) printf("%lg ", A[N*i+eigenvalue]);
          myatom_id = id2index[int(i/3)+1];
          if(eigenvalue == 11 && i == 11) printf("\n--------------Row: %d 's index is: %d; 's id is: %d--------------\n", i, myatom_id, myatom[myatom_id]->id);
          atomtype = myatom[myatom_id]->type;
          if(atomtype==1) {M_I = 15.9994;} else {M_I = 28.0855;}
          //eigenvector[i] = 1.0;
          eigenvector[i] = A[N*i+eigenvalue] / sqrt(M_I) * V_sqrt;
      }
      //printf("\n");
      w_n_cm = sqrt(w_n_sqr) * unit2_Rad_per_Second * 33.35640 / 2.0 / Pi; // w_n is in unit of cm^-1 (*98.21 * 33.35640)
      round_w = int(w_n_cm);

      R_11 = cblas_ddot(NN,eigenvector, IND, D_alpha[0][0], IND);
      R_22 = cblas_ddot(NN,eigenvector, IND, D_alpha[1][1], IND);
      R_33 = cblas_ddot(NN,eigenvector, IND, D_alpha[2][2], IND);
      R_12 = cblas_ddot(NN,eigenvector, IND, D_alpha[0][1], IND);
      R_13 = cblas_ddot(NN,eigenvector, IND, D_alpha[0][2], IND);
      R_23 = cblas_ddot(NN,eigenvector, IND, D_alpha[1][2], IND);

      if(eigenvalue==1){printf("\nfirst_wn: R11 %lg  R22 %lg R33 %lg R12 %lg\n", R_11, R_22, R_33, R_12);}
      if(round_w <= 0 || round_w >= 1500) continue;
      I_P_HH_count[round_w]+=1;
      w_n_THz = sqrt(w_n_sqr) * unit2_Rad_per_Second / 2.0 / Pi; // sqrt(w_n_sqr) * unit2_Rad_per_Second is in unit of Rad/s
      a_n = 1.0 / 3.0 * (R_11 + R_22 + R_33);
      b_n_square = 0.5 * ((R_11-R_22)*(R_11-R_22)+(R_11-R_33)*(R_11-R_33)+(R_22-R_33)*(R_22-R_33))
             + 3.0 * (R_12*R_12 + R_13*R_13 + R_23*R_23);
      I_n_HH = a_n * a_n + 4.0 / 45.0 * b_n_square;
      I_n_HV = 3.0 / 45.0 * b_n_square;
      //I_n_HH_sum[round_w] += I_n_HH;
      //I_P_a_n_square[round_w] += a_n * a_n;
      //I_P_b_n_square[round_w] += b_n_square;
      //w_L = 654.7 * 33.35640; // cm^-1, 457.9 nm laser
      w_L = 563.5 * 33.35640; // cm^-1, 532 nm laser from Madoka
      //n_w = 1.0/(exp(h/k_B/T*w_n_THz*1.0e12)-1.0);
      n_w = 1.0/(exp(w_n_THz/TEMPER*47.994)-1.0);
      //double h = 4.135668e-15; // eV*s; Planck constant
      //double k_B = 8.617e-5; // boltzmann constant  eV/K;
      //double Pi= 3.1415926;
      //double C = 2.998e8; // light speed, 299800000 m/s
      //I_P_factor[round_w] +=  1.0/w_n_cm*(n_w+1.0)*pow(abs(w_L - w_n_cm), 4.0);
      I_P_HH[round_w] += 1.0/w_n_cm*(n_w+1.0)*pow(abs(w_L - w_n_cm), 4.0) * I_n_HH * V;
      I_P_HV[round_w] += 1.0/w_n_cm*(n_w+1.0)*pow(abs(w_L - w_n_cm), 4.0) * I_n_HV * V;
}


  FILE * raman;
  char filename[100];
  sprintf(filename,"%s.RAMAN.%d", output_prefix, timestep);
  raman = fopen(filename, "w");
  fprintf(raman, "#Frequency(cm^-1) Intensity(I_P_HH) Intensity(I_P_HV)\n");
    //fprintf(raman, "#Frequency(cm^-1) Intensity(I_P_HH) Intensity(I_P_HV)  Intensity(I_n_HH) Intensity(constant_related_to w) a_n^2  b_n^2  I_P_HH_count\n");
  for (int i=0; i < 1500; i++) {
      if (I_P_HH_count[i] > 0)
       fprintf(raman, "%d   %lg  %lg\n", i, I_P_HH[i], I_P_HV[i]);
       //fprintf(raman, "%d   %lg  %lg  %lg  %lg %lg %lg %d\n", i, I_P[i], I_P_HH[i], I_P_HV[i], I_P_factor[i], I_P_a_n_square[i], I_P_b_n_square[i], I_P_HH_count[i]);
  }
  fclose(raman);
}

void calculate_DOS(double *E, int N, char * output_prefix, int timestep, double unit2_Rad_per_Second){
    std::vector<int> DOS = std::vector<int>(1500, 0);
    double w_n_THz;
    int w_n_cm;
    printf("Total Number of Eigenvalue=%d\n",N);
    for (int eigenvalue = 0; eigenvalue < N; eigenvalue++) {
        //printf("Eigenvalue=%lg\n",E[eigenvalue]);
        if (E[eigenvalue]>0) {
            w_n_THz = sqrt(E[eigenvalue]) * unit2_Rad_per_Second / 2.0 / 3.1415926;
            w_n_cm = int(w_n_THz * 33.3564);
            if (w_n_cm <= 1500) DOS[w_n_cm]++;
        }
    }
    FILE * dos;
    char filename[100];
    sprintf(filename,"%s.DOS.%d", output_prefix, timestep);
    dos = fopen(filename, "w");
    fprintf(dos, "#Frequency(cm^-1) frequency Frequency(THz)\n");
    // print DOS
    for (int i=0; i < 1500; i++) {
        fprintf(dos, "%f   %d  %f\n", float(i), DOS[i], float(i)/33.3564);
    }
    fclose(dos);
}

void readdx_upper(double * A, char * filename, int number_atom){

    FILE * fp;
    fp = fopen(filename,"r");
    if(fp==NULL) {
        printf("File %s does not exist!\n", filename);
        exit(0);
    }
    int n2 = number_atom * number_atom * 9;
    int n1 = number_atom * 3;
    for (int i=0;i< n2;i++){A[i]=0.0;}
    double f1, f2, f3;
    double avef1, avef2, avef3;
    char templine[ONELINEMAX];
    // ith atom  and jth atom
    printf("\n\n>>>Please check if the dynamical matrix is symmetrical!\n"
                   ">>>The following are [i(0,1,2)][i(0,1,2)] example terms "
                   "along the diagonal:\n\n");
    for (int i=0; i < number_atom; i++) {
        for (int k=0;k<3;k++){
            for (int j = 0; j < number_atom; j++) {
                if(fgets(templine, ONELINEMAX, fp) == NULL) {
                    printf("Dynamical matrix reading error!\n");
                    printf("file: %s\n",filename);
                    fclose(fp); exit(1);
                }
                sscanf(templine, "%lg %lg %lg", &f1, &f2, &f3);
                if(i<=2 && i==j) {std::cout << f1 << " " << f2 << " " << f3 <<  std::endl; if((k+1)%3==0) std::cout<< std::endl;}
                if (i * 3 + k >= number_atom*3 || j*3 +2 >= number_atom*3){
                    printf("Dynamical matrix index error: i*3+k: %d, j*3+2:%d\n",i * 3 + k, j*3 +2);
                    fclose(fp); exit(1);
                }
                if(i*3+k == j*3+0){A[(i*3+k) * n1 + j*3+0] = f1;}
                if(i*3+k == j*3+1){A[(i*3+k) * n1 + j*3+1] = f2;}
                if(i*3+k == j*3+2){A[(i*3+k) * n1 + j*3+2] = f3;}
                if(i*3+k <  j*3+0){A[(i*3+k) * n1 + j*3+0] = f1;}
                if(i*3+k <  j*3+1){A[(i*3+k) * n1 + j*3+1] = f2;}
                if(i*3+k <  j*3+2){A[(i*3+k) * n1 + j*3+2] = f3;}
                if(i*3+k >  j*3+0){avef1=0.5*(f1+A[(j*3+0)*n1 + i*3+k]); A[(j*3+0)*n1+i*3+k]=avef1;}
                if(i*3+k >  j*3+1){avef2=0.5*(f2+A[(j*3+1)*n1 + i*3+k]); A[(j*3+1)*n1+i*3+k]=avef2;}
                if(i*3+k >  j*3+2){avef3=0.5*(f3+A[(j*3+2)*n1 + i*3+k]); A[(j*3+2)*n1+i*3+k]=avef3;}
            }
        }
    }
    fclose(fp);
    //for(int i=0;i<n2;i++){printf("i, A[i] = %d %f  ", i, A[i]);}
    printf("Successfully imported the Dynamical matrix!\n");
} // end of readdx


// calculate distance between atoms
float bondlen(const Atom * A, const Atom * B, const float &lx, const float &ly, const float &lz)
{
    float dist=0;
    float x1x,x1y,x1z;
    x1x=A->x[0]-B->x[0];
    x1y=A->x[1]-B->x[1];
    x1z=A->x[2]-B->x[2];
    while (x1x>lx/2.0) {x1x=x1x-lx;}
    while (x1x<-lx/2.0){x1x=x1x+lx;}
    while (x1y>ly/2.0) {x1y=x1y-ly;}
    while (x1y<-ly/2.0){x1y=x1y+ly;}
    while (x1z>lz/2.0) {x1z=x1z-lz;}
    while (x1z<-lz/2.0){x1z=x1z+lz;}
    dist=sqrt(x1x*x1x+x1y*x1y+x1z*x1z);
    return dist;
}

// calculate angle between A-B-C
//float bondangle( Atom& A,  Atom& B,  Atom& C)
float bondangle(const Atom * A, const Atom * B, const Atom * C, const float &lx, const float &ly, const float &lz)
{
    float angle=0;
    float v1x,v1y,v1z;
    float v2x,v2y,v2z;
    v1x=A->x[0]-B->x[0];
    v1y=A->x[1]-B->x[1];
    v1z=A->x[2]-B->x[2];
    v2x=C->x[0]-B->x[0];
    v2y=C->x[1]-B->x[1];
    v2z=C->x[2]-B->x[2];
    while (v1x>lx/2.0) {v1x=v1x-lx;}
    while (v1x<-lx/2.0){v1x=v1x+lx;}
    while (v1y>ly/2.0) {v1y=v1y-ly;}
    while (v1y<-ly/2.0){v1y=v1y+ly;}
    while (v1z>lz/2.0) {v1z=v1z-lz;}
    while (v1z<-lz/2.0){v1z=v1z+lz;}
    while (v2x>lx/2.0) {v2x=v2x-lx;}
    while (v2x<-lx/2.0){v2x=v2x+lx;}
    while (v2y>ly/2.0) {v2y=v2y-ly;}
    while (v2y<-ly/2.0){v2y=v2y+ly;}
    while (v2z>lz/2.0) {v2z=v2z-lz;}
    while (v2z<-lz/2.0){v2z=v2z+lz;}
    angle=acos((v1x*v2x+v1y*v2y+v1z*v2z)/(sqrt(v1x*v1x+v1y*v1y+v1z*v1z)*sqrt(v2x*v2x+v2y*v2y+v2z*v2z)))*180.0/PI;
    if(angle<0) angle+=180;
    if(angle>180) angle-=180;
    return angle;
}

class Cell {   // this is about the cell formed by dividing the box
public:
    int member[MaxCellAtoms];
    int nmember;
    int nbr[MaxCellNbr];
    int nnbr;
    float origin[3];
    float size[3];
    Cell() {nmember=0;nnbr=0;};
};

class SimulationBox {
public:
    Atom **myatom;
    Cell ****mycell;
    int ncell[3]; //dimension of the cells
    float L[3];
    float origin[3]; //origin of the simulation box, default will be (0,0,0)
    int nmols;
    void SetupCell(float cellsize) {
#if VERBOSE ==1
        printf("Initiating setupcell\n");
#endif
        int i,j,k,p,q,r,pi,qi,ri,l;
        if(nmols==0) {
            printf("Simulation Box is empty!\n");
            exit(0);
        }
        for(int d=0;d<3;d++) {
            ncell[d] = (int)(L[d]/cellsize)-1; //make it a little larger than necessary
            //if(ncell[d]<4) {ncell[d]=1;} //if the cells are too few, there are no advantages of doing link-cell.
        }
#if VERBOSE ==1
        printf("cell[%d,%d,%d]\n", ncell[0],ncell[1],ncell[2]);
#endif
        //if(ncell[0]==1||ncell[1]==1||ncell[2]==1){printf ("Warning!!! Cells are too few to use linkcell!!! Please decrese cellsize.\n"); exit(0);}
        mycell = new Cell ***[ncell[0]];
        for(i=0;i<ncell[0];i++) {
            mycell[i] = new Cell **[ncell[1]];
            for(j=0;j<ncell[1];j++) {
                mycell[i][j] = new Cell *[ncell[2]];
                for(k=0;k<ncell[2];k++) {
                    mycell[i][j][k] = new Cell();
                }
            }
        }

        //Setup NBR cells, Coding here
///////////////////////////////////////////
///////////////////////////////////////////
///////Yongjian code following/////////////
///////////////////////////////////////////
///////////////////////////////////////////
        //loop of cells
        for(i=0;i<ncell[0];i++) {
            for(j=0;j<ncell[1];j++) {
                for(k=0;k<ncell[2];k++) {
//printf("k has a value: %d\n",k);

                    //loop of neighbour cells
                    mycell[i][j][k]->nnbr=27;l=0;
                    for (pi=i-1;pi<i+2;pi++){
                        for (qi=j-1;qi<j+2;qi++){
                            for (ri=k-1;ri<k+2;ri++){
//printf("ri has a value: %d\n",ri);
                                //convert neighbour cells outside the range (periodic boundary condition)
                                p=pi;q=qi;r=ri;l=(pi-i+1)*9+(qi-j+1)*3+(ri-k+1);
                                if (l<0 || l>26) {printf( "l:%d is wrong",l); exit(0);}
                                if (p==-1){p=ncell[0]-1;}
                                if (p==ncell[0]){p=0;}  //if (p==-1||p==ncell[0]){p=abs(abs(p)-ncell[0]);}
                                if (q==-1){q=ncell[1]-1;}
                                if (q==ncell[1]){q=0;}  //if (q==-1||q==ncell[1]){q=abs(abs(q)-ncell[1]);}
                                if (r==-1){r=ncell[2]-1;}
                                if (r==ncell[2]){r=0;}  //if (r==-1||r==ncell[2]){r=abs(abs(r)-ncell[2]);}
                                //convert every neighbour's vector index into a scalor index
                                mycell[i][j][k]->nbr[l]=p*ncell[1]*ncell[2]+q*ncell[2]+r;
//printf("l has a value: %d\nnbr[l] has a value:%d\ni=%d j=%d k=%d\n",l,mycell[i][j][k]->nbr[l],i,j,k);

                            } } }
//printf("mycell[i][j][k] has %d neighbour\ni j k are:%d %d %d\n", mycell[i][j][k]->nnbr,i,j,k);
                }
            }
        }
        printf("i j k are  %d %d %d\n", i,j,k);
    };
//////////////////////////////////////////
//////////////////////////////////////////
////////Yongjian code ends////////////////
//////////////////////////////////////////
//////////////////////////////////////////

    SimulationBox() {ncell[0]=ncell[1]=ncell[2]=0; nmols=0;mycell=NULL;};

    ~SimulationBox() {
        int i,j,k;
        if(mycell!=NULL) {
#if VERBOSE == 1
            printf("Now destroy mycell.\n");
#endif
            for(i=0;i<ncell[0];i++) {
                for(j=0;j<ncell[1];j++) {
                    for(k=0;k<ncell[2];k++) {
                        delete mycell[i][j][k];
                    }
                    delete [] mycell[i][j];
                }
                delete [] mycell[i];
            }
            delete [] mycell;
        }
    };
    void PrintInfo() {
        printf("SimulationBox: ncell(%d,%d,%d), L(%f,%f,%f), origin(%f,%f,%f),nmols(%d)\n",
               ncell[0],ncell[1],ncell[2],
               L[0],L[1],L[2],
               origin[0],origin[1],origin[2],
               nmols);
    };
};

/////////////////
// Subroutines //
/////////////////

////////////////////
// Main functions //
////////////////////

int main(int argc, char **argv)
{

    int i=0,j=0,k=0;
    int il=0,jl=0,kl=0;
    char inputfilename[100], DMfilename[100], outputfilename[100];
    int example_atom;
    int charge;
    double unit2_Rad_per_Second;
    double TEMPER;
    FILE *ifp, *ofp;
    char templine[ONELINEMAX], str1[100], str2[100];
    long current_timestep, nmols,nmols_speciesA,nmols_speciesB;
    long initial_timestep = -1;
    float xlo, xhi, ylo, yhi, zlo, zhi;
    int id, type;
    float iq, ix, iy, iz, ivx, ivy, ivz, ifx, ify, ifz;
    SimulationBox **mysystem;
    Atom **myatom; //convenient pointer
    int iframe=0; //input number of frames
    int counter = 1; // for running average of frames
    int num_frames = 0; // for running average of frames
    int totalframes;
    int speciesA;
    int speciesB;
    float cellsize;
    float bondlength;

    //************* Command line parameter processing ************/

    if(argc != 12) {
        // Provide the correct syntax of this analysis code
        printf("Correct syntax: md.analysis.cxx 1.inputfile1[lata4olivia] 2.inputfile2[dynamical_matrix_file] "
                       "3.outputfile[prefix_to_Raman/DOS_dat] \n4.example_atom(print atomid) 5.cellsize "
                       "6.speciesA_type 7.speciesB_type 8.if_lata4olivia_contain_charge[1/0] 9.[bond-cutoff] "
                       "10.[unit conversion to Rad/s] 11.Temperature(K)\n");
        exit(0);
    }
    else {
        //Note that argv[0] is "md.analysis.cxx"
        sscanf(argv[1], "%s", inputfilename);
        sscanf(argv[2], "%s", DMfilename);
        sscanf(argv[3], "%s", outputfilename);
        sscanf(argv[4], "%d", &example_atom);
        sscanf(argv[5], "%f", &cellsize);
        sscanf(argv[6], "%d", &speciesA);
        sscanf(argv[7], "%d", &speciesB);
        sscanf(argv[8], "%d", &charge);
        sscanf(argv[9], "%f", &bondlength);
        sscanf(argv[10], "%lg", &unit2_Rad_per_Second);
        sscanf(argv[11], "%lg", &TEMPER);
    }

    /************** Import the data file ********************/
    ifp = fopen(inputfilename, "r");
    //ofp = fopen(outputfilename, "w");

    if(ifp==NULL) {
        printf("File %s does not exist!\n", inputfilename);
        exit(0);
    }




    //initialize Simulation Boxes
    mysystem = new SimulationBox *[MAXDUMP];
    for(i=0;i<MAXDUMP;i++) {mysystem[i] = new SimulationBox();}
    printf("Initialization of Simulation Boxes has been finished! %d\n",i);

    //initialize bins
    int nbins=181;
    long anglenumber=0;
    int bin[nbins];
    int bbin[10][nbins];
    float g[nbins];
    float gg[10][nbins];
    //initialize the array!
    for(i=0;i<nbins;i++) {bin[i]=0;g[i]=0.0;
        for(int j=0;j<10;j++){bbin[j][i]=0;gg[j][i]=0.0;}
    }

    //Read in the entire file

    while(!feof(ifp)) {
        fgets(templine, ONELINEMAX, ifp); // TIMESTEP line
        sscanf(templine, "%s %s", str1, str2);

        if(strcmp(str1, "ITEM:")==0) {
            fgets(templine, ONELINEMAX, ifp); //actual timesteps
            sscanf(templine, "%ld", &current_timestep);
            if(initial_timestep==-1){initial_timestep=current_timestep;}
            fgets(templine, ONELINEMAX, ifp); // ITEM: NUMBER OF ATOMS
            fgets(templine, ONELINEMAX, ifp); // nmols
            sscanf(templine, "%ld", &nmols);  //number of atoms

            fgets(templine, ONELINEMAX, ifp); // ITEM: BOX BOUND
            fgets(templine, ONELINEMAX, ifp); // xlo xhi
            sscanf(templine, "%f %f", &xlo, &xhi);
            fgets(templine, ONELINEMAX, ifp); // ylo yhi
            sscanf(templine, "%f %f", &ylo, &yhi);
            fgets(templine, ONELINEMAX, ifp); // zlo zhi
            sscanf(templine, "%f %f", &zlo, &zhi);
            fgets(templine, ONELINEMAX, ifp); //empty line

            mysystem[iframe]->origin[0] = xlo;
            mysystem[iframe]->origin[1] = ylo;
            mysystem[iframe]->origin[2] = zlo;

            mysystem[iframe]->L[0] = xhi - xlo;
            mysystem[iframe]->L[1] = yhi - ylo;
            mysystem[iframe]->L[2] = zhi - zlo;

            double V = (xhi-xlo) * (yhi-ylo) * (zhi-zlo);
            double V_sqrt = sqrt(V);

            printf("current_timestep %d\n",current_timestep);


            std::vector<int> id2index = std::vector<int>(nmols + 1, 0); // from atom_ID to index in  myatom[]

            //initialize myatom
            if(nmols>0) {
                mysystem[iframe]->nmols = nmols;
                mysystem[iframe]->myatom = new Atom *[nmols];
                myatom = mysystem[iframe]->myatom;
                for(i=0;i<nmols;i++) {
                    myatom[i] = new Atom();
                }

                printf("Initialization of myatom has been successfully finished! %d \n",i);
            }
            else {
                printf("Error in nmols (%d) \n", nmols);
                exit(0);
            }

            //Be aware that id may start from 1
            //id is the TRUE identification of atoms
            //your i index is NOT the identification of atoms

#if VERBOSE == 1
            printf("Importing %d atoms...\n", nmols);
#endif
            nmols_speciesA=0;
            nmols_speciesB=0;
            for(i=0;i<nmols;i++) {
                fgets(templine, ONELINEMAX, ifp);

                if(charge){
                    sscanf(templine, "%d %d %f %f %f %f %f %f %f %f %f %f",
                           &id, &type, &iq, &ix, &iy, &iz, &ivx, &ivy, &ivz, &ifx, &ify, &ifz);
                }
                else {
                    sscanf(templine, "%d %d %f %f %f %f %f %f %f %f %f",
                           &id, &type, &ix, &iy, &iz, &ivx, &ivy, &ivz, &ifx, &ify, &ifz);
                }
                myatom[i]->id = id; if (id==4000){printf("%d\n",id);}
                myatom[i]->type = type;
                id2index[id] = i; // very important line to bookkeeping the atom_ID to index;

                if(type==speciesA) {nmols_speciesA++;}
                if(type==speciesB || speciesB==-1) {nmols_speciesB++;}
                myatom[i]->x[0] = ix;
                myatom[i]->x[1] = iy;
                myatom[i]->x[2] = iz;
            }
            //Setup the LinkCell
            mysystem[iframe]->SetupCell(cellsize);
            printf("%d\n",iframe);
            //printf("mycell-111 nmember:%d nnbr:%d\n",mysystem[iframe]->mycell[1][1][1]->nmember,mysystem[iframe]->mycell[1][1][1]->nnbr);

            //Any analysis routine

            //readdx(A, myfile, 27, 2000);  // just read one frame
//////////////////////////////////////////////////////
//////Yongjian Coded the following part///////////////
//////////////////////////////////////////////////////

            //Throw in atoms, Coding here
            float x1=0,x2=0,x3=0,originx,originy,originz,cellsizex,cellsizey,cellsizez;
            int p1,p2,p3,memberid;
            originx= mysystem[iframe]->origin[0]; printf("originx=%f\n",originx);
            originy= mysystem[iframe]->origin[1];
            originz= mysystem[iframe]->origin[2];
            cellsizex=(float)(mysystem[iframe]->L[0]/mysystem[iframe]->ncell[0]); printf("cellsizex=%f\n",cellsizex);
            cellsizey=(float)(mysystem[iframe]->L[1]/mysystem[iframe]->ncell[1]);
            cellsizez=(float)(mysystem[iframe]->L[2]/mysystem[iframe]->ncell[2]);
            //loop of every atoms of one frame in dump file
            for (i=0;i<nmols;i++) {
                x1=myatom[i]->x[0]; // if (i==400){printf("%d\n",i);}
                x2=myatom[i]->x[1];
                x3=myatom[i]->x[2];
                //according to the atom's position, allocate them to the right cell. Since cellsizex is float number, to make it safe, when the index from atom's position is calculated from cellsizex, bring it back to within the cell.
                p1=(int)((x1-originx)/cellsizex);
                if (p1==mysystem[iframe]->ncell[0]){p1=p1-1;}
                p2=(int)((x2-originy)/cellsizey);
                if (p2==mysystem[iframe]->ncell[1]){p2=p2-1;}
                p3=(int)((x3-originz)/cellsizez);
                if (p3==mysystem[iframe]->ncell[2]){p3=p3-1;}
                memberid=mysystem[iframe]->mycell[p1][p2][p3]->nmember;
                //allocate the atom's index to the right cell
                mysystem[iframe]->mycell[p1][p2][p3]->member[memberid]=i;
                //mysystem[iframe]->mycell[p1][p2][p3]->member[memberid]=myatom[i]->id;
                //increase the nmember of cell when a new atom is added into it.
                mysystem[iframe]->mycell[p1][p2][p3]->nmember=mysystem[iframe]->mycell[p1][p2][p3]->nmember+1;
            }

//for (i=0;i<7;i++) {printf("mycell[0][%d][0] has nmember:%d\n",i,mysystem[iframe]->mycell[0][i][0]->nmember);}

/////////////////////////////////////////////
            //Do Derivative_Polarizability calculations, Coding here
/////////////////////////////////////////////

            //initialize alpha derivative in  x  y  z directions
            // Derivative of local bond polarizability:
            int num_of_atom;
            num_of_atom = nmols;
            double ***D_alpha;
            D_alpha = new double **[3];
            for(int i=0;i<3;i++) {
                D_alpha[i] = new double *[3];
                for (int j = 0; j < 3; j++) {
                    D_alpha[i][j] = new double[num_of_atom * 3];
                    for (int k = 0; k < num_of_atom * 3; k++) {
                        D_alpha[i][j][k] = 0.0;
                    }
                }
            }


            float lx=mysystem[iframe]->L[0], ly=mysystem[iframe]->L[1],lz=mysystem[iframe]->L[2];
            float dr=1,dx,dy,dz,dist,area,angle;
            float Rc[3];
            int index;
            double delta_ij,delta_ik,delta_jk;

            //for (iframe=0, iframe<MAXDUMP,iframe++) //loop of different atoms in different cell;
            int cellx,celly,cellz,nmember,imember; //cell index
            int ncellx,ncelly,ncellz; //neighbour cell index
            int ncell[3],iii,ii,c,inmember;
            int comparei,comparej;
            int itype,jtype;

            for(iii=0;iii<3;iii++) {ncell[iii]=mysystem[iframe]->ncell[iii];} //initialize array ncell[3]
            //loop among different cell
            for (cellx=0;cellx<ncell[0];cellx++){
                for (celly=0;celly<ncell[1];celly++){
                    for (cellz=0;cellz<ncell[2];cellz++){
                        nmember=mysystem[iframe]->mycell[cellx][celly][cellz]->nmember;

//loop among all the atoms contained in one cell, this CENTER atom
//##################################################################
                        //std::vector<int> neighbours;
                        //std::priority_queue<std::pair<float,int>, std::vector<std::pair<float,int> >,CompareDist> q;
                        for (imember=0;imember<nmember;imember++){
                            //loop among all neighbour cells
                            for (ii=0;ii<27;ii++){
                                c=mysystem[iframe]->mycell[cellx][celly][cellz]->nbr[ii];

                                ncellx=(int)(c/(ncell[1]*ncell[2]));
                                ncelly=(int)(fmod((int)(c/ncell[2]),ncell[1]));
                                ncellz=(int)(fmod(c,ncell[2]));
                                inmember=mysystem[iframe]->mycell[ncellx][ncelly][ncellz]->nmember;
                                //loop among all the neighbour's atoms
                                for (j=0;j<inmember;j++){
                                    comparei=mysystem[iframe]->mycell[cellx][celly][cellz]->member[imember]; // get atom index
                                    comparej=mysystem[iframe]->mycell[ncellx][ncelly][ncellz]->member[j]; // get atom index
                                    itype=myatom[comparei]->type; // get type
                                    jtype=myatom[comparej]->type; // get type
                                    //Avoid double counting of pair comparei-comparej
                                    if (
                                            comparei!=comparej
                                            && ((itype==speciesA && jtype==speciesB) || (itype=speciesB && jtype==speciesA))
                                            )
                                      {
                                        float xi,yi,zi,xj,yj,zj;
                                        xi=myatom[comparei]->x[0];xj=myatom[comparej]->x[0];
                                        yi=myatom[comparei]->x[1];yj=myatom[comparej]->x[1];
                                        zi=myatom[comparei]->x[2];zj=myatom[comparej]->x[2];
                                        dx=xi-xj; //printf("dx=%f\n",dx);
                                        dy=yi-yj;
                                        dz=zi-zj;
                                        while (dx>lx/2.0) {dx=dx-lx;}
                                        while (dx<-lx/2.0){dx=dx+lx;}
                                        while (dy>ly/2.0) {dy=dy-ly;}
                                        while (dy<-ly/2.0){dy=dy+ly;}
                                        while (dz>lz/2.0) {dz=dz-lz;}
                                        while (dz<-lz/2.0){dz=dz+lz;}
                                        dist=sqrt(dx*dx+dy*dy+dz*dz);
                                        Rc[0]=dx/dist;Rc[1]=dy/dist;Rc[2]=dz/dist;
                                        if(dist<=bondlength){
                                            //if(dist<=2.0){
                                            //if(dist==0) std::cout<<"dx dy dz" << dx <<":"<<dy<<";"<<dz<< std::endl;
                                            //std::cout << "compute atom I:J " << comparei<<":" << comparej << std::endl;
                                            for (int i = 0; i < 3; i++)
                                                for (int j = 0; j < 3; j++)
                                                    for (int k = 0; k < 3; k++){
                                                        if (i==j) delta_ij = 1.0;
                                                        else delta_ij = 0.0;
                                                        if (i==k) delta_ik = 1.0;
                                                        else delta_ik = 0.0;
                                                        if (j==k) delta_jk = 1.0;
                                                        else delta_jk = 0.0;
                                                        //std::cout << "i:j:k;dist " << i <<":" << j <<":"<<k<<":"<<dist << std::endl;
                                                        if((myatom[comparei]->id - 1)*3 + k >= num_of_atom * 3) {
                                                            printf("D_alpha index go beyond range!\n"); exit(1);
                                                        }
                                                        //printf("Here 100! i j k index %d %d %d %d\n", i, j, k, (myatom[comparei]->id -1)*3 + k);
                                                        D_alpha[i][j][(myatom[comparei]->id - 1) * 3 + k] += 1.0/ V * (
                                                                             1.0/3.0 * 0.771 * delta_ij * Rc[k]
                                                                           + 0.196 * (Rc[i] * Rc[j] - 1.0/3.0 * delta_ij) * Rc[k]
                                                                           + 0.056 * (delta_ik * Rc[j] + delta_jk * Rc[i] - 2.0 * Rc[i]*Rc[j]*Rc[k])
                                                                                                                      );

                                                    }
                                        }
                                    }
                                }
                            }

                        }
                    }
                }
            }

           // // print D_alpha
           // double sum = 0.0;
           // printf("D_alpha[0][1]:\n");
           // for (int i = 0; i < num_of_atom * 3;i++) {
           //     printf("%lg  ", D_alpha[0][1][i]);
           //     sum += D_alpha[0][1][i];
           // }
           // printf("\nsum = %lg\n", sum);
         //#################### do DOS calculation here:
            int MM=0;
            double *A;
            double *E, *X; //E store eigenvalue; X store eigenvector;
            int DimN = num_of_atom * 3;
            E = new double[DimN];
            //X = new double[DimN * DimN];
            A = new double[DimN * DimN];
            // using lapack LAPACK_dsyevr
            readdx_upper(A, DMfilename, num_of_atom);  // read Dynamical Matrix from output by Dynamical_Matrix
            LAPACK_dsyev(DimN, A, E);

            std::cout << "Diagonalization of the matrix is finished!!" << std::endl;
            //std::cout << "first eigenvalue and eigenvector are:" << std::endl << es.eigenvalues()[0].real() << std::endl << es.eigenvectors().col(0) << std::endl;
           calculate_DOS(E, DimN, outputfilename, current_timestep, unit2_Rad_per_Second);
           calculate_RAMAN(E, A, D_alpha, DimN, V_sqrt, myatom, id2index,
                           outputfilename, current_timestep, unit2_Rad_per_Second, TEMPER);
            //deallocate memory
            for(int i=0;i<3;i++) {
                for(int j=0;j<3;j++) {
                    delete [] D_alpha[i][j];
                }
                delete [] D_alpha[i];
            }
            delete [] D_alpha;

            delete [] A;
            delete [] E;

            printf( "This is frame:(%d) \n", iframe);
            printf("speciesA:%d, speciesB:%d\n",speciesA,speciesB);
            printf("nmols_speciesA:%d, nmols_speciesB:%d\n",nmols_speciesA,nmols_speciesB);

            //The simulation box information is in mysystem[iframe]. For instance mysystem[iframe]->L[0] will be X-dimension of the system
            //The atom information is in myatom. For instance myatom[i]->x[0] is the x-coordinate of myatom
////////////////////////////////////////////////////
/////////// End of Yongjian's coding////////////////
////////////////////////////////////////////////////
            //Increment
            iframe++;
            if(iframe>=MAXDUMP) {
                printf("Too many dumps (%d), increase MAXDUMP!\n", iframe);
                exit(0);
            }

        }
        else {
            if(!feof(ifp))
                printf("Not the right syntax, end of file.\n");
            else
                printf("Finish reading in input file\n");
        }
    }
    totalframes = iframe;
    printf("Total number of frames is %d\n", totalframes);

    //output

    for(j=0;j<totalframes;j++) {
        for(i=0;i<mysystem[j]->nmols;i++) {
            if(mysystem[j]->myatom[i]->id==example_atom)
                fprintf(ofp, "%d %f %f %f\n", j,
                        mysystem[j]->myatom[i]->x[0],
                        mysystem[j]->myatom[i]->x[1],
                        mysystem[j]->myatom[i]->x[2]
                );
        }
    }

    fclose(ifp);
    //fclose(ofp);

    //deallocate memory
    for(int j=0;j<totalframes;j++) {
        myatom = mysystem[j]->myatom;
        for(int i=0;i<mysystem[j]->nmols;i++) {
            delete myatom[i];
        }
        delete [] myatom;
    }
    for(j=0;j<MAXDUMP;j++) {
        delete mysystem[j];
    }
    delete [] mysystem;

}  // main
