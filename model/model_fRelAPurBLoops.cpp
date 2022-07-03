#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>
#include <deque>
#include <math.h>

#include <stdlib.h>
#include <stdio.h>

# include <sstream>
 using namespace std;
void ErrorhandlerExit(const char* p1,const char* p2="")
{
  cerr << p1 << "  " << p2 << endl;
  std::exit(1);
}

int main(int argc,char **argv){

//TIME
double t=0.0;
double dt=0.000001, time =100;

//MODEL PARAMETERS
double k0=1.0,k1=1.0, k2=1;
double d0= 1.0, d1= 1.0, d2= 1.0;
double mRNA=0.0, ptf=0.0, target=0.0;
double frequency=0.000001, amplitude = 1;
double frequency_array[] {0.001, 0.01, 0.1, 0.5, 0.0139, 0.0417, 0.01042, 0.02038, 1, 2, 4, 10};

//FILE & PROGRAM PARAMETERS
long counter=0;
double max_target=0.0; 
int reading_rate = 5;

int i, j;
string TF_Prot="", double_as_string, outfile_name;
stringstream double_to_string;
ofstream file_out, file_out_2;

file_out_2.open("parameters.txt");
file_out_2 << "## \tfrequency" << "\t"<< "d0" << "\t"<< "d1" << "\t"<< "d2" << "\t"<< "k0" << "\t"<< "k1" << "\t"<< "k2" << "\t"<< "dt" <<"\t"<< "time" <<endl;

//for loop to run different TF & protein parameter sets
for(j=0;j<2;j++)
{
    if(j==0)
    {
        TF_Prot="RelA_VCAM1_";
        k1=45.0;
        d0=0.08741;
        d1=0.02981;
        d2=0.06945;
        amplitude=1.0;
        dt=0.000037;
        time =1000;
    }
    else
    {
        TF_Prot="PurB_ACTA2_";
        k1=592.0;
        d0=0.2280;
        d1=0.0044;
        d2=0.03775;
        amplitude=1.0;
        dt=0.000037;
        time =2500;
    }

    //for loop to run a set of frequencies
    for(i=18; i<27; i++)
    {

        //setting parameters for new round of model calculation
        mRNA=0;
        ptf=0;
        target=0;
        max_target=0;
        t=0;
        counter=0;
        frequency = 0.000001*pow(2,i);

        //opening and naming of the output file
        double_to_string.str("");
        double_to_string << frequency;
        double_as_string = double_to_string.str();
        outfile_name = "Cpp_model_" + TF_Prot + "f_" + double_as_string + ".txt";
        file_out.open(outfile_name.c_str());

        //writing current model parameters to the output file
        file_out_2 << "##" << "\t" << frequency << "\t"<< d0 << "\t"<< d1 << "\t"<< d2 << "\t"<< k0 << "\t"<< k1 << "\t"<< k2 << "\t" << dt << "\t" << time << endl;
        file_out << /*"##"<<*/"t" << "\t" << "k0" << "\t" << "TF_mRNA" <<  "\t" << "TF_protein" << "\t" << "target_mRNA" << "\t" << "max_target" << "\t" << "f" << "\t" << "par_set" << endl;

        //core while loop for numeric solution of ODEs
        while(t<time)
        {
            //definition of sinus input
            k0=amplitude*sin(2*3.1415*frequency*t);

            //ODEs
            mRNA+= (k0 - mRNA*d0)*dt;
            ptf+= (k1*mRNA - ptf*d1)*dt;
            target+= (k2*ptf - target*d2)*dt;

            //reading of max target
            if((target>max_target)&&(t>(0.75*time))){
                max_target=target;
                }

            //output
            if(counter%reading_rate==0){
                file_out << t << "\t" << k0 << "\t" << mRNA <<  "\t" << ptf << "\t" << target << "\t" << max_target << "\t" << frequency << "\t"<<  TF_Prot << endl;
                }
            counter++;
            t+=dt;

            } //while loop: ODEs
        file_out.close();
        } //for loop: frequencies
    } //for loop: TF and Protein
file_out_2.close();
} //main




