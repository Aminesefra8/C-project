#include <iostream> 
#include <cmath>
using namespace std; 
#include<cstring>
#include "BIS.h"



int main(){

    //Vérification

    //1 Givens

    float c,s;
    givens(1,2,c,s);
    cout<<"c:"<<c<<endl;
    cout<<"s:"<<s<<endl;

    //2 Householder

    float beta_1;
    float tab[2]={-1,0};
    Vecteur V_1(tab,2);
    Vecteur V_2;
    V_2=householder(V_1,beta_1);
    V_2.affiche();
    cout<<"beta= "<<beta_1<<endl;

    //

    float beta_2;
    float tab_bis[2]={1/sqrt(2),1/sqrt(2)};
    Vecteur V_3(tab_bis,2);
    Vecteur V_4=householder(V_3,beta_2);
    V_4.affiche();
    cout<<"beta= "<<beta_2<<endl;

    //
    
    float beta_3;
    float tab_tilde[1]={-4};
    Vecteur V_5(tab_tilde,1);
    Vecteur V_6=householder(V_5,beta_3);
    V_6.affiche();
    cout<<"beta= "<<beta_3<<endl;
    
    //qrsym

    float tabl_1[2]={-20,-6};
    float tabl_2[2]={-6,10};
    Vecteur Vect_1(tabl_1,2);
    Vecteur Vect_2(tabl_2,2);
    Vecteur tabVect[2]={Vect_1, Vect_2};
    int dim[2]={2,2};
    Matrice M(tabVect,dim);
    Matrice Q(dim[0],dim[1]);
    qrsym(M,Q);
    Q.affiche();

    //qrpivot

    float T1[3]={1,1/2,1/3};
    float T2[3]={1/2,1/3,1/4};
    float T3[3]={1/3,1/4,1/5};
    Vecteur Vec1(T1,3);
    Vecteur Vec2(T2,3);
    Vecteur Vec3(T3,3);
    Vecteur tabvec[3]={Vec1, Vec2, Vec3};
    int dime[3]={3,3};
    Matrice M_tilde(tabvec,dime);
    Matrice Q_tilde(3,3);
    Matrice PI=qrpivot(M_tilde, Q_tilde);
    cout<<"R: "<<endl;
    M_tilde.affiche();
    cout<<"Q: "<<endl;
    Q_tilde.affiche();
    cout<<"PI: "<<endl;
    PI.affiche();
    


    return 0;
}