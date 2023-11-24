#include <iostream> 
#include <cmath>
using namespace std; 
#include<cstring>
#include "BIS.h"



int main(){

    //Vérification pour la partie 2

    //1

    //Création de la matrice A

    float T_1[3]={1,1,0};
    Vecteur vect_1(T_1,3);
    float T_2[3]={-0.5,2,-1};
    Vecteur vect_2(T_2,3);
    float T_3[3]={0,-1,1};
    Vecteur vect_3(T_3,3);
    Vecteur V[3]={vect_1,vect_2,vect_3};
    int dim[2]={3,3};
    Matrice A(V,dim);

    //Création de la matrice B

    float T_1_bis[2]={-2,0};
    Vecteur vect_1_bis(T_1_bis,2);
    float T_2_bis[2]={3,1};
    Vecteur vect_2_bis(T_2_bis,2);
    Vecteur V_bis[2]={vect_1_bis, vect_2_bis};
    int dim_bis[2]={2,2};
    Matrice B(V_bis,dim_bis);

    //Affichage

    A.affiche();
    B.affiche();

    //2

    //Copie de B dans C

    Matrice C=B;
    B[1][0]=0;
    B.affiche();
    C.affiche();

    //3

    Matrice D=A.submat(1,3,1,2);
    D.affiche();

    //4

    float tab[3]={3,2,1};
    Vecteur vector(tab,3);
    Matrice E(vector);
    E.affiche();
   
   //5

   (B+C).affiche();
   (C-B).affiche();
   (D*C).affiche();

   //6

   cout<<"Norme de Frobenius de C: "<<C.norm()<<endl;

   //7

   (0.5*(B+B.transpose())).affiche();




    return 0;
}