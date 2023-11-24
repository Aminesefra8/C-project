#include <iostream> 
#include <cmath>
using namespace std; 
#include<cstring>
#include "Vecteur_bis.h"





int main(){

   //1.2 Validation
   
   //1
   float tab_1[]={1,1,1}; 
   float tab_2[]={3,4,0,0};
   Vecteur u(tab_1,3); //Premier vecteur
   Vecteur v(tab_2,4); //Deuxième vecteur
   u.affiche();
   v.affiche();
   //2
   Vecteur t=u; //On utilise la fonction de l'opérateur d'affectation
   //3
   float tab_3[]={1,1,0};
   Vecteur t_tilde(tab_3,3); //On crée le nouveau vecteur (1,1,0), et on l'affecte à u
   u=t_tilde;
   t.affiche();
   u.affiche();
   //4
   cout<<"Le produit scalaire de v est de "<<dot(v,v)<<endl;
   cout<<"La norme de v est de "<<norm(v)<<endl;
   //5
   (norm(v)*v).affiche(); //On utilise la fonction opérateur *
   //6
    Vecteur w(3);
    w=v.subvec(2,4);
    w.affiche(); //On affiche le vecteur w=(4,0,0) en utilisant la fonction subvec()
    v.affiche();

   //7
   (u+w).affiche();
   (u-w).affiche();

   //2.2 Validation pour les Matrices





   return 0;
}