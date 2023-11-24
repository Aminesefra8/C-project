#include <iostream> 
#include <cmath>
using namespace std; 
#include<cstring>


class Vecteur{

    float *tab;
    int dim;
    
    public:

       Vecteur():tab(nullptr), dim(0){}

       Vecteur(int i){

           dim=i;
           tab= new float[dim];
           for(int e=0;e<i;e++){
               tab[e]=0;
           }
        }

        Vecteur(float *t, int n){

            dim=n;
            tab= new float[dim];
            for(int i=0;i<n;i++){
                tab[i]=t[i];
            }
        }

        Vecteur(const Vecteur& v){

            dim=v.dim;
            tab= new float[dim];
            for(int i=0;i<dim;i++){
                tab[i]=v.tab[i];
            }
        }

        ~Vecteur() { delete[] tab; }

         void affiche() const {
            std::cout << "Vecteur de taille " << dim << " : ";
            for(int i=0; i<dim; i++){
                 std::cout << tab[i] << " ";
            }
            std::cout << std::endl;
        }

        Vecteur& operator=(const Vecteur& v) {
            if(this != &v) {
               delete[] tab;
               dim = v.dim;
               tab = new float[dim];
               for(int e=0;e<dim;e++){
                   tab[e]=v.tab[e];
               }
            }
            return *this;
        }

        friend Vecteur operator+(const Vecteur& v1, const Vecteur& v2) {
             Vecteur res(v1.dim);
            for(int i=0; i<v1.dim; i++){
                 res[i] = v1[i] + v2[i];
            }
             return res;
        }

        friend Vecteur operator-(const Vecteur& v1, const Vecteur& v2) {
             Vecteur res(v1.dim);
            for(int i=0; i<v1.dim; i++){
                 res[i] = v1[i] - v2[i];
            }
            return res;
        }

        float& operator[](int i){
            return tab[i];
        }

        const float& operator[](int i) const{
            return tab[i];
        }

        Vecteur subvec(int i, int j) const {
            Vecteur res(j-i+1);
            for(int k=0; k<j-i+1; k++){
               res[k] = tab[i+k-1];
            }
            return res;
        }
        
        friend float dot(const Vecteur& v1, const Vecteur& v2) {
           float res = 0.0;
            for(int i=0; i<v1.dim; i++){
               res += v1[i] * v2[i];
            }
            return res;
        }

        friend float norm(const Vecteur& v) {
           float res = 0.0;
           for(int i=0; i<v.dim; i++){
               res += v[i] * v[i];
           }
           return sqrt(res);
        }


        friend Vecteur operator*(float f, const Vecteur& v) {
            Vecteur res(v.dim);
            for(int i=0; i<v.dim; i++){
               res[i] = f * v[i];
            }
            return res;
        }
        
        int Getdim(){
            return dim;
        }

        float Gettab(int i){
            return tab[i];
        }



    
};