#include <iostream> 
#include <cmath>
using namespace std; 
#include<cstring>
#include<stdexcept>
#include "Vecteur_bis.h"



class Matrice{ //Création de la classe Matrice
    Vecteur *mat;
    int dims[2];
    public:

        void affiche();
        
        Matrice(int i, int j){

           this->dims[0]=i;
           this->dims[1]=j;
           this->mat= new Vecteur[j];
           for(int e=0;e<j;e++){
              this->mat[e]= Vecteur(i);
            }

        }

        Matrice(Vecteur v){

           this->dims[0]= v.Getdim();
           this->dims[1]= v.Getdim();
           this->mat= new Vecteur[v.Getdim()];
           for(int i=0;i<v.Getdim();i++){
              this->mat[i]=Vecteur(v.Getdim());
              (this->mat[i])[i]=v.Gettab(i);
           }
        }

        Matrice(Vecteur *V, int dims[2]){
     
           this->dims[0]=dims[0];
           this->dims[1]= dims[1];
           this->mat= new Vecteur[dims[1]];
           for(int e=0;e<dims[1];e++){
              if(V[e].Getdim() != dims[0]){
                 throw std::invalid_argument("Les vecteurs n'ont pas la même dimension.");
              }
              this->mat[e]=V[e];
           }

        }

        ~Matrice(){

           delete[] this->mat;
        }

        Matrice(const Matrice &M){

          this->dims[0]=M.dims[0];
          this->dims[1]=M.dims[1];
          this->mat= new Vecteur[M.dims[1]];
          for(int e=0;e<M.dims[1];e++){
             this->mat[e]=M.mat[e];
          }
        }

        Matrice& operator= (const Matrice &M);
        Vecteur& operator[] (int i);
        Matrice operator+ (const Matrice &M);
        Matrice operator- (const Matrice &M);
        Matrice operator* (const Matrice &M);
        Vecteur mvprod(Vecteur V);
        Matrice transpose();
        Matrice submat(int l1,int l2,int c1,int c2);
        float norm();
        Matrice operator* (float r);
        friend Matrice operator* (float r, Matrice M);
        Matrice outer(Vecteur v1, Vecteur v2);
        friend void givens(float x, float Z, float& c, float& s);
        friend Vecteur householder(Vecteur x, float &beta);
        friend Matrice reductridiag(Matrice &D);
        float Getmat(int i,int j);
        int* dimat();
        friend float signe(float d);
        friend void qrsym(Matrice &A, Matrice &Q);
        bool isdiag();
        bool istridiag();
        void find_p_q(Matrice T, int &p, int &q);
        Matrice qrpivot(Matrice &A, Matrice &Q);
        void svd(Matrice &A, Matrice* &tabM);
        
};

void Matrice::affiche(){

    cout<<"Les dimensions de la matrice: "<<"("<<this->dims[0]<<","<<this->dims[1]<<")"<<endl;
    for(int i=0;i<this->dims[0];i++){
        for(int j=0;j<this->dims[1];j++){
             std::cout<<(this->mat[j])[i]<<" ";
        }
        std::cout<< std::endl;
    }
}



Matrice& Matrice::operator= (const Matrice &M){

    if(this != &M){
        delete[] this->mat;
        this->dims[0]=M.dims[0];
        this->dims[1]=M.dims[1];
        this->mat= new Vecteur[M.dims[1]];
        for(int e=0;e<M.dims[1];e++){
             this->mat[e]=M.mat[e];
        }
    }
    return *this;
}

Vecteur& Matrice::operator[] (int i){

    return this->mat[i];
}

Matrice Matrice::operator+ (const Matrice &M){

    if(this->dims[0] != M.dims[0] || this->dims[1] != M.dims[1]){
         throw std::invalid_argument("Erreur dans les dimension.");
    }
    else{
        Matrice A(M.dims[0],M.dims[1]);
        for(int e=0;e<M.dims[1];e++){
            A.mat[e]=this->mat[e]+M.mat[e];
        }
        return A;
    }

}


Matrice Matrice::operator- (const Matrice &M){

    if(this->dims[0] != M.dims[0] || this->dims[1] != M.dims[1]){
         throw std::invalid_argument("Erreur dans les dimension.");
    }
    else{
        Matrice A(M.dims[0],M.dims[1]);
        for(int e=0;e<M.dims[1];e++){
            A.mat[e]=this->mat[e]-M.mat[e];
        }
        return A;
    }

}


Matrice Matrice::operator* (const Matrice &M){

    if(this->dims[1]!=M.dims[0]){
        throw std::invalid_argument("Erreur dans les dimension.");
    }
    else{
        
        int dim_0= this->dims[0];
        int dim_1= M.dims[1];
        Matrice A(dim_0,dim_1);
        for(int j=0;j<dim_1;j++){
            for(int i=0;i<dim_0;i++){
                float S=0;
                for(int k=0;k<this->dims[1];k++){
                    S=S+(this->mat[k][i]*M.mat[j][k]);
                }
                A.mat[j][i]=S;
            }
        }
        return A;

    }
}


Vecteur Matrice::mvprod(Vecteur V){
     
     if(this->dims[1] != V.Getdim()){
         throw std::invalid_argument("Erreur dans les dimension.");
     }
     else{
         Vecteur vect(this->dims[0]);
         for(int e=0;e<this->dims[1];e++){
             vect= vect+(V.Gettab(e)*this->mat[e]);
         }
         return vect;
     }

}

Matrice Matrice::transpose(){

    Matrice A(this->dims[1], this->dims[0]);
    for(int i=0;i<this->dims[1];i++){
        for(int j=0;j<this->dims[0];j++){
            A.mat[j][i]=this->mat[i][j];
        }
    }
    return A;
}

Matrice Matrice::submat(int l1, int l2, int c1, int c2){


    Matrice A(l2-l1+1,c2-c1+1);
    for(int e=0; e<l2-l1+1;e++){
        for(int k=0;k<c2-c1+1;k++){
            A.mat[k][e]=this->mat[c1-1+k][l1-1+e];
        }
    }
    return A;
}


float Matrice::norm(){

    float S=0;
    for(int i=0;i<this->dims[0];i++){
        for(int j=0;j<this->dims[1];j++){
              S=S+(this->mat[j][i]*this->mat[j][i]);
        }
        
    }
    return sqrt(S);
}

Matrice Matrice::operator* (float r){

    Matrice A(this->dims[0],this->dims[1]);
    for(int e=0;e<this->dims[1];e++){
        A.mat[e]=r*this->mat[e];
    }
    return A;
}

Matrice operator* (float r, Matrice M){

    return M*r;
}

Matrice outer(Vecteur v1, Vecteur v2){

    if(v1.Getdim() != v2.Getdim()){
        throw std::invalid_argument("Erreur dans les dimension.");
    }
    else{

        int dimension=v1.Getdim();
        Vecteur *tab;
        tab= new Vecteur[dimension];
        int tab_dim[2];
        for(int i=0;i<dimension;i++){
            tab[i]= v2.Gettab(i)*v1;
        }
        tab_dim[0]=dimension;
        tab_dim[1]=dimension;
        Matrice A(tab,tab_dim);

        return A;

    }
}

void givens(float x, float Z, float &c, float &s){

    if(Z==0){
        c=1;
        s=0;
    }
    else{
        if(abs(Z)>abs(x)){
            float tau=-x/Z;
            s=1/sqrt(1+tau*tau);
            c=s*tau;
        }
        else{
            float tau=-Z/x;
            c=1/sqrt(1+tau*tau);
            s=c*tau;
        }
    }
}


Vecteur householder(Vecteur x, float &beta){

    int n=x.Getdim();
    float sigma=0;
    for(int i=1;i<n;i++){
        sigma=sigma+x.Gettab(i)*x.Gettab(i); //NON
    }
    float tab[n];
    tab[0]=1;
    for(int i=1;i<n;i++){
        tab[i]=x.Gettab(i);
    }
    Vecteur v(tab,n);

    if(sigma==0 & x.Gettab(0)>=0){
        beta=0;
    }
    else if(sigma==0 & x.Gettab(0)<0){
        beta=2;
    }
    else{
        float mu=sqrt(x.Gettab(0)*x.Gettab(0)+sigma);
        if(x.Gettab(0)<=0){
            v[0]=x.Gettab(0)-mu;
        }
        else{
            v[0]=-sigma/(x.Gettab(0)+mu);
        }
        beta= (2*v.Gettab(0)*v.Gettab(0))/(sigma+v.Gettab(0)*v.Gettab(0));
        v=(1/v.Gettab(0))*v;
    }

    return v;
    
}

float Matrice::Getmat(int i, int j){

    return this->mat[j][i];
}

int* Matrice::dimat(){

    int *dim;
    dim= new int[2];
    dim[0]= this->dims[0];
    dim[1]= this->dims[1];
    return dim;
}

float signe(float d){

    if(d>=0){
        return 1;
    }
    else{
        return -1;
    }
}


bool Matrice::isdiag(){

    int n=this->dims[0];
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(this->mat[j][i] != 0 & i!=j){
                return false;
            }
            if(this->mat[j][i]==0 & i==j){
                return false;
            }
        }
    }
    return true;
}

bool Matrice::istridiag(){

    int n=this->dims[0];
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(this->mat[j][i]==0 & i==j){
                return false;
            }
            if(this->mat[j][i] ==0 & i==j+1){
                return false;
            } 
            if(this->mat[j][i]==0 & j==i+1){
                return false;
            }
            if(this->mat[j][i]!=0 & j>i+1){
                return false;
            }
            if(this->mat[j][i]!=0 & i>j+1){
                return false;
            }
        }
    }

    return true;
}


void find_p_q(Matrice T, int &p, int &q){

    int n= T.dimat()[0];
    for(int i=1;i<n;i++){
        for(int j=1;i<n;j++){
            if(i+j<n){
                Matrice T1=T.submat(1,i,1,i);
                Matrice T2=T.submat(i+1,n-j,i+1,n-j);
                Matrice T3=T.submat(n-j+1,n,n-j+1,n);
                if(T1.isdiag()==true & T2.istridiag()==true & T3.isdiag()==true){
                      p=i;
                      q=j;
                      break;
                }
            }  
        
        }
    }
}




Matrice reductridiag(Matrice &D){

   int l= D.dimat()[0];
   float d= (D.Getmat(l-2,l-2)-D.Getmat(l-1,l-1))/2.0;
   float mu= (D.Getmat(l-1,l-1)-D.Getmat(l-1,l-2)*D.Getmat(l-1,l-2))/(d+signe(d)*sqrt(d*d+D.Getmat(l-1,l-2)*D.Getmat(l-1,l-2)));
   float x=D.Getmat(0,0)-mu;
   float Zm=D.Getmat(1,0);
   float tab[l];
   for(int i=0;i<l;i++){
       tab[i]=1;
   }
   Vecteur V(tab,l);
   Matrice Z(V);
   
   for(int k=0;k<l-1;k++){
       float c;
       float s;
       givens(x,Zm,c,s);
       //Application de la rotation de Givens à droite
       for(int j=0;j<l;j++){
           float tau_1=D.Getmat(j,k);
           float tau_2=D.Getmat(j,k+1);
           D[k][j]=c*tau_1-s*tau_2;
           D[k+1][j]=s*tau_1+c*tau_2;
           tau_1=Z[k][j];
           tau_2=Z[k+1][j];
           Z[k][j]=c*tau_1-s*tau_2;
           Z[k+1][j]=s*tau_1+c*tau_2;
        }

       //Application de la rotation de Givens à gauche

       for(int j=0;j<l;j++){
           float tau_1=D.Getmat(k,j);
           float tau_2=D.Getmat(k+1,j);
           D[j][k]=c*tau_1-s*tau_2;
           D[j][k+1]=s*tau_1+c*tau_2;
        }

       if(k<l-2){
           x=D.Getmat(k+1,k);
           Zm=D.Getmat(k+2,k);
        }
   }

   return Z;

}

void qrsym(Matrice &A, Matrice &Q){

    int n= A.dimat()[0];
    float tab[n];
    for(int e=0;e<n;e++){
        tab[e]=1;
    }
    Vecteur V(tab,n);
    Q=Matrice(V);

    //Phase 1: Tridiagonalisation de la matrice A
    
    for(int k=0;k<n-2;k++){
        
        float beta;
        Vecteur v=householder(A[k].subvec(k+1,n),beta);
        Vecteur p=beta*((A.submat(k+1,n,k+1,n)).mvprod(v));
        float tran=0;
        for(int e=0;e<n-k;e++){
            tran=tran+p.Gettab(e)*v.Gettab(e);
        }
        Vecteur w=p-((beta/2)*tran)*v;
        A[k][k+1]= (A.submat(k+1,n,k,k)).norm();
        A[k+1][k]=A.Getmat(k+1,k);
        int tab[4]={k+1,n,k+1,n};
        A.submat(k+1,n,k+1,n)=A.submat(k+1,n,k+1,n)-outer(v,w)-outer(w,v);
        Q.submat(k+1,n,k+1,n)=Q.submat(k+1,n,k+1,n)-beta*(outer(v,v)*Q.submat(k+1,n,k+1,n));

    }

    Matrice T(n,n);
    for(int j=n-1;j>=0;j--){
        T[j][j]=A[j][j];
        if(j>1){
            T[j][j-1]=A.Getmat(j,j-1);
            T[j-1][j]=T[j][j-1];
        }
    }
    //Phase 2: Diagonalisation de T et mise à jour de Q

    while(T.isdiag()==false){
        for(int i=0;i<n-1;i++){
             //Gestion des erreurs numériques
            if(abs(T.Getmat(i+1,i))+abs(T.Getmat(i,i+1))<= pow(10,-9)*(abs(T.Getmat(i,i))+abs(T.Getmat(i+1,i+1)))){
                T[i+1][i]=0;
                T[i][i+1]=0;
            }
             
        }
         
        int p;
        int q;
        find_p_q(T,p,q);
        if(p+q<n){
            Matrice T2=T.submat(p+1,n-q,p+1,n-q);
            Matrice Z=reductridiag(T2);
            Matrice T_hat=T;
            T.submat(p+1,n-q,p+1,n-q)=T2;
            T=0.5*(T_hat+T_hat.transpose());

            Matrice Q_bis(n,n);
            float i_tab[p];
            for(int e=0;e<p;e++){
                i_tab[e]=1;
            }
            Vecteur V(i_tab,p);
            Matrice I_1(V);
            float i_tab_bis[q];
            for(int e=0;e<q;e++){
                i_tab_bis[e]=1;
            }
            Vecteur V_bis(i_tab_bis,q);
            Matrice I_2(V_bis);
            Q_bis.submat(1,p,1,p)=I_1;
            Q_bis.submat(p+1,n-q,p+1,n-q)=Z;
            Q_bis.submat(n-q+1,n,n-q+1,n)=I_2;
            Q=Q*Q_bis;
        }

    }

}



Matrice qrpivot(Matrice &A, Matrice &Q){

    //Initialisation de la matrice PI, on lui affecte la matrice identité

    int n= A.dimat()[1];
    int m=A.dimat()[0];
    float tab[n];
    for (int i=0;i<n;i++){
        tab[i]=1;
    }  
    Vecteur V(tab,n);
    Matrice PI(V);

    //

    float c[n];
    for(int j=0;j<n;j++){
        c[j]=dot(A[j].subvec(1,m),A[j].subvec(1,m));
    }
    
    int r=-1;
    float tau;
    float max=c[0];
    for(int i=1;i<n;i++){
        if(max<=c[i]){
            max=c[i];
        }
    }
    tau=max;
    while(tau>0 && r<n-1){
        
        r=r+1;
        int k=r;
        double max_value=c[r];
        for(int i=r+1;i<n;i++){
            if(c[i]>max_value){
                max_value=c[i];
                k=i;
            }
        }
        
        Matrice temp=A.submat(1,m,r+1,r+1);
        A.submat(1,m,r+1,r+1)=A.submat(1,m,k+1,k+1);
        A.submat(1,m,k+1,k+1)=temp;
        float temp_1=c[k];
        c[r]=c[k];
        c[k]=temp_1;
        Matrice temp_2=PI.submat(1,n,r+1,r+1);
        A.submat(1,n,r+1,r+1)=A.submat(1,n,k+1,k+1);
        A.submat(1,n,k+1,k+1)=temp_2;
        Vecteur v;
        float beta;
        v=householder(A[r].subvec(r+1,m),beta);
        A.submat(r+1,m,r+1,n)=A.submat(r+1,m,r+1,n)-beta*outer(v,v)*A.submat(r+1,m,r+1,n);
        A[r].subvec(r+2,m)=v.subvec(2,m-r+2);
        for(int i=r+1;i<n;i++){
            c[i]=c[i]*A.Getmat(r,i)*A.Getmat(r,i);

        }

        float tau;
        if(r<n){
            float max=c[r+1];
            for(int i=r+2;i<n;i++){
               if(max<=c[i]){
                   max=c[i];
                }
            }
           tau=max;
        }

        else{
            tau=0;
        }


    }

    //Calcul de Q

    int n_1= A.dimat()[0];
    float tab_p[n_1];
    for (int i=0;i<n_1;i++){
        tab_p[i]=1;
    }  
    Vecteur V_bis(tab_p,n_1);
    Matrice Q_bis(V_bis);
    Q=Q_bis;

    Vecteur v(m);
    for(int j=n_1-1;j>=0;j--){
        v[j]=1;
        v.subvec(j+2,m)=A[j].subvec(j+2,m);
        float beta=2/(1+norm(A[j].subvec(j+2,m))*norm(A[j].subvec(j+2,m)));
        Q.submat(j+1,m,j+1,m)=Q.submat(j+1,m,j+1,m)-beta*outer(v.subvec(j+1,m),v.subvec(j+1,m))*Q.submat(j+1,m,j+1,m);

    }

    return PI;


}






