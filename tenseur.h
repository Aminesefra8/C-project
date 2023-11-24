
#include "BIS.h"

class Tenseur{

    int ordre;
    int* dimensions;
    int nbelts;
    Vecteur vectTenseur;

    public:

        //Constructeur

        Tenseur(int* pDimsTab, int pSiZeTab);
        Tenseur(int* pDimsTab, int pSiZeTab, Vecteur pVect);
        Tenseur(int* pDimsTab, int pSiZeTab, int k, Matrice A);

        //Destructeur

        ~Tenseur();

        //Constructeur de recopie
        
        Tenseur(const Tenseur& pTenseur);


        const Tenseur& operator=(const Tenseur& pTenseur);
	    float& operator[](int pIndice);
	    Tenseur operator+(const Tenseur& pTenseur);
	    Tenseur operator-(const Tenseur& pTenseur);

        
        Matrice mode(int k);
        void affiche_Tenseur();
        
        static Tenseur pmod(Tenseur& S, int k, Matrice& M);
        void CalculeTenseur(int i, int k, int& ligne, int& colonne);
        



};



