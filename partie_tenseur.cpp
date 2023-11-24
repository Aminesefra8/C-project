#include "tenseur.h"
#include <iostream>

using namespace std;


Tenseur::Tenseur(int* pDimsTab, int pSizeTab)
{
	if (pDimsTab == nullptr) {
		cerr << "\n" << __PRETTY_FUNCTION__ << "Mauvaise adresse" << endl;
		exit(EXIT_FAILURE);
	}

	dimensions = new int[pSizeTab];
	ordre = pSizeTab;

	nbelts = 1;

	for (int i = 0; i < pSizeTab; i++) {
		nbelts *= pDimsTab[i];
		dimensions[i] = pDimsTab[i];
	}

	vectTenseur = Vecteur(nbelts);
}


Tenseur::Tenseur(int* pDimsTab, int pSizeTab, Vecteur pVect)
{
	if (pDimsTab == nullptr) {
		cerr << "\n" << __PRETTY_FUNCTION__ << "Mauvaise adresse" << endl;
		exit(EXIT_FAILURE);
	}

	dimensions = new int[pSizeTab];
	ordre = pSizeTab;

	nbelts = 1;

	for (int i = 0; i < pSizeTab; i++) {
		nbelts *= pDimsTab[i];
		dimensions[i] = pDimsTab[i];
	}

	if (pVect.Getdim() != nbelts) {
		cerr << "\n" << __PRETTY_FUNCTION__ << " : La taille de pVect doit etre egale au produit des dimensions." << endl;
		exit(EXIT_FAILURE);
	}


	vectTenseur = pVect;
}


Tenseur::Tenseur(int* pDimsTab, int pSizeTab, int k, Matrice A)
{
	nbelts = 1;

	for (int i = 0; i < pSizeTab; i++) nbelts *= pDimsTab[i];

	if (pDimsTab == nullptr) {
		cerr << "\n" << __PRETTY_FUNCTION__ << "Mauvaise adresse" << endl;
		exit(EXIT_FAILURE);
	}

	if (k < 1 || k > pSizeTab) {
		cerr << "\n" << __PRETTY_FUNCTION__ << " : k doit etre compris dans l'intervalle [ 1, " << pSizeTab << " ] ( k = " << k << " )" << endl;
		exit(EXIT_FAILURE);
	}

	if (A.dimat()[0] != pDimsTab[k-1]) {
		cerr << "\n" << __PRETTY_FUNCTION__ << " : Le nombre de ligne de la matrice A doit etre de " << pDimsTab[k - 1]  << endl;
		exit(EXIT_FAILURE);
	}

	if (A.dimat()[1] != (nbelts/pDimsTab[k - 1])) {
		cerr << "\n" << __PRETTY_FUNCTION__ << " : Le nombre de colonne de la matrice A doit etre de " << nbelts / pDimsTab[k - 1] << endl;
		exit(EXIT_FAILURE);
	}

	dimensions = new int[pSizeTab];
	ordre = pSizeTab;

	for (int i = 0; i < pSizeTab; i++) {
		dimensions[i] = pDimsTab[i];
	}

	vectTenseur = Vecteur(nbelts);

	int colonne;
	int ligne;

	for (int i = 0; i < nbelts; i++) {
		
		CalculeTenseur(i, k, ligne, colonne);
		vectTenseur[i] = A[colonne][ligne];
	}

	
}


Tenseur::~Tenseur()
{
	if (dimensions != nullptr) {
		delete[] dimensions;
		dimensions = nullptr;
	}
}


Tenseur::Tenseur(const Tenseur& pTenseur)
{
	ordre = pTenseur.ordre;

	dimensions = new int[ordre];

	for (int i = 0; i < ordre; i++) dimensions[i] = pTenseur.dimensions[i];

	nbelts = pTenseur.nbelts;
	vectTenseur = pTenseur.vectTenseur;
}

const Tenseur& Tenseur::operator=(const Tenseur& pTenseur)
{
	ordre = pTenseur.ordre;

	dimensions = new int [ordre];
	for (int i = 0; i < ordre; i++) dimensions[i] = pTenseur.dimensions[i];

	nbelts = pTenseur.nbelts;
	vectTenseur = pTenseur.vectTenseur;

	return *this;
}


float& Tenseur::operator[](int pIndice)
{
	if (pIndice >= nbelts || pIndice < 0) {
		cerr << "\n" << __PRETTY_FUNCTION__ << " : L'indice du vecteur doit etre compris entre 0 et " << nbelts - 1 << "." << endl;
		exit(EXIT_FAILURE);
	}

	return vectTenseur[pIndice];
}


Tenseur Tenseur::operator+(const Tenseur& pTenseur)
{
	if (this->ordre != pTenseur.ordre) {
		cerr << "\n" << __PRETTY_FUNCTION__ << " : Les tenseurs a additionnes doivent avoir le meme ordre." << endl;
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < ordre; i++) {
		if (this->dimensions[i] != pTenseur.dimensions[i]) {
			cerr << "\n" << __PRETTY_FUNCTION__ << " : Les dimensions des tenseurs a additionnes doivent etre identiques." << endl;
			exit(EXIT_FAILURE);
		}
	}

	Tenseur tens = *this;
	tens.vectTenseur = tens.vectTenseur + pTenseur.vectTenseur;

	return tens;
}



Tenseur Tenseur::operator-(const Tenseur& pTenseur)
{
	if (this->ordre != pTenseur.ordre) {
		cerr << "\n" << __PRETTY_FUNCTION__ << " : Les tenseurs a additionnes doivent avoir le meme ordre." << endl;
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < ordre; i++) {
		if (this->dimensions[i] != pTenseur.dimensions[i]) {
			cerr << "\n" << __PRETTY_FUNCTION__ << " : Les dimensions des tenseurs a additionnes doivent etre identiques." << endl;
			exit(EXIT_FAILURE);
		}
	}

	Tenseur tens = *this;
	tens.vectTenseur = tens.vectTenseur - pTenseur.vectTenseur;
	return tens;
}


Matrice Tenseur::mode(int k)
{
	if (k < 1 || k > ordre) {
		cerr << "\n" << __PRETTY_FUNCTION__ << " : k doit etre compris dans l'intervalle [ 1, " << ordre << " ] ( k = " << k << " )" << endl;
		exit(EXIT_FAILURE);
	}

	int nbLine = this->dimensions[k - 1];
	int nbColumn = this->nbelts / nbLine;

	Matrice res(nbLine, nbColumn);

	int colonne;
	int ligne;

	for (int i = 0; i < nbelts; i++) {
		CalculeTenseur(i, k, ligne, colonne);
		res[colonne][ligne] = vectTenseur[i];
	}

	return res;
}

void Tenseur::affiche_Tenseur()
{
	vectTenseur.affiche();
}


Tenseur Tenseur::pmod(Tenseur& S, int k, Matrice& M)
{
	if (S.dimensions[k-1] != M.dimat()[1]) {
		cerr << "\n" << __PRETTY_FUNCTION__ << " : Le nombre de colonne de la matrice M doit etre egal a la dimension " << k << " du tenseur S." << endl;
		exit(EXIT_FAILURE);
	}

	if (k < 1 || k > S.ordre) {
		cerr << "\n" << __PRETTY_FUNCTION__ << " : k doit etre compris dans l'intervalle [ 1, " << S.ordre << " ] ( k = " << k << " )" << endl;
		exit(EXIT_FAILURE);
	}

	Matrice A = S.mode(k);
	Matrice B = M * A;

	
	Tenseur temp(S.dimensions, S.ordre);
	temp.dimensions[k - 1] = M.dimat()[0];
	
	return Tenseur(temp.dimensions, temp.ordre, k, B);
}


void Tenseur::CalculeTenseur(int i, int k, int& ligne, int& colonne)
{
	int f;
	int f_t;
	int f_t_next;
	int indice;

	int* indices = new int[ordre];
	int* newIndicesTab = nullptr;
	int* newDimensionTab = nullptr;

	f = i + 1;
	indice = f % dimensions[ordre - 1];
	indices[ordre - 1] = (indice == 0) ? dimensions[ordre - 1] : indice;
	f_t_next = f;

	for (int t = ordre - 1; t > 0; t--) {
		f_t = ((f_t_next - indices[t]) / dimensions[t]) + 1;
		indice = f_t % dimensions[t - 1];
		indices[t - 1] = (indice == 0) ? dimensions[t - 1] : indice;
		f_t_next = f_t;
	}


	ligne = indices[k - 1] - 1;


	int sz = ordre - 1;
	newIndicesTab = new int[sz];	
	newDimensionTab = new int[sz];	

	for (int j = 1, l = 0; j <= ordre; j++) {
		if (j == k) continue;

		newIndicesTab[l] = indices[j - 1];
		newDimensionTab[l++] = dimensions[j - 1];
	}

	colonne = newIndicesTab[sz - 1];
	int prodDim = 1;

	for (int l = sz - 1; l > 0; l--) {
		prodDim *= newDimensionTab[l];
		colonne += prodDim * (newIndicesTab[l - 1] - 1);
	}

	colonne--;
	
	delete[] indices;
	delete[] newIndicesTab;
	delete[] newDimensionTab;

	indices = nullptr;
	newIndicesTab = nullptr;
	newDimensionTab = nullptr;
}



int main(){


    int dims[3] = { 2, 2, 2 };
	Tenseur T(dims, 3);

	T.affiche_Tenseur();

	// Création et affichage du tenseur U:
	float tab1[8] = { 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f };
	Vecteur v(tab1,8);
	Tenseur U(dims, 3, v);

	U.affiche_Tenseur();

	// Création et affichage du tenseur V = U + T:
	Tenseur V = U + T;
	V.affiche_Tenseur();

	// Création et affichage du tenseur W = U - T:
	Tenseur W = U - T;
	W.affiche_Tenseur();

	// Affichage de U(2,2,2):
	cout << "\nU(2,2,2) = " << U[7] << endl;

	// Modification de U(2,2,2) en -1 et affichage de U et V:
	U[7] = -1.0f;
	U.affiche_Tenseur();
	V.affiche_Tenseur();


	// Calcul et affichage du mode 1 du tenseur T:
	Matrice A = T.mode(1);
	A.affiche();

	// Modification des coefficients du tenseur T:
	float tab2[2] = { 1.f, 4.f };
	float tab3[2] = { 3.f, 1.f / 3 };
	float tab4[2] = { 0.f, 1.5f };
	float tab5[2] = { -1.f, 2.f };

	Vecteur v1(tab2,2);
	Vecteur v2(tab3,2);
	Vecteur v3(tab4,2);
	Vecteur v4(tab5,2);
	

	Vecteur tv1[4] = { v1, v2, v3, v4 };
    int dim[2]={2,4};

	Matrice M(tv1, dim);
	int mode = 2;
	int ordre = 3;
	Tenseur temp(dims, ordre, mode, M);
	T = temp;

	T.affiche_Tenseur();

	M = T.mode(2);
	M.affiche();


	// Calcul et affichage du tenseur S = T *3 A:
	float tab6[3] = { 3.f, 0.f, 0.f };
	float tab7[3] = { -1.f, 6.f, -3.f };

	Vecteur v5(tab6,3);
	Vecteur v6(tab7,3);

	Vecteur tab_bis_2[2] = {v5,v6};
    int dim_1[2]={3,2};

	A = Matrice(tab_bis_2, dim_1); 
    A.affiche();
	Tenseur S = Tenseur::pmod(T, 3, A); 
	S.affiche_Tenseur();

	// Calul et affichage du tenseur R = S + S:
	Tenseur R = S + S;
	R.affiche_Tenseur();





    return 0;
}