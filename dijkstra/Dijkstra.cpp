#include <set>
#include <utility>
#include <iostream>
#include <fstream>
#include <math.h>
#include <SFML/Graphics.hpp>

const double infini = 10000;
const int larg = 15;
const int temps_max = 27;
const int pas_temps = 6;
const int pas_x = 1;
const int pas_y = 1;
const int Nx = 71;
const int Ny = 39;

double square (double x){return x*x;};

struct noeud {
	int i,j, temps;
	double cout;
	bool marque;
	noeud(): i(0), j(0), marque(false) {};
	noeud(int x, int y, int t ): i(x), j(y), temps(t), marque(false){};
	struct noeud * parent = NULL;
	bool operator== (const noeud& a ){return (i == a.i && j == a.j && temps == a.temps); };
};

struct comp_vent{
	double u,v;
	bool obs;
	comp_vent(): u(0), v(0), obs(false) {};
	comp_vent(double x, double y, bool z): u(x), v(y), obs(z) {};
};

struct noeud_compare{
	bool operator()(const noeud& a, const noeud& b){return a.cout < b.cout;};
};

comp_vent grille[Nx][Ny][temps_max];

noeud G[Nx][Ny][temps_max];

double vect_to_angle (double u, double v)
{
	double th;
	if ( u == 0)
	{
		if (v > 0) th = 1.57; else th = 4.71;
	}
	if (u > 0) 
	{
		if (v > 0 ) th = atan (v / u); else th = 6.28 + atan (v / u);
	}
	else th = 3.14 + atan (v / u);
	return th;
}

double poly (double x)
{
	double a5 = -0.0992;
	double a4 = 0.7048;
	double a3 = -1.71;
	double a2 = 1.519;
	double a1 = 0.072;
	double a0 = 0.01;
	return ( a5 * pow (x,5) + a4*pow(x,4) + a3*pow(x,3) + a2*pow(x,2) + a1*x + a0 );
}

double cout_arc (noeud co, noeud vo)
{
	int t = co.temps;
	double u = 0.5 * (grille[co.i][co.j][t].u + grille[vo.i][vo.j][t].u);
	double v = 0.5 * (grille[co.i][co.j][t].v + grille[vo.i][vo.j][t].v);
	double a = vect_to_angle (double (vo.i - co.i), double(vo.j - co.j));
	double b = vect_to_angle (-u, -v);
	//pour ramener l'angle dans [0,pi]
	double theta;
	if (abs(b-a) > 3.14) theta = abs(b-a) - 3.14;
	else theta = abs(b-a);
	double vi =   0.04 * sqrt(0.1 + (u*u + v*v))*poly(theta);
	double cout = sqrt (square(pas_x * (vo.i -co.i)) + square( pas_x * (vo.j -co.j)))/vi;
	return std::min(infini , cout);
}

void pop_word (std::string &str)
{
	if (not str.empty())
	{
		int k = 0;
		while (str[k] != ' ' ) k++;
		while (str[k] == ' ') k++;
		str = str.substr(k, str.length() - 1);
	}
}

void file_to_grid (comp_vent grille[Nx][Ny][temps_max])
{
	std::ifstream file;
	file.open ("all_data.txt");
	std::string line;
	std::string str;
	std::getline(file,line);
	std::getline(file, line);
	int latt_max =0;
	int longi_max =0;
	int t_max =0;
	while (line[0] != 'F')
	{
		if (line[0] == ' ')
		{
			int n = line.length();
			int latt = std::stoi (line.substr (3,8));
			double longi = std::stod (line.substr (11, 17));
			str = line.substr (19, n-1);
			std::string strr = str;
			pop_word (strr);
			pop_word (strr);
			int temps = std::stoi (strr);
			std::cout << temps << ' '<< latt << ' '<< longi <<' ';
			if (line.find ('u') != -1)
			{
					std::cout << 'u' << ' '<< int (4*(longi - 180)) / 4 << ' ' ;
				if ( str[0] == 'm' )
				{
					grille [int (4*(longi - 180)) / 4][latt - 22][temps/6 - 1].u = 0;
					grille [int (4*(longi - 180)) / 4][latt - 22][temps/6 - 1].obs = true;
				}
				else 
				{
					grille [int (4*(longi - 180)) / 4][latt - 22][temps/6 - 1].u = std::stod (str);
					grille [int (4*(longi - 180)) / 4][latt - 22][temps/6 - 1].obs = false;
				}
					std::cout << grille [int (4*(longi - 180)) / 4][latt - 22][temps/6-1].u << '\n';
			}
			else
			{
					std::cout << 'v' << ' ' ;
				if ( str[0] == 'm' )
				{
					grille [int (4*(longi - 180)) / 4][latt - 22][temps/6 - 1].v = 0;
					grille [int (4*(longi - 180)) / 4][latt - 22][temps/6 - 1].obs  = true;
				}
				else 
				{
					grille [int (4*(longi - 180)) / 4][latt - 22][temps/6 - 1].v = std::stod (str);
					grille [int (4*(longi - 180)) / 4][latt - 22][temps/6 - 1].obs = false ;
				}
					std::cout <<  grille [int (4*(longi - 180)) / 4][latt - 22][temps/6-1].v << '\n';
			}
			t_max = std::max( t_max, temps/6 - 1);
			longi_max = std::max (longi_max, int (4*(longi -180)) / 4);
			latt_max = std::max ( latt_max, latt - 22);
		}
		std::getline(file,line);
	}
	std::cout << t_max<< ' ' << longi_max << ' ' << latt_max << ' ';
}

void bordures (comp_vent grille[Nx][Ny][temps_max])
{		
		for (int i=0; i < Nx; i++)
		{
			for (int t=0; t < temps_max; t++)
		  {
				grille[i][0][t].obs = true;
				grille[i][Ny-1][t].obs = true;
				grille[i][1][t].obs = true;
				grille[i][Ny-2][t].obs = true;
		  }
		}
		for (int j=0; j < Ny; j++)
		{
			for (int t=0; t < temps_max; t++)
		  {
				grille[0][j][t].obs = true;
				grille[Nx-1][j][t].obs = true;
				grille[1][j][t].obs = true;
				grille[Nx-2][j][t].obs = true;
		  }
		}
}
				

void initialise (noeud G[Nx][Ny][temps_max], noeud deb)
{		
		for (int i=0; i < Nx; i++)
		  {for (int j=0; j < Ny; j++) 
		    {for (int t=0; t < temps_max; t++)
			  {
				(G[i][j][t]).cout = infini;		
			  }
			}
		  }
		G[deb.i][deb.j][deb.temps].cout = 0;
}

typedef std::multiset<noeud, noeud_compare>::iterator It;  

void supprimer (std::multiset<noeud, noeud_compare> &V, noeud v)
{
	std::pair<It, It> range = V.equal_range(v);
	It k = range.first;
	while ( k != range.second){
		if (v == (*k))
			V.erase (k++);
		else k++;}
}

void visuel ( noeud deb, noeud fin)
{
	int l, n;
	noeud b = G[fin.i][fin.j][fin.temps];
	l=0;
	while (not(b == deb) && b.parent != NULL)
	{
		b = *b.parent;
		l=l+1;
	}
	b = G[fin.i][fin.j][fin.temps];
	noeud chemin[l];
	for (int k = 0; k <l; k++)
	{
		noeud a = *b.parent;
		chemin[l-1-k] = a;
		b=a;
	}
	sf::Image img;
	for (int i = 1 ; i <= 6; i++)
	{
	img.loadFromFile("grille" + std::to_string(i) + ".jpg");
	int tx=img.getSize().x - 160;
	int ty=img.getSize().y - 80;
	int k = 0;
		while (k < l and chemin[k].cout < 30*i)
		{
			img.setPixel(80 + floorf(chemin[k].i * (double (tx) /Nx)), 40 + floorf(chemin[k].j * (double (ty) /Ny)), sf::Color::Red);
			img.setPixel(80 + floorf(chemin[k].i * (double (tx) /Nx)) + 1, 40 + floorf(chemin[k].j * (double (ty) /Ny)), sf::Color::Red);
			img.setPixel(80 + floorf(chemin[k].i * (double (tx) /Nx)) - 1, 40 + floorf(chemin[k].j * (double (ty) /Ny)), sf::Color::Red);
			img.setPixel(80 + floorf(chemin[k].i * (double (tx) /Nx)), 40 + floorf(chemin[k].j * (double (ty) /Ny)) + 1, sf::Color::Red);
			img.setPixel(80 + floorf(chemin[k].i * (double (tx) /Nx)), 40 + floorf(chemin[k].j * (double (ty) /Ny)) - 1, sf::Color::Red);
			img.setPixel(80 + floorf(chemin[k].i * (double (tx) /Nx)) + 1, 40 + floorf(chemin[k].j * (double (ty) /Ny)) + 1, sf::Color::Red);
			img.setPixel(80 + floorf(chemin[k].i * (double (tx) /Nx)) + 1, 40 + floorf(chemin[k].j * (double (ty) /Ny)) - 1, sf::Color::Red);
			img.setPixel(80 + floorf(chemin[k].i * (double (tx) /Nx)) - 1, 40 + floorf(chemin[k].j * (double (ty) /Ny)) + 1, sf::Color::Red);
			img.setPixel(80 + floorf(chemin[k].i * (double (tx) /Nx)) - 1, 40 + floorf(chemin[k].j * (double (ty) /Ny)) - 1, sf::Color::Red);
			k++;
		}
		img.saveToFile("result" + std::to_string(i) + ".jpg");			
	}
}


void dijkstra ( noeud deb, noeud fin)
{
	initialise (G, deb);
	// liste triée selon la distance à deb des voisins du sous-graphe minimal
	std::multiset <noeud, noeud_compare> V;
	V.insert(deb);
	noeud courant;
	while (not V.empty() and not (courant.i == fin.i and courant.j == fin.j)) 
	{
		courant = *(V.begin());
		V.erase(V.begin());
		G[courant.i][courant.j][courant.temps].marque = true;
		//mise a jour de V pour les voisins de courant
		for ( int i = -2; i <= 2; i++)
		{
			for (int j = -2; j <= 2; j++)
			{
				if ((i != 0 or j != 0) and (abs (i) != 2 or j != 0) and (i != 0 or abs (j) != 0) and (abs (i) != 2 or abs(j) != 2))
				{
					noeud v;
					v.i=courant.i + i;
					v.j=courant.j + j;
					//on calcule le cout s'il n'y a pas d'obstacle et sinon rien n'est a mettre à jour
					if (not grille[v.i][v.j][courant.temps].obs )
					{
						v.cout = courant.cout + cout_arc( courant, v);
						v.temps = int (courant.cout) / pas_temps;
						std :: cout << v.i << ' ' << v.j << ' ' << v.temps<< ' ' << cout_arc ( courant,v) << '\n';
						//le sommet est mis à jour si il n'est pas marqué, et le cout est plus bas
						if ( G[v.i][v.j][v.temps].marque == false && v.cout < G[v.i][v.j][v.temps].cout) 
						{
							 std:: cout << "ajouté!"<< '\n';
							 v.parent = &(G[courant.i][courant.j][courant.temps]);
							 G[v.i][v.j][v.temps] = v;
							 supprimer(V, G[v.i][v.j][v.temps]); 
							 V.insert(v);
						}
			    }
				}
			}
		}
	}
	visuel( deb, courant);
}



int main(){
	file_to_grid(grille);
	bordures (grille);
	noeud fin (15, 10,1);
	noeud deb (50,27,0);
	dijkstra (fin, deb);
	return 0;
}
