#include <utility>
#include <iostream>
#include <math.h>
#include <SFML/Graphics.hpp>
#include <random>

//creation de la fenêtre sfml
sf::RenderWindow window(sf::VideoMode(700, 700, 32), "SFML Graphics");

const double dt = 0.04;
const double m = 80;
const double J = 2000;
const double pi = 3.14159;
const int Nv = 34;

double coeff[Nv] = {};


double sgn (double x)
{if (x >= 0) return 1;
else return -1;
}

//définition de la structure de vecteur

struct vect{
    double x,y,z;
    vect (){};
    vect (double i, double j): x(i), y(j), z(0) {};
    vect (double i, double j, double k): x(i), y(j), z(k){};
};
 
 
double norme (const vect u){return (sqrt(u.x*u.x+u.y*u.y));};
double N2 (const vect u){return (u.x*u.x+u.y*u.y);};
 
vect operator + (const vect u, const vect v)
{
    return vect(u.x + v.x, u.y + v.y, u.z + v.z);
}

 
vect operator * (const double k,const vect u)
{
	 return vect (k*u.x, k*u.y, k*u.z);
}
 
double operator * (const vect u, const vect v)
{
    return (u.x*v.x + u.y*v.y + u.z*v.z);
}


vect operator - (const vect u, const vect v)
{
    return (u + (-1)*v);
}
 
vect operator ^ (const vect u, const vect v)
{
    return vect (u.y*v.z-v.y*u.z, u.z*v.x-v.z*u.x, u.x*v.y-v.x*u.y);
}


//structure de réseau de Petri


struct place
{
		double jetons;
		double seuil;
		place (double s): jetons(0), seuil(s) {};
		place (double s, double j): jetons(j), seuil(s) {};
};


bool valide (place q)
{
		return (q.jetons >= q.seuil);
}

struct transition
{
		bool * valide;
		transition () {valide = NULL;}
		transition (bool b):valide(&b) {}
		void tirer();
};


void tirer (int t, std::vector< std::vector<int> > C, std::vector <place> p)
{	
		int n = C[0].size();
		for (int k = 0; k<n; k++)
		{
			p[k].jetons = p[k].jetons + C[t][k];
		}
}


void maj_rp (std::vector< std::vector<int> > C, std::vector <place> p, std::vector <transition> t)
{
		int n = C.size();
		for (int j = 0; j < n; j++)
		{
				if (*t[j].valide)
					tirer (j, C, p);
		}
}



//modélisation des fluctuations du vent

double poly (std::vector <double> l, double t)
{
	double s = 0;
	double tk = 1;
	for (int k =0; k < l.size(); k++)
	{
		s += l[k]*tk;
		tk *= t;
	}
	return (s);
}

double vect_to_angle (vect u, vect v)
{
	double angle = sgn (u.x*v.y - v.x*u.y) * acos (u*v/(norme(u)*norme(v)));
	return angle;	
}


double mod_2pi (double theta)
{ 
	while (theta > 2*pi)
		theta -= 2*pi;
	while (theta <= 0)
		theta += 2*pi;
	return theta;
}


double spectre (double f)
{
	return (-4000*(f-0.003)*(f-0.3) + 5);
}


void coeff_spectre (double coeff[Nv])
{
  double df = 0.009;
  double f = 0.003;
  for (int k = 0; k < Nv; k++)
  {
    coeff[k] = spectre(f)*df;
    f += df;
  }
}



double maj_vent (  double t)
{  
	std::default_random_engine generator;
  std::normal_distribution<double> distribution (0, 0.5);
  double v = 10, f = 0.003, df = 0.009;
  for (int k =0; k < Nv; k++)
  {
    double phi = distribution (generator);
    v += coeff[k]*cos(2*pi*f*t + phi);
    f += df;
  }
  return v;
}


vect force_voile (double beta, double delta, double kv, vect va)
{
	vect Fv(0, 0);
	std::vector <double> l = {0.615636, 1.80597, 6.23573, 7.23569, 3.98324, 1.04682, 0.105118};
	if (abs(beta) < abs(delta) && delta > 28*pi/180)
	{
	    double F = kv*pow(norme(va), 2)*poly(l, abs(delta-abs(beta)));
		if (beta >= 0)
		{
			Fv.x = F*sin(beta);
			Fv.y = -F*cos(beta);
		}
		if (beta < 0)
		{
			Fv.x = -F*sin(beta);
			Fv.y = F*cos(beta);
		}
	}
	return Fv;
}


vect rot (double theta, vect e)
{
	vect u (e.x*cos(theta) - e.y*sin(theta), e.x*sin(theta) + e.y*cos(theta));
	return u;
}


void maj_voilier (vect etat[5], const vect w)
{
	//x, y, vx, vy, theta, alpha, beta
	double theta = etat[3].x;
	double omega = etat[3].y;
	double alpha = etat[4].x;
	double beta = etat[4].y;
	double kg =2, klat = 40, klong = 1, kv = 0.8;
	vect v = etat[1];
	vect vr = rot (-theta, v);
	vect a = etat[2];
	vect ez (0, 0, 1);
	vect er (cos(theta), sin(theta));
	vect Gg (-2,0);
	vect var = rot(-theta, w - v);
	double delta = acos ((-1)*er*var/norme(var));
	vect Fg (-cos(alpha)*kg*pow(vr.x,2), -sin(alpha)*kg*pow(vr.x,2));
	vect Fva ( -klong*pow(vr.x,2), -klat*abs(vr.y)*vr.y);
	vect Fv = force_voile (beta, delta, kv, var);
	a = 1/m*rot(theta, (Fg + Fv + Fva));
	vect dom =dt*v + 0.5*dt*dt*a;
	v = v + dt*a;
	etat[0] = etat[0] + dom;
	double domega = 1/J*(Gg^Fg)*ez;
	theta += dt*omega + 0.5*dt*dt*domega;
	omega += dt*domega;
	etat[1] = v;
	etat[2] = a;
	etat[3] = vect (theta, omega, domega);
}


void commande(vect etat[5], vect w, double cap)
{
	double kp = 0.18 , kd = -20;
	double theta = etat[3].x;
	vect dir (cos(theta), sin(theta));
	vect dir_cap (cos(cap), sin(cap));
	double ecart = vect_to_angle (dir, dir_cap);

	etat[4].x = kp*ecart + kd*etat[3].y ;

	//reglage de beta
	vect v = etat[1];
	vect var = rot(-theta, w - v);
	vect er (cos(theta), sin(theta));
	double delta = vect_to_angle(er, -1*var);
	etat[4].y = delta/2;
}





void simulation (int n )
{
	//affichage fenêtre
	sf::Image grille;
	sf::Texture rendu;
	sf::Sprite affichage;
	grille.loadFromFile("grille.png");
	rendu.loadFromImage(grille);
	affichage.setTexture(rendu);
	window.clear();
	window.draw(affichage);
	window.display();
	sf::Texture img_voilier;
	img_voilier.loadFromFile("voilier.png");
	sf::Sprite voilier;
	voilier.setTexture(img_voilier);
	sf::Texture img_fleche;
	img_fleche.loadFromFile("fleche.jpg");
	sf::Sprite fleche;
	fleche.setTexture(img_fleche);
	voilier.scale(0.1, 0.1);
	voilier.setOrigin(voilier.getLocalBounds().width /2, voilier.getLocalBounds().width/2);
	fleche.setOrigin(fleche.getLocalBounds().width /2, fleche.getLocalBounds().height/2);
	fleche.scale(0.2, 0.2);

    int tx=grille.getSize().x;
    int ty=grille.getSize().y;
	fleche.setPosition(tx-100,70);
	double cap = 0.5, phi = 1.5, phi_p = 0;
    double temps = 0;
	//vecteur (x,y,vx,vy,theta,alpha,beta):alpha = gouvernail, beta = voile
    vect etat [5];
	etat[0].x = 300;
	etat[0].y = 300;
	etat[3].x = -0;
	etat[4].y = 1;
    vect w (0,0);

	//création du réseau de Petri
	std::vector <place> q;
	std::vector <transition> t;
	std::vector< std::vector<int> > pre;
	std::vector< std::vector<int> > post;
	std::vector< std::vector<int> > C;
	int n_p = q.size();
	int n_t = t.size();
	for (int i =0; i < n_t; i++)
	{for (int j = 0; j < n_p; j++) C[i][j] = pre[i][j] - post[i][j];}

	//angle aléatoire
	std::default_random_engine generator;
    std::normal_distribution<double> distribution (phi, 0.5);

    for (int k = 0; k < n; k++)
    {
	if (k % 500 == 0) phi_p = distribution(generator);
	phi += (phi_p - phi)/100;
	double nw = maj_vent(temps);
    w.x = nw*cos(phi);
    w.y = nw*sin(phi);
    maj_voilier (etat, w);
    commande (etat, w, cap);
    temps += dt;
	int x = floorf(etat[0].x);
	int y = floorf(etat[0].y);
	if (etat[0].x > tx - 5) {etat[0].x -= tx-5; grille.loadFromFile("grille.png");}
	if (etat[0].y > ty - 5) {etat[0].y -= ty-5; grille.loadFromFile("grille.png");}
	if (etat[0].x <  5) {etat[0].x += tx-5; grille.loadFromFile("grille.png");}
	if (etat[0].y <  5) {etat[0].y += ty-5; grille.loadFromFile("grille.png");}
	std::cout << nw << "  "<< norme(etat[1]) << "\n";

		//affichage
		if ( x>5 && y > 5 && x <tx-5 && y< ty-5)
		{
			grille.setPixel(x, y, sf::Color::Blue);
			grille.setPixel(x+1, y, sf::Color::Blue);
			grille.setPixel(x-1, y, sf::Color::Blue);
			grille.setPixel(x, y+1, sf::Color::Blue);
			grille.setPixel(x, y-1, sf::Color::Blue);
			voilier.setPosition(x,y);
			voilier.setRotation(etat[3].x*180/pi + 90);
			fleche.setRotation(phi*180/pi);

			rendu.update(grille);
			affichage.setTexture(rendu);
			window.clear();
			window.draw(affichage);
			window.draw(voilier);
			window.draw(fleche);
			window.display();
		}
	}
    sf::Texture texture;
    texture.create(window.getSize().x, window.getSize().y);
    texture.update(window);
    sf::Image screenshot = texture.copyToImage();
    screenshot.saveToFile("result.png");
}


void simulation_route (int n, int route[], int N )
{
	//affichage fenêtre
	sf::Image grille;
	sf::Texture rendu;
	sf::Sprite affichage;
	grille.loadFromFile("grille.png");
	rendu.loadFromImage(grille);
	affichage.setTexture(rendu);
	window.clear();
	window.draw(affichage);
	window.display();
	sf::Image copie_grille = grille;

    int tx=grille.getSize().x;
    int ty=grille.getSize().y;
	double cap = -1.2, phi = -1.5;
    double t = 0;
	//vecteur (x,y,vx,vy,theta,alpha,beta):alpha = gouvernail, beta = voile
    vect etat [5];
	etat[0].x = 10;
	etat[0].y = 10;
	etat[3].x = 0.8;
	etat[4].y = 2;
    vect w (-5,7);

	//angle aléatoire
	std::default_random_engine generator;
    std::normal_distribution<double> distribution (0, 2.5);

	int k=0, l=0;
	while (k < n && l <N)
  {
		vect u (route[l]-etat[0].x, route[l+1]-etat[0].y);
		cap = 0;
		phi += distribution(generator);
		double nw = maj_vent(t);
    w.x = nw*cos(phi);
    w.y = nw*sin(phi);
    maj_voilier (etat, w);
    //commande (etat, w, cap);
    t += dt;
	k++;
	if (norme(u) <= 10)	l +=2;

	//affichage
	if ( floorf(etat[0].x)>5 && floorf(etat[0].y) > 5 && floorf(etat[0].x) <tx-5 && floorf(etat[0].y)< ty-5)
	{
		grille.setPixel(floorf(etat[0].x), floorf(etat[0].y), sf::Color::Blue);
		grille.setPixel(floorf(etat[0].x)+1, floorf(etat[0].y), sf::Color::Blue);
		grille.setPixel(floorf(etat[0].x)-1, floorf(etat[0].y), sf::Color::Blue);
		grille.setPixel(floorf(etat[0].x), floorf(etat[0].y)+1, sf::Color::Blue);
		grille.setPixel(floorf(etat[0].x), floorf(etat[0].y)-1, sf::Color::Blue);
    }
		rendu.update(grille);
		affichage.setTexture(rendu);
		window.clear();
		window.draw(affichage);
		window.display();
	}
  grille.saveToFile("result.png");
}


int main ()
{
	coeff_spectre(coeff);
	simulation(20000);
	int route[12] = {350,350, 600,130, 600,600, 130,130, 130,600, 350,350};
	//simulation_route(12000,route,12);
	window.close();
	return 1;
}	
    
    
