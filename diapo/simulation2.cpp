#include <utility>
#include <iostream>
#include <math.h>
#include <SFML/Graphics.hpp>
#include <random>

//creation de la fenêtre sfml
sf::RenderWindow window(sf::VideoMode(700, 700, 32), "SFML Graphics");

const double dt = 0.04;
const double m = 700;
const double J = 500;
const double pi = 3.14159;
const int Nv = 30;

double coeff[Nv] = {};


double sgn (double x)
{if (x >= 0) return 1;
else return -1;
}


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
    return (u.x*v.x + u.y*v.y + u.z*v.z); }


vect operator - (const vect u, const vect v)
{
    return (u + (-1)*v);
}
 
vect operator ^ (const vect u, const vect v)
{
    return vect (u.y*v.z-v.y*u.z, u.z*v.x-v.z*u.x, u.x*v.y-v.x*u.y);
}

double poly (std::vector <double> l, double t)
{
	double s = 0;
	double tk = 1;
	for (int k = 0; k < l.size(); k++)
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
	return (-300000*(f-0.003)*(f-0.06));
}


void coeff_spectre (double coeff[Nv])
{
  double df = 0.0020;
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
  double v = 10, f = 0.003, df = 0.0020;
  for (int k =0; k < Nv; k++)
  {
    double phi = distribution (generator);
    v += coeff[k]*cos(2*pi*f*t + phi);
    f += df;
  }
  return v;
}


vect rot (double theta, vect e)
{
	vect u (e.x*cos(theta) - e.y*sin(theta), e.x*sin(theta) + e.y*cos(theta));
	return u;
}


vect force_voile (double beta, double delta, double kv, vect va)
{
	vect Fv(0, 0);
	std::vector <double> cl = {-3.420129, 32.100585, -84.563164, 117.32685, -91.27235, 36.880786, -5.9792068};
	std::vector <double> cd = {0.130438, 0.085364, 0.097265,0.193015};
	if (abs(beta) < abs(delta) && abs(delta) >= 20*pi/180 && abs(beta) <= 85*pi/180)
	{
	    Fv.y = kv*pow(norme(va), 2)*poly(cl, abs(abs(delta)-abs(beta)));
	    Fv.x = kv*pow(norme(va), 2)*poly(cd, abs(abs(delta)-abs(beta)));
		if (beta >= 0)
		{
			Fv = rot(-(pi/2 - beta), Fv);
		}
		if (beta < 0)
		{
			Fv.y = -Fv.y;
			Fv = rot(pi/2 - beta, Fv);
		}
	}
	return Fv;
}



void maj_voilier (vect etat[5], const vect w)
{
	//x, y, vx, vy, theta, alpha, beta
	double theta = etat[3].x;
	double omega = etat[3].y;
	double alpha = etat[4].x;
	double beta = etat[4].y;
	double kg = 2, klat = 300, klong = 18, kv = 5;
	vect v = etat[1];
	vect vr = rot (-theta, v);
	vect a = etat[2];
	vect ez (0, 0, 1);
	vect er (cos(theta), sin(theta));
	vect ex (1, 0, 0);
	vect Gg (-2,0);
	vect var = rot(-theta, w - v);
	vect Fg (-cos(alpha)*kg*pow(vr.x,2), -sin(alpha)*kg*pow(vr.x,2));
	vect Fva ( -klong*pow(vr.x,2), -klat*abs(vr.y)*vr.y);
	double delta = vect_to_angle(-1*ex, var);
	vect Fv = force_voile (beta, delta, kv, var);
	a = 1/m*rot(theta, (Fg + Fv + Fva));
	vect dom =dt*v + 0.5*dt*dt*a;
	v = v + dt*a;
	etat[0] = etat[0] + dom;
	double domega = 1/J*((Gg^Fg)*ez - 0.000001*abs(omega)*omega);
	theta += dt*omega + 0.5*dt*dt*domega;
	omega += dt*domega;
	etat[1] = v;
	etat[2] = a;
	etat[3] = vect (theta, omega, domega);
}


void commande(vect etat[5], vect w, double cap)
{
	double k = 0.40 , kd = -2.0 ;
	double theta = etat[3].x;
	vect dir (cos(theta), sin(theta));
	vect dir_cap (cos(cap), sin(cap));
	double ecart = vect_to_angle (dir, dir_cap);

	double x = k*ecart + kd*etat[3].y ;
	if (not isnan(x)) etat[4].x = x;

	//reglage de beta
	vect v = etat[1];
	vect var = rot(-theta, w - v);
	vect ex (1, 0, 0);
	double delta = vect_to_angle(-1*ex, var);
	if (abs(delta) >= 27*pi/180 && abs(delta) < 60*pi/180)
		etat[4].y = delta - 22*pi/180*sgn(delta);
	if (abs(delta) >= 60*pi/180 )
		etat[4].y = delta - 27*pi/180*sgn(delta);
	if (abs(delta) >= 90*pi/180 )
		etat[4].y = delta - 32*pi/180*sgn(delta);
	if (abs(delta) >= 120*pi/180 )
		etat[4].y = 75;
	if (abs(delta) >= 150*pi/180 )
		etat[4].y = 85;
}



void simulation_cap (int n )
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
	double cap = -0.3, phi = 90*pi/180, phi_p = 90*pi/180;
    double t = 0;
	//vecteur (x,y,vx,vy,theta,alpha,beta):alpha = gouvernail, beta = voile
    vect etat [5];
	etat[0].x = 300;
	etat[0].y = 300;
	etat[1].x = 0;
	etat[1].y = 0;
	etat[3].x = -0.3;
	etat[4].y = 1;
    vect w (0,0);
	double tnw [5000];

	//angle aléatoire
	std::default_random_engine generator;
    std::normal_distribution<double> distribution (phi, 0.5);

    for (int k = 0; k < n; k++)
    {
	if (k % 500 == 0) phi_p = distribution(generator);
	//phi += (phi_p - phi)/100;
	//double nw = 10;
	double nw = maj_vent(t);
    w.x = nw*cos(phi);
    w.y = nw*sin(phi);
	//if (nw < 9.7) cap = -0.54;
	if (norme(etat[1]) > 8) cap = -0.3;
    maj_voilier (etat, w);
    commande (etat, w, cap);
    t += dt;
	int x = floorf(etat[0].x);
	int y = floorf(etat[0].y);
	if (etat[0].x > tx - 5) {etat[0].x -= tx-5; grille.loadFromFile("grille.png");}
	if (etat[0].y > ty - 5) {etat[0].y -= ty-5; grille.loadFromFile("grille.png");}
	if (etat[0].x <  5) {etat[0].x += tx-5; grille.loadFromFile("grille.png");}
	if (etat[0].y <  5) {etat[0].y += ty-5; grille.loadFromFile("grille.png");}
	if (k >= 800 ){
	  std::cout << t <<"             "<< nw << "              " <<  norme(etat[1]) << "    " << etat[3].x << "   "<< etat[0].y <<"\n";
	  //if ((tnw[k] - tnw[k-100]) < -0.15) cap = -0.54;
	}
	//if (k % 500 == 0) cap += 0.1;
	//if (k % 1000 == 0) cap += 1*pi/180;
	//std::cout << etat[3].x*180/pi + 90 << "              " << norme(etat[1]) << "\n";
	//cap += 0.0003;
	if (k >= 500) tnw[k] = nw;

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


int main ()
{
	coeff_spectre(coeff);
	simulation_cap(4000);
	int route[12] = {350,350, 600,130, 600,600, 130,130, 130,600, 350,350};
	//simulation_route(12000,route,12);
	window.close();
	return 1;
}	
    
    

