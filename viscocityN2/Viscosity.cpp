#include <iostream>
#include <cmath>
#include <array>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
/*
	Viscosity and Thermal Conductivity Equations for Nitrogen, Oxigen
	Argon and Air. Lemon and Jacobsen 2003
*/

namespace constants //Table 1
{
	double Tc{126.192}; 	//K
	double rhoc{11.1839};	//mol.dm^-3
	double Pc{3.3958};		//Mpa
	double M{28.01348};		//g.mol^-1
	double e_k{98.94};		//K
	double sigma{0.3656};	//nm
	double xi0{0.17};		//nm
	double BigGamma{0.055};
	double qd{0.40};		//nm
	double Tref{252.384};	//K
	
}

namespace coeffCollisionIntegral //Table 2
{
	std::vector<double> b{0.431,-0.4623,0.08406,0.005341,-0.00331};
}

namespace coeffResidualVisco //Table 3
{
	std::vector<double> N{10.72,0.03989,0.001208,-7.402,4.620};
	std::vector<double> t{0.1,0.25,3.2,0.9,0.3};
	std::vector<double> d{2.,10.,12.,2.,1.};
	std::vector<double> l{0.,1.,1.,2.,3.};
	
	/*
	In the article they say that gamma[i] is zero when l[i] is zero,
	and one when l[i] is not zero
	*/
	std::vector<double> gamma{0.,1.,1.,1.,1.};
}

namespace coeffResidualCond //Table 4
{
	std::vector<double> N{1.511, 2.117,-3.332, 8.862, 31.11, -73.13, 20.03, -0.7096, 0.2672};
	std::vector<double> t{0., -1., -0.7, 0., 0.03, 0.2, 0.8, 0.6, 1.9};
	std::vector<double> d{0, 0, 0, 1,2,3,4,8,10};
	std::vector<double> l{0, 0, 0, 0,0,1,2,2,2};
	/*
	In the article they say that gamma[i] is zero when l[i] is zero,
	and one when l[i] is not zero
	*/
	std::vector<double> gamma{0, 0, 0, 0,0,1,1,1,1};
}


double collisionIntegral(double t)
{
	double collision{0.};
	for(int i{0};i<5;i++)
	{
		collision+=coeffCollisionIntegral::b[i]*pow(log(t/constants::e_k),i);
	}
	return exp(collision);
}

double diluteViscosity(double t)
{
	return 0.0266958*std::sqrt(constants::M*t)/(pow(constants::sigma,2)*collisionIntegral(t));
}

double residualVisco(double rho, double t)
{
	double tau = constants::Tc/t;
	double delta = rho/constants::rhoc;
	double resVisc = 0.;
	
	for(int i{0};i<5;i++)
	{
		resVisc += coeffResidualVisco::N[i]*pow(tau,coeffResidualVisco::t[i])*pow(delta,coeffResidualVisco::d[i])*exp(-coeffResidualVisco::gamma[i]*pow(delta,coeffResidualVisco::l[i]));
	}
	
	return resVisc;
}

double viscosity(double rho, double t)
{
	return diluteViscosity(t)+residualVisco(rho,t);
}
///////////////////////////////////////////////////////////////
//////////////////////Conductivity/////////////////////////////
///////////////////////////////////////////////////////////////

double diluteCond(double t)
{
	double tau = constants::Tc/t;
	return coeffResidualCond::N[0]*(diluteViscosity(t)) + coeffResidualCond::N[1]*pow(tau,coeffResidualCond::t[1]) + coeffResidualCond::N[2]*pow(tau,coeffResidualCond::t[2]);
}

double residualCond(double rho, double t)
{
	double tau = constants::Tc/t;
	double delta = rho/constants::rhoc;
	double resCond = 0.;
	for(int i{3};i<9;i++)
	{
		resCond += coeffResidualCond::N[i]*pow(tau,coeffResidualCond::t[i])*pow(delta, coeffResidualCond::d[i])*exp(-coeffResidualCond::gamma[i]*pow(delta, coeffResidualCond::l[i]));
	}
	return resCond;
}

double conductivity(double rho, double t)
{
	return diluteCond(t)+residualCond(rho,t);
}

int main()
{
	std::cout<<coeffCollisionIntegral::b[0]<<std::endl;
	std::cout<<viscosity(0.,100)<<std::endl;
	std::cout<<viscosity(0.,300)<<std::endl;
	std::cout<<viscosity(25,100)<<std::endl;
	std::cout<<viscosity(10,200)<<std::endl;
	std::cout<<viscosity(5,300)<<std::endl;
	std::cout<<viscosity(11.18,126.195 )<<std::endl;
	std::cout<<"//////////////////"<<std::endl;
	std::cout<<conductivity(0.,100)<<std::endl;
	std::cout<<conductivity(0,300)<<std::endl;
	std::cout<<conductivity(25.,100)<<std::endl;
	std::cout<<conductivity(10,200)<<std::endl;
	std::cout<<conductivity(5,300)<<std::endl;
	std::cout<<conductivity(11.18,126.195)<<std::endl;
	return 0;
}