#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

int main()
{
	// Defining physical properties all with units
    int rhoReel=1;
    double csReel=8;
    double kinematicViscosity=1e-5;// 0.001//0.000001
    double height= 0.01;
    double deltaXReel=height/21.0;
    double Velocity= 3.125;//3 //3000  //velo 30 fuerza 0.25 visco1e-7///
    double Fxreel=2.5;//240
	
	//Defining lattice unites
    // I'm going to use the classical Deltax, deltat, and rho0 equal to 1
    int rho0=1;
    double cs=sqrt(1.0/3.0);
	
	//Defining interesant parameters
    //double CFL=0.9; 
    double deltaTreel= deltaXReel*cs/csReel;//pow((Velocity+csReel)/(deltaXReel*CFL),-1);

    double Cx= deltaXReel;
    double Ct= deltaTreel;
    double Cu=Cx/Ct;
    double Cv=pow(Cx,2)/Ct;
    double Crho=rhoReel;
    double Cf=Crho*Cx/pow(Ct,2);
	
	double tau = 3*kinematicViscosity/Cv; //3*kinematicViscosity/Cv+0.5;
	double sigma = 1.;//1.0001;//1.02 //1.00105;
	//double k = (sigma - 1.)/(2.*tau);
	double k = (sqrt(sigma)*sqrt(sigma + 8*sigma*tau + 16*pow(tau,2))-sigma -4*tau)/(4*tau);
	double tauBar = tau/sigma*(1.+k) + 0.5;
	
	
	double Fx=Fxreel/Cf;
    double Fy=0;
	
	// Defining velocities using D2Q9
    double w[9]={4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
    int cx[9]={0,1,0,-1,0,1,-1,-1,1};
    int cy[9]={0,0,1,0,-1,1,1,-1,-1};
    double e[9][2];

    for(int i=0; i<9; i++)
    {
        e[i][0]=cx[i];
        e[i][1]=cy[i];
    }

    //Defining of lattice nodes 
    int ny = height/Cx;
    int nx=1;
	
	//Definition de la vitesse initial
    double u0=Velocity/Cu;
    double ux[nx][ny];
    double uy[nx][ny];
    //double rhoFoisUx[nx][ny];
    //double rhoFoisUy[nx][ny];

    for(int a=0;a<nx;a++)
    {
        for(int b=0;b<ny;b++)
        {
            ux[a][b]=u0;
            uy[a][b]=0;
        }
    }
    
    // Definition of equilibrium fonction, distribution fonction and post-collision function
    double  f[9][nx][ny];
    double feq[9][nx][ny];
    double fcollision[9][nx][ny];

    //Initialisation of distribution function
    for(int a=0; a<nx;a++)
    {
        for(int b=0;b<ny;b++)
        {
            for(int i=0; i<9; i++)
            {
                double A=(ux[a][b]*e[i][0]+uy[a][b]*e[i][1])/pow(cs,2);
                double B=pow(A,2)/2;
                double C=-(ux[a][b]*ux[a][b]+uy[a][b]*uy[a][b])/(2*pow(cs,2));
                feq[i][a][b]=w[i]*rho0*(1+A+B+C);
                f[i][a][b]=feq[i][a][b];
            }
        }
    } 
	

    int t=0;
    int tmax=5000000;
    double S[9][nx][ny];
    double rho[nx][ny];

    while(t<tmax)
    {
        //Definition of source term with collision
        for(int a=0; a<nx;a++)
        {
            for(int b=0;b<ny;b++)
            {
                for(int i=0; i<9; i++)
                {
                    double A=((e[i][0]-ux[a][b])*Fx+(e[i][1]-uy[a][b])*Fy)/pow(cs,2);
                    double B=(e[i][0]*ux[a][b]+e[i][1]*uy[a][b])/pow(cs,4);
                    double C=e[i][0]*Fx+e[i][1]*Fy;
                    S[i][a][b]=(1-1/(2*tauBar))*w[i]*(A+B*C);
                    fcollision[i][a][b]=f[i][a][b]+(-1./tauBar*(f[i][a][b]-feq[i][a][b])+S[i][a][b]);
                }
            }
        } 

        // Streaming function
        int ia=0, ja=0; // this are artificial nodes that I use to make the peridoic function when is necessary
        for(int a=0; a<nx;a++)
        {
            for(int b=0; b<ny;b++)
            {
                for(int i=0; i<9;i++)
                {
                    ia=a-e[i][0];
                    ja=b-e[i][1];
                    if(ia==nx){ia=0;}
                    if(ia==-1){ia=nx-1;}
                    if(ja==ny){ia=ny-1;}
                    if(ja==-1){ia=0;}
                    f[i][a][b]=fcollision[i][ia][ja];
                }
            }
        }

        // Boundary function
        for (int a =0;a<nx;a++)
        {
            f[7][a][ny-1]=fcollision[5][a][ny-1];
            f[4][a][ny-1]=fcollision[2][a][ny-1];
            f[8][a][ny-1]=fcollision[6][a][ny-1];
            f[5][a][0]=fcollision[7][a][0];
            f[2][a][0]=fcollision[4][a][0];
            f[6][a][0]=fcollision[8][a][0];
        }
        
        t+=1;

        // Update properties
        for(int a=0;a<nx;a++)
        {
            for(int b=0;b<ny;b++)
            {
                double rhoi=0;
                for(int i=0;i<9;i++)
                {
                    rhoi+=f[i][a][b];
                }
                rho[a][b]=rhoi;
                ux[a][b]=((f[1][a][b]+f[5][a][b]+f[8][a][b])-(f[3][a][b]+f[6][a][b]+f[7][a][b])+Fx/2)/rho[a][b];
                uy[a][b]=((f[2][a][b]+f[5][a][b]+f[6][a][b])-(f[4][a][b]+f[7][a][b]+f[8][a][b])+Fy/2)/rho[a][b];
            }
        }

        // update equilibrium distribution function
        for(int a=0; a<nx; a++)
        {
            for(int b=0;b<ny;b++)
            {
                for(int i=0; i<9; i++)
                {
                    double A=(ux[a][b]*e[i][0]+uy[a][b]*e[i][1])/pow(cs,2);
                    double B=pow(A,2)/2;
                    double C=-(ux[a][b]*ux[a][b]+uy[a][b]*uy[a][b])/(2*pow(cs,2));
                    feq[i][a][b]=w[i]*rho[a][b]*(1+A+B+C);
                }
            }
        }


    }
	std::cout<<"Tau = "<<tau<<'\n';
    std::cout<<"TauBar = "<<tauBar<<'\n';
    std::cout<<"deltaTreel = "<< deltaTreel<<'\n';
	std::cout<<"sigma = "<< sigma<<'\n';
	std::cout<<"k = "<<k<<std::endl;
    for (int b=0;b<ny;b++)
    {
        std::cout<<ux[(nx-1)/2][b]*Cx/Ct<<'\n';
    }
	
	std::string fileName = "visco1e-7BGK.txt";
	std::ofstream outputFile(fileName);

    if (outputFile.is_open()) {
		outputFile << "b\tvelocity\n";
        for (int b = 0; b < ny; b++) {
            // Escribe la posición y la velocidad en el archivo
            outputFile << b*deltaXReel + deltaXReel/2 << "\t" << ux[(nx-1)/2][b] * Cx / Ct << "\n";
        }
        outputFile.close();
        std::cout << "Resultados guardados en " << fileName<<std::endl;
    } else {
        std::cerr << "No se pudo abrir el archivo para escribir.\n";
    }

    return 0;
}
