#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

/////////////////////////////////////////////
//2nd-order finite-volume implementation of linear advection
//with piecewise linear slope reconstruction
//a_t+u*a_x=0


/////////////////////////////////////////////
//limiter is to reduce the slope near extreme
//in case of overshoot or undershoot see Page76 of Zingale
double minmod(double a, double b)
{
    if(abs(a) < abs(b) && a*b > 0.0)
    {
        return a;
    }
    else if(abs(b) < abs(a) && a*b > 0.0)
    {
        return b;
    }
    else
    {
        return 0.0;
    }
}

/////////////////////////////////////////
//fill the boundary
void fill_ghostcells(double a[],int indlo,int indhi,int ng)
{
///////////////fill the boundary///////////////
    for(int i=0;i<ng;i++)
    {
        a[indlo-1-i] = a[indhi-i]; //left
        a[indhi+1+i] = a[indlo+i];
    }
}

////////////////////////////////////////////
//calculate the left and right interface states
void cal_states_update(double a[], double anew[], double u, double dx, double dt, int indlo, int indhi,int nx, int ng)
{
    double slope[nx+2*ng];
    double al[nx+2*ng];
    double ar[nx+2*ng];
    double flux[nx+2*ng];
    int limitertype = 1;//0 = centered difference 1 = minmod
    for(int i=0;i<nx+2*ng;i++)
    {
        slope[i] = 0.0;
        al[i] = 0.0;
        ar[i] = 0.0;
        flux[i] = 0.0;
    }
    //unlimited centered difference slopes
    if(limitertype==0)
    {
        for(int i=indlo-1;i<indhi+2;i++)
        {
            slope[i] = 0.5*(a[i+1]-a[i-1])/dx;
        }
    }

    if(limitertype==1)
    {
        for(int i=indlo-1;i<indhi+2;i++)
        {
            slope[i] = minmod((a[i]-a[i-1])/dx,(a[i+1]-a[i])/dx);
        }
    }

    //left and right state on the current interface comes from zone i-1 and zone i
    for(int i=indlo; i<indhi+2; i++)
    {
        al[i] = a[i-1] + 0.5*dx*(1.0-u*dt/dx)*slope[i-1];
        ar[i] = a[i] - 0.5*dx*(1.0+u*dt/dx)*slope[i];
    }
//riemann solver for advection -- this is simply upwinding
    if(u>0.0)
    {
        for(int i=0;i<nx+2*ng;i++)
        {
            flux[i] = u*al[i];
        }
    }
    else
    {
        for(int i=0;i<nx+2*ng;i++)
        {
            flux[i] = u*ar[i];
        }
    }
    //update the solution conservatively
    
///////////////////////////////////////////
    for(int i=indlo; i<=indhi; i++)
    {
        anew[i] = a[i] + dt/dx*(flux[i]-flux[i+1]);
    }
}

///////////////////////////////////

int main()
{
    //assign the initial conditions and grid 
    int nx = 64;//nums of cells not include the ghost cells
    int ng = 2;//nums of ghost cells each boundary
    double xmin = 0.0;//left boundary of domian
    double xmax = 1.0;//right boundary of domian
    double u = 1.0;//veloctiy of the physical quantity >0 is go right
    double cfl = 0.7;//cfl num
    double time = 0.0;//begin time
    int difftype = 0;//0==tophat 1=sine 2=gaussian
    int index_low = ng;//the first domain cell index
    int index_high = ng+nx-1;//the last domian cell index
    double dx = (xmax-xmin)/nx;//the cell width
    double x[nx+2*ng];//cell center position
    for(int i=0;i<nx+2*ng;i++)
    {
        x[i] = xmin + (i-ng+0.5)*dx;//use illustration to calculate
        //cout << (i-ng) << '\t' <<x[i] << endl;
    }
    double timestep = cfl*dx/u;
    double period = (xmax-xmin)/u;
    double tmax = 5.0 * period;// 5 period
    double xl[nx+2*ng];//cell left boundary
    for(int i=0;i<nx+2*ng;i++)
    {
        xl[i] = xmin + (i-ng)*dx;
        //cout << (i-ng) << '\t' <<xl[i] << endl;
    }
    double xr[nx+2*ng];//cell right boundary
    for(int i=0;i<nx+2*ng;i++)
    {
        xr[i] = xmin + (i-ng+1.0)*dx;
        //cout << (i-ng) << '\t' << xr[i] << endl;
    }
    double a[nx+2*ng];//the solution
    double anew[nx+2*ng];//the update solution
    for(int i=0;i<nx+2*ng;i++)
    {
        a[i] = 0.0;
        anew[i] = 0.0;
    }

    if(difftype==0) //init the profile
    {
        for(int i=0;i<nx+2*ng;i++)
        {
            if(x[i] > 0.333 && x[i] < 0.666)
            {
                a[i] = 1.0;
            }
        }
    }

//main evoulution loop
    while(time < tmax)
    {
        fill_ghostcells(a,index_low,index_high,ng);
        if(time+timestep>tmax)
        {
            timestep = tmax -time;
        }
        //get the interface states //cal riemann problem and update the solution
        cal_states_update(a, anew, u, dx, timestep, index_low, index_high, nx, ng);
        for(int i=0;i<nx+2*ng;i++)
        {
            a[i] = anew[i];
        }
        time = time + timestep;
    }

    ofstream myfile("advcetion_minmod.dat");
    myfile.precision(10);
    for(int i=ng;i<=ng+nx-1;i++)
    {
        myfile << x[i] << '\t' << a[i] << endl;
    }





    myfile.close();


    return 0;
}






































