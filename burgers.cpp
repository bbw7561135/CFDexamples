#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

/////////////////////////////////////////////
//2nd-order finite-volume implementation of inviscid Burger's
//with piecewise linear slope reconstruction
//u_t + u*u_x = 0 with outflow conditions


/////////////////////////////////////////////
//limiter is to reduce the slope near extreme
//in case of overshoot or undershoot see Page76 of Zingale
double max_u(double u[], int indlo,int indhi)
{
    double temp = 0.0;
    for(int i=indlo; i<=indhi; i++)
    {
        if(fabs(u[i])>temp)
        {
            temp = fabs(u[i]);
        }
    }
    return temp;
}

/////////////////////////////////////////
//fill the boundary
void fill_ghostcells(double u[],int indlo,int indhi,int ng)
{
///////////////fill the boundary///////////////
    for(int i=0;i<ng;i++)
    { //outflow
        u[indlo-1-i] = u[indlo]; //left
        u[indhi+1+i] = u[indhi];
    }
}

////////////////////////////////////////////
//calculate the left and right interface states
void cal_states_update(double u[], double unew[], double dx, double dt, int indlo, int indhi,int nx, int ng)
{
    double slope[nx+2*ng];
    double ul[nx+2*ng];
    double ur[nx+2*ng];
    double flux[nx+2*ng];
//MC limiter
    int ibegin = indlo -1;
    int iend = indhi + 1;
    double dc[nx+2*ng];
    double dl[nx+2*ng];
    double dr[nx+2*ng];
    double d1[nx+2*ng];
    double d2[nx+2*ng];
    double ldeltau[nx+2*ng];
    double shock_speed[nx+2*ng];
    double ushock[nx+2*ng];
    double urare[nx+2*ng];
    double us[nx+2*ng];
    for(int i=0;i<nx+2*ng;i++)
    {
        slope[i]=0.0;
        ul[i]=0.0;
        ur[i]=0.0;
        flux[i]=0.0;
        dc[i]=0.0;
        dl[i]=0.0;
        dr[i]=0.0;
        d1[i]=0.0;
        d2[i]=0.0;
        ldeltau[i]=0.0;
        shock_speed[i]=0.0;
        ushock[i]=0.0;
        urare[i]=0.0;
        us[i]=0.0;
    }
//////////////////////////////////////MC limiter
    for(int i=ibegin; i<=iend; i++)
    {
        dc[i] = 0.5*(u[i+1]-u[i-1]);
        dl[i] = u[i+1]-u[i];
        dr[i] = u[i] - u[i-1];
    }

    for(int i=ibegin; i<=iend; i++)
    {
        if(fabs(dl[i])<fabs(dr[i]))
        {
            d1[i] = dl[i];
        }
        else
        {
            d1[i] = dr[i];
        }
    }

    for(int i=ibegin; i<=iend; i++)
    {
        if(fabs(dc[i])<fabs(d1[i]))
        {
            d2[i] = dc[i];
        }
        else
        {
            d2[i] = d1[i];
        }
    }

    for(int i=ibegin; i<=iend; i++)
    {
        if(dl[i]*dr[i]>0.0)
        {
            ldeltau[i] = d2[i];
        }
        else
        {
            ldeltau[i] = 0.0;
        }
    }
/////////////////////////////////////////////
    for(int i=ibegin; i<=iend; i++)
    {
        ul[i+1] = u[i] + 0.5*(1.0-u[i]*dt/dx)*ldeltau[i];
        ur[i] = u[i] - 0.5*(1.0+u[i]*dt/dx)*ldeltau[i];
    }

/////Riemann solver
    for(int i=0;i<nx+2*ng;i++)
    {
        shock_speed[i]= 0.5*(ul[i]+ur[i]);
    }

    for(int i=0;i<nx+2*ng;i++)
    {
        if(shock_speed[i]>0.0)
        {
            ushock[i]=ul[i];
        }
        else
        {
            ushock[i]=ur[i];
        }
        if(shock_speed[i]<1.0e-5)
        {
            ushock[i] = 0.0;
        }
    }

    for(int i=0;i<nx+2*ng;i++)
    {
        if(ur[i] < 0.0)
        {
            urare[i]=ur[i];
        }
        else
        {
            urare[i]=0.0;
        }
    }

    for(int i=0;i<nx+2*ng;i++)
    {
        if(ul[i] > 0.0)
        {
            urare[i]=ul[i];
        }
        else
        {
            urare[i]=urare[i];
        }
    }

    for(int i=0;i<nx+2*ng;i++)
    {
        if(ul[i]>ur[i])
        {
            us[i]=ushock[i];
        }
        else
        {
            us[i]=urare[i];
        }
    }

    for(int i=indlo;i<=indhi;i++)
    {
        unew[i] = u[i] + dt/dx*(0.5*us[i]*us[i]-0.5*us[i+1]*us[i+1]);//u^n+0.5_i+0.5 in eq 6.16
    }
}

///////////////////////////////////

int main()
{
    //assign the initial conditions and grid 
    int nx = 256;//nums of cells not include the ghost cells
    int ng = 2;//nums of ghost cells each boundary
    double xmin = 0.0;//left boundary of domian
    double xmax = 1.0;//right boundary of domian
    double cfl = 0.8;//cfl num
    double time = 0.0;//begin time
    double tmax = 0.22;//(xmax-xmin)/1.0;//tmax based on period for unit velocity
    int index_low = ng;//the first domain cell index
    int index_high = ng+nx-1;//the last domian cell index
    double dx = (xmax-xmin)/nx;//the cell width
    double x[nx+2*ng];//cell center position
    for(int i=0;i<nx+2*ng;i++)
    {
        x[i] = xmin + (i-ng+0.5)*dx;//use illustration to calculate
        //cout << (i-ng) << '\t' <<x[i] << endl;
    }
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
    double u[nx+2*ng];//the solution
    double unew[nx+2*ng];//the update solution
    for(int i=0;i<nx+2*ng;i++)
    {
        u[i] = 0.0;
        unew[i] = 0.0;
    }

//init the profile
    for(int i=0;i<nx+2*ng;i++)
    {
        u[i] = 1.0;
        if(x[i] > 0.5)
        {
            u[i] = 2.0;
        }
    }
    double timestep = cfl*dx/max_u(u,index_low,index_high);
//main evoulution loop
    while(time < tmax)
    {
        fill_ghostcells(u,index_low,index_high,ng);
        timestep=cfl*dx/max_u(u,index_low,index_high);
        if(time+timestep>tmax)
        {
            timestep = tmax - time;
        }
        //get the interface states //cal riemann problem and update the solution
        cal_states_update(u, unew, dx, timestep, index_low, index_high, nx, ng);
        for(int i=0;i<nx+2*ng;i++)
        {
            u[i] = unew[i];
        }
        time = time + timestep;
    }

    ofstream myfile("burgers022.dat");
    myfile.precision(10);
    for(int i=ng;i<=ng+nx-1;i++)
    {
        myfile << x[i] << '\t' << u[i] << endl;
    }





    myfile.close();


    return 0;
}






































