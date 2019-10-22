#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

//finite-difference implementation of upwind for linear
//advection i.e. a_t+u*a_x=0
int main()
{
    int i;//loop index
    double u = 1.0;// the velocity of physical quantity >0 is go right
    double C = 0.9;//dt<=C*dx/u C is the CFL numbers
    int ncells = 65;//total 65 cells in computational domain
    //here cell in fact means point in the domian not the interval
    int nghost = 1;//each boundary owns 1 //in total 2 ghost cells
    double xmin = 0.0;//left boundary
    double xmax = 1.0;//right boundary
    double dx = (xmax-xmin)/(ncells-1.0);//65 cells = 64 intervals
    //dx the interval between two cells or two points
    double dt = C*dx/u;//time step
    double t_zhouqi = (xmax-xmin)/u;
    double t_tot = 1.0 * t_zhouqi; //one zhou qi
    double t = 0.0;//time
    int total_cells = ncells+2*nghost;
    double cell_pos[total_cells];
    double cell_val[total_cells];
    double cell_val_update[total_cells]={0.0};
    //index of the left ghost cell is from 0 to nghost-1
    //index of the most left cell is nghost
    //index of the most right cell is nghost+ncells-1
    //index of the right ghost cell is from nghost+ncells to total_cells-1
    
    //assign the domain i.e. not include the ghost cells
    for(i=0;i<ncells;i++)
    {
        cell_pos[i+nghost] = xmin + dx*i; //position of each cell
    }
    for(i=nghost;i<=nghost+ncells-1;i++)
    { //each cell's physical value
        if(cell_pos[i]<1.0/3.0)
        {
            cell_val[i] = 0.0;
        }
        else if(cell_pos[i]<2.0/3.0)
        {
            cell_val[i] = 1.0;
        }
        else
        {
            cell_val[i] = 0.0;
        }
    }

    //assign the ghost cells //here use periodic boundary conditions
    for(i=0;i<=nghost-1;i++) //left ghost cells
    {
        cell_pos[i] = -1.0; //no need for pos of ghost cells
        cell_val[i] = cell_val[nghost+ncells-1-(nghost-i)];//in case nghost large than 1
    }//use illustrtion to get the expression in []
    for(i=nghost+ncells;i<=total_cells-1;i++)
    {
        cell_pos[i] = -1.0;
        cell_val[i] = cell_val[i-ncells+1];//use illustrtion to get the expression in []
    }
    for(i=0;i<total_cells;i++) //check for initial profile
    {
        cout << cell_pos[i] << '\t' << cell_val[i] <<endl;
    }

    while(t<t_tot)
    {
        //update the BCs from last timestep
            //assign the ghost cells //here use periodic boundary conditions
        for(i=0;i<=nghost-1;i++) //left ghost cells
        {
            cell_pos[i] = -1.0; //no need for pos of ghost cells
            cell_val[i] = cell_val[nghost+ncells-1-(nghost-i)];//in case nghost large than 1
        }//use illustrtion to get the expression in []
        for(i=nghost+ncells;i<=total_cells-1;i++)
        {
        cell_pos[i] = -1.0;
        cell_val[i] = cell_val[i-ncells+1];//use illustrtion to get the expression in []
        }

        //update domain using upwind method
        for(i=nghost;i<=nghost+ncells-1;i++)
        { //eq 4.2 -4.5 in Zingale book page 61
            //cell_val[i] = cell_val[i]- C*(cell_val[i]-cell_val[i-1]);
            //the above expression is wrong as we use the newly calculated value in current timestep
            cell_val_update[i] = cell_val[i]- C*(cell_val[i]-cell_val[i-1]);
        }
        for(i=nghost;i<=nghost+ncells-1;i++)
        {
            cell_val[i] = cell_val_update[i];
        }
        t = t + dt;
    }

    ofstream myfile("up1dlinearadvectionc09.dat");
    myfile.precision(10);

    for(i=nghost;i<=nghost+ncells-1;i++)
    {
        myfile << cell_pos[i] << '\t' << cell_val[i] << endl;
    }


    myfile.close();
    return 0;
}

