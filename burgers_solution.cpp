//#include <iostream>
#include <fstream>
#include <print>
#include <cmath>
#include <string>

using namespace std;



// global variables
const int idim = 1000; // Max iterations
const float timeto = 1.0e-7; // numerical zero

const string input_file_path = "params.ini";
const string output_file_path = "output/output.txt";


struct SimulationParams {
    bool first_output = true;

    int test_number = 0, n_cells = 0, output_freq = 0, max_iter = 0;
    float cflcoe = 0.0f, domlen = 0.0f;
    float time_max = 0.0f, time_now = 0.0f;
    float dx = 0.0f, dt = 0.0f;
    float smax = 0.0f; // another magic number

    float flux[idim + 2];
    float u[idim + 2];
};


// prototypes
int reader (SimulationParams& params);
int output (SimulationParams& params);
int initia (SimulationParams& params);
int set_boundary_conditions (SimulationParams& params);
int cflcon (SimulationParams& params);
int compute_fluxes (SimulationParams& params);
int update (SimulationParams& params);
int riemann (float ul, float ur, float& ustar);



int main (void)
{
    SimulationParams params;

    // read parameters
    reader(params);

    println("{}", params.test_number);

    // initallizing conditions
    initia(params);

    // main loop -> called "time_now marching procedure" in the original code
    params.time_now = 0.0f;
    for (int n=1; n <= params.max_iter; n++)
    {
        // setting boundary conditions 
        // only periodic boundaries are considered at the moment
        set_boundary_conditions(params);
        
        // Courant-Friedrichs-Lewy (CFL) condition imposed
        cflcon(params);
        params.time_now += params.dt;

        // intercell numerical fluxes are computed
        compute_fluxes(params);

        // updating solution with conservative formula
        update(params);

        if (n % params.output_freq == 0)
        {
            println("{:<4d} {:>2.2f}", n, params.time_now);
            output(params);
        }

        float delta = params.time_now - params.time_max;
        if (fabs(delta) <= timeto) // weird definition, but okay
        {
            output(params);
            println("---------------------------------------------------");
            println("  Current time step {:5d}, Simulation time {:>3.2f}",
                n, params.time_now);
            break;
        }
    }

    return 0;
}



/// Reads the parameter file
/// expects a file containing the simulation parameters in the order
///           Courant number coefficient
///           Domain length
///           Test problem
///           Number of n_cells in domain
///           Output frequency to screen
///           Maximum number of time steps
///           Output time
/// with one value at a time.
/// Test problem chooses the actual test to run.
int reader (SimulationParams& params)
{

    string buffer_string;
    ifstream input_file(input_file_path);

    if (!input_file.is_open())
    {
        println("Error opening file");
        return 1;
    }

    input_file >> params.cflcoe;
    input_file >> params.domlen;
    input_file >> params.test_number;
    input_file >> params.n_cells;
    input_file >> params.output_freq;
    input_file >> params.max_iter;
    input_file >> params.time_max;
    
    input_file.close();

    println("--------------------------------");
    println("Data read in is echoed to screen");
    println("--------------------------------");
    println("CFLCOE  = {}", params.cflcoe);     // Courant number coefficient
    println("DOMLEN  = {}", params.domlen);     // Domain length
    println("ITEST   = {}", params.test_number);// Test problem
    println("CELLS   = {}", params.n_cells);    // Number of n_cells in domain
    println("NFREQ   = {}", params.output_freq);// Output frequency to screen
    println("NTMAXI  = {}", params.max_iter);   // Maximum number of time steps
    println("TIMEOUT = {}", params.time_max);   // Output time
    println("--------------------------------");

    return 0;
}



/// Outputs the u[i] field.
/// in the first time it runs, a new output.txt is created.
int output (SimulationParams& params)
{
    ofstream output_file;
    if (params.first_output)
    {
        output_file.open(output_file_path, ios::trunc);
        params.first_output = false;
    }
    else
    {
        output_file.open(output_file_path, ios::app);
    }

    if (!output_file.is_open())
    {
        println("Error while creating output file");
        return 1;
    }

    for(int i = 1; i < params.n_cells; i++)
        if (i < params.n_cells -1)
            output_file << params.u[i] << ",";
        else
            (output_file << params.u[i]);

    output_file << "\n";

    output_file.close();
    
    return 0;
}




int initia (SimulationParams& params)
{
    // Purpose: to set initial conditions for solution U
    // and initialise other variables. There are
    // two choices of initial conditions,
    // determined by ITEST
    // Local variables:
    // Name                    Description
    // ====                      ===========
    // DX                        Spatial mesh size
    // I                         Variable in do loop
    // ITEST                     Defines test problem
    // FLUX                      Array for intercell fluxes
    // U                         Array for numerical solution
    // XPOS                      Position along x-axis
    // XLEFT                     Left diaphragm
    // XMIDDL                    Middle diaphragm  
    // XRIGHT                    Right diaphragm

    //int i, test_number, n_cells
    // float domlen
    float xpos;

    params.dx = params.domlen / float(params.n_cells);

    // zeroing vectors
    for (int i = 0; i < idim + 2; i++)
    {
        params.flux[i] = 0.0f;
        params.u[i] = 0.0f;
    }

    /*
    if (test_number == 1)
    {
        // Test 1: smooth 
        xpos = -1.0;
        for (int i = 1; i < n_cells; i++)
        {
            xpos = xpos + 2.0 / float(n_cells);
            u[i] = expf(-8.0 * xpos * xpos);
        }
    }
    else
    {
        // Test 2: square waves

        xleft = 0.1 * domlen;
        xmiddl = 0.5 * domlen;
        xright = 0.9 * domlen;

        for (int i = 1; i < n_cells; i++)
        {
            xpos = (float(i) - 1.0) * dx;

            if (xpos < xleft)
                u[i] = -1.0f;
            else if ((xpos >= xleft) && (xpos <= xmiddl))
                u[i] = 1.0f;
            else if ((xpos > xmiddl) && (xpos <= xright))
                u[i] = 0.0f;
            else //if (xpos > xright)
                u[i] = -1.0f;
        }
    }
    */


    for (int i = 0; i <= params.n_cells + 1; i++)
    {
        xpos = i * params.dx;

        switch (params.test_number) {
            case 1:
                // Gaussian curve
                //xpos = i * dx / domlen;
                params.u[i] = expf(-100 * (xpos - 0.5) * (xpos - 0.5));
                break;

            case 2: {
                // book square waves
                const float xleft = 0.1f * params.domlen;
                const float xmiddl = 0.5f * params.domlen;
                const float xright = 0.9f * params.domlen;
                if (xpos < xleft)
                    params.u[i] = -1.0f;
                else if ((xpos >= xleft) && (xpos <= xmiddl))
                    params.u[i] = 1.0f;
                else if ((xpos > xmiddl) && (xpos <= xright))
                    params.u[i] = 0.0f;
                else
                    params.u[i] = -1.0f;
                break;
            }

            case 3:
                // smooth sin curve
                params.u[i] = sin(2.0f * M_PI * xpos);
                break;
    
            case 4:
                // simetric square wave
                params.u[i] = (xpos > 0.25f && xpos < 0.75f) ? 1.0f : 0.0f;
                break;
    
            case 5:
                // shockwave (Riemann scalar)
                params.u[i] = (xpos < 0.5f) ? 1.0f : 0.0f;
                break;
    
            case 6:
                // centered triangular wave
                if (xpos < 0.5f)
                    params.u[i] = 2.0f * xpos;
                else
                    params.u[i] = 2.0f * (1.0f - xpos);
                break;
    
            case 7:
                // narrow delta pulse
                params.u[i] = (fabs(xpos - 0.5f) < 0.0f) ? 0.0f : 1.0f;
                break;

            case 8:
                // step function centered at 0.5 domlen
                params.u[i] = xpos - 0.5f < 0.0f ? -1.0f : 0.0f;
                println("{}", xpos);
                break;
            
            case 9:
                params.u[i] = 1 - tanhf(10.0f * (xpos - 1.0f / 2.0f));
                break;
        }
    }

    return 0;
}



/// sets periodic boundaries
int set_boundary_conditions (SimulationParams& params)
{
    params.u[0] = params.u[params.n_cells];
    params.u[params.n_cells + 1] = params.u[1];

    return 0;
}


/// get dt based on current max u
int cflcon (SimulationParams& params)
{
    //float smax;
    params.smax = 0.0f;
       
    // find maximum characteristic speed
    for (int i = 0; i <= params.n_cells + 1; i++)
        if (fabs(params.u[i]) > params.smax)
            params.smax = fabs(params.u[i]);
    
    params.dt = params.cflcoe * params.dx / params.smax;

    // avoid exceding output time
    if ((params.time_now + params.dt) > params.time_max)
        params.dt = params.time_max - params.time_now;

    return 0;
}



/// updates the u[i] states based on current fluxes
int update (SimulationParams& params)
{
    float dtodx = params.dt / params.dx;

    for (int i = 1; i <= params.n_cells; i++)
        params.u[i] = params.u[i] + dtodx * (params.flux[i-1] - params.flux[i]);

    return 0;
}



/// computes intercell fluxes according to the Godunov first-order upwind
/// method, in conjunction with the exact Riemman solver
/// gets flux and u
int compute_fluxes (SimulationParams& params)
{
    float ul, ur, ustar;

    // compute intercell flux flux[i], 0 <= i <= n_cells
    for (int i = 0; i <= params.n_cells; i++)
    {
        // defining states ul and ur for local riemann problem RP(UL, UR)
        ul = params.u[i];
        ur = params.u[i+1];

        riemann(ul, ur, ustar);

        params.flux[i] = 0.5 * ustar * ustar;
    }
    return 0;
}



/// Calculates the Riemann problem across ul and ur
int riemann (float ul, float ur, float& ustar)
{
    float s;

    if (ul > ur)
    {
        // shockwave solution
        // computing shock speed s
        s = 0.5 * (ul + ur);

        // sample the state along the t-axis
        if (s >= 0.0f)
            ustar = ul;
        else
            ustar = ur;
    }
    else
    {
        // rarefaction wave solution with 3 cases
        if (ul >= 0.0f)
            ustar = ul; // right supersonic rarefaction

        if (ur <= 0.0f)
            ustar = ur; // left supersonic rarefaction

        if((ul <= 0.0f) && (ur >= 0.0f))
            ustar = 0.0; // transonic rarefaction
    }

    return 0;
}
