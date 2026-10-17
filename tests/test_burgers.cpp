#include <cstdlib>
#include <iostream>
#include <numeric>

// Exercise the original single-file solver without changing its CLI entry point.
#define main burgers_example_main
#include "../burgers_solution.cpp"
#undef main

void require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "FAILED: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

int main() {
    SimulationParams params;
    params.test_number = 2;
    params.n_cells = 100;
    params.cflcoe = 0.8f;
    params.time_max = 0.01f;

    for (const float length : {1.0f, 3.0f}) {
        params.domlen = length;
        initia(params);
        require(params.u[4] == -1.0f, "square-wave left plateau");
        require(params.u[25] == 1.0f, "square-wave positive plateau");
        require(params.u[75] == 0.0f, "square-wave zero plateau");
        require(params.u[95] == -1.0f, "square-wave right plateau");
        const double before = std::accumulate(params.u + 1, params.u + params.n_cells + 1, 0.0);
        set_boundary_conditions(params);
        cflcon(params);
        compute_fluxes(params);
        update(params);
        const double after = std::accumulate(params.u + 1, params.u + params.n_cells + 1, 0.0);
        require(std::abs(after - before) < 1e-5, "periodic Godunov update conserves the cell sum");
        for (int i = 1; i <= params.n_cells; ++i) {
            require(std::isfinite(params.u[i]), "updated solution remains finite");
            require(params.u[i] >= -1.0f && params.u[i] <= 1.0f,
                    "Godunov update introduces no new extrema");
        }
    }

    float star = 0.0f;
    riemann(1.0f, -0.5f, star);
    require(star == 1.0f, "right-moving shock samples the left state");
    riemann(-1.0f, 1.0f, star);
    require(star == 0.0f, "transonic rarefaction samples zero");
    riemann(-2.0f, -1.0f, star);
    require(star == -1.0f, "left-moving rarefaction samples the right state");
    std::cout << "Burgers square-wave and conservation checks passed\n";
    return EXIT_SUCCESS;
}
