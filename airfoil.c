#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <mpi.h>

#include "data.h"
#include "vtk.h"
#include "setup.h"
#include "boundary.h"
#include "args.h"

/**
 * @brief Computation of tentative velocity field (f, g)
 * 
 */
void compute_tentative_velocity(MPI_Domain* p_mpi_domain) {
    int starting_X_index = p_mpi_domain->starting_X_index;
    int ending_X_index = p_mpi_domain->ending_X_index;

    for (int i = starting_X_index; i <= ending_X_index; i++) {
        for (int j = 1; j < jmax+1; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i+1][j] & C_F)) {
                double du2dx = ((u[i][j] + u[i+1][j]) * (u[i][j] + u[i+1][j]) +
                                y * fabs(u[i][j] + u[i+1][j]) * (u[i][j] - u[i+1][j]) -
                                (u[i-1][j] + u[i][j]) * (u[i-1][j] + u[i][j]) -
                                y * fabs(u[i-1][j] + u[i][j]) * (u[i-1][j]-u[i][j]))
                                / (4.0 * delx);
                double duvdy = ((v[i][j] + v[i+1][j]) * (u[i][j] + u[i][j+1]) +
                                y * fabs(v[i][j] + v[i+1][j]) * (u[i][j] - u[i][j+1]) -
                                (v[i][j-1] + v[i+1][j-1]) * (u[i][j-1] + u[i][j]) -
                                y * fabs(v[i][j-1] + v[i+1][j-1]) * (u[i][j-1] - u[i][j]))
                                / (4.0 * dely);
                double laplu = (u[i+1][j] - 2.0 * u[i][j] + u[i-1][j]) / delx / delx +
                                (u[i][j+1] - 2.0 * u[i][j] + u[i][j-1]) / dely / dely;
   
                f[i][j] = u[i][j] + del_t * (laplu / Re - du2dx - duvdy);
            } else {
                f[i][j] = u[i][j];
            }
        }
    }
    for (int i = starting_X_index; i <= ending_X_index; i++) {
        for (int j = 1; j < jmax; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i][j+1] & C_F)) {
                double duvdx = ((u[i][j] + u[i][j+1]) * (v[i][j] + v[i+1][j]) +
                                y * fabs(u[i][j] + u[i][j+1]) * (v[i][j] - v[i+1][j]) -
                                (u[i-1][j] + u[i-1][j+1]) * (v[i-1][j] + v[i][j]) -
                                y * fabs(u[i-1][j] + u[i-1][j+1]) * (v[i-1][j]-v[i][j]))
                                / (4.0 * delx);
                double dv2dy = ((v[i][j] + v[i][j+1]) * (v[i][j] + v[i][j+1]) +
                                y * fabs(v[i][j] + v[i][j+1]) * (v[i][j] - v[i][j+1]) -
                                (v[i][j-1] + v[i][j]) * (v[i][j-1] + v[i][j]) -
                                y * fabs(v[i][j-1] + v[i][j]) * (v[i][j-1] - v[i][j]))
                                / (4.0 * dely);
                double laplv = (v[i+1][j] - 2.0 * v[i][j] + v[i-1][j]) / delx / delx +
                                (v[i][j+1] - 2.0 * v[i][j] + v[i][j-1]) / dely / dely;

                g[i][j] = v[i][j] + del_t * (laplv / Re - duvdx - dv2dy);
            } else {
                g[i][j] = v[i][j];
            }
        }
    }

    /* f & g at external boundaries */
    for (int j = 1; j < jmax+1; j++) {
        f[0][j]    = u[0][j];
        f[imax][j] = u[imax][j];
    }
    for (int i = starting_X_index; i <= ending_X_index; i++) {
        g[i][0]    = v[i][0];
        g[i][jmax] = v[i][jmax];
    }
}


/**
 * @brief Calculate the right hand side of the pressure equation 
 * 
 */
void compute_rhs(MPI_Domain* p_mpi_domain) {
    int starting_X_index = p_mpi_domain->starting_X_index;
    int ending_X_index = p_mpi_domain->ending_X_index;
    
    for (int i = starting_X_index; i <= ending_X_index; i++) {
        for (int j = 1;j < jmax+1; j++) {
            if (flag[i][j] & C_F) {
                /* only for fluid and non-surface cells */
                rhs[i][j] = ((f[i][j] - f[i-1][j]) / delx + 
                             (g[i][j] - g[i][j-1]) / dely)
                             / del_t;
            }
        }
    }
}


/**
 * @brief Red/Black SOR to solve the poisson equation.
 * 
 * @return Calculated residual of the computation
 * 
 */
double poisson(MPI_Domain* p_mpi_domain) {
    int starting_X_index = p_mpi_domain->starting_X_index;
    int ending_X_index = p_mpi_domain->ending_X_index;
    
    double rdx2 = 1.0 / (delx * delx);
    double rdy2 = 1.0 / (dely * dely);
    double beta_2 = -omega / (2.0 * (rdx2 + rdy2));

    double local_p0 = 0.0;
    /* Calculate sum of squares */
    for (int i = starting_X_index; i <= ending_X_index; i++) {
        for (int j = 1; j < jmax+1; j++) {
            if (flag[i][j] & C_F) { local_p0 += p[i][j] * p[i][j]; }
        }
    }
    
    double global_p0;
    MPI_Allreduce(&local_p0, &global_p0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    global_p0 = sqrt(global_p0 / fluid_cells); 
    if (global_p0 < 0.0001) { global_p0 = 1.0; }

    exchange_array(p, p_mpi_domain);

    /* Red/Black SOR-iteration */
    int iter;
    double local_res = 0.0;
    double global_res;
    for (iter = 0; iter < itermax; iter++) {
        for (int rb = 0; rb < 2; rb++) {
            for (int i = starting_X_index; i <= ending_X_index; i++) {
                for (int j = 1; j < jmax+1; j++) {
                    if ((i + j) % 2 != rb) { continue; }
                    if (flag[i][j] == (C_F | B_NSEW)) {
                        /* five point star for interior fluid cells */
                        p[i][j] = (1.0 - omega) * p[i][j] - 
                              beta_2 * ((p[i+1][j] + p[i-1][j] ) *rdx2
                                  + (p[i][j+1] + p[i][j-1]) * rdy2
                                  - rhs[i][j]);
                    } else if (flag[i][j] & C_F) { 
                        /* modified star near boundary */

                        double eps_E = ((flag[i+1][j] & C_F) ? 1.0 : 0.0);
                        double eps_W = ((flag[i-1][j] & C_F) ? 1.0 : 0.0);
                        double eps_N = ((flag[i][j+1] & C_F) ? 1.0 : 0.0);
                        double eps_S = ((flag[i][j-1] & C_F) ? 1.0 : 0.0);

                        double beta_mod = -omega / ((eps_E + eps_W) * rdx2 + (eps_N + eps_S) * rdy2);
                        p[i][j] = (1.0 - omega) * p[i][j] -
                            beta_mod * ((eps_E * p[i+1][j] + eps_W * p[i-1][j]) * rdx2
                                + (eps_N * p[i][j+1] + eps_S * p[i][j-1]) * rdy2
                                - rhs[i][j]);
                    }
                }
            }

            exchange_array(p, p_mpi_domain);
        }
        
        /* computation of residual */
        for (int i = starting_X_index; i <= ending_X_index; i++) {
            for (int j = 1; j < jmax+1; j++) {
                if (flag[i][j] & C_F) {
                    double eps_E = ((flag[i+1][j] & C_F) ? 1.0 : 0.0);
                    double eps_W = ((flag[i-1][j] & C_F) ? 1.0 : 0.0);
                    double eps_N = ((flag[i][j+1] & C_F) ? 1.0 : 0.0);
                    double eps_S = ((flag[i][j-1] & C_F) ? 1.0 : 0.0);

                    /* only fluid cells */
                    double add = (eps_E * (p[i+1][j] - p[i][j]) - 
                        eps_W * (p[i][j] - p[i-1][j])) * rdx2  +
                        (eps_N * (p[i][j+1] - p[i][j]) -
                        eps_S * (p[i][j] - p[i][j-1])) * rdy2  -  rhs[i][j];
                    local_res += add * add;
                }
            }
        }

        MPI_Allreduce(&local_res, &global_res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        global_res = sqrt(global_res / fluid_cells) / global_p0;
        
        /* convergence? */
        if (global_res < eps) break;
    }

    return global_res;
}


/**
 * @brief Update the velocity values based on the tentative
 * velocity values and the new pressure matrix
 */
void update_velocity(MPI_Domain* p_mpi_domain) {
    int starting_X_index = p_mpi_domain->starting_X_index;
    int ending_X_index = p_mpi_domain->ending_X_index;

    for (int i = starting_X_index; i <= ending_X_index; i++) {
        for (int j = 1; j < jmax-1; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i+1][j] & C_F)) {
                u[i][j] = f[i][j] - (p[i+1][j] - p[i][j]) * del_t / delx;
            }
        }
    }
    
    for (int i = starting_X_index; i <= ending_X_index; i++) {
        for (int j = 1; j < jmax-2; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i][j+1] & C_F)) {
                v[i][j] = g[i][j] - (p[i][j+1] - p[i][j]) * del_t / dely;
            }
        }
    }
}


/**
 * @brief Set the timestep size so that we satisfy the Courant-Friedrichs-Lewy
 * conditions. Otherwise the simulation becomes unstable.
 */
void set_timestep_interval(MPI_Domain* p_mpi_domain) {
    /* del_t satisfying CFL conditions */
    if (tau >= 1.0e-10) { /* else no time stepsize control */
        double local_umax = 1.0e-10;
        double local_vmax = 1.0e-10; 
        
        for (int i = p_mpi_domain->starting_X_index; i <= p_mpi_domain->ending_X_index; i++) {
            for (int j = 1; j <= jmax; j++) {
                local_umax = fmax(fabs(u[i][j]), local_umax);
            }
        }

        for (int i = p_mpi_domain->starting_X_index; i <= p_mpi_domain->ending_X_index; i++) {
            for (int j = 1; j <= jmax; j++) {
                local_vmax = fmax(fabs(v[i][j]), local_vmax);
            }
        }

        // Use a reduction to find the maximum global umax
        double global_umax;
        MPI_Allreduce(&local_umax, &global_umax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        // Use a reduction to find the maximum global vmax
        double global_vmax;
        MPI_Allreduce(&local_vmax, &global_vmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        double deltu = delx / global_umax;
        double deltv = dely / global_vmax; 
        double deltRe = 1.0 / (1.0 / (delx * delx) + 1 / (dely * dely)) * Re / 2.0;

        if (deltu < deltv) {
            del_t = fmin(deltu, deltRe);
        } else {
            del_t = fmin(deltv, deltRe);
        }
        del_t = tau * del_t; /* multiply by safety factor */
    }
}

/**
 * @brief The main routine that sets up the problem and executes the solving routines routines
 * 
 * @param argc The number of arguments passed to the program
 * @param argv An array of the arguments passed to the program
 * @return int The return value of the application
 */
int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    set_defaults();
    parse_args(argc, argv);
    setup();

    allocate_arrays();
    problem_set_up();

    // Setup domain information for this process.
    int min_domain_X_span = imax / size;
    int remaining_X_cells = imax % size; // Add remaining cells to first few ranks.
    MPI_Domain process_domain;
    process_domain.rank = rank;
    process_domain.size = size;
    process_domain.owned_X_span = min_domain_X_span + (rank < remaining_X_cells ? 1 : 0);
    process_domain.starting_X_index = 1 + rank * min_domain_X_span + (rank < remaining_X_cells ? rank : remaining_X_cells);
    process_domain.ending_X_index = process_domain.starting_X_index + process_domain.owned_X_span - 1;

    if (rank == 0 && verbose) print_opts();

    double res;

    /* Main loop */
    int iters = 0;
    double t;
    for (t = 0.0; t < t_end; t += del_t, iters++) {
        if (!fixed_dt)
            set_timestep_interval(&process_domain);

        exchange_array(u, &process_domain);
        exchange_array(v, &process_domain);

        compute_tentative_velocity(&process_domain);

        exchange_array(f, &process_domain);
        exchange_array(g, &process_domain);

        compute_rhs(&process_domain);

        res = poisson(&process_domain);

        // Halo exchange p array.
        exchange_array(p, &process_domain);

        update_velocity(&process_domain);

        exchange_array(u, &process_domain);
        exchange_array(v, &process_domain);

        apply_boundary_conditions();

        if (rank == 0 && (iters % output_freq == 0)) {
            printf("Step %8d, Time: %14.8e (del_t: %14.8e), Residual: %14.8e\n", iters, t+del_t, del_t, res);
 
            if ((!no_output) && (enable_checkpoints))
                write_checkpoint(iters, t+del_t);
        }
    } /* End of main loop */

    if (rank == 0) {
        printf("Step %8d, Time: %14.8e, Residual: %14.8e\n", iters, t, res);
        printf("Simulation complete.\n");

        if (!no_output)
            write_result(iters, t);
    }

    free_arrays();

    MPI_Finalize();
}

