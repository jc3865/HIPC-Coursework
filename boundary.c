#include "data.h"
#include "boundary.h"

/**
 * @brief Given the boundary conditions defined by the flag matrix, update
 * the u and v velocities. Also enforce the boundary conditions at the
 * edges of the matrix.
 */
void apply_boundary_conditions(MPI_Domain* p_mpi_domain) {
    int rank = p_mpi_domain->rank;
    int size = p_mpi_domain->size;
    int starting_X_index = p_mpi_domain->starting_X_index;
    int ending_X_index = p_mpi_domain->ending_X_index;
    
    if (rank == 0) {
        for (int j = 0; j < jmax+2; j++) {
            /* Fluid freely flows in from the west */
            u[0][j] = u[1][j];
            v[0][j] = v[1][j];
        }
    }
    
    if (rank == p_mpi_domain->size - 1) {
        for (int j = 0; j < jmax+2; j++) {
            /* Fluid freely flows out to the east */
            u[imax][j] = u[imax-1][j];
            v[imax+1][j] = v[imax][j];
        }
    }

    for (int i = starting_X_index; i <= ending_X_index; i++) {
        /* The vertical velocity approaches 0 at the north and south
         * boundaries, but fluid flows freely in the horizontal direction */
        v[i][jmax] = 0.0;
        u[i][jmax+1] = u[i][jmax];

        v[i][0] = 0.0;
        u[i][0] = u[i][1];
    }

    /* Apply no-slip boundary conditions to cells that are adjacent to
     * internal obstacle cells. This forces the u and v velocity to
     * tend towards zero in these cells.
     */
    for (int i = starting_X_index; i <= ending_X_index; i++) {
        for (int j = 1; j < jmax+1; j++) {
            if (flag[i][j] & B_NSEW) {
                switch (flag[i][j]) {
                    case B_N: 
                        v[i][j]   = 0.0;
                        u[i][j]   = -u[i][j+1];
                        u[i-1][j] = -u[i-1][j+1];
                        break;
                    case B_E: 
                        u[i][j]   = 0.0;
                        v[i][j]   = -v[i+1][j];
                        v[i][j-1] = -v[i+1][j-1];
                        break;
                    case B_S:
                        v[i][j-1] = 0.0;
                        u[i][j]   = -u[i][j-1];
                        u[i-1][j] = -u[i-1][j-1];
                        break;
                    case B_W: 
                        u[i-1][j] = 0.0;
                        v[i][j]   = -v[i-1][j];
                        v[i][j-1] = -v[i-1][j-1];
                        break;
                    case B_NE:
                        v[i][j]   = 0.0;
                        u[i][j]   = 0.0;
                        v[i][j-1] = -v[i+1][j-1];
                        u[i-1][j] = -u[i-1][j+1];
                        break;
                    case B_SE:
                        v[i][j-1] = 0.0;
                        u[i][j]   = 0.0;
                        v[i][j]   = -v[i+1][j];
                        u[i-1][j] = -u[i-1][j-1];
                        break;
                    case B_SW:
                        v[i][j-1] = 0.0;
                        u[i-1][j] = 0.0;
                        v[i][j]   = -v[i-1][j];
                        u[i][j]   = -u[i][j-1];
                        break;
                    case B_NW:
                        v[i][j]   = 0.0;
                        u[i-1][j] = 0.0;
                        v[i][j-1] = -v[i-1][j-1];
                        u[i][j]   = -u[i][j+1];
                        break;
                }
            }
        }
    }

    if (rank == 0) {
        /* Finally, fix the horizontal velocity at the  western edge to have
        * a continual flow of fluid into the simulation.
        */
        v[0][0] = 2 * vi-v[1][0];
        for (int j = 1; j < jmax+1; j++) {
            u[0][j] = ui;
            v[0][j] = 2 * vi - v[1][j];
        }
    }
}