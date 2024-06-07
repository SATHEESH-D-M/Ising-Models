#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<unistd.h> // to get pid 
#include<mpi.h>

#define MASTER_RANK 0

/***** RNG between 0 and 1 *****/
double random_number(){
    double r;
    pid_t pid = getpid();
    srand48(time(NULL) * pid);
    r = (double)rand() / RAND_MAX;
    return r;
}

/***** Linear-space function *****/
double* temperature_points(double start, double end, int num_points) {
    double* array = (double*)malloc(num_points * sizeof(double));
    if (array == NULL) {
        printf("temp_array memory allocation failed.\n");
        exit(1);
    }

    double step = (end - start) / (num_points - 1);
    for (int i = 0; i < num_points; i++) {
        array[i] = start + step * i;
    }
    return array;
}

/***** function to create a state of lattice *****/
int*** create3DLattice(int L) {
    int*** arr = (int***)malloc(L * sizeof(int**));
    for (int i = 0; i < L; i++) {
        arr[i] = (int**)malloc(L * sizeof(int*));
        for (int j = 0; j < L; j++){
        arr[i][j] = (int*)malloc(L * sizeof(int));
        }  
    }
    return arr;
}

/***** frees up the dynamic memory allocated 
    TO AVOID MEMORY LEAKAGE *****/
void free3DLattice(int*** arr, int L) {
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            free(arr[i][j]);
        }
        free(arr[i]);
    }
    free(arr);
}

/***** Initialise the lattice *****/
void initialize_random_state(int*** state, int L){
    pid_t pid = getpid();
    srand(time(NULL) * pid);         // different seed values 
    for(int i = 0; i < L; i++){
        for(int j = 0; j < L; j++){
            for (int k = 0; k < L; k++){
                state[i][j][k] = 2 * (rand() % 2) - 1;
                // printf("%d\t", state[i][j][k]); 
            }  
            // printf("\n");  
        }   
        // printf("\n\n"); 
    }
}

/***** METROPOLIS ALGORITHM for spin flipping *****/
void metropolis(int*** state, int L, double J, float beta){
    int a = 0, b = 0;
    int up, down, left, right, front, back, spin, NN;
    float delta_E, w, r;

    for (int i = 0; i < L; i++){
        for (int j = 0; j < L; j++){
            for (int k = 0; k < L; k++){
                up = i - 1;
                down = i + 1;
                left = j - 1;
                right = j + 1;
                front = k - 1;
                back = k + 1;

                /* Periodic boundary condition */
                up = (up == -1) ? (L - 1) : up;
                left = (left == -1) ? (L - 1) : left;
                front = (front == -1) ? (L - 1) : front;
                right = (right == L) ? 0 : right;
                down = (down == L) ? 0 : down;
                back = (back == L) ? 0 : back;

                //calculating delta_E, r, w
                spin = state[i][j][k];
                NN = J * (state[down][j][k] + state[up][j][k] + state[i][right][k] +\
                             state[i][left][k] + state[i][j][front] + state[i][j][back]);
                delta_E = 2 * spin * NN;
                r = random_number();
                w = exp(-1 * delta_E * beta);

                // Metropolis selection process
                if (delta_E <= 0){
                    spin = -1 * spin;
                } else if (delta_E > 0 && r <= w){
                    spin = -1 * spin;
                } else if (delta_E > 0 && r > w){
                    spin = 1 * spin;
                }

                //replace the spin 
                state[i][j][k] = spin;

                up = down = left = right = front = back = 0;
             }  
        }
    }
}

/***** Energy calculation *****/ 
float lattice_energy(int*** state, int L, double J){
    int up, down, left, right, front, back, spin;
    float eng = 0.0;
    float NN;
    for (int i = 0; i < L; i++){
        for (int j = 0; j < L; j++){
             for (int k = 0; k < L; k++){
                up = i - 1;
                down = i + 1;
                left = j - 1;
                right = j + 1;
                front = k - 1;
                back = k + 1;

                /* Periodic boundary condition */
                up = (up == -1) ? (L - 1) : up;
                left = (left == -1) ? (L - 1) : left;
                front = (front == -1) ? (L - 1) : front;
                right = (right == L) ? 0 : right;
                down = (down == L) ? 0 : down;
                back = (back == L) ? 0 : back;
                
                spin = state[i][j][k];
                NN = J * (state[down][j][k] + state[up][j][k] + state[i][right][k] +\
                             state[i][left][k] + state[i][j][front] + state[i][j][back]);
                eng += (-1 * spin * NN);
                /*printf("eng[%d][%d][%d] = %f \n",i,j,k,eng);*/
                up = down = left = right = front = back = 0;
             }   
        }
    } 
    return eng / 2;
}

/***** Magnetisation calculation *****/
float lattice_magnetisation(int*** state, int L){
    float mag;
    mag = 0.0;
    for (int i = 0; i < L; i++){
        for (int j = 0; j < L; j++){
            for (int k = 0; k < L; k++){
                mag = mag + state[i][j][k];
            }   
        }
    }
    return mag;
}


int main(int argc, char** argv){
    time_t start_time;
    time(&start_time);

    int L = 5;                  // Lattice dimension (L x L)
    int eq_steps =   pow(10,4);
    int avg_steps =  pow(10,4);
    int t_points = 100;
    float J = 1;

    // Temperature array initiation
    double t_start = 2;
    double t_end = 7;
    double* T = temperature_points(t_start, t_end, t_points);  // Temperature points

    // Initialise the MPI
    int processes, rank;
        // the "processes" entered below here should be used in the terminal as well 
        // and it should be chosen such that "count" remains an integer
    processes = 6;
    int count = t_points / (processes - 1);
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Master process 
    if (rank == MASTER_RANK){ 
        // Send the temperature points to slave processes
        int remaining_points = t_points;
        for (int slave_rank = 1; slave_rank < processes; slave_rank++){
            MPI_Send(&T[t_points - remaining_points], count, MPI_DOUBLE, slave_rank, 0, MPI_COMM_WORLD);
            remaining_points = remaining_points - count;
        }
    } else { //(rank != 0) slave processes
        // Recieve the temperature points send from master
        double T_buffer[count];
        MPI_Recv(T_buffer, count, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        //Initialise buffer arrays in slave to hold results temporarily
        long double E_buffer[count];
        long double M_buffer[count];
        long double C_buffer[count];
        long double X_buffer[count];
        long double U_buffer[count];
    
        // Begin the Metropolis algotithm for T_buffer values
        for (int i = 0; i < count; i++) {
            double beta = 1.0/ T_buffer[i];       // Beta defenition
            
            // Create memory for dynamic array 
            int*** state = create3DLattice(L);

            // initialize the state 
            initialize_random_state(state , L);

            // Initiate the markov chain and let it reach target distribution
            for (int i = 0; i < eq_steps; i++){
                metropolis(state, L, J, beta);
            }
            
            // Quantities calculation
            long double energy, mag, Energy_sum, Mag_sum, Energy_sq_sum, Mag_sq_sum;
            long double Mag_four_sum, Mag_2_sum, Mag_4_sum;
            Energy_sum = Energy_sq_sum = Mag_sq_sum = Mag_sum =  0.0;
            Mag_four_sum = Mag_2_sum = Mag_4_sum = 0.0;
            
            for (int i = 0; i < avg_steps; i++){
                metropolis(state, L, J, beta);
                energy = lattice_energy(state, L, J);
                mag = lattice_magnetisation(state, L);

                Energy_sum = Energy_sum + energy;
                Mag_sum = Mag_sum + mag;
                Energy_sq_sum = Energy_sq_sum + (energy * energy);
                Mag_sq_sum = Mag_sq_sum + (mag * mag);
                
                Mag_2_sum = Mag_sq_sum / avg_steps;
                Mag_4_sum = Mag_4_sum + pow(mag, 4);
                Mag_four_sum = Mag_4_sum / avg_steps;
            }

            E_buffer[i] = Energy_sum / (avg_steps);
            M_buffer[i] = Mag_sum / (avg_steps);
            C_buffer[i] = ((Energy_sq_sum/avg_steps)-pow(E_buffer[i],2)) * pow(beta,2);
            X_buffer[i] = ((Mag_sq_sum/avg_steps)-pow(M_buffer[i],2)) * beta;
            U_buffer[i] = 1 - Mag_four_sum/(3 * pow(Mag_2_sum, 2));

            E_buffer[i] = E_buffer[i]/(L*L*L);
            M_buffer[i] = M_buffer[i]/(L*L*L);
            C_buffer[i] = C_buffer[i]/(L*L*L);
            X_buffer[i] = X_buffer[i]/(L*L*L);

            // free the dynamic memory alloted
            free3DLattice(state, L);
        }

        // Send the buffer results to the master
        MPI_Send(E_buffer, count, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        MPI_Send(M_buffer, count, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        MPI_Send(C_buffer, count, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
        MPI_Send(X_buffer, count, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
        MPI_Send(U_buffer, count, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
        MPI_Send(T_buffer, count, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
    }
    
    // To synchronise the outputs
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Back to master process again
    if (rank == MASTER_RANK){
        // Array initiation for final results
        double E_array[t_points];
        double M_array[t_points];
        double C_array[t_points];
        double X_array[t_points];
        double U_array[t_points];

        // Recieve the results from the slave processes
        for (int slave_rank = 1; slave_rank < processes; slave_rank++){
            MPI_Recv(&E_array[count * (slave_rank - 1)], count, MPI_DOUBLE, slave_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&M_array[count * (slave_rank - 1)], count, MPI_DOUBLE, slave_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&C_array[count * (slave_rank - 1)], count, MPI_DOUBLE, slave_rank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&X_array[count * (slave_rank - 1)], count, MPI_DOUBLE, slave_rank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&U_array[count * (slave_rank - 1)], count, MPI_DOUBLE, slave_rank, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&T[count * (slave_rank - 1)], count, MPI_DOUBLE, slave_rank, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // Print the result to the output terminal
        printf("T \t\t\t\t E \t\t\t\t M \t\t\t\t C \t\t\t\t X \t\t\t\t U\n");
        for (int i = 0; i < t_points; i++){
            printf("%f   \t %f   \t %f   \t %f   \t %f   \t %f  \n", T[i], E_array[i], M_array[i], C_array[i], X_array[i] , U_array[i]);  
        }

        time_t end_time;
        time(&end_time);
        double timeDiff = difftime(end_time, start_time);

        /***** store data as a txt file *****/ 
        FILE* file;
        char filename[] = "3d_mpi.txt"; 

        // Get the current time
        time_t current_time;
        struct tm* time_info;
        time(&current_time);
        time_info = localtime(&current_time);

        // Open the file 
        file = fopen(filename, "a");
        if (file == NULL) {
            printf("Error opening the file.\n");
            return 1;
        }

        // Write data to the file
        for (int i = 0; i < 5; i++) {
            fprintf(file, "*#*#* NEW DATA BEGINS *#*#*");
        }
        fprintf(file, "\n");
        fprintf(file, "Current Date and Time: %04d-%02d-%02d %02d:%02d:%02d\n", time_info->tm_year + 1900, time_info->tm_mon + 1, time_info->tm_mday, time_info->tm_hour, time_info->tm_min, time_info->tm_sec);
        fprintf(file, "PARAMETERS\n");
        fprintf(file, "L = %d ;  Lattice of size (L x L)\n", L);
        fprintf(file, "Eq_steps = %d\n", eq_steps);
        fprintf(file, "Avg_steps = %d\n", avg_steps);
        fprintf(file, "J = %f\n", J);
        fprintf(file, "Total points computed = %d\n", t_points);
        fprintf(file, "Runtime of the code  = %0.2f seconds\n", timeDiff);
        fprintf(file, "\n\nOUTPUT TABLE: \n");
        fprintf(file, "T\t\tE\t\tM\t\tC\t\tX\t\tU\n");
        for (int i = 0; i < t_points; i++){
            fprintf(file, "%f\t%f\t%f\t%f\t%f\t%f\n", T[i], E_array[i], M_array[i], C_array[i], X_array[i] , U_array[i]);  
        }
        fprintf(file, "\n");
        for (int i = 0; i < 5; i++) {
            fprintf(file, "-x-x-x- DATA ENDS -x-x-x-");
        }

        fprintf(file, "\n\n\n\n");

        // Close the file
        fclose(file);   
    }
   
    MPI_Finalize();

    return 0;
}