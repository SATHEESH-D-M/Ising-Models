#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<unistd.h>          // to get pid 
#include<mpi.h>

#define MASTER_RANK 0
// Number of arrays
#define NUM_NEIGHBOURS 6    // No. of neighbours in the model

// Global variables - Parameters
#define L  5              // System size (L x L x L)
#define J  1.0            // Interaction strength
#define T_start  0.5
#define T_end  7   // Starting and ending temperature values
#define T_points  60            // No. of temperature points
#define eq_steps  10000
#define avg_steps  10000

/***** function to create a state of lattice *****/
int*** create3DLattice() {
    int*** arr = (int***)malloc(L * sizeof(int**));
    for (int i = 0; i < L; i++) {
        arr[i] = (int**)malloc(L * sizeof(int*));
        for (int j = 0; j < L; j++){
        arr[i][j] = (int*)malloc(L * sizeof(int));
        }  
    }
    return arr;
}

/***** frees up the dynamic memory allocated - TO AVOID MEMORY LEAKAGE *****/
void free3DLattice(int*** arr) {
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            free(arr[i][j]);
        }
        free(arr[i]);
    }
    free(arr);
}

// Global definition of Neighbour table
int *neighbour_table[NUM_NEIGHBOURS];       // Matrix with all neighbour indices

// Process the neighbour table
void Process_neighbour_table(){
    // Dynamically allocate memory for each array of neighbour table
    for (int i = 0; i < NUM_NEIGHBOURS; i++) {
        neighbour_table[i] = (int *)malloc(pow(L, 3) * sizeof(int));
    }

    // Create memory for dynamic array 
    int*** lattice_table = create3DLattice();

    // initialise the index notation (map from 3D to 1D for ease of computation)
    for(int k = 0; k < L; k++){
        for(int j = 0; j < L; j++){
            for (int i = 0; i < L; i++){
                int m = (k * pow(L, 2)) + (j * L) + i;
                lattice_table[i][j][k] = m;
                //printf("%d\t", m); 
            }  
            //printf("\n");  
        }   
        //printf("\n"); 
    }

    // Assign index values for the neighbour table in the mapped form 3D -> 1D 
    for(int k = 0; k < L; k++){
        for(int j = 0; j < L; j++){
            for (int i = 0; i < L; i++){
                int up, down, left, right, front, back;
                // coordinate definition to locate the spins in lattice
                up = j - 1;
                down = j + 1;
                left = i - 1;
                right = i + 1;
                front = k - 1;
                back = k + 1;

                /* Periodic boundary condition */
                up = (up == -1) ? (L - 1) : up;
                left = (left == -1) ? (L - 1) : left;
                front = (front == -1) ? (L - 1) : front;
                right = (right == L) ? 0 : right;
                back = (back == L) ? 0 : back;
                down = (down == L) ? 0 : down;
                
                int m = (k * pow(L, 2)) + (j * L) + i;
                // Nearest neighbours index
                neighbour_table[0][m] = lattice_table[left][j][k]; 
                neighbour_table[1][m] = lattice_table[right][j][k]; 
                neighbour_table[2][m] = lattice_table[i][up][k];
                neighbour_table[3][m] = lattice_table[i][down][k];
                neighbour_table[4][m] = lattice_table[i][j][front];
                neighbour_table[5][m] = lattice_table[i][j][back];
            }    
        }    
    }

    // free the dynamic memory alloted for 3D array
    free3DLattice(lattice_table);
}

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

/***** Initialise the lattice *****/
void initialize_random_state(int state_array[]){
    // Create memory for dynamic array 
    int*** state_3D = create3DLattice();
    
    // assign spins to 3D array
    pid_t pid = getpid();
    srand(time(NULL) * pid);         // different seed values 
    for(int i = 0; i < L; i++){
        for(int j = 0; j < L; j++){
            for (int k = 0; k < L; k++){
                state_3D[i][j][k] = 2 * (rand() % 2) - 1;
                // printf("%d\t", state[i][j][k]); 
            }  
            // printf("\n");  
        }   
        // printf("\n\n"); 
    }

    // map 3D -> 1D
    for(int k = 0; k < L; k++){
        for(int j = 0; j < L; j++){
            for (int i = 0; i < L; i++){
                int m = (k * pow(L, 2)) + (j * L) + i;
                state_array[m] = state_3D[i][j][k]; 
                //printf("%d\t", state_array[m]);
            } 
            //printf("\n");   
        }
        //printf("\n");    
    }

    // free the dynamic memory alloted for 3D array
    free3DLattice(state_3D);
}

// Function to add an element to the array
int* add_spin_index(int* list, size_t* size, int value) {
    // Increase the size of the array
    (*size)++;

    // Resize the array using realloc
    list = (int*)realloc(list, (*size) * sizeof(int));

    // Check if realloc was successful
    if (list == NULL) {
        perror("Memory allocation failed");
        exit(EXIT_FAILURE);
    }

    // Add the new element at the end
    list[*size - 1] = value;

    return list;
}

// Function to remove the first element from an array and resize it 
int* remove_spin_index(int* arr, size_t* size) {
    if (*size == 0) {
        printf("Array is empty, cannot remove the first element.\n");
        return arr;
    }

    // Create a new array with one less element
    int* newArr = malloc((*size - 1) * sizeof(int));

    if (newArr == NULL) {
        // Handle malloc failure
        perror("Memory allocation failed");
        exit(EXIT_FAILURE);
    }

    // Copy the remaining elements from the original array to the new array
    for (size_t i = 1; i < *size; ++i) {
        newArr[i - 1] = arr[i];
    }

    // Decrease the size of the array
    (*size)--;

    // Free the memory of the original array
    free(arr);

    return newArr;
}

// Function to check if a number is present in the list
int isInList(int list[], size_t size, int target) {
    for (size_t i = 0; i < size; ++i) {
        if (list[i] == target) {
            return 1;  // Number is found
        }
    }
    return 0;  // Number is not found
}

/***** Magnetisation calculation *****/
float state_magnetisation(int state_array[]){
    float mag;
    mag = 0.0;
    for (int i = 0; i < pow(L, 3); i++){    
        mag = mag + state_array[i];
    }   
    return fabs(mag);
}

/***** Energy calculation *****/
float state_energy(int state_array[]){
    float eng = 0.0;
    for (int i = 0; i < pow(L, 3); i++){
        int NN = state_array[neighbour_table[0][i]] + state_array[neighbour_table[1][i]] + state_array[neighbour_table[2][i]] +\
                state_array[neighbour_table[3][i]] + state_array[neighbour_table[4][i]] + state_array[neighbour_table[5][i]];
        eng += (-J * state_array[i] * NN)/2;
        /*printf("eng[%d][%d][%d] = %f \n",i,j,k,eng);*/
    } 
    return eng;
}

/********** WOLFF CLUSTER UPDATE **********/
void metropolis_cluster_update(int state_array[], float beta){
    float delta_E, w, r, spin, NN;
    int* cluster = NULL;
    int* buffer = NULL;
    size_t cluster_size = 0;
    size_t buffer_size = 0;

    // Step 1 - select random spin and add to buffer and cluster
    int spin_index = rand() % (L*L*L);
    buffer = add_spin_index(buffer, &buffer_size, spin_index);
    cluster = add_spin_index(cluster, &cluster_size, spin_index);

    // Step 2 - checking the buffer array spins too to grow cluster
    do{ 
        for (int i = 0; i < NUM_NEIGHBOURS; i++)
        {
            if (state_array[buffer[0]] == state_array[neighbour_table[i][buffer[0]]] && isInList(cluster, cluster_size, neighbour_table[i][buffer[0]]) == 0){
                r = random_number();
                w = 1 - exp(-2 * J * beta);
                if (r < w){
                    buffer = add_spin_index(buffer, &buffer_size, neighbour_table[i][buffer[0]]);
                    cluster = add_spin_index(cluster, &cluster_size, neighbour_table[i][buffer[0]]);
                }else {}  
            }else {}
        }
        buffer = remove_spin_index(buffer, &buffer_size); 
    } while (buffer_size != 0);

    // step 3 - FLIP SPINS FROM CLUSTER AND CLEAR CLUSTER ARRAY
    do{
        state_array[cluster[0]] = -1 * state_array[cluster[0]];     // Flip the spin cluster
        cluster = remove_spin_index(cluster, &cluster_size);        // clear the spin from cluster
    } while (cluster_size != 0);

}

int main(int argc, char** argv){
    time_t start_time;
    time(&start_time);

    // Temperature array initiation
    double* T = temperature_points(T_start, T_end, T_points);  

    // Initialise the MPI
    int processes, rank;
        // the "processes" entered below here should be used in the terminal as well 
        // and it should be chosen such that "count" remains an integer
    processes = 6;
    int count = T_points / (processes - 1);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Master process 
    if (rank == MASTER_RANK){ 
        // Send the temperature points to slave processes
        int remaining_points = T_points;
        for (int slave_rank = 1; slave_rank < processes; slave_rank++){
            MPI_Send(&T[T_points - remaining_points], count, MPI_DOUBLE, slave_rank, 0, MPI_COMM_WORLD);
            remaining_points = remaining_points - count;
        }
        }else { //(rank != 0) slave processes
        // Recieve the temperature points send from master
        double T_buffer[count];
        MPI_Recv(T_buffer, count, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        //Initialise buffer arrays in slave to hold results temporarily
        long double E_buffer[count];
        long double M_buffer[count];
        long double C_buffer[count];
        long double X_buffer[count];
        long double U_buffer[count];

        // Prepare the neighbour table
        Process_neighbour_table();
        
        // Begin the Metropolis algotithm for T_buffer values
        for (int i = 0; i < count; i++){
            double beta = 1.0/ T_buffer[i];       // Beta defenition

            // 1D array of present state for easy manipulation 
            int state_array[L*L*L];

            // initialize the state 
            initialize_random_state(state_array);

            // Initiate the markov chain and let it reach target distribution
            for (int i = 0; i < eq_steps; i++){
                metropolis_cluster_update(state_array, beta);
            }
            
            // Quantities calculation
            long double energy, mag, Energy_sum, Mag_sum, Energy_sq_sum, Mag_sq_sum;
            long double Mag_four_sum, Mag_2_sum, Mag_4_sum;
            Energy_sum = Energy_sq_sum = Mag_sq_sum = Mag_sum =  0.0;
            Mag_four_sum = Mag_2_sum = Mag_4_sum = 0.0;

            for (int i = 0; i < avg_steps; i++){
                metropolis_cluster_update(state_array, beta);
                energy = state_energy(state_array);
                mag = state_magnetisation(state_array);

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
        double E_array[T_points];
        double M_array[T_points];
        double C_array[T_points];
        double X_array[T_points];
        double U_array[T_points];

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
        for (int i = 0; i < T_points; i++){
            printf("%f   \t %f   \t %f   \t %f   \t %f   \t %f  \n", T[i], E_array[i], M_array[i], C_array[i], X_array[i] , U_array[i]);  
        }

        time_t end_time;
        time(&end_time);
        double timeDiff = difftime(end_time, start_time);

        /***** store data as a txt file *****/ 
        FILE* file;
        char filename[] = "/Users/satheeshdm/Desktop/3d_mpi_cluster.txt"; 

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
        fprintf(file, "Total points computed = %d\n", T_points);
        fprintf(file, "Runtime of the code  = %0.2f seconds\n", timeDiff);
        fprintf(file, "\n\nOUTPUT TABLE: \n");
        fprintf(file, "T\t\tE\t\tM\t\tC\t\tX\t\tU\n");
        for (int i = 0; i < T_points; i++){
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


