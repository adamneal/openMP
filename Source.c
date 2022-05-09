#include<stdio.h>
#include<omp.h>
#include <time.h>
int main()
{
    //variables for backwards sub
    float ratio;
    float values[1000];
    float backwardsSum = 0.0;

    //begin clock for performance analysis
    double time_spent = 0.0;
    clock_t begin = clock();

    //set multithreading
    omp_set_num_threads(4);
    
   
    
    // Matrix in form A|B
    float augmentedMatrix[3][4] = {
        {2,1,1,10},
        {3,2,3,18},
        {1,4,9,16}


    };
    int sizeOfArray = (sizeof(augmentedMatrix) / sizeof(augmentedMatrix[0])); //order of matrix (amount of rows)

    int j;
    int i;
    int k;
  
    
    //forward elimination
    for (k = 0; k <= sizeOfArray-1; k++) {
        //init parallel for loop
        #pragma omp parallel shared(k, i, sizeOfArray, ratio, augmentedMatrix) private(j)
        #pragma omp for
        for (i = 0; i <= sizeOfArray -1; i++) {
            if (i > k) { // as we are only looking for upper triangle
                    ratio = augmentedMatrix[i][k] / augmentedMatrix[k][k]; //ratio for multiplier
                        
                        //init parallel for loop
                        #pragma omp parallel shared(k,i, sizeOfArray, ratio, augmentedMatrix) private(j)
                        #pragma omp for 
                        for (j = 0; j <= sizeOfArray; j++) { //does forward elim calcs with the ratio
                          augmentedMatrix[i][j] = augmentedMatrix[i][j] - ratio * augmentedMatrix[k][j];
                        } 
            }
        }
    }


    
    //backwards sub
    values[sizeOfArray -1] = augmentedMatrix[sizeOfArray -1][sizeOfArray] / augmentedMatrix[sizeOfArray -1][sizeOfArray -1];

    for (int i = sizeOfArray - 2; i >= 0; i--) {
        backwardsSum = 0;
        for (int j = i; j <= sizeOfArray -1; j++) {
            backwardsSum = backwardsSum + augmentedMatrix[i][j] * values[j];
        }
        values[i] = (augmentedMatrix[i][sizeOfArray] - backwardsSum) / augmentedMatrix[i][i];
    }


    //prints out the unknown values that were solved
    printf("\nThe solution for unknowns are: \n");
    for (int i = 0; i <= sizeOfArray -1; i++)
    {
        printf("\n%f\t", values[i]);
    }

    //calculates time taken and prints it
    clock_t end = clock();
    time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
    printf("\n\nThe elapsed time is %f seconds", time_spent);

    return(0);
}