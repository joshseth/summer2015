#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_blas.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
/***************************************************************** List of Functions **********************************************************/


/* g := length of proteins and regulatory regions.
 * N := number of generations. 
 * pop := population size. 
 */
static const int isolate = 2000000;
static const int g = 20;
static const int N = 500000;
static const int pop = 100;
static const int mutate = 1000000;
/* weight_rand := Smirnov Transform Algorithm to weight random number to be proportional to individual fitness score.
 * The random number generated below determines the probability an individual will be chosen to reproduce and is essential to simulate Natural Selection.
 */


/* double weight_rand (double *fitnesses) {

 double sum_of_weight = 0;
 double cum_fitness[pop];
 int i; 
  for (int i = 0; i < pop; i++) {
     sum_of_weight += fitnesses[i];
     cum_fitness[i] = sum_of_weight;
  }
    double rnd = ((double)rand()/(double)RAND_MAX)*sum_of_weight;
  for (i = 0; i < pop; i++) {
    if (rnd <= cum_fitness[i])
    {
      return i;
      break;
    }
      rnd -= cum_fitness[i];
  }
  return i;
}
*/

double weight_rand (double *fitnesses)
{
  double rnd = (double) rand()/(double) RAND_MAX;
  double totalfitness = 0;
  double probability[pop];
  double weights = 0;
  // weights MUST be set to 0, otherwise function cannot properly be called more than once! 
  int survivor;
  for (int i = 0; i < pop; i++)
  {
    totalfitness += fitnesses[i]; 
  }
  for (int j = 0; j < pop; j++)
  {
    probability[j] = fitnesses[j]/totalfitness; 
  }
  for (int k = 0; k < pop; k++)
  {
    weights += probability[k];
    if (rnd < weights)
    {
      survivor = k;
      break;
    }
  }
  return survivor;
}

double geography_rand (double *fitnesses, int *geog)
{
  double rnd = (double) rand()/(double) RAND_MAX;
  double totalfitness = 0;
  double probability[pop];
  double mateprob[pop];
  double weights = 0;
  int survivor;
  for (int l = 0; l < pop; l++)
  {
    mateprob[l] = fitnesses[l] * exp(-geog[l]);
  }
  for (int i = 0; i < pop; i++)
  {
    totalfitness += mateprob[i]; 
  }
  for (int j = 0; j < pop; j++)
  {
    probability[j] = mateprob[j]/totalfitness; 
  }
  for (int k = 0; k < pop; k++)
  {
    weights += probability[k];
    if (rnd < weights)
    {
      survivor = k;
      break;
    }
  }
  return survivor;
}
/*double geography_rand (double *fitnesses, int *geog) {

 double sum_of_weight = 0;
 int i; 
  for (int i = 0; i < pop; i++) {
     sum_of_weight += (fitnesses[i]*exp(-geog[i]/10));
  }
    double rnd = fmod((double)rand()/((double)RAND_MAX), sum_of_weight);
  for (i = 0; i < pop; i++) {
    if (rnd < (fitnesses[i]*exp(-geog[i/10]))){
      return i;
      break;
      }
      rnd -= (fitnesses[i]*exp(-geog[i]/10));
  }
  return i;
}*/

/* vectors_dot_prod := dot product of vectors x and y of dimension size n. 
 */

double vectors_dot_prod(double *x, double *y, int n)
{
    double res = 0.0;
    int i;
    for (i = 0; i < n; i++) {
        res += x[i] * y[i];
    
	}

    return res;
}

/* matrix_vector_mult :=
 * Arguments: 5 x 5 matrix, a vector, and row and column lengths. 
 */

double mult_result[5] = {0, 0, 0, 0, 0};
double* matrix_vector_mult(double mat[5][5], double *vec, int rows, int cols)
{ // in matrix form: result = mat * vec;
	int i;
		for (i = 0; i < rows; i++) {
        	mult_result[i] = vectors_dot_prod(mat[i], vec, cols);
		}
	return mult_result;
} 

/* vec_add := vector addition and squashing. 
 * Adds vectors, HOWEVER, this also constrains the value of any element to be between [0, 1].
 */

double add_result[5] = {0, 0, 0, 0, 0};
double* vec_add (double vec1[5], double vec2[5]) {
    for (int i = 0; i < 5; i++) {
    add_result[i] = vec1[i] + vec2[i];
   //   if (add_result[i] < 0) {
   //    add_result[i] = 0; 
   //   }
   //   if (add_result[i] > 1) {
   //     add_result[i] = 1; 
   //   } 
    }
  return add_result;
}

/* mat_power := takes a 5x5 matrix, vector, and row/col lentgh arguments. 
 * Enables iterative matrix multiplication. 
 */

double* mat_power(double mat[5][5], double *vec, int rows, int cols, int power) {
  if (power == 0)
  {
    return vec;
  } 
  else
  {
    vec =  vec_add(vec, matrix_vector_mult(mat, vec, rows, cols)); 
    double *solution2 = mat_power (mat, vec, rows, cols, (power - 1));
   /* for (int i = 0; i < 5; i++) {
      solution2[i] = solution2[i] + solution[i];
        if (solution2[i] < 0) {
          solution2[i] = 0; 
        }
        if (solution2[i] > 10) {
          solution2[i] = 10; 
        } 
    }*/
    return vec; 
  }

}

////////

////////

/*double* mat_power2(double mat[5][5], double *vec, int rows, int cols, int power)
{

}
*/
/* gene := a "gene" struct contains 2 arrays denotes "regulator" and "protein",
 * and integers for protein and regulator directions (the influence these exert when interacting in a network.
 *
 * Note: only one direction should be used.
 * Note2: remove one of the direction from code.
 */

/*********************************************************** List of Structs ******************************************************************/

typedef struct {
//	int gene_size;
	int regulator[g];
	int protein[g];
	int reg_direction;
  int pro_direction;

	} gene;

/* genome := a "genome" struct contains 5 "gene" structs (defined above),
 * a 5x5 GRN (Genetic Regulatory Network),
 * a "phenotype", and a "fitness" score.
 *
 * allel[5] is an array containing 5 "gene" structs.
 * Note: in the future, to create organisms with larger genomes the "5s" should be converted to a variable.
 */

typedef struct {

      gene allele[5];
      double GRN[5][5];
      double *phenotype;
      double fitness;
      int location;

	} genome;

/* population := a "population" struct is an array of 100 organism structs.
 *
 * Note: to simulate populations larger (or smaller) than 100, the "100" should be changed to a variable in the code.
 */

typedef struct {

  genome organism[pop];

} population;


/* lineage := a "lineage" struct is an array of multiple populations.
 * Every single generation (represented as a population struct) in evolutionary time is an element in the lineage array.
 *
 * Originally the struct was defined as
 *                                     {
 *                                      population generation[1000];
 *                                     }
 *                                     however, stack overflow errors prevented assigning lineage arrays that were very large,
 *                                     so it is now cast as a pointer and the size of the array is determined in int main ().
 */

typedef struct {

// population generation[1000];
 population *generation;
} lineage;


/* **************************************************** BEGINING OF CODE EXECUTION ************************************************************
 *
 * int main().
 *
 */

int  main () {

  double meanfit;
  
  double meanfit1;
  double meanfit2;

  FILE *fit; 
  fit = fopen ("fitness.dat", "w");

  FILE *flocation;
  flocation = fopen("location.dat", "w");

  FILE *hybridfit;
  hybridfit = fopen ("Hybrid_Fitnesses.txt", "w");

  FILE *allo1fit;
  allo1fit = fopen("Allo1_Fitness.txt", "w");

  FILE *allo2fit;
  allo2fit = fopen("Allo2_Fitness.txt", "w");

  time_t time;
  srand ((unsigned) (&time));
  /* "time_t" above ensures better random number generation since it uses time.
   */

  //lineage lin;
  //lin.generation = (population*) malloc (N * sizeof(population));
  /* Assigning a lineage struct to to "lin"
   * lin.generation is an array of 1,000,000 populations. 
   */
  population currpop;
  population prevpop;

  population allo1curr;
  population allo1prev;
  population allo2curr;
  population allo2prev;

  int count;

   double y[5]; // = {1, 1, 1, 1, 1};  
  double optimal_phenotype[5]; // = {2, 2, 2, 2, 2}; 
  for (int ph = 0; ph < 5; ph++)
  {
    y[ph] = ((double) rand()/(double) RAND_MAX)*(rand()%10);
    optimal_phenotype[ph] = ((double) rand()/(double) RAND_MAX)*(rand()%10);
  }
  /* y[5] is a 5 dimensional initial vector.
   * Interpreted as the "environment."
   */

  /* the below for loop generates ONLY the first generation (currpop) of the simulation. 
   * All individuals in this generation are identical and random.
   */

    currpop.organism[0].location = 0;   

    for (count = 0; count < 5; count++) {

      currpop.organism[0].allele[count].reg_direction = rand()%2;
      currpop.organism[0].allele[count].pro_direction = rand()%2;  

      int i;
      int p_rand[g];
      int r_rand[g]; 
      printf ("GENE %d\n", count); 

      for (i = 0; i < g; i++) {
      p_rand[i] = rand()%4;
      r_rand[i] = rand()%4;
      currpop.organism[0].allele[count].regulator[i] = r_rand[i];
      currpop.organism[0].allele[count].protein[i] = p_rand[i];	
      printf ("%d\t%d\n", currpop.organism[0].allele[count].regulator[i], currpop.organism[0].allele[count].protein[i] );
      }
  }
  int h;
  int j;
  int k;
  for (h = 0; h < 5; h++) {
    for (j = 0; j < 5; j++) {
      for (k = 0; k < g; k++) {
        if (currpop.organism[0].allele[h].regulator[k] == currpop.organism[0].allele[j].protein[k]) {
         currpop.organism[0].GRN[h][j] = (currpop.organism[0].GRN[h][j]) + 1;
        /*  if (currpop.organism[org].allele[h].reg_direction == 1) {
            currpop.organism[org].GRN[h][j] = (currpop.organism[org].GRN[h][j] * (-1));
            }
         */ }
      }
      if (currpop.organism[0].allele[h].reg_direction == currpop.organism[0].allele[j].pro_direction) {
         currpop.organism[0].GRN[h][j] = (currpop.organism[0].GRN[h][j] * (-1)); 
      }
    }
	}

  for (h = 0; h < 5; h++)
  {
    for (j = 0; j < 5; j++)
    {
    //  currpop.organism[org].GRN[h][j] = pow(0.5 + (exp(6 - currpop.organism[org].GRN[h][j]))/5, -1);

      if (currpop.organism[0].allele[h].reg_direction == 1)
      {
        currpop.organism[0].GRN[h][j] = (currpop.organism[0].GRN[h][j] * (-1)); 
      }
    }
  }
	
  int q, r; 
  for (q = 0; q < 5; q++) {
    for (r = 0; r < 5; r++) {
    printf ("%f ", currpop.organism[0].GRN[q][r]);
    }
  printf ("\n"); 
  }


 // double *Mat[] = {currpop.organism[org].GRN[0], currpop.organism[org].GRN[1], currpop.organism[org].GRN[2], currpop.organism[org].GRN[3], currpop.organism[org].GRN[4]};


  int i; 
/* double *ans =  matrix_vector_mult (Mat, y, 5, 5);
  printf ("\n"); 
  for (i = 0; i < 5; i++) {
    printf ("%f\n", ans[i]); 
  }
*/

  printf("\nMatrix iteration:\n");
  currpop.organism[0].phenotype = mat_power(currpop.organism[0].GRN, y, 5, 5, 10);
  for (i = 0; i < 5; i++) {
    printf ("%f\n", currpop.organism[0].phenotype[i]);
  }

  double euc_dist = sqrt (vectors_dot_prod (currpop.organism[0].phenotype, currpop.organism[0].phenotype, 5) + vectors_dot_prod (optimal_phenotype, optimal_phenotype, 5) - (2 * vectors_dot_prod (optimal_phenotype, currpop.organism[0].phenotype, 5) ) ); 
  printf ("norm = %f\n", euc_dist); 
  
  currpop.organism[0].fitness = exp ((-1) * (pow(euc_dist, 2)));

  printf ("fitness-%d = %f\n", 0, currpop.organism[0].fitness);    

 // free(ans); 
  

for (int tst = 1; tst < pop; tst++) {
currpop.organism[tst] = currpop.organism[0];
}
for (count = 0; count < 5; count ++){
  printf( "copy gene %d\n", count); 
  for (int i = 0; i < g; i++) {
    printf ("%d\t%d\n", currpop.organism[1].allele[count].regulator[i], currpop.organism[0].allele[count].protein[i] );
  }
}

/***************************** Simulation after initial (0th) generation *************************
 *
 * Generations gen = 1 through gen = 100,000.
 */

  int gen;
  int hap; 

  for (gen = 1; gen < N; gen ++) {
    prevpop = currpop; 
    int org;
    meanfit = 0;
    double fitnesses[pop];
     // if (gen % 20000 == 0)
     // {
     //   for (int ph = 0; ph < 5; ph++)
     //   {
     //     optimal_phenotype[ph] = ((double) rand()/(double) RAND_MAX) * (rand()%10);
     //   }
     // }
      for (org = 0; org < pop; org++) {
        fitnesses[org] = prevpop.organism[org].fitness;
       // printf ("org %d fitness - %f\n ", org, fitnesses[org]);
      }
      for (org = 0; org < pop; org++) {
        
        int parent1 = weight_rand(fitnesses);
        //printf("gen-%d  org-%d parent-%d\n", gen, org, parent1);
      //  int geog[pop];
      //    for (int loc = 0; loc < pop; loc++)
      //    {
      //       geog[loc] = abs(prevpop.organism[parent1].location - prevpop.organism[loc].location);
      //    }
        
       //  int parent2 = geography_rand(fitnesses, geog);  
        int parent2 = weight_rand(fitnesses);
  //      printf ("organism %d has parent 1 = %d and parent 2 = %d\n" , org, parent1, parent2);     
        for (hap = 0; hap < 5; hap++){
          int rec = rand()%2;
   //       if (org < (pop/2))
   //       {
   //         currpop.organism[org].allele[hap] = prevpop.organism[parent1].allele[hap]; 
   //       }
          if (/* org > (pop/2) && */rec == 0) {
            currpop.organism[org].allele[hap] = prevpop.organism[parent1].allele[hap];
            currpop.organism[org].allele[hap].reg_direction = prevpop.organism[parent1].allele[hap].reg_direction;
            currpop.organism[org].allele[hap].pro_direction = prevpop.organism[parent1].allele[hap].pro_direction; 
        }
          if (/* org > (pop/2) &&  */rec == 1) {
            currpop.organism[org].allele[hap] = prevpop.organism[parent2].allele[hap]; 
            currpop.organism[org].allele[hap].reg_direction = prevpop.organism[parent2].allele[hap].reg_direction;
            currpop.organism[org].allele[hap].pro_direction = prevpop.organism[parent2].allele[hap].pro_direction; 
          }
          
          for (int nuc = 0; nuc < g; nuc++) {
            int mut_reg1 = rand()%mutate;
            int mut_reg2 = rand()%4; 
            int mut_pro1 = rand()%mutate;
            int mut_pro2 = rand()%4; 
            
            if (mut_reg1 == 1 && mut_reg2 == 0) {  
              currpop.organism[org].allele[hap].regulator[nuc] = 0; 
            } 
            if (mut_reg1 == 1 && mut_reg2 == 1) {
              currpop.organism[org].allele[hap].regulator[nuc] = 1; 
            }
            if (mut_reg1 == 1 && mut_reg2 == 2) {
              currpop.organism[org].allele[hap].regulator[nuc] = 2; 
            }
            if (mut_reg1 == 1 && mut_reg2 == 3) {
              currpop.organism[org].allele[hap].regulator[nuc] = 3;
            }

            if (mut_pro1 == 1 && mut_pro2 == 0) {
              currpop.organism[org].allele[hap].protein[nuc] = 0; 
            }
            if (mut_pro1 == 1 && mut_pro2 == 1) {
              currpop.organism[org].allele[hap].protein[nuc] = 1; 
            }
            if (mut_pro1 == 1 && mut_pro2 == 2) {
              currpop.organism[org].allele[hap].protein[nuc] = 2;
            }
            if (mut_pro1 == 1 && mut_pro2 == 3) {
              currpop.organism[org].allele[hap].protein[nuc] = 3;
            }

//      printf ("%d\t%d\n", currpop.organism[0].allele[hap].protein[nuc], currpop.organism[org].allele[hap].protein[nuc] );
          }
//      printf ("\n"); 
            int mut_dir  = rand()%mutate;
            int mut_dir2 = rand()%mutate;  
          if (mut_dir == 1 && currpop.organism[org].allele[hap].reg_direction == 1) {
            currpop.organism[org].allele[hap].reg_direction = 0; 
          }
          else if (mut_dir == 1 && currpop.organism[org].allele[hap].reg_direction == 0) {
            currpop.organism[org].allele[hap].reg_direction = 1;
          }
          if (mut_dir2 == 1 && currpop.organism[org].allele[hap].pro_direction  == 1) {
            currpop.organism[org].allele[hap].pro_direction = 0;
          }
          else if (mut_dir2 == 1 && currpop.organism[org].allele[hap].pro_direction  == 0) {
            currpop.organism[org].allele[hap].pro_direction = 1;
          }
   

        }
         // int mig_rand = rand()%3;
         // int mig_dist;

         // if (mig_rand == 0)
         // {
         //   mig_dist = 0;
         // }
         // if (mig_rand == 1)
         // {
         //   mig_dist = 1;
         // }
         // if (mig_rand == 2)
         // {
         //   mig_dist = -1;
         // }
        //    currpop.organism[org].location = prevpop.organism[parent1].location + mig_dist;
        //  if (gen % 100 == 0)
        //  {
        //    fprintf (flocation, "%d,", currpop.organism[org].location);
        //    if (org == (pop-1))
        //    {
        //      fprintf(flocation, "\n"); 
        //    }
        //  }
          

        int h, j, k;

        for (h = 0; h < 5; h++)
        {
          for (j = 0; j < 5; j++)
          {
            currpop.organism[org].GRN[h][j] = 0;
          }
        }

        for (h = 0; h < 5; h++)
        {
          for (j = 0; j < 5; j++)
          {
            for (k = 0; k < g; k++)
            {
              if (currpop.organism[org].allele[h].regulator[k] == currpop.organism[org].allele[j].protein[k])
              {
                currpop.organism[org].GRN[h][j] = (currpop.organism[org].GRN[h][j] + 0.001);
              }
            }

            // if (currpop.organism[org].allele[h].reg_direction == 0)
            // {
            //   currpop.organism[org].GRN[h][j] = (currpop.organism[org].GRN[h][j] * (-1));
            // }
          }
         }


  for (h = 0; h < 5; h++)
  {
    for (j = 0; j < 5; j++)
    {
        // currpop.organism[org].GRN[h][j] = pow(0.5 + (exp(6 - currpop.organism[org].GRN[h][j]))/5, -1);
      if (currpop.organism[org].allele[h].reg_direction != currpop.organism[org].allele[j].pro_direction)
      {
        currpop.organism[org].GRN[h][j] = (currpop.organism[org].GRN[h][j] * (-1));
      }
    }
  }
  

  /* don't need pmat now that matrix_vec_mult has mat[5][5] instead **mat now */
  //       double *pmat[] = {currpop.organism[org].GRN[0], currpop.organism[org].GRN[1], currpop.organism[org].GRN[2], currpop.organism[org].GRN[3], currpop.organism[org].GRN[4]};
         currpop.organism[org].phenotype = mat_power (currpop.organism[org].GRN, y, 5, 5, 100); 

       double euc_dist = sqrt (vectors_dot_prod (currpop.organism[org].phenotype, currpop.organism[org].phenotype, 5) + vectors_dot_prod (optimal_phenotype, optimal_phenotype, 5) - (2 * vectors_dot_prod (optimal_phenotype, currpop.organism[org].phenotype, 5) ) );
//  printf ("norm = %f\n", euc_dist);

       currpop.organism[org].fitness = exp (((-1) * (pow ((euc_dist), 2))));
       if (currpop.organism[org].fitness <= 0.0000001)
       {
         currpop.organism[org].fitness = 0.0000001; 
       }

    /*   if (currpop.organism[org].fitness <= 0.000001)
       {
         currpop.organism[org].fitness = 0.000001;
       }*/
       if ( (gen == 1 || gen % 1000 == 0) && org == 1 ){
       printf ("GENERATION %d fitness-%d = %f\n", gen, org, currpop.organism[org].fitness);  
         int q, r;
          for (q = 0; q < 5; q++) {
            for (r = 0; r < 5; r++) {
              printf ("%f ", currpop.organism[org].GRN[q][r]);
            }
            printf ("\n");
          }
    /*  for (hap = 0; hap < 5; hap++) {
        printf ("Gene %d\n", hap); 
        for (int i = 0; i < g; i++) {

          printf ("%d\t%d\n", currpop.organism[org].allele[hap].regulator[i], currpop.organism[org].allele[hap].protein[i] );

          }
       } */
     }
         meanfit += currpop.organism[org].fitness*0.01;
     if (gen % 1 == 0)
     {
         if (org == (pop - 1))
         {
           fprintf (fit, "%f\n,", meanfit);
         }
     }

      }
  }

  lineage allo1;
  lineage allo2;
    allo1.generation = malloc (isolate * sizeof(population)); 
    allo2.generation = malloc (isolate * sizeof(population)); 
  allo1curr = currpop;
  allo2curr = currpop;
  
  
  for (int gee = 0; gee < isolate; gee ++) {
    meanfit1 = 0;
    allo1prev = allo1curr; 
    int org;
    double fitnesses[pop];
      for (int org = 0; org < pop; org++) {
        fitnesses[org] = allo1prev.organism[org].fitness;
       // printf ("org %d fitness - %f\n ", org, fitnesses[org]);
      }
      for (org = 0; org < pop; org++) {
        
        int parent1 = weight_rand(fitnesses);
        
        int parent2 = weight_rand(fitnesses);
        for (hap = 0; hap < 5; hap++){
          int rec = rand()%2;
          if (/* org > (pop/2) && */rec == 0) {
            allo1curr.organism[org].allele[hap] = allo1prev.organism[parent1].allele[hap];
        }
          if (/* org > (pop/2) &&  */rec == 1) {
            allo1curr.organism[org].allele[hap] = allo1prev.organism[parent2].allele[hap]; 
          }
          
          for (int nuc = 0; nuc < g; nuc++) {
            int mut_reg1 = rand()%mutate;
            int mut_reg2 = rand()%4; 
            int mut_pro1 = rand()%mutate;
            int mut_pro2 = rand()%4; 
            
            if (mut_reg1 == 1 && mut_reg2 == 0) {  
              allo1curr.organism[org].allele[hap].regulator[nuc] = 0; 
            } 
            if (mut_reg1 == 1 && mut_reg2 == 1) {
              allo1curr.organism[org].allele[hap].regulator[nuc] = 1; 
            }
            if (mut_reg1 == 1 && mut_reg2 == 2) {
              allo1curr.organism[org].allele[hap].regulator[nuc] = 2; 
            }
            if (mut_reg1 == 1 && mut_reg2 == 3) {
              allo1curr.organism[org].allele[hap].regulator[nuc] = 3;
            }

            if (mut_pro1 == 1 && mut_pro2 == 0) {
              allo1curr.organism[org].allele[hap].protein[nuc] = 0; 
            }
            if (mut_pro1 == 1 && mut_pro2 == 1) {
              allo1curr.organism[org].allele[hap].protein[nuc] = 1; 
            }
            if (mut_pro1 == 1 && mut_pro2 == 2) {
              allo1curr.organism[org].allele[hap].protein[nuc] = 2;
            }
            if (mut_pro1 == 1 && mut_pro2 == 3) {
              allo1curr.organism[org].allele[hap].protein[nuc] = 3;
            }

//      printf ("%d\t%d\n", allo1curr.organism[0].allele[hap].protein[nuc], allo1curr.organism[org].allele[hap].protein[nuc] );
          }
//      printf ("\n"); 
            int mut_dir  = rand()%mutate;
            int mut_dir2 = rand()%mutate;  
          if (mut_dir == 1 && allo1curr.organism[org].allele[hap].reg_direction == 1) {
            allo1curr.organism[org].allele[hap].reg_direction = 0; 
          }
          else if (mut_dir == 1 && allo1curr.organism[org].allele[hap].reg_direction == 0) {
            allo1curr.organism[org].allele[hap].reg_direction = 1;
          }
          if (mut_dir2 == 1 && allo1curr.organism[org].allele[hap].pro_direction  == 1) {
            allo1curr.organism[org].allele[hap].pro_direction = 0;
          }
          else if (mut_dir2 == 1 && allo1curr.organism[org].allele[hap].pro_direction  == 0) {
            allo1curr.organism[org].allele[hap].pro_direction = 1;
          }
   

        }

        int h, j, k;

        for (h = 0; h < 5; h++)
        {
          for (j = 0; j < 5; j++)
          {
            allo1curr.organism[org].GRN[h][j] = 0;
          }
        }

        for (h = 0; h < 5; h++)
        {
          for (j = 0; j < 5; j++)
          {
            for (k = 0; k < g; k++)
            {
              if (allo1curr.organism[org].allele[h].regulator[k] == allo1curr.organism[org].allele[j].protein[k])
              {
                allo1curr.organism[org].GRN[h][j] = (allo1curr.organism[org].GRN[h][j] + 0.001);
              }
            }

            // if (allo1curr.organism[org].allele[h].reg_direction == 0)
            // {
            //   allo1curr.organism[org].GRN[h][j] = (allo1curr.organism[org].GRN[h][j] * (-1));
            // }
          }
         }


  for (h = 0; h < 5; h++)
  {
    for (j = 0; j < 5; j++)
    {
        // allo1curr.organism[org].GRN[h][j] = pow(0.5 + (exp(6 - allo1curr.organism[org].GRN[h][j]))/5, -1);
      if (allo1curr.organism[org].allele[h].reg_direction != allo1curr.organism[org].allele[j].pro_direction)
      {
        allo1curr.organism[org].GRN[h][j] = (allo1curr.organism[org].GRN[h][j] * (-1));
      }
    }
  }
  

  /* don't need pmat now that matrix_vec_mult has mat[5][5] instead **mat now */
  //       double *pmat[] = {allo1curr.organism[org].GRN[0], allo1curr.organism[org].GRN[1], allo1curr.organism[org].GRN[2], allo1curr.organism[org].GRN[3], allo1curr.organism[org].GRN[4]};
         allo1curr.organism[org].phenotype = mat_power (allo1curr.organism[org].GRN, y, 5, 5, 100); 

       double euc_dist = sqrt (vectors_dot_prod (allo1curr.organism[org].phenotype, allo1curr.organism[org].phenotype, 5) + vectors_dot_prod (optimal_phenotype, optimal_phenotype, 5) - (2 * vectors_dot_prod (optimal_phenotype, allo1curr.organism[org].phenotype, 5) ) );

       allo1curr.organism[org].fitness = exp (((-1) * (pow ((euc_dist), 2))));
       if (allo1curr.organism[org].fitness <= 0.0000001)
       {
         allo1curr.organism[org].fitness = 0.0000001; 
       }

         meanfit1 += allo1curr.organism[org].fitness*0.01;
     if (gee % 1 == 0)
     {
         if (org == (pop - 1))
         {
           fprintf (allo1fit, "%f\n,", meanfit1);
         }
     }
       if ( (gee == 1 || gee % 1000 == 0) && org == 1 ){
       printf ("geeERATION %d fitness-%d = %f\n", gee, org, allo1curr.organism[org].fitness);  
         int q, r;
          for (q = 0; q < 5; q++) {
            for (r = 0; r < 5; r++) {
              printf ("%f ", allo1curr.organism[org].GRN[q][r]);
            }
            printf ("\n");
          }
     }

      }
        allo1.generation[gee].organism[0] = allo1curr.organism[0];
  }
  
  for (int eeg = 0; eeg < isolate; eeg ++) {
    meanfit2=0;
    allo2prev = allo2curr;
    int org;
    meanfit = 0;
    double fitnesses[pop];
      for (org = 0; org < pop; org++) {
        fitnesses[org] = allo2prev.organism[org].fitness;
      }
      for (org = 0; org < pop; org++) {
        
        int parent1 = weight_rand(fitnesses);
        int parent2 = weight_rand(fitnesses);
        for (hap = 0; hap < 5; hap++){
          int rec = rand()%2;
          if (rec == 0) {
            allo2curr.organism[org].allele[hap] = allo2prev.organism[parent1].allele[hap];
            allo2curr.organism[org].allele[hap].reg_direction = allo2prev.organism[parent1].allele[hap].reg_direction;
            allo2curr.organism[org].allele[hap].pro_direction = allo2prev.organism[parent1].allele[hap].pro_direction; 
        }
          if (rec == 1) {
            allo2curr.organism[org].allele[hap] = allo2prev.organism[parent2].allele[hap]; 
            allo2curr.organism[org].allele[hap].reg_direction = allo2prev.organism[parent2].allele[hap].reg_direction;
            allo2curr.organism[org].allele[hap].pro_direction = allo2prev.organism[parent2].allele[hap].pro_direction; 
          }
          
          for (int nuc = 0; nuc < g; nuc++) {
            int mut_reg1 = rand()%mutate;
            int mut_reg2 = rand()%4; 
            int mut_pro1 = rand()%mutate;
            int mut_pro2 = rand()%4; 
            
            if (mut_reg1 == 1 && mut_reg2 == 0) {  
              allo2curr.organism[org].allele[hap].regulator[nuc] = 0; 
            } 
            if (mut_reg1 == 1 && mut_reg2 == 1) {
              allo2curr.organism[org].allele[hap].regulator[nuc] = 1; 
            }
            if (mut_reg1 == 1 && mut_reg2 == 2) {
              allo2curr.organism[org].allele[hap].regulator[nuc] = 2; 
            }
            if (mut_reg1 == 1 && mut_reg2 == 3) {
              allo2curr.organism[org].allele[hap].regulator[nuc] = 3;
            }

            if (mut_pro1 == 1 && mut_pro2 == 0) {
              allo2curr.organism[org].allele[hap].protein[nuc] = 0; 
            }
            if (mut_pro1 == 1 && mut_pro2 == 1) {
              allo2curr.organism[org].allele[hap].protein[nuc] = 1; 
            }
            if (mut_pro1 == 1 && mut_pro2 == 2) {
              allo2curr.organism[org].allele[hap].protein[nuc] = 2;
            }
            if (mut_pro1 == 1 && mut_pro2 == 3) {
              allo2curr.organism[org].allele[hap].protein[nuc] = 3;
            }

          }
            int mut_dir  = rand()%mutate;
            int mut_dir2 = rand()%mutate;  
          if (mut_dir == 1 && allo2curr.organism[org].allele[hap].reg_direction == 1) {
            allo2curr.organism[org].allele[hap].reg_direction = 0; 
          }
          else if (mut_dir == 1 && allo2curr.organism[org].allele[hap].reg_direction == 0) {
            allo2curr.organism[org].allele[hap].reg_direction = 1;
          }
          if (mut_dir2 == 1 && allo2curr.organism[org].allele[hap].pro_direction  == 1) {
            allo2curr.organism[org].allele[hap].pro_direction = 0;
          }
          else if (mut_dir2 == 1 && allo2curr.organism[org].allele[hap].pro_direction  == 0) {
            allo2curr.organism[org].allele[hap].pro_direction = 1;
          }
   

        }

        int h, j, k;

        for (h = 0; h < 5; h++)
        {
          for (j = 0; j < 5; j++)
          {
            allo2curr.organism[org].GRN[h][j] = 0;
          }
        }

        for (h = 0; h < 5; h++)
        {
          for (j = 0; j < 5; j++)
          {
            for (k = 0; k < g; k++)
            {
              if (allo2curr.organism[org].allele[h].regulator[k] == allo2curr.organism[org].allele[j].protein[k])
              {
                allo2curr.organism[org].GRN[h][j] = (allo2curr.organism[org].GRN[h][j] + 0.001);
              }
            }

          }
         }


  for (h = 0; h < 5; h++)
  {
    for (j = 0; j < 5; j++)
    {
      if (allo2curr.organism[org].allele[h].reg_direction != allo2curr.organism[org].allele[j].pro_direction)
      {
        allo2curr.organism[org].GRN[h][j] = (allo2curr.organism[org].GRN[h][j] * (-1));
      }
    }
  }
  

         allo2curr.organism[org].phenotype = mat_power (allo2curr.organism[org].GRN, y, 5, 5, 100); 

       double euc_dist = sqrt (vectors_dot_prod (allo2curr.organism[org].phenotype, allo2curr.organism[org].phenotype, 5) + vectors_dot_prod (optimal_phenotype, optimal_phenotype, 5) - (2 * vectors_dot_prod (optimal_phenotype, allo2curr.organism[org].phenotype, 5) ) );

       allo2curr.organism[org].fitness = exp (((-1) * (pow ((euc_dist), 2))));
       if (allo2curr.organism[org].fitness <= 0.0000001)
       {
         allo2curr.organism[org].fitness = 0.0000001; 
       }

         meanfit2 += allo2curr.organism[org].fitness*0.01;
     if (eeg % 1 == 0)
     {
         if (org == (pop - 1))
         {
           fprintf (allo2fit, "%f\n,", meanfit2);
         }
     }
       if ( (eeg == 1 || eeg % 1000 == 0) && org == 1 ){
       printf ("eegERATION %d fitness-%d = %f\n", eeg, org, allo2curr.organism[org].fitness);  
         int q, r;
          for (q = 0; q < 5; q++) {
            for (r = 0; r < 5; r++) {
              printf ("%f ", allo2curr.organism[org].GRN[q][r]);
            }
            printf ("\n");
          }
     }

      }
        allo2.generation[eeg].organism[0] = allo2curr.organism[0];
  }
  
  lineage hybrid;
    hybrid.generation = malloc (isolate * sizeof(population)); 

  for (int hy = 0; hy < isolate; hy++)
  {
    for (int a = 0; a < 5; a++)
    {
      int mate = rand()%2;
      if (mate == 0)
      {
        hybrid.generation[hy].organism[0].allele[a] = allo1.generation[hy].organism[0].allele[a];
      }
      if (mate == 1)
      {
        hybrid.generation[hy].organism[0].allele[a] = allo2.generation[hy].organism[0].allele[a];
      }
    }
        int h, j, k;
//
//        for (h = 0; h < 5; h++)
//        {
//          for (j = 0; j < 5; j++)
//          {
//            hybrid.generation[hy].organism[0].GRN[h][j] = 0;
//          }
//        }

        for (h = 0; h < 5; h++)
        {
          for (j = 0; j < 5; j++)
          {
            for (k = 0; k < g; k++)
            {
              if (hybrid.generation[hy].organism[0].allele[h].regulator[k] == hybrid.generation[hy].organism[0].allele[j].protein[k])
              {
                hybrid.generation[hy].organism[0].GRN[h][j] = (hybrid.generation[hy].organism[0].GRN[h][j] + 0.001);
              }
            }
          }
        }
  for (h = 0; h < 5; h++)
  {
    for (j = 0; j < 5; j++)
    {
      if (hybrid.generation[hy].organism[0].allele[h].reg_direction != hybrid.generation[hy].organism[0].allele[j].pro_direction)
      {
        hybrid.generation[hy].organism[0].GRN[h][j] = (hybrid.generation[hy].organism[0].GRN[h][j] * (-1));
      }
    }
  }
  

         hybrid.generation[hy].organism[0].phenotype = mat_power (hybrid.generation[hy].organism[0].GRN, y, 5, 5, 100); 

       double euc_dist = sqrt (vectors_dot_prod (hybrid.generation[hy].organism[0].phenotype, hybrid.generation[hy].organism[0].phenotype, 5) + vectors_dot_prod (optimal_phenotype, optimal_phenotype, 5) - (2 * vectors_dot_prod (optimal_phenotype, hybrid.generation[hy].organism[0].phenotype, 5) ) );

       hybrid.generation[hy].organism[0].fitness = exp (((-1) * (pow ((euc_dist), 2))));
       fprintf (hybridfit, "%f\n", hybrid.generation[hy].organism[0].fitness);
  
        }
  

  fclose(fit);
  fclose(flocation);
  fclose(hybridfit);
  fclose(allo1fit);
  fclose(allo2fit);

  free(allo1.generation);
  free(allo2.generation);
  free(hybrid.generation);
/*  FILE *f; 
 *  f = fopen ("long-run-lin1.dat", "wb"); 
 *  fwrite (lin.generation, sizeof(population), 10000000, f); 
 *  fclose(f);  
 */ //  free (lin.generation); 
}  



	
