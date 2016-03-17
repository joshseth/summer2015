#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_blas.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <libiomp/omp.h>

/***************************************************************** List of Functions **********************************************************/


/* g := length of proteins and regulatory regions.
 * N := number of generations. 
 * pop := population size. 
 */
static const int g = 20;
static const int N = 10000;
static const int pop = 100;
static const int MUTATE = 10000;
static const int CHROMOSOMES = 5; 
/* select_mate := Smirnov Transform Algorithm to weight random number to be proportional to individual fitness score.
 * The random number generated below determines the probability an individual will be chosen to reproduce and is essential to simulate Natural Selection.
 */

double select_mate (double *fitnesses)
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
      if (add_result[i] < 0) {
       add_result[i] = 0; 
      }
      if (add_result[i] > 10) {
        add_result[i] = 10; 
      } 
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



 int point_mutate(int seq[], int mutation_rate) 
    {
      int *genetic_sequence;
      genetic_sequence = &seq[0]; 
          for (int nuc = 0; nuc < g; nuc++) 
            {
              int mutation_probability = rand()%mutation_rate;
              int mutation_target = rand()%4; 

              for (int tar = 0; mutation_probability == 1 && mutation_target == tar && tar < 4; tar++)
              {
                    genetic_sequence[nuc] = tar;
              }
            }
          return (genetic_sequence[0]); 
    }




  int mutate_regulatory_direction(int expression_direction, int direction_change_rate)
    {
      int mut_dir  = rand()%direction_change_rate;
      if (mut_dir == 1 && expression_direction == 1)
        {
          expression_direction = 0; 
        }
      else if (mut_dir == 1 && expression_direction == 0)
        {
          expression_direction = 1;
        }
       return expression_direction; 
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
  //genome *organism; 
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



/* **************************************************** BEGINING OF CODE EXECUTION ************************************************************
 *
 *
 */

  int  main ()
  {

    double meanfit;
    FILE *fit; 
    fit = fopen ("evo_fitness1e3.dat", "w");

    FILE *flocation;
    flocation = fopen("evo_location1e3.dat", "w");

    time_t time;
    srand ((unsigned) (&time));
    population *currpop, *prevpop;
    currpop = (population*) malloc(sizeof(population));

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

      currpop->organism[0].location = 0;   

      for (count = 0; count < 5; count++) {

        currpop->organism[0].allele[count].reg_direction = rand()%2;
        currpop->organism[0].allele[count].pro_direction = rand()%2;  

        int i;
        int p_rand[g];
        int r_rand[g]; 
        printf ("GENE %d\n", count); 

        for (i = 0; i < g; i++) {
        p_rand[i] = rand()%4;
        r_rand[i] = rand()%4;
        currpop->organism[0].allele[count].regulator[i] = r_rand[i];
        currpop->organism[0].allele[count].protein[i] = p_rand[i];	
        printf ("%d\t%d\n", currpop->organism[0].allele[count].regulator[i], currpop->organism[0].allele[count].protein[i] );
        }
    }
    int h;
    int j;
    int k;
    for (h = 0; h < 5; h++) {
      for (j = 0; j < 5; j++) {
        for (k = 0; k < g; k++) {
          if (currpop->organism[0].allele[h].regulator[k] == currpop->organism[0].allele[j].protein[k]) {
           currpop->organism[0].GRN[h][j] = (currpop->organism[0].GRN[h][j]) + 1;
            }
        }
        if (currpop->organism[0].allele[h].reg_direction == currpop->organism[0].allele[j].pro_direction) {
           currpop->organism[0].GRN[h][j] = (currpop->organism[0].GRN[h][j] * (-1)); 
        }
      }
    }

    for (h = 0; h < 5; h++)
    {
      for (j = 0; j < 5; j++)
      {

        if (currpop->organism[0].allele[h].reg_direction == 1)
        {
          currpop->organism[0].GRN[h][j] = (currpop->organism[0].GRN[h][j] * (-1)); 
        }
      }
    }
    
    int q, r; 
    for (q = 0; q < 5; q++) {
      for (r = 0; r < 5; r++) {
      printf ("%f ", currpop->organism[0].GRN[q][r]);
      }
    printf ("\n"); 
    }


    int i; 

    currpop->organism[0].phenotype = mat_power(currpop->organism[0].GRN, y, 5, 5, 10);

    double euc_dist = sqrt (vectors_dot_prod (currpop->organism[0].phenotype, currpop->organism[0].phenotype, 5) + vectors_dot_prod (optimal_phenotype, optimal_phenotype, 5) - (2 * vectors_dot_prod (optimal_phenotype, currpop->organism[0].phenotype, 5) ) ); 
    
    currpop->organism[0].fitness = exp ((-1) * (pow(euc_dist, 2)));

    for (int tst = 1; tst < pop; tst++)
    {
      currpop->organism[tst] = currpop->organism[0];
    }

/***************************** Simulation after initial (0th) generation *************************
 *
 */


  for ( int gen = 1; gen < N; gen++ )

  {

    prevpop = currpop;
    int org;
    meanfit = 0;
    double fitnesses[pop];

      for ( org = 0; org < pop; org++ )

      {

        fitnesses[org] = prevpop->organism[org].fitness;

      }


#pragma omp parallel for shared(currpop, prevpop) private(org, add_result, mult_result) 

      for (org = 0; org < pop; org++)

      {
        int h, j, k;
/////////////////////////////////////////////////////// Choose Two Parents to Make a new Organism ////////////////////////////////////////////////////////////////////////

        int parent1 = select_mate(fitnesses);
        int parent2 = select_mate(fitnesses);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        for (int hap = 0; hap < CHROMOSOMES; hap++)

        {
//////////////////////////////////////////////////////////////////// Recombine the Two Parental Genomes to Form a New Organism ////////////////////////////////////////////

          int rec = rand() % 2 == 0 ? parent1 : parent2;
          currpop->organism[org].allele[hap] = prevpop->organism[rec].allele[hap]; 

///////////////////////////////// Point Mutations of Both Regulatory and Protein Coding Regions for all genes  ////////////////////////////////////////////////////////////

           currpop->organism[org].allele[hap].regulator[0] = point_mutate ( currpop->organism[org].allele[hap].regulator, MUTATE );
           currpop->organism[org].allele[hap].protein[0] = point_mutate ( currpop->organism[org].allele[hap].protein, MUTATE );

///////////////////////////////// Expression Direction Mutation for Both Regulatory and Protein Interaction terms /////////////////////////////////////////////////////////

           currpop->organism[org].allele[hap].reg_direction = mutate_regulatory_direction ( currpop->organism[org].allele[hap].reg_direction, MUTATE );
           currpop->organism[org].allele[hap].pro_direction = mutate_regulatory_direction ( currpop->organism[org].allele[hap].pro_direction, MUTATE );

        }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Below for loop ensures entries of the GRN Matrix are Zeroed Out. It is necessary ////////////////////////////////////////////////////

        for ( h = 0; h < CHROMOSOMES; h++ )

        {

          for ( j = 0; j < CHROMOSOMES; j++ )

          {

            currpop->organism[org].GRN[h][j] = 0;

          }

        }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        for ( h = 0; h < CHROMOSOMES; h++ )

        {

          for ( j = 0; j < CHROMOSOMES; j++ )

          {

            for ( k = 0; k < g; k++ )

            {

              if ( currpop->organism[org].allele[h].regulator[k] == currpop->organism[org].allele[j].protein[k] )

              {

                currpop->organism[org].GRN[h][j] = ( currpop->organism[org].GRN[h][j] + 0.001 );
              }

            }

          }

        }


        for ( h = 0; h < CHROMOSOMES; h++ )

        {

          for ( j = 0; j < CHROMOSOMES; j++ )

          {

            if ( currpop->organism[org].allele[h].reg_direction != currpop->organism[org].allele[j].pro_direction )

            {

              currpop->organism[org].GRN[h][j] = ( currpop->organism[org].GRN[h][j] * ( -1 ) );

            }
          }
        }


         currpop->organism[org].phenotype = mat_power ( currpop->organism[org].GRN, y, CHROMOSOMES, CHROMOSOMES, 100 );


       double euc_dist =
         sqrt ( vectors_dot_prod ( currpop->organism[org].phenotype, currpop->organism[org].phenotype, 5 )
           + vectors_dot_prod ( optimal_phenotype, optimal_phenotype, 5 )
           - ( 2 * vectors_dot_prod ( optimal_phenotype, currpop->organism[org].phenotype, 5 ) ) );


       currpop->organism[org].fitness = exp (((-1) * (pow ((euc_dist), 2))));
       if (currpop->organism[org].fitness <= 0.0000001)
       {
         currpop->organism[org].fitness = 0.0000001; 
       }

       if ( ( gen == 1 || gen % 1000 == 0 ) && org == 1 )

       {

         printf ( "GENERATION %d fitness-%d = %f\n", gen, org, currpop->organism[org].fitness );

         int q, r;

            for ( q = 0; q < CHROMOSOMES; q++ )

            {

              for ( r = 0; r < CHROMOSOMES; r++ )

              {

                printf ( "%f ", currpop->organism[org].GRN[q][r] );

              }

              printf ( "\n" );

            }
       }
//         meanfit += currpop->organism[org].fitness*0.01;
//     if (gen % 1 == 0)
//     {
//         if (org == (pop - 1))
//         {
//           fprintf (fit, "%f\n,", meanfit);
//         }
//     }

      } /* <<-- END of SINGLE GENERATION (THE "org" FOR LOOP) */

  } /* <<-- END of ALL GENERATIONS (THE "gen" FOR LOOP) */

  free(currpop);
  fclose(fit);
  fclose(flocation);



/*  FILE *f;
 *  f = fopen ("long-run-lin1.dat", "wb");
 *  fwrite (lin.generation, sizeof(population), 10000000, f);
 *  fclose(f);
 */ //  free (lin.generation);


}
