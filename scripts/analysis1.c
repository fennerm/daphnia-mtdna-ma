#include "/home/eang33/libhdr"
#include "/home/eang33/nrdir/nrhdr.h"
#include "/home/eang33/nrdir/nrutil.h"
#include "/home/eang33/tmatrix/bctmat01_haploid.hd"

/* Compile with:
gcc -O2 -o analysis1 analysis1.c ~/genlib.o ~/tmatrix/bctmat01_haploid.c ~/nrdir/*.o -lm -lgsl -lgslcblas
*/

#define max_segregating 100000
#define undefined -99999999999.99
#define max_pop_size         1000
#define maxt 1000
#define max_amoeba_evaluations 1000
#define max_n_evals 1000
#define max_data_file 3
#define max_results_vec 7
#define madrid_generations 262
#define florida_generations 187
#define c_elegans_generations 363
#define max_n_profile_ve 1000000
#define florida 0
#define madrid 1
#define florida_madrid 2
#define c_elegans 3

int n[max_data_file], t[max_data_file], eval_num, n_data_files, silent_mode, dataset;
int segregating[max_data_file], nonsegregating[max_data_file],
   n_n_evals1, n_n_evals2, same_u, fixed_mut_rate_flag, n_profile_ve,
   analyse_real_data;
double fixed_mut_rate;
double segregating_freqs[maxn + 1][max_data_file], f0[2][max_data_file],
   error_var[2][max_data_file], fac, criterion, f0_start[max_data_file],
   error_var_start[max_data_file];
double expected_freq_vec[maxn + 1][max_data_file],
   expected_freq_vec_unscaled[maxn + 1][max_data_file],
   results_vec[max_n_evals][max_n_evals][max_results_vec],
   profile_ve[2][max_n_profile_ve];
FILE *outfile,  *fopen();


read_profile_ve_file()
{
   FILE *profile_ve_dat, *fopen();
   int stat;
   profile_ve_dat = openforread("profile_ve.dat");
   n_profile_ve = 0;
   for (;;)
   {
      if (n_profile_ve >= max_n_profile_ve) 
         gabort("read_profile_ve_file: n_profile_ve >= max_n_profile_ve", 0);
      stat = fscanf(profile_ve_dat, "%lf %lf\n", &(profile_ve[0][n_profile_ve]), 
         &(profile_ve[1][n_profile_ve]));
      if (stat == EOF) break;
      if (stat != 2) gabort("read_profile_ve_file: read error 1", 1);
      n_profile_ve++;
   }
   printf("read_profile_ve_file: n_profile_ve %d\n",n_profile_ve);
//   monitorinput();
   fclose(profile_ve_dat);
}

read_data_file(char *data_file, int data_file_ind, int *nonsegregating,
   int *segregating, double segregating_freqs[maxn + 1][max_data_file], int start)
{
   FILE *infile, *fopen();
   int stat, idum, i, temp_int;
   struct acc freq_acc;
   char line_name[10];
   initacc(&freq_acc);
   if (data_file_ind >=max_data_file)
      gabort("data_file_ind >=max_data_file", data_file_ind);
   infile = openforread(data_file);
   stat = fscanf(infile, "%d %d", &idum, &(nonsegregating[data_file_ind]));
   if (stat!=2) gabort("read_data_file: read error 1", 1);
   if (idum != 0) gabort("read_data_file: read error 2", 2);
   stat = fscanf(infile, "%d %d", &idum, &temp_int);
   if (stat!=2) gabort("read_data_file: read error 1", 1);
   if (idum != 1) gabort("read_data_file: read error 3", 3);
   segregating[data_file_ind] += temp_int;
   if (segregating[data_file_ind] > max_segregating)
      gabort("read_data_file: segregating > max_segregating", 0);
   for (i=start; i<=segregating[data_file_ind]; i++)
   {
      stat = fscanf(infile, "%s %lf", line_name, &(segregating_freqs[i][data_file_ind]));
//      printf("i %d segregating_freqs[%d][%d] %lf\n", i, i, data_file_ind,
//           segregating_freqs[i][data_file_ind]);
      if (stat!=2) gabort("read_data_file: read error 4", 4);
      printf("segregating_freqs[%d][%d] %lf\n", i, data_file_ind,
           segregating_freqs[i][data_file_ind]);
      if (segregating_freqs[i][data_file_ind] != 1)
      {
         accum(&freq_acc, segregating_freqs[i][data_file_ind]);
      }
   }
   fclose(infile);
   printf("nonsegregating[%d] %d segregating[%d] %d\n", data_file_ind, nonsegregating[data_file_ind], data_file_ind, segregating[data_file_ind]);
   printmse("Mean frequency of segregating sites ", &freq_acc);
//   monitorinput();
}


generate_expected_frequency(double expected_freq_vec[maxn + 1][max_data_file], int file)
{
   int i;
   double tot;
   static double a[maxn+1][maxn+1], startingfreq[maxn + 1],
      steadystatefreq[4][maxn+1];
   mutq_haploid(a, n[file], t[file], 0.0, startingfreq, steadystatefreq, true);
   tot = 0;
   for (i=0; i<=n[file]; i++)
   {
      expected_freq_vec[i][file] = steadystatefreq[3][i];
      tot += expected_freq_vec[i][file];
   }
   for (i=0; i<=n[file]; i++)
   {
      expected_freq_vec[i][file]/=tot;
//      printf("expected_freq_vec[%d][%d] %lf\n", i, file, expected_freq_vec[i][file]);
   }
}

scale_vector(int n, double expected_freq_vec_unscaled[maxn + 1][max_data_file],
   double expected_freq_vec[maxn + 1][max_data_file], double f0[2][max_data_file],
   int file)
{
   int i;
   for (i=0; i<=n; i++)
   {
      expected_freq_vec[i][file] = expected_freq_vec_unscaled[i][file];
      expected_freq_vec[i][file] *= (1.0 - f0[1][file]);
   }
   expected_freq_vec[0][file] += f0[1][file];
}


print_results_string(FILE *outfile, double log_p)
{
   int i;
   double mu;
   for (i=1; i<=n_data_files; i++)
   {
      fprintf(outfile, "f0_%d %6.5lf ", i, f0[1][i]);
      fprintf(outfile, "Ve%d %lf ", i, error_var[1][i]);
      fprintf(outfile, "N%d %d ", i, n[i]);
   }
   fprintf(outfile, "L %6.5lf ", log_p);
   if (fixed_mut_rate_flag)
   {
      mu = (1.0 - f0[1][1])/(double)(t[1]*n[1]);
      fprintf(outfile, "mu %12.12lf ", mu);
   }
   fprintf(outfile, "\n");
   fflush(outfile);
}

double loss_ve_function(double ve)
{
   int i;
   double num, denom, ratio, res;
   res = undefined;
   for (i=0; i<n_profile_ve; i++)
   {
      if (ve <= profile_ve[0][i])
      {
         if (i==0)
         {
//            printf("loss_ve_function: i= 0, returning undefined\n");
//            monitorinput();
            return undefined;
         }
         if (i==n_profile_ve - 1) return profile_ve[1][i];
//         printf("ve %2.12lf before %2.12lf after %2.12lf\n", ve, profile_ve[0][i-1], profile_ve[0][i]);
//         printf("before %2.12lf after %2.12lf\n", profile_ve[1][i-1], profile_ve[1][i]);
         num = ve - profile_ve[0][i-1];
         denom = profile_ve[0][i] - profile_ve[0][i-1];
         ratio = num/denom;
//         printf("ratio %lf\n", ratio);
         res = profile_ve[1][i-1];
         res += (profile_ve[1][i] - profile_ve[1][i-1])*ratio;
//         printf("res %lf\n", res);
//         monitorinput();
         break;
      }
   }
   return (res);
}

double compute_likelihood_func()
{
   int i, j, file;
   double q, p, p_tot, logl_file, logl = 0, lik_non_seg, exp_freq, var,
      recip_root_two_pi_var, x_upper, x_lower, cdf0, cdf1, logl_ve;
   static double egf_vec[maxn + 1], cdf[maxn + 1];
   if (fixed_mut_rate_flag)
   {
      f0[1][1] = 1 - (double)(t[1]*n[1])*fixed_mut_rate;
   }
   if (same_u)
   {
      f0[1][2] = 1 - (double)(t[2]*n[2])*(1 - f0[1][1])/(double)(t[1]*n[1]);
   }
   for (file = 1; file<=n_data_files; file++)
   {
      logl_file = 0;
//       for (i=0; i<=n[file]; i++)
//       {
//          egf_vec[i] = expected_freq_vec_unscaled[i][file];
//       }
//       dumpvector(egf_vec, 0, n[file], "scaled expected_freq_vec:");
//    monitorinput();
//    printf("compute_likelihood_func f0[1][%d] %lf error_var[1][%d] %lf\n",
//         file, f0[1][file], file, error_var[1][file]);
//      monitorinput();
      scale_vector(n[file], expected_freq_vec_unscaled, expected_freq_vec, f0, file);
      for (i=0; i<=n[file]; i++)
      {
         egf_vec[i] = expected_freq_vec[i][file];
      }
//      dumpvector(egf_vec, 0, n[file], "expected_freq_vec:");
      var = error_var[1][file];
      recip_root_two_pi_var = 1.0/sqrt(2.0*pi*var);
      for (j=1; j<=n[file]; j++)         // Omit the zero class
      {
         exp_freq = (double)j/(double)n[file];
         x_upper = (1 - exp_freq)/sqrt(var);
         x_lower = (0 - exp_freq)/sqrt(var);
         cdf0 = normal_cumulative_probability_function(x_lower);
         cdf1 = normal_cumulative_probability_function(x_upper);
         cdf[j] = (cdf1 - cdf0);
      }
      for (i=1; i<=segregating[file]; i++)
      {
         q = segregating_freqs[i][file];
//        printf("i %d q %lf\n", i, q);
//          monitorinput();

//         if (q == 1.0)
//         {
//          printf("q == 1.0\n");
//            p_tot = expected_freq_vec[n[file]][file];
//         }
//         else
//         {
            p_tot = 0;
            for (j=1; j<=n[file]; j++)         // Omit the zero class
            {
               exp_freq = (double)j/(double)n[file];
               p = recip_root_two_pi_var*exp(-0.5*(exp_freq - q)*(exp_freq - q)/var);
               p /= cdf[j];      // Divide by cumulative normal distribution
//             printf("var %6.10lf exp_freq %lf j %d p %2.16lf\n", var, exp_freq, j, p);
//             monitorinput();
//               if (j==0)
//                  p_tot += p*expected_freq_vec[j][file]*(1 - f0[1][file]);
//               else
                  p_tot += p*expected_freq_vec[j][file];
            }
//         }
         logl_file += log(p_tot);
//         printf("q %lf p_tot %2.16lf logl_file %lf\n", q, p_tot, logl_file);
//         monitorinput();
      }
//    printf("nonsegregating %d expected_freq_vec[0][file] %lf\n", nonsegregating[file], expected_freq_vec[0][file]);
      lik_non_seg = (double)nonsegregating[file]*log(expected_freq_vec[0][file]);
//     printf("Logl segregating %lf Logl nonsegregating %lf\n", logl_file, lik_non_seg);
      logl_file += lik_non_seg;
//      printf("File %d logl_file %lf\n", file, logl_file);
      logl += logl_file;
      logl_ve = 0;
      if (analyse_real_data)
      {
         if (dataset != c_elegans)
         {
            logl_ve = loss_ve_function(error_var[1][file]);
         }
      }
      logl += logl_ve;
//    monitorinput();
//      print_results_string(stdout, logl_file);   
   }
   return(logl);
}


unloadparam2(double paramvec[2][max_data_file], int ind, int *curpos, double x[])
{
   if (paramvec[0][ind]!=0.0)
   {
      *curpos = *curpos + 1;
      paramvec[1][ind] = x[*curpos];
   }
}

double compute_likelihood(double x[])
{
   int curpos, i, undefined_flag = 0;
   double logl;
   curpos = 0;
   eval_num++;
//   printf("compute_likelihood: about to unloadparam2\n"); monitorinput();
   for (i=1; i<=n_data_files; i++)
   {
      unloadparam2(f0, i, &curpos, x);
      unloadparam2(error_var, i, &curpos, x);
   }
//   printf("compute_likelihood: done unloadparam2\n"); monitorinput();
   for (i=1; i<=n_data_files; i++)
   {
//      printf("compute_likelihood: f0[1][%d] %lf error_var[1][%d] %lf\n", i,
//           f0[1][i], i, error_var[1][i]);
//      monitorinput();
      if (error_var[1][i] <= 0)
      {
         logl = undefined;
         undefined_flag = 1;
      }
      if ((f0[1][i] < 0)||(f0[1][i] > 1))
      {
          logl = undefined;
          undefined_flag = 1;
      }
   }
   if (!undefined_flag)
   {
      logl =  compute_likelihood_func();
   }
   if (!silent_mode) print_results_string(stdout, logl);
   return (-logl);
}


computevariableparamnum(int *n)
{
   int i;
   *n = 0;
   for (i=1; i<=n_data_files; i++)
   {
      if (f0[0][i] == 1.0) *n = *n + 1;
      if (error_var[0][i] == 1.0) *n = *n + 1;
//      printf("computevariableparamnum: i %d f0[0][i] %lf error_var[0][i] %lf\n", 
//         i, f0[0][i], error_var[0][i]);
   }
//   if (*n == 0) gabort("computevariableparamnum: n = 0", 0);
}


loadparam2(double paramvec[2][max_data_file], int ind, int *curpos, double **p)
{
   if (paramvec[0][ind]!=0.0)
   {
      *curpos = *curpos + 1;
      p[1][*curpos] = paramvec[1][ind];
   }
}


double dolik()
{
   double l;
   double **p, *y, min, logl;
   int i, j, nfunc, n_parm, curpos;
   double totfreq;
   static double egf_vec[maxn + 1];

   for (i=1; i<=n_data_files; i++)
   {
      generate_expected_frequency(expected_freq_vec_unscaled, i);
//      for (j=0; j<=n[i]; j++)
//      {
//         egf_vec[j] = expected_freq_vec_unscaled[j][i];
//         printf("expected_freq_vec_unscaled[%d][%d] %lf egf_vec[%d] %lf\n",
//            j, i, expected_freq_vec_unscaled[j][i], j, egf_vec[j]);
//      }
//      dumpvector(egf_vec, 0, n[i], "expected_freq_vec_unscaled:"); monitorinput();
   }

   eval_num = 0;
   computevariableparamnum(&n_parm);
//   printf("Done computevariableparamnum: n_parm %d\n", n_parm); monitorinput();
   if (n_parm!=0)
   {
      p=dmatrix(1, n_parm+1, 1, n_parm);
      y=dvector(1, n_parm+1);
      curpos = 0;
      for (i=1; i<=n_data_files; i++)
      {
         loadparam2(f0, i, &curpos, p);
         loadparam2(error_var, i, &curpos, p);
      }
      if (curpos!=n_parm) gabort("dolik: curpos!=n_parm\n", curpos);
   }
   else
   {
      p=dmatrix(1, 2, 1, 1);
      y=dvector(1, 2);
   }
//   printf("About to compute_likelihood\n");   
//   monitorinput();
   y[1] = compute_likelihood(p[1]);

//   printf("Done compute_likelihood\n");   
//   monitorinput();

   if (n_parm==0)
   {
      n_parm = 1;
      goto skip1;
   }

//   exit(0);				// TEMPORARY

   for(i=2;i<=(n_parm+1);i++)
   {
      for(j=1; j<=n_parm; j++) p[i][j] = p[1][j];
      p[i][i-1]/=fac;
//      printf("p[%d][%d] %lf\n", i, i-1, p[i][i-1]);
      y[i] = compute_likelihood(p[i]);
   }

//   exit (0);				// TEMPORARY

   amoeba(p,y,n_parm,criterion, compute_likelihood, &nfunc,
      max_amoeba_evaluations);
skip1: min=y[1];
//   printf("No of evaluations - %d\n", nfunc);
//   printf("MLEs: ");
//   for(i=1; i<=n_parm; i++) printf("%lf ", p[1][i]);
//   printf("LogL %lf\n", -min);
   print_results_string(stdout, -min);
   if (n_parm!=0)
   {
      free_dmatrix(p, 1, n_parm+1, 1, n_parm);
      free_dvector(y, 1, n_parm+1);
   }
   if (!fixed_mut_rate_flag) print_results_string(outfile, -min);
   return(-min);
}


load_results_vec(int ind1, int ind2, double logl)
{
   if (max_results_vec<=6)
      gabort("load_results_vec: max_results_vec too small", max_results_vec);
   results_vec[ind1][ind2][0] = logl;
   results_vec[ind1][ind2][1] = (double)n[1];
   results_vec[ind1][ind2][2] = (double)n[2];
   results_vec[ind1][ind2][3] = f0[1][1];
   results_vec[ind1][ind2][4] = f0[1][2];
   results_vec[ind1][ind2][5] = error_var[1][1];
   results_vec[ind1][ind2][6] = error_var[1][2];
}


search_over_n(double (*f)(), int n_lower, int n_upper, int file, double ve)
{
   int i, j, k, n_lower2, n_upper2;
   double logl, prev_log_l, diff_log_l;
   if (n_upper - n_lower + 1 >=max_n_evals) gabort("Too many n evaluations", 0);
   if (n_data_files==1)
   {
      n_lower2 = 1;
      n_upper2 = 1;      // Dummy values;
   }
   else
   {
      n_lower2 = n_lower;
      n_upper2 = n_upper;
   }
   n_n_evals1 = 0;
//   n_lower = 30;               // TEMPORARY
//   n_upper = 30;
//   n_lower2 = 30;
//   n_upper2 = 30;
//   printf("n_lower2 %d n_upper2 %d\n", n_lower2, n_upper2); monitorinput();
   for (i = n_lower; i<=n_upper; i++)
   {
      n_n_evals2 = 0;
      n_n_evals1++;
      for (j = n_lower2; j<=n_upper2; j++)
      {
         n[1] = i;
         n[2] = j;
         n_n_evals2++;
         for (k=1; k<=n_data_files; k++)
         {
            error_var[1][k] = error_var_start[k];
            f0[1][k] = f0_start[k];
         }
         if (file != -1) error_var[1][file] = ve;
//         error_var[0][1] = 0;                    // TEMPORARY
//         error_var[1][1] = 0.000975;
//         f0[0][1] = 0;
//         f0[1][1] = 0.999859;
         prev_log_l = -9999999999999999999.0;
         for (;;)
         {
            printf("error_var[0][1] %lf error_var[1][1] %lf\n", error_var[0][1], error_var[1][1]);
//            monitorinput();
            logl = (*f)();
            diff_log_l = logl - prev_log_l;
            prev_log_l = logl;
            if (diff_log_l<0.01)
            {
//             printf("prev_log_l %lf logl %lf - breaking\n", prev_log_l, logl);
               break;
            }
         }
//         print_results_string(outfile, logl);
//      printf("n[1] %d n[2] %d logl %lf\n", n[1], n[2], logl); monitorinput();
         load_results_vec(n_n_evals1, n_n_evals2, logl);
//        printf("results_vec[n_n_evals1][n_n_evals2][3] %lf results_vec[n_n_evals1][n_n_evals2][4] %lf\n", results_vec[n_n_evals1][n_n_evals2][3], results_vec[n_n_evals1][n_n_evals2][4]);
//         monitorinput();
      }
   }
}


unload_results_vec(double *ml, int ind1, int ind2)
{
   *ml = results_vec[ind1][ind2][0];
   n[1] = (int)results_vec[ind1][ind2][1];
   n[2] = (int)results_vec[ind1][ind2][2];
   f0[1][1] = results_vec[ind1][ind2][3];
   f0[1][2] = results_vec[ind1][ind2][4];
   error_var[1][1] = results_vec[ind1][ind2][5];
   error_var[1][2] = results_vec[ind1][ind2][6];
}

double find_ml_val_cis()
{
   int i, j;
   double ml;
   ml = undefined;
   for (i=1; i<=n_n_evals1; i++)
   {
      for (j=1; j<=n_n_evals2; j++)
      {
         if ((ml==undefined)||(results_vec[i][j][0] > ml))
         {
            unload_results_vec(&ml, i, j);
         }
      }
   }
   if (ml == undefined) gabort("find_ml_val_cis: ml undefined", 0);
   print_results_string(stdout, ml);
   return(ml);
}


double dolik_search_ve(int file, int do_search, int n_lower, int n_upper)
{
   double ve_start, ve_end, ve_inc, ve, logl, ml;
   int steps, i, j;
   steps = 100;
   if (steps >=max_n_evals) gabort("dolik_search_ve: too many steps");
   ve_start = 0.00005;
   ve_end = 0.005;
   ve = ve_start;
   ve_inc = (ve_end - ve_start)/(double)(steps - 1);
   error_var[0][file] = 0;
   steps = 1;
   for (i=1; i<=steps; i++)
   {
      if (do_search)
      {
         error_var[0][file] = 1;
         error_var[1][file] = 0.001;
         search_over_n(dolik, n_lower, n_upper, file, ve);
//         printf("Done search over N \n"); monitorinput();
         logl = find_ml_val_cis();
      }
      else
      {
         for (j=1; j<=n_data_files; j++)
         {
            error_var[1][j] = error_var_start[j];
            f0[1][j] = f0_start[j];
         }
         error_var[1][file] = ve;
         logl = dolik();
      }
      load_results_vec(i, 0, logl);
//      monitorinput();
      ve += ve_inc;
   }
   ml = undefined;
   for (i=1; i<=steps; i++)
   {
      if ((ml==undefined)||(results_vec[i][0][0] > ml))
      {
         unload_results_vec(&ml, i, 0);
      }
   }
//   print_results_string(outfile, ml);
   if (ml == undefined) gabort("dolik_search_ve: ml undefined", 0);
   printf("dolik_search_ve: ml %lf\n", ml);
//   monitorinput();
   error_var[0][file] = 1;
   logl = dolik();
   return (logl);
}


#define point 0
#define indel 1
#define point_indel 2
#define point_g_to_a 3
#define point_not_g_to_a 4
#define sim_data_file "segregating_freqs.out"
#define madrid_point_data_file "madrid_point.dat"
#define madrid_indel_data_file "madrid_indel.dat"
#define madrid_point_g_to_a_data_file "madrid_g_to_a.dat"
#define madrid_point_not_g_to_a_data_file "madrid_not_g_to_a.dat"
#define florida_point_data_file "florida_point.dat"
#define florida_indel_data_file "florida_indel.dat"
#define florida_point_g_to_a_data_file "florida_g_to_a.dat"
#define florida_point_not_g_to_a_data_file "florida_not_g_to_a.dat"
#define c_elegans_all_data_file "c-elegans-all.dat"

main(int argc, char *argv[])
{
   double u, u_approx, tot,  logl, tot_seg_freq, mu;

   int i, j, test_run, do_search, n_lower, n_upper,
      mutation_type, point_indel_flag, temp_int, stat;
   char *data_file[3], *data_file_indel[3];

   get_silent_mode(argc, argv, &silent_mode);

   point_indel_flag = 0;
   same_u = 0;   
   outfile = openforwrite("analysis1.out", "w");
//   monitorinput();

/*
   getint("Enter N ", &n, 1, infinity);
   getint("Generations ", &t, 3, maxt);
*/
   test_run = 1;
   if (test_run)
   {
      printf("TEST RUN - *************************************************\n");
      analyse_real_data = 1;
      do_search = 1;
      if (analyse_real_data)
      {
         same_u = 1;
         read_profile_ve_file();
         dataset = c_elegans;
         mutation_type=point;
         if (dataset==florida)
         {
            n_data_files = 1;
            t[1] = florida_generations;
            if (mutation_type==point)
               data_file[1] = florida_point_data_file;
            else if (mutation_type==point_g_to_a)
               data_file[1] = florida_point_g_to_a_data_file;
            else if (mutation_type==point_not_g_to_a)
               data_file[1] = florida_point_not_g_to_a_data_file;
            else if (mutation_type==indel)
               data_file[1] = florida_indel_data_file;
            else if (mutation_type==point_indel)
            {
               data_file[1] = florida_point_data_file;
               data_file_indel[1] = florida_indel_data_file;
               point_indel_flag = 1;
            }
            else
               gabort("Invalid value for mutation_type", mutation_type);
         }
         else if (dataset==madrid)
         {
            n_data_files = 1;
            t[1] = madrid_generations;
            if (mutation_type==point)
               data_file[1] = madrid_point_data_file;
            else if (mutation_type==point_g_to_a)
               data_file[1] = madrid_point_g_to_a_data_file;
            else if (mutation_type==point_not_g_to_a)
               data_file[1] = madrid_point_not_g_to_a_data_file;
            else if (mutation_type==indel)
               data_file[1] = madrid_indel_data_file;
            else if (mutation_type==point_indel)
            {
               data_file[1] = madrid_point_data_file;
               data_file_indel[1] = madrid_indel_data_file;
               point_indel_flag = 1;
            }
            else
               gabort("Invalid value for mutation_type", mutation_type);
         }
         else if (dataset==florida_madrid)
         {
            n_data_files = 2;
            t[1] = madrid_generations;
            if (mutation_type==point)
               data_file[1] = madrid_point_data_file;
            else if (mutation_type==indel)
               data_file[1] = madrid_indel_data_file;
            else if (mutation_type==point_indel)
            {
               data_file[1] = madrid_point_data_file;
               data_file_indel[1] = madrid_indel_data_file;
               point_indel_flag = 1;
            }
            else
               gabort("Invalid value for mutation_type", mutation_type);
            t[2] = florida_generations;
            if (mutation_type==point)
               data_file[2] = florida_point_data_file;
            else if (mutation_type==indel)
               data_file[2] = florida_indel_data_file;
            else if (mutation_type==point_indel)
            {
               data_file[2] = florida_point_data_file;
               data_file_indel[2] = florida_indel_data_file;
               point_indel_flag = 1;
            }
            else
               gabort("Invalid value for mutation_type", mutation_type);
            getint("Different [0] or same [1] mutation rates ", &same_u, 0, 1);
            if (same_u)
            {
               f0[0][2] = 0;
            }
            else
            {
               f0[0][2] = 1.0;
            }
         }
         else if (dataset==c_elegans)
         {
            n_data_files = 1;
            t[1] = c_elegans_generations;
            data_file[1] = c_elegans_all_data_file;
         }
         else
         {
            gabort("Unknown data set ", dataset);
         }
      }
      else
      {
         n_data_files = 1;
         t[1] = 200;
         data_file[1] = sim_data_file;
      }
      n_lower = 10;
      n_upper = 200;
      for (i=1; i<=n_data_files; i++)
      {
         f0[0][i] = 1.0;
         f0[1][i] = 0.9996;
// CHANGE THESE 
         error_var[0][i] = 1.0;
         error_var[1][i] = 0.001;
//         printf("Enter error_var[1][%d] ", i); scanf("%lf", &(error_var[1][i]));
      }
   }
   else
   {
      getint("Analyse simulated [0] or real [1] data ", &analyse_real_data, 0, 1);
      if (analyse_real_data) read_profile_ve_file();
      getint("Carry out search? over n [0/1] ", &do_search, 0, 1);
      if (analyse_real_data)
      {
         if (do_search)
         {
            n_lower = 5;
            n_upper = 200;
            getint("Is the mutation rate fixed? ", &fixed_mut_rate_flag, 0, 1);
            if (fixed_mut_rate_flag)
            {
               printf("Enter fixed mutation rate ");
               stat = scanf("%lf", &fixed_mut_rate);
               f0[0][1] = 0;
            }
         }
         else
         {
            fixed_mut_rate_flag = 0;
            printf("Enter N "); stat = scanf("%d", &(n[1]));
         }
         n_data_files = 1;                    // default
         printf("Analyse Florida [%d], Madrid [%d] or combined Madrid-Florida [%d] data ", florida, madrid, florida_madrid);
         stat = scanf("%d", &dataset);
         printf("Analyse point [%d], indel [%d], combined [%d], point G->A [%d] not point G->A [%d] ",
            point, indel, point_indel, point_g_to_a, point_not_g_to_a);
         stat = scanf("%d", &mutation_type);
         if (dataset==florida)
         {
            t[1] = florida_generations;
            if (mutation_type==point)
               data_file[1] = florida_point_data_file;
            else if (mutation_type==point_g_to_a)
               data_file[1] = florida_point_g_to_a_data_file;
            else if (mutation_type==point_not_g_to_a)
               data_file[1] = florida_point_not_g_to_a_data_file;
            else if (mutation_type==indel)
               data_file[1] = florida_indel_data_file;
            else if (mutation_type==point_indel)
            {
               data_file[1] = florida_point_data_file;
               data_file_indel[1] = florida_indel_data_file;
               point_indel_flag = 1;
            }
            else
               gabort("Invalid value for mutation_type", mutation_type);
         }
         else if (dataset==madrid)
         {
            t[1] = madrid_generations;
            if (mutation_type==point)
               data_file[1] = madrid_point_data_file;
            else if (mutation_type==point_g_to_a)
               data_file[1] = madrid_point_g_to_a_data_file;
            else if (mutation_type==point_not_g_to_a)
               data_file[1] = madrid_point_not_g_to_a_data_file;
            else if (mutation_type==indel)
               data_file[1] = madrid_indel_data_file;
            else if (mutation_type==point_indel)
            {
               data_file[1] = madrid_point_data_file;
               data_file_indel[1] = madrid_indel_data_file;
               point_indel_flag = 1;
            }
            else
               gabort("Invalid value for mutation_type", mutation_type);
         }
         else if (dataset==florida_madrid)
         {
            n_data_files = 2;
            t[1] = madrid_generations;
            if (mutation_type==point)
               data_file[1] = madrid_point_data_file;
            else if (mutation_type==indel)
               data_file[1] = madrid_indel_data_file;
            else if (mutation_type==point_indel)
            {
               data_file[1] = madrid_point_data_file;
               data_file_indel[1] = madrid_indel_data_file; 
               point_indel_flag = 1;
            }
            else
               gabort("Invalid value for mutation_type", mutation_type);
            t[2] = florida_generations;
            if (mutation_type==point)
               data_file[2] = florida_point_data_file;
            else if (mutation_type==indel)
               data_file[2] = florida_indel_data_file;
            else if (mutation_type==point_indel)
            {
               data_file[2] = florida_point_data_file;
               data_file_indel[2] = florida_indel_data_file;
               point_indel_flag = 1;
            }
            else
               gabort("Invalid value for mutation_type", mutation_type);
            getint("Different [0] or same [1] mutation rates ", &same_u, 0, 1);
            if (same_u)
            {
               f0[0][2] = 0;        
            }
            else
            {
               f0[0][2] = 1.0;        
            }
         }
         else
         {
            gabort("Unknown data set ", dataset);
         }
         for (i=1; i<=n_data_files; i++)
         {
            if (!fixed_mut_rate_flag) f0[1][i] = 0.99;
            error_var[1][i] = 0.02;
         }
      }
      else
      {
         n_data_files = 1;                    // default
         if (do_search)
         {
            getint("Enter N-lower ", &n_lower, 1, infinity);
            getint("Enter N-upper ", &n_upper, 1, infinity);
         }
         else
         {
            for (i=1; i<=n_data_files; i++)
            {
               getint("Enter N ", &(n[i]), 1, infinity);
            }
         }
         getint("Generations ", &(t[1]), 3, maxt);
         for (i=1; i<=n_data_files; i++)
         {
            printf("Enter starting value for f0 %d ", i);
            stat = scanf("%lf", &(f0[1][i]));
            printf("Enter starting value for error variance %d ", i);
            stat = scanf("%lf", &(error_var[1][i]));
         }
         data_file[1] = sim_data_file;
      }
      if (!fixed_mut_rate_flag) f0[0][1] = 1.0;        
      for (i=1; i<=n_data_files; i++)
      {
        error_var[0][i] = 1.0;        
      }
   }
   for (i=1; i<=n_data_files; i++)
   {
      segregating[i] = 0;
      if (point_indel_flag)
      {
         read_data_file(data_file[i], i, nonsegregating, segregating,
            segregating_freqs, 1);
         temp_int = segregating[i];
         read_data_file(data_file_indel[i], i, nonsegregating, segregating,
            segregating_freqs, segregating[i]+1);
         nonsegregating[i] -= temp_int;
         printf("nonsegregating[%d] %d\n", i, nonsegregating[i]);
//         monitorinput();
      }
      else
      {
//         printf("About to read_data_file data_file[%d] %s\n", i, data_file[i]);
//         monitorinput();
         read_data_file(data_file[i], i, nonsegregating, segregating,
            segregating_freqs, 1);
      }
      error_var_start[i] = error_var[1][i];        
      f0_start[i] = f0[1][i];
   }
   fac = 1.1;
   criterion = 0.0000000001;
   if ((n_data_files==2)||((n_data_files==1)&&((t[1]==madrid_generations)||(t[1]==florida_generations))))
   {
      printf("About to dolik_search_ve n_data_files %d t[1] %d\n", n_data_files, t[1]);
      monitorinput();
      for (i=1; i<=n_data_files; i++)
      {
         logl = dolik_search_ve(i, do_search, n_lower, n_upper);
      }
   }
   else
   {
      if (do_search)
      { 
//         printf("About to search_over_n\n"); monitorinput();
         search_over_n(dolik, n_lower, n_upper, -1, 0);
         logl = find_ml_val_cis();
      }
      else
      {
//         printf("About to dolik\n"); monitorinput();
         logl = dolik();
      }
   }
   for (i=1; i<=n_data_files; i++)
   {
      tot_seg_freq = 0;
      for (j=1; j<=segregating[i]; j++)
      {
         tot_seg_freq += segregating_freqs[j][i]/(double)t[i];
      }
      mu = (1.0 - f0[1][i])/(double)(t[i]*n[i]);
//      printf("segregating[%d] %d nonsegregating[%d] %d tot_seg_freq %lf\n",
//         i, segregating[i], i, nonsegregating[i], tot_seg_freq);
      u_approx = tot_seg_freq/(double)(nonsegregating[i] + segregating[i]);
      printf("Mutation rate estimate %d: ML %10.12lf approx %10.12lf\n",
         i, mu, u_approx);
   }
   print_results_string(outfile, logl);
   print_results_string(stdout, logl);
   fclose (outfile);
   exit(0);
}
