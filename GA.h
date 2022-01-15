#ifndef GA_H_
#define GA_H_

#include <iostream>
#include <cmath>

constexpr int N_bit_parameter_beta_1 { 8 };
constexpr int N_bit_parameter_beta_2 { 8 };
constexpr int N_bit_parameter_h_prime_1 { 8 };
constexpr int N_bit_parameter_h_prime_2 { 8 };
constexpr int N_bit_total { N_bit_parameter_beta_1 + N_bit_parameter_beta_2 + N_bit_parameter_h_prime_1 + N_bit_parameter_h_prime_2 };
constexpr int N_bit_parameters_beta { N_bit_parameter_beta_1 + N_bit_parameter_beta_2 };
constexpr int N_bit_parameters_h_prime { N_bit_parameter_h_prime_1 + N_bit_parameter_h_prime_2 };
// constexpr double t[0] { 6.16667 }; /*6:10*/
// constexpr double t[1] { 6.33333 }; 
// constexpr double t[2] { 6.5 }; 


constexpr double i32 { 4294967296.0 }; /* 2^32 */
constexpr double MUTATION { 0.03 }; /* 突然変異の確率 */
constexpr int Number_of_Individual { 48 };  /*n個体*/
constexpr int Number_of_Generation { 50 };  /*n世代*/

//constexpr int MAX_Generation { 120 };

constexpr int step_km { 10 }; /*100km~900kmでstep_km間隔で観測点を置く*/
constexpr int GA_Nr { (800 / step_km) + 1 };
constexpr int M { 2 };
constexpr double p_beta { 1e3 };
constexpr double p_h_prime { 1.0 };

void create_ind(class Agent *agent);
void cal_ind(double **S, double **Ei_tm, int i, double **parameter, double *score);
int bin2dec(const int N_bit_initial, const int N_bit_end, bool *binary);
void compose_roulette(const int N, class Agent *agent, double *roulette, double *score_average, int n_generation);
void crossover(int head, class Agent *p, class Agent *c, int *s);
void selection_crossover(double *roulette, Agent *p, Agent *c);
double fitting(double parameter_beta_1, double parameter_beta_2, 
            double parameter_h_prime_1, double parameter_h_prime_2, double **S, double **Ei_tm);
void sort_ind(class Agent *p);
void mutate_ind(Agent *c);
void final_cal_ind(double **S,double **Ei_tm, int i,double **parameter, double *score);
double *cal_fdtd(double beta, double h_prime, int t, double **Ei_tm);
void input(double **S, int time, std::string filename);
//double **allocate_memory2d(int m, int n, double ini_v);
#endif



