#include <random>
#include <fstream>

#include "GA.h"
#include "agent.h"

#include "memory_allocate.h"


double fitting(double parameter_beta_1, double parameter_beta_2, 
            double parameter_h_prime_1, double parameter_h_prime_2, 
            double **S, double **Ei_tm, double *time, double *beta, double *h_prime){
    double v = 0.0; /*score*/
    double **s = allocate_memory2d(3, 801, 0.0);  /*探索する個体の電界強度の変化率*/
    // double beta[3], h_prime[3];
    // double time[3];
    // //double compare_1, compare_2, compare_3;  /*評価関数の3つの部分を比較する変数*/
    // time[0] = 8.16667;
    // time[1] = 8.33333;
    // time[2] = 8.5; 
    // beta[0] = 0.493661;  /*パラメタを規格化 beta_0 = 0.7 beta[0] = beta0 / beta_0*/
    // h_prime[0] = 77.69128;  /*パラメタを規格化 h_prime0 = 75*/

    for(int t = 1; t <= 2; t++){
        beta[t] = parameter_beta_2 * pow((time[t] - time[0]), 2) + parameter_beta_1 * (time[t] - time[0]) + beta[0];
        h_prime[t] = parameter_h_prime_2 * pow((time[t] - time[0]), 2) + parameter_h_prime_1 * (time[t] - time[0]) + h_prime[0];
        cal_fdtd(beta[t], h_prime[t] * 1e3, t, Ei_tm); /*betaとh'を代入して電界を返す*/
    }
    
    //double sum;
    // ////////////////*100km~900kmのstep_km毎に観測点を設置*/////////////////////
    for(int t_m = 1; t_m <= M; t_m++){
        for(int i = 0; i < 801; i = i + step_km){ 
            s[t_m][i] = (Ei_tm[t_m][i] - Ei_tm[t_m - 1][i]) / (time[t_m] - time[t_m -1]);

            v += std::pow(std::abs( S[t_m][i] - s[t_m][i] ), 2);
            //sum += std::pow(std::abs( S[t_m][i] - s[t_m][i] ), 2);
        }
    }
    //v += p_beta * std::pow(std::abs( beta[M] - beta[0]), 2) 
                // + p_h_prime * std::pow(std::abs( h_prime[M] - h_prime[0]), 2);

    // std::ofstream ofs("../data/compare" + std::to_string(Number_of_Generation) + "," + std::to_string(Number_of_Individual) + ".dat");
    // ofs << "compare_1= " << v << std::endl 
    //     << "std::pow(std::abs( S[t_m][i] - s[t_m][i] ), 2)= " << sum << std::endl 
    //     << "compare_3= " << std::pow(std::abs( h_prime[M] - h_prime[0]), 2) << std::endl;
    // ofs.close();

    deallocate_memory2d(s);
    

    v = GA_Nr * M / v;

    return v;
}

void cal_ind(double **S, double **Ei_tm, int i, double **parameter,
             double *score, double *time, double *beta, double *h_prime){
    score[i] = fitting(parameter[0][i], parameter[1][i], parameter[2][i], parameter[3][i], 
                       S, Ei_tm, time, beta, h_prime); 
}

void create_ind(Agent *agent){
    std::random_device rnd;
    std::mt19937 mt(rnd());
    //std::mt19937 rnd(1); 
    for(int i = 0; i < Number_of_Individual; i++){
        for(int n = 0; n < N_bit_total; n++){
            agent[i].Gene[n] = rnd() % 2;
        }
    }
}

void compose_roulette(const int N, Agent *agent, double *roulette, double *score_average, int n_generation){/*ルーレット作成*/
    double sum = 0.0; 
    for(int i = 0; i < Number_of_Individual; i++){
        sum += agent[i].score;
    }
    score_average[n_generation] = sum / Number_of_Individual;
    roulette[0] = agent[0].score / sum;
    for(int i = 1; i < Number_of_Individual ; i++){
        roulette[i] = roulette[i-1] + agent[i].score / sum;
    }
}

void crossover(int head, Agent *p, Agent *c, int *s){ /*交叉*/
    //std::mt19937 rnd(1); 
    std::random_device rnd;
    std::mt19937 mt(rnd());
  
    for(int i = 0; i < N_bit_total; i++){
        if(rnd() / i32 < 0.5){  /*入れ替えなし*/  
            c[head].Gene[i] = p[s[0]].Gene[i];  
            c[head+1].Gene[i] = p[s[1]].Gene[i];
        }
        else{                   /*入れ替えあり*/
            c[head].Gene[i] = p[s[1]].Gene[i];
            c[head+1].Gene[i] = p[s[0]].Gene[i];
        }
    }
}

void selection_crossover(double *roulette, Agent *p, Agent *c){
    std::random_device rnd;
    std::mt19937 mt(rnd());
    //std::mt19937 rnd(1);
    for(int i = 2; i < Number_of_Individual; i+=2){
            int sict[2];
            for(int j = 0; j <2 ; j++){
                double rnd_num = rnd() / i32;
                int k = 0;
                while( roulette[k] < rnd_num){  /*親を2体選ぶルーレット*/
                    k++;
                }
                sict[j] = k;
            }
            crossover(i, p, c, sict); /*交叉*/
    }
    /*エリート戦略*/
    for(int i = 0; i < 2; i++){ /*スコア上位2体を無条件に選択*/
        for(int n = 0; n < N_bit_total; n++){
            c[i].Gene[n] = p[i].Gene[n];
        }
    }
}

void mutate_ind(Agent *c){
    std::random_device rnd;
    std::mt19937 mt(rnd());
    //std::mt19937 rnd(1);
    for(int i = 2; i < Number_of_Individual; i++){
        for (int j = 0; j < N_bit_total; j++){
            if(rnd() / i32 < MUTATION ){
                c[i].Gene[j] = !c[i].Gene[j];  /*0と1を反転*/
            }
        }
    }
}

void final_cal_ind(double **S, double **Ei_tm, int i,double **parameter, 
                   double *score, double *time, double *beta, double *h_prime){
    score[i] = fitting(parameter[0][i], parameter[1][i], parameter[2][i], parameter[3][i], 
                        S, Ei_tm, time, beta, h_prime); 
}

void sort_ind(Agent *p){
    Agent tmp;
    for(int i = 0; i < Number_of_Individual - 1; i++){
        for(int j = i + 1 ; j < Number_of_Individual; j++){
            if( p[i].score < p[j].score ){
                tmp = p[i];
                p[i] = p[j];
                p[j] = tmp;
            } 
        }
    }
}


/*2進数→10進数変換*/
int bin2dec(const int N_bit_initial, const int N_bit_end, bool *binary){
    int v { 0 };
    int base { 1 };

    for(int i = N_bit_end - 1; i >= N_bit_initial; i--){
        v += int(binary[i]) * base;
        base *= 2;
    }
    return v;
}

// double **allocate_memory2d(int m, int n, double ini_v){
//   double **v = new double* [m];
//   v[0] = new double [m*n];
//   for(int i = 0; i < m; i++){
//     v[i] = v[0] + i*n;
//     for(int j = 0; j < n; j++){
//       v[i][j] = ini_v;
//     }
//   }
//   return v;
// }

