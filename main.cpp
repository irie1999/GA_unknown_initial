#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <mpi.h>

#include "GA.h"
#include "agent.h"
#include "memory_allocate.h"

int main(int argc, char **argv){
    std::random_device rnd;
    std::mt19937 mt(rnd());
    //std::mt19937 rnd(1); 
    
    MPI::Init(argc, argv);

    const int rank = MPI::COMM_WORLD.Get_rank();
    const int size = MPI::COMM_WORLD.Get_size();
    
    /*MPI parameter*/
    const int Range { Number_of_Individual / size };

    std::cout << "rank= " << rank << ", size= " << size << std::endl;
    
    double **S = allocate_memory2d(3, 801, 0.0);  /*観測した電界強度の変化率*/
    double **Ei_tm = allocate_memory2d(3, 801, 0.0); /*個体のβとh'を与えた得られた電界強度"*/
    double **parameter = allocate_memory2d(4, Number_of_Individual, 0.0);  /*各プロセスに渡すパラメタ*/
    double score[Number_of_Individual];     /*受信するときに用いるスコア*/
    double MAX[Number_of_Generation + 1];   /*スコア最大値を格納*/
    double **MAX_parameter = allocate_memory2d(4, Number_of_Generation, 0.0);  /*スコアの最も高いパラメタを格納*/
    double score_average[Number_of_Generation + 1]; /*平均値を格納*/
    double roulette[Number_of_Individual];  /*ルーレット*/
    Agent agent[2][Number_of_Individual];   /*個体*/
    double beta[3], h_prime[3];
    double time[3];


    ////////////////////////////////////////
    int searched_time = 6;
    time[0] = 6.16667;
    time[1] = 6.33333;
    time[2] = 6.5; 
    beta[0] = 0.493661;  
    h_prime[0] = 77.69128;  
    ///////////////////////////////////////////
    ////////////////////////////////////////
    // int searched_time = 8;
    // time[0] = 8.16667;
    // time[1] = 8.33333;
    // time[2] = 8.5; 
    // beta[0] = 0.867372;  
    // h_prime[0] = 68.451;  
    ///////////////////////////////////////////
    ////////////////////////////////////////
    // int searched_time = 10;
    // time[0] = 10.16667;
    // time[1] = 10.33333;
    // time[2] = 10.5; 
    // beta[0] = 0.933215;  
    // h_prime[0] = 66.2318;  
    ///////////////////////////////////////////


    create_ind(agent[0]); /*初期ランダム遺伝子の作成*/
    
    input(S, 1, "Si_tm_1.dat");  /*t_1の時の観測した電界強度*/
    input(S, 2, "Si_tm_2.dat");  /*t_2の時の観測した電界強度*/
    input(Ei_tm, 0, "E_i_0.dat"); /*最初の時刻の電界強度*/
    
    if(rank == 0){
        std::cout << "世代= " << Number_of_Generation << std::endl
                  << "個体= " << Number_of_Individual << std::endl
                  << "観測間隔= " << step_km << "km" << std::endl
                  << "探索時間= " << searched_time << ":00" << std::endl << std::endl;
    }
    
    for(int n_generation = 0; n_generation < Number_of_Generation - 1; n_generation++){
        if(rank ==0){
            std::cout << "世代= " << n_generation << std::endl;
        }
        const int PARENT { n_generation % 2 };
        const int CHILD { (n_generation + 1) % 2};

        /*rank0がパラメタを求める*/
        if(rank == 0){
            for(int i = 0; i < Number_of_Individual; i++){ /*Range=3*/
                agent[PARENT][i].set_parameter();
                parameter[0][i] = agent[PARENT][i].parameter_beta_1;
                parameter[1][i] = agent[PARENT][i].parameter_beta_2;
                parameter[2][i] = agent[PARENT][i].parameter_h_prime_1;
                parameter[3][i] = agent[PARENT][i].parameter_h_prime_2;
            }   
        }

        /*rank0から各rankへパラメタを送信*/
        if(rank == 0){
            for(int r =1; r < size; r++){
                for(int k = 0; k < 4; k++){
                MPI::COMM_WORLD.Send( parameter[k] + r*Range, Range, MPI::DOUBLE, r, k);
                }
            }
        }
        else{   /*rank0以外はパラメタを受信*/
            for(int k = 0; k < 4; k++){
            MPI::COMM_WORLD.Recv( parameter[k] + rank*Range, Range, MPI::DOUBLE, 0, k);
            }
        }
        
        /*スコアを計算*/
        for(int i = rank * Range; i < (rank + 1) * Range; i++){ /*全個体を各コアに分割*/
        cal_ind(S, Ei_tm, i, parameter, score, time, beta, h_prime);
        }
        
        if(rank != 0){  /*rank0にスコアを送信*/
            MPI::COMM_WORLD.Send(score + rank*Range, Range, MPI::DOUBLE, 0, 0);
        }
        else{          /*rank0以外はスコアを受信*/
            for(int r = 1; r < size; r++){
                MPI::COMM_WORLD.Recv(score + r*Range, Range, MPI::DOUBLE, r, 0);
            }
        }

        if(rank == 0){  /*受け取ったスコアを個体のスコアに代入*/
            for(int n_individual = 0; n_individual < Number_of_Individual; n_individual++){
                agent[PARENT][n_individual].score = score[n_individual];
            }

            /*スコア順にソート*/
            sort_ind(agent[PARENT]);

            /*スコアが最も高いパラメタを格納*/
            MAX_parameter[0][n_generation] = agent[PARENT][0].parameter_beta_1;
            MAX_parameter[1][n_generation] = agent[PARENT][0].parameter_beta_2;
            MAX_parameter[2][n_generation] = agent[PARENT][0].parameter_h_prime_1;
            MAX_parameter[3][n_generation] = agent[PARENT][0].parameter_h_prime_2;

            /*各世代の最大値を格納*/
            MAX[n_generation] = agent[PARENT][0].score;

            /*ルーレットと平均値作成*/
            compose_roulette(Number_of_Individual, agent[PARENT], roulette, score_average, n_generation);    
        
            /*選択と交叉*/
            selection_crossover(roulette, agent[PARENT], agent[CHILD]);
        
            /*突然変異*/
            mutate_ind(agent[CHILD]);
        }
        
    }
    
    ///////*最終世代の最もスコアが高いものを判断*/////////
    const int PARENT { (Number_of_Generation - 1) % 2 };
    if(rank ==0){
        std::cout << "世代= " << Number_of_Generation - 1 << std::endl;
    }

    /*rank0がパラメタを求める*/
    if(rank == 0){
        for(int i = 0; i < Number_of_Individual; i++){ /*Range=3*/
            agent[PARENT][i].set_parameter();
            parameter[0][i] = agent[PARENT][i].parameter_beta_1;
            parameter[1][i] = agent[PARENT][i].parameter_beta_2;
            parameter[2][i] = agent[PARENT][i].parameter_h_prime_1;
            parameter[3][i] = agent[PARENT][i].parameter_h_prime_2;
        }   
    }

    /*rank0から各rankへパラメタを送信*/
    if(rank == 0){
        for(int r =1; r < size; r++){
            for(int k = 0; k < 4; k++){
            MPI::COMM_WORLD.Send( parameter[k] + r*Range, Range, MPI::DOUBLE, r, k);
            }
        }
    }
    else{
        for(int k = 0; k < 4; k++){
        MPI::COMM_WORLD.Recv( parameter[k] + rank*Range, Range, MPI::DOUBLE, 0, k);
        }
    }

    /*スコアを計算*/
    for(int i = rank * Range; i < (rank +1) * Range; i++){
        final_cal_ind(S, Ei_tm, i, parameter, score, time, beta, h_prime);
    }

    if(rank != 0){  /*rank0にスコアを送信*/
        MPI::COMM_WORLD.Send(score + rank*Range, Range, MPI::DOUBLE, 0, 0);
    }
    else{ 
        for(int r = 1; r < size; r++){
            MPI::COMM_WORLD.Recv(score + r*Range, Range, MPI::DOUBLE, r, 0);
        }
    }

    if(rank == 0){
        double sum = 0.0;
        for(int n_individual = 0; n_individual < Number_of_Individual; n_individual++){
            agent[PARENT][n_individual].score = score[n_individual];
            sum += score[n_individual];
        }
        
        /*最終世代の平均値を格納*/
        score_average[Number_of_Generation - 1] = sum / Number_of_Individual;

        /*スコア順にソート*/
        sort_ind(agent[PARENT]);

        /*スコアが最も高いパラメタを格納*/
        MAX_parameter[0][Number_of_Generation - 1] = agent[PARENT][0].parameter_beta_1;
        MAX_parameter[1][Number_of_Generation - 1] = agent[PARENT][0].parameter_beta_2;
        MAX_parameter[2][Number_of_Generation - 1] = agent[PARENT][0].parameter_h_prime_1;
        MAX_parameter[3][Number_of_Generation - 1] = agent[PARENT][0].parameter_h_prime_2;

        /*各世代の最大値を格納*/
        MAX[Number_of_Generation - 1] = agent[PARENT][0].score;

        std::ofstream ofs("../data/score,time" + std::to_string(searched_time) + ",step " + std::to_string(step_km) + "km,gene " + std::to_string(Number_of_Generation) + ",ind " + std::to_string(Number_of_Individual) + ".dat");
        for(int n_generation = 0; n_generation < Number_of_Generation; n_generation++){
            ofs << n_generation << " " << MAX[n_generation] << " " << score_average[n_generation] << std::endl;
        }
        ofs.close();

        std::cout << "beta_1= " << MAX_parameter[0][Number_of_Generation - 1] << " beta_2= " << MAX_parameter[1][Number_of_Generation - 1]  << " h_prime_1= " << MAX_parameter[2][Number_of_Generation - 1]  
                << " h_prime_2= " << MAX_parameter[3][Number_of_Generation - 1] << " max= " << MAX[Number_of_Generation - 1] << std::endl;

        // double beta[3], h_prime[3];
        // double time[3];
        // time[0] = 6.16667;
        // time[1] = 6.33333;
        // time[2] = 6.5; 
        // beta[0] = 0.493661;  
        // h_prime[0] = 77.69128;

        for(int t = 1; t < 3; t++){
        std::cout << "beta_" + std::to_string(t) << "= " << ( MAX_parameter[1][Number_of_Generation - 1] * pow((time[t] - time[0]),2) + MAX_parameter[0][Number_of_Generation - 1] * (time[t] - time[0]) + beta[0] ) << std::endl
                  << "h_prime" + std::to_string(t)  << "= " << (MAX_parameter[3][Number_of_Generation - 1] * pow((time[t] - time[0]),2) + MAX_parameter[2][Number_of_Generation - 1] * (time[t] - time[0]) + h_prime[0] ) << std::endl;
        }

        std::ofstream ofs_1("../data/Answer,time" + std::to_string(searched_time) + ",step " + std::to_string(step_km) + "km,gene " + std::to_string(Number_of_Generation) + ",ind " + std::to_string(Number_of_Individual) + ".dat");
        for(int t = 1; t < 3; t++){
            ofs_1   << "Answer beta_" + std::to_string(t) << "= " << ( MAX_parameter[1][Number_of_Generation - 1] * pow((time[t] - time[0]),2) + MAX_parameter[0][Number_of_Generation - 1] * (time[t] - time[0]) + beta[0] ) << std::endl
                    << "Answer h_prime" + std::to_string(t)  << "= " << ( MAX_parameter[3][Number_of_Generation -1] * pow((time[t] - time[0]),2) + MAX_parameter[2][Number_of_Generation -1] * (time[t] - time[0]) + h_prime[0] ) << std::endl << std::endl;
        }

        std::cout << std::endl << std::endl;

        for(int i = 0; i < Number_of_Generation; i++){
            for(int t = 1; t < 3; t++){
            ofs_1   << "beta_" + std::to_string(t) << "= " << ( MAX_parameter[1][i] * pow((time[t] - time[0]),2) + MAX_parameter[0][i] * (time[t] - time[0]) + beta[0] ) << std::endl
                    << "h_prime" + std::to_string(t)  << "= " << ( MAX_parameter[3][i] * pow((time[t] - time[0]),2) + MAX_parameter[2][i] * (time[t] - time[0]) + h_prime[0] ) << std::endl << std::endl;
            }
        }
        ofs_1.close();
    
    }

    deallocate_memory2d(S);
    deallocate_memory2d(parameter);
    deallocate_memory2d(MAX_parameter);

    MPI::Finalize();

    return 0;

}
