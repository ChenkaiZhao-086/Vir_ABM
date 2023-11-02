#include <RcppArmadillo.h> // 如果使用这个库，则不需要再调用Rcpp
#include <cmath>
#include <omp.h>

using namespace Rcpp;
using namespace std;
using namespace arma;

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

// 计算每一列累计病例数
vec ComputeInfeCase(const mat &agent_status)
{
    umat binaryMatrix = (agent_status == 1);         // Binary matrix with 1s where agent_status is 1
    return conv_to<vec>::from(sum(binaryMatrix, 0)); // Sum each column and convert, 0 is stand for column
    /*
    conv_to<T>::from是Armadillo提供的一个模板函数，用于将数据从一种类型转换为另一种类型。
    在这里，它将上面的求和结果转换（或强制转换）为rowvec类型（即行向量类型）。
    虽然sum函数的结果在大多数情况下已经是行向量，但这个转换确保了类型的一致性，避免了编译器的类型不匹配错误。
    */
}

// [[Rcpp::export]]
// given iniitial bool value of NPI at the begining of the function, it will be changed to true after the function

mat ModelSimCpp(List Parm, int ncores, bool NPI = false, bool BaseImmu = false)
{
    size_t years = Parm["years"];
    float dt = Parm["dt"];
    size_t steps = (52 * years) / dt; // 如果明确不会有负数，可以使用这个
    size_t initial_seeds = Parm["initial_seeds"];
    size_t added_cases = Parm["added_cases"];
    size_t N = Parm["num_of_agent"];
    float beta_seasonal = Parm["beta_seasonal"];
    float beta_phi = Parm["phi"];
    // float gamma = Parm["gamma"];
    // float gamma_prob = 1 - exp(-gamma * dt);

    // beta of each virus
    vec beta_value = {Parm["beta_IFVA"], Parm["beta_IFVB"], Parm["beta_RSV"], Parm["beta_HPIV"],
                      Parm["beta_HMPV"], Parm["beta_HCoV"], Parm["beta_HRV"], Parm["beta_HAdV"]};
    // gamma of each virus
    vec gamma_value = {Parm["gamma_IFVA"], Parm["gamma_IFVB"], Parm["gamma_RSV"], Parm["gamma_HPIV"],
                       Parm["gamma_HMPV"], Parm["gamma_HCoV"], Parm["gamma_HRV"], Parm["gamma_HAdV"]};

    // competition value
    vec comp_value = {Parm["comp_IFVA"], Parm["comp_IFVB"], Parm["comp_RSV"], Parm["comp_HPIV"],
                      Parm["comp_HMPV"], Parm["comp_HCoV"], Parm["comp_HRV"], Parm["comp_HAdV"]};
    // omega value
    vec omega_value = {Parm["omega_IFVA"], Parm["omega_IFVB"], Parm["omega_RSV"], Parm["omega_HPIV"],
                       Parm["omega_HMPV"], Parm["omega_HCoV"], Parm["omega_HRV"], Parm["omega_HAdV"]};

    // Create a series of matrix
    // 使用Armadillo创建一个零矩阵，mat创建的是双精度浮点数，umat创建的是无符号整数
    mat agent_status(N, 8, arma::fill::zeros);           // Store status for each agent
    agent_status.rows(0, initial_seeds - 1).fill(1);     // Put initial infectious seeds
    mat State_probability(N, 8, fill::zeros);            // Create state transition probability matrix
    mat Comparing_matrix(N, 8, fill::zeros);             // Create comparing matrix
    mat results(steps + 1, 9, fill::zeros);              // Create storage for simulation results
    results.col(0) = linspace(0, steps * dt, steps + 1); // 使用linspace来创建一个等距序列，操作与R的(0:30)一样

    // NPI parameters
    size_t NPI_start; // 注意即使不使用这些变量，也需要在这里声明，否则会报错
    size_t NPI_end;
    float NPI_value;
    float decay_coef;
    if (NPI == true)
    {
        NPI_start = Parm["NPI_start"];
        NPI_end = Parm["NPI_end"];
        NPI_value = Parm["NPI_value"];
        decay_coef = Parm["decay_coef"];
    }

    // Base immune parameters
    size_t base_immune_IFVA;
    size_t base_immune_IFVB;
    size_t base_immune_RSV;
    size_t base_immune_HPIV;
    size_t base_immune_HMPV;
    size_t base_immune_HCoV;
    size_t base_immune_HRV;
    size_t base_immune_HAdV;
    if (BaseImmu == true)
    {
        base_immune_IFVA = Parm["base_immune_IFVA"];
        base_immune_IFVB = Parm["base_immune_IFVB"];
        base_immune_RSV = Parm["base_immune_RSV"];
        base_immune_HPIV = Parm["base_immune_HPIV"];
        base_immune_HMPV = Parm["base_immune_HMPV"];
        base_immune_HCoV = Parm["base_immune_HCoV"];
        base_immune_HRV = Parm["base_immune_HRV"];
        base_immune_HAdV = Parm["base_immune_HAdV"];

        agent_status.submat(initial_seeds, 0, initial_seeds + base_immune_IFVA, 0).fill(-1); // 起始行、起始列、结束行和结束列
        agent_status.submat(initial_seeds, 1, initial_seeds + base_immune_IFVB, 1).fill(-1);
        agent_status.submat(initial_seeds, 2, initial_seeds + base_immune_RSV, 2).fill(-1);
        agent_status.submat(initial_seeds, 3, initial_seeds + base_immune_HPIV, 3).fill(-1);
        agent_status.submat(initial_seeds, 4, initial_seeds + base_immune_HMPV, 4).fill(-1);
        agent_status.submat(initial_seeds, 5, initial_seeds + base_immune_HCoV, 5).fill(-1);
        agent_status.submat(initial_seeds, 6, initial_seeds + base_immune_HRV, 6).fill(-1);
        agent_status.submat(initial_seeds, 7, initial_seeds + base_immune_HAdV, 7).fill(-1);
    }

    // Store the initial conditions
    vec infe_cases = ComputeInfeCase(agent_status);
    results.row(0).subvec(1, 8) = infe_cases.t();

    size_t i, j;

    for (size_t ts = 0; ts < steps; ts++) // time
    {
        infe_cases = ComputeInfeCase(agent_status); // Calculate the number of currently infected agents
        results.row(ts + 1).subvec(1, 8) = infe_cases.t();

        int time = floor(results(ts, 0));
        float seasonal_force = (1 + beta_seasonal * cos((2 * M_PI * time / 52) - beta_phi));
        vec beta = beta_value * seasonal_force;

        if (NPI && (time >= NPI_start && time <= NPI_end))
        {
            beta = beta * (1 - (NPI_value * exp(-decay_coef * (time - NPI_start))));
        }

#pragma omp parallel for private(i, j) num_threads(ncores)
        for (size_t i = 0; i < N; i++) // row
        {
            for (size_t j = 0; j < 8; j++) // column
            {
                switch ((int)agent_status(i, j))
                {
                case 0:
                {
                    // pow是元素级别的求指数，与R中的(comp_j / comp_k)^inf_states功能一致，prod是连乘，与R一致
                    vec inf_states = conv_to<vec>::from(agent_status.row(i) == 1);
                    float delta = prod(pow((comp_value(j) / comp_value), inf_states));
                    float lambda = beta(j) * delta * infe_cases(j) / N;
                    State_probability(i, j) = 1 - exp(-lambda * dt);
                }
                break;

                case 1:
                    // State_probability(i, j) = gamma_prob;
                    State_probability(i, j) = 1 - exp(-gamma_value(j) * dt);
                    break;

                case -1:
                    State_probability(i, j) = 1 - exp(-omega_value(j) * dt);
                    break;
                }
            }
        }

        // Calculate state transition matrix based on random number (if State_probability bigger than random number, the agent should transition the state)
        //  arma_rng::set_seed(seeds); //ChatGPT建议这里加上随机数设置，但是似乎在R中使用set.seed()就可以控制
        umat Comparing_matrix = (State_probability > randu<mat>(State_probability.n_rows, State_probability.n_cols));

        // Locate the state that needs to be transformed
        umat condition_S = (agent_status == 0) % Comparing_matrix;
        umat condition_I = (agent_status == 1) % Comparing_matrix;
        umat condition_R = (agent_status == -1) % Comparing_matrix;
        // State transition
        agent_status.elem(find(condition_S)).fill(1);
        agent_status.elem(find(condition_I)).fill(-1);
        agent_status.elem(find(condition_R)).fill(0);

        // Add new cases per dt
        for (size_t col = 0; col < 8; col++)
        {
            uvec S_Locate = find(agent_status.col(col) == 0); // S_Locate是一个标记了哪里是0的向量，如agent_status中第1，3，5个为0，那么S_Locate则返回1，3，5

            if (S_Locate.n_elem >= added_cases)
            {
                S_Locate = shuffle(S_Locate); // 这里将S_Locate的顺序打乱，这样在下一步找到10个填入病例的时候才能实现随机
                for (size_t k = 0; k < added_cases; k++)
                {
                    agent_status(S_Locate(k), col) = 1;
                }
            }
        }
    }
    return results;
}
