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
rowvec ComputeInfeCase(const mat &agent_status)
{
    umat binaryMatrix = (agent_status == 1);            // Binary matrix with 1s where agent_status is 1
    return conv_to<rowvec>::from(sum(binaryMatrix, 0)); // Sum each column and convert
    /*
    conv_to<T>::from是Armadillo提供的一个模板函数，用于将数据从一种类型转换为另一种类型。
    在这里，它将上面的求和结果转换（或强制转换）为rowvec类型（即行向量类型）。
    虽然sum函数的结果在大多数情况下已经是行向量，但这个转换确保了类型的一致性，避免了编译器的类型不匹配错误。
    */
}

// 计算状态转换概率
mat ComputeStateTrans(int time, float dt, float beta_seasonal, float beta_phi, size_t N, float gamma,
                      const mat &agent_status, rowvec infe_cases, rowvec beta_value,
                      rowvec comp_value, rowvec omega_value)
{
    mat state_probability(N, 8, fill::zeros);
    float seasonal_force = (1 + beta_seasonal * cos(2 * M_PI * time / 52) - beta_phi);
    float gamma_prob = 1 - exp(-gamma * dt);

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < 8; j++)
        {
            switch ((int)agent_status(i, j))
            {
            case 0:
            {
                float beta = beta_value(j) * seasonal_force;
                // pow是元素级别的求指数，与R中的(comp_j / comp_k)^inf_states功能一致，prod是连乘，与R一致
                vec inf_states = conv_to<vec>::from(agent_status.row(i) == 1);
                float delta = prod(pow((comp_value(j) / comp_value), inf_states));
                float lambda = beta * delta * infe_cases(j) / N;
                state_probability(i, j) = 1 - exp(-lambda * dt);
            }
            break;

            case 1:
                state_probability(i, j) = gamma_prob;
                break;
            case -1:
                state_probability(i, j) = 1 - exp(-omega_value(j) * dt);
                break;
            }
            /*
             if (agent_status(i, j) == 0)
                        {
                            float beta = beta_value(j) * seasonal_force;
                            float delta = 1;
                            float comp_j = comp_value(j);
                            // pow是元素级别的求指数，与R中的(comp_j / comp_k)^inf_states功能一致，prod是连乘，与R一致
                            vec inf_states = conv_to<vec>::from(agent_status.row(i) == 1);
                            delta = prod(pow((comp_j / comp_value), inf_states));

                            float lambda = beta * delta * infe_cases(j) / N;
                            state_probability(i, j) = 1 - exp(-lambda * dt);
                        }
                        else if (agent_status(i, j) == 1)
                        {
                            state_probability(i, j) = gamma_prob;
                        }
                        else if (agent_status(i, j) == -1)
                        {
                            state_probability(i, j) = 1 - exp(-omega_value(j) * dt);
                        }
            */
        }
    }

    return state_probability;
}

// 状态转换
mat stateTransition(mat &State_transition, mat &agent_status)
{
    // Create the random matrix and comparing matrix
    umat Comparing_matrix = (State_transition > randu<mat>(State_transition.n_rows, State_transition.n_cols));

    /*
        mat agent_status2 = agent_status;
        //在mat中使用等号复制一个矩阵会产生一个新的副本，不会对原始的agent_status产生影响
        //似乎这里不需要复制一个agent_status2了，因为下面的三个condition已经将要修改的位置限定死了
    */

    // State transition
    umat condition_S = (agent_status == 0) % Comparing_matrix;
    umat condition_I = (agent_status == 1) % Comparing_matrix;
    umat condition_R = (agent_status == -1) % Comparing_matrix;

    agent_status.elem(find(condition_S)).fill(1);
    agent_status.elem(find(condition_I)).fill(-1);
    agent_status.elem(find(condition_R)).fill(0);

    return agent_status;
}

// 增加新的易感者
mat AddCases(mat agent_status, int added_cases)
{
    int num_cols = agent_status.n_cols;

    for (size_t col = 0; col < num_cols; col++)
    {
        uvec zero_indicater = find(agent_status.col(col) == 0); // zero_indicater是一个标记了哪里是0的向量，如agent_status中第1，3，5个为0，那么zero_indicater则返回1，3，5

        if (zero_indicater.n_elem >= added_cases)
        {
            zero_indicater = shuffle(zero_indicater); // 这里将zero_indicater的顺序打乱，这样在下一步找到10个填入病例的时候才能实现随机
            for (size_t k = 0; k < added_cases; k++)
            {
                agent_status(zero_indicater(k), col) = 1;
            }
        }
    }

    return agent_status;
}

// [[Rcpp::export]]
mat ModelSimCpp(List Parm)
{
    size_t years = Parm["years"];
    float dt = Parm["dt"];
    size_t steps = (52 * years) / dt; // 如果明确不会有负数，可以使用这个
    size_t initial_seeds = Parm["initial_seeds"];
    size_t added_cases = Parm["added_cases"];
    size_t N = Parm["num_of_agent"];
    float beta_seasonal = Parm["beta_seasonal"];
    float beta_phi = Parm["phi"];
    float gamma = Parm["gamma"];
    float gamma_prob = 1 - exp(-gamma * dt);

    // beta of each virus
    vec beta_value = {Parm["beta_IFVA"], Parm["beta_IFVB"], Parm["beta_RSV"], Parm["beta_HPIV"],
                      Parm["beta_HMPV"], Parm["beta_HCoV"], Parm["beta_HRV"], Parm["beta_HAdV"]};
    // competition value
    vec comp_value = {Parm["comp_IFVA"], Parm["comp_IFVB"], Parm["comp_RSV"], Parm["comp_HPIV"],
                      Parm["comp_HMPV"], Parm["comp_HCoV"], Parm["comp_HRV"], Parm["comp_HAdV"]};
    // omega value
    vec omega_value = {Parm["omega_IFVA"], Parm["omega_IFVB"], Parm["omega_RSV"], Parm["omega_HPIV"],
                       Parm["omega_HMPV"], Parm["omega_HCoV"], Parm["omega_HRV"], Parm["omega_HAdV"]};

    mat agent_status(N, 8, arma::fill::zeros);
    // Put initial seeds
    agent_status.rows(0, initial_seeds - 1).fill(1);

    // agent_status2 = copy(agent_status);

    // Create state transition probability matrix
    mat State_probability(N, 8, fill::zeros);

    // Create comparing matrix
    mat Comparing_matrix(N, 8, fill::zeros);

    // Create storage for simulation results
    mat results(steps + 1, 9, fill::zeros);              // 使用Armadillo创建一个零矩阵，mat创建的是双精度浮点数，umat创建的是无符号整数
    results.col(0) = linspace(0, steps * dt, steps + 1); // 使用linspace来创建一个等距序列，操作与R的(0:30)一样

    // Store the initial conditions
    rowvec infe_cases = ComputeInfeCase(agent_status);
    results.row(0).subvec(1, 8) = infe_cases;

    /*
    下面这一段是使用Rcpp标准库完成，可能效率低，但是不会出错
        NumericVector infe_cases(8);
        int nrows = agent_status.nrow();
        // #pragma omp parallel for schedule(static)
        for (int j = 0; j < 8; j++)
        {
            int count = 0;
            for (int i = 0; i < nrows; i++)
            {
                if (agent_status(i, j) == 1)
                {
                    count++;
                }
            }
            infe_cases[j] = count;
        }
        results.row(0).subvec(1, 8) = infe_cases.t(); // .t() is for transpose since vec is a column vector 这一行是使用Armadillo实现

    */
    size_t ts, i, j;

    for (size_t ts = 0; ts < steps * dt; ts++)
    {
        infe_cases = ComputeInfeCase(agent_status);
        results.row(ts + 1).subvec(1, 8) = infe_cases;

        int time = floor(results(ts, 0));
        float seasonal_force = (1 + beta_seasonal * cos((2 * M_PI * time / 52) - beta_phi));

#pragma omp parallel for private(i, j)
        for (size_t i = 0; i < N; i++)
        {
            for (size_t j = 0; j < 8; j++)
            {
                switch ((int)agent_status(i, j))
                {
                case 0:
                {
                    float beta = beta_value(j) * seasonal_force;
                    // pow是元素级别的求指数，与R中的(comp_j / comp_k)^inf_states功能一致，prod是连乘，与R一致
                    vec inf_states = conv_to<vec>::from(agent_status.row(i) == 1);
                    float delta = prod(pow((comp_value(j) / comp_value), inf_states));
                    float lambda = beta * delta * infe_cases(j) / N;
                    State_probability(i, j) = 1 - exp(-lambda * dt);
                }
                break;

                case 1:
                    State_probability(i, j) = gamma_prob;
                    break;

                case -1:
                    State_probability(i, j) = 1 - exp(-omega_value(j) * dt);
                    break;
                }
            }
        }

        // State_probability = ComputeStateTrans(time, dt, beta_seasonal, beta_phi, N, gamma, agent_status,
        //                                      infe_cases, beta_value, comp_value, omega_value);

        umat Comparing_matrix = (State_probability > randu<mat>(State_probability.n_rows, State_probability.n_cols));

        /*
            mat agent_status2 = agent_status;
            //在mat中使用等号复制一个矩阵会产生一个新的副本，不会对原始的agent_status产生影响
            //似乎这里不需要复制一个agent_status2了，因为下面的三个condition已经将要修改的位置限定死了
        */

        // State transition
        umat condition_S = (agent_status == 0) % Comparing_matrix;
        umat condition_I = (agent_status == 1) % Comparing_matrix;
        umat condition_R = (agent_status == -1) % Comparing_matrix;

        agent_status.elem(find(condition_S)).fill(1);
        agent_status.elem(find(condition_I)).fill(-1);
        agent_status.elem(find(condition_R)).fill(0);

        // agent_status = stateTransition(State_probability, agent_status);

        int num_cols = agent_status.n_cols;

        for (size_t col = 0; col < num_cols; col++)
        {
            uvec zero_indicater = find(agent_status.col(col) == 0); // zero_indicater是一个标记了哪里是0的向量，如agent_status中第1，3，5个为0，那么zero_indicater则返回1，3，5

            if (zero_indicater.n_elem >= added_cases)
            {
                zero_indicater = shuffle(zero_indicater); // 这里将zero_indicater的顺序打乱，这样在下一步找到10个填入病例的时候才能实现随机
                for (size_t k = 0; k < added_cases; k++)
                {
                    agent_status(zero_indicater(k), col) = 1;
                }
            }
        }
    }
    return results;
}
