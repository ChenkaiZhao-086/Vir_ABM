#include <RcppArmadillo.h>
#include <cmath>
#include <omp.h>

using namespace Rcpp;
using namespace std;
using namespace arma;

// [[Rcpp::plugins("cpp17")]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat ModelSimCpp(List Parm, int ncores, bool NPI = false, bool BaseImmu = false, double virus4_start_time = -1)
{
    // Base parameters
    const size_t NUM_VIRUSES = 4;
    size_t N = Parm["num_of_agent"];
    double beta_seasonal = Parm["beta_seasonal"];
    double beta_phi = Parm["phi"];
    size_t years = Parm["years"];
    double dt = Parm["dt"];
    size_t steps = (365 * years) / dt;
    size_t initial_seeds = Parm["initial_seeds"];
    size_t added_cases = Parm["added_cases"];

    double Penal = Parm["Penal"];

    vec beta_value = {Parm["beta_IFVA"], Parm["beta_IFVB"],
                      Parm["beta_RSV"], Parm["beta_RV"]};
    vec gamma_value = {Parm["gamma_IFVA"], Parm["gamma_IFVB"],
                       Parm["gamma_RSV"], Parm["gamma_RV"]};
    vec comp_value = {Parm["comp_IFVA"], Parm["comp_IFVB"],
                      Parm["comp_RSV"], Parm["comp_RV"]};
    vec omega_value = {Parm["omega_IFVA"], Parm["omega_IFVB"],
                       Parm["omega_RSV"], Parm["omega_RV"]};

    // gamma and omega probability
    vec gamma_prob = 1 - exp(-gamma_value * dt);
    vec omega_prob = 1 - exp(-omega_value * dt);
    const double N_inv = 1.0 / N;
    const double dt_seasonal = 2 * M_PI / 365;

    bool virus4_delayed = (virus4_start_time > 0);

    // -----------------------------
    // 初始化状态
    // -----------------------------
    // 状态编码：
    //    0: S (易感)
    //    1: I (感染)
    //   -1: R (康复/免疫)
    // -----------------------------
    mat agent_status(N, NUM_VIRUSES, fill::zeros);
    // 对于每个病毒，如果是病毒4且为延时启动，则跳过，否则设置初始感染（seeds）
    if (initial_seeds > 0)
    {
        for (size_t virus = 0; virus < NUM_VIRUSES; virus++)
        {
            if (virus == 3 && virus4_delayed)
                continue;
            agent_status.submat(0, virus, initial_seeds - 1, virus).fill(1);
        }
    }
    mat State_probability(N, NUM_VIRUSES, fill::zeros);
    mat results(steps + 1, NUM_VIRUSES + 1, fill::zeros);
    results.col(0) = linspace(0, steps * dt, steps + 1);

    // 非药物干预 (NPI)
    double NPI_value = 0;
    double decay_coef = 0;
    if (NPI)
    {
        NPI_value = Parm["NPI_value"];
        decay_coef = Parm["decay_coef"];
    }

    // Base immunity：只处理非延迟启动的病毒4
    if (BaseImmu)
    {
        double base_immune = Parm["base_immune"];
        int ImmuPop = floor(N * base_immune);
        if (ImmuPop + initial_seeds > N)
        {
            ImmuPop = N - initial_seeds;
        }
        for (size_t col = 0; col < NUM_VIRUSES; col++)
        {
            // 若是病毒4且为延迟启动则跳过
            if (col == 3 && virus4_delayed)
                continue;
            uvec S_Locate = find(agent_status.col(col) == 0);
            S_Locate = shuffle(S_Locate);
            for (size_t k = 0; k < ImmuPop; k++)
            {
                agent_status(S_Locate(k), col) = -1;
            }
        }
    }

    // 记录初始状态
    vec infe_cases = conv_to<vec>::from(sum(agent_status == 1, 0));
    results.row(0).subvec(1, NUM_VIRUSES) = infe_cases.t();

    // Penal lookup table
    // 最大同时感染数为 4 则 PenalNum 的取值范围 ∈ {0,1,2,3,4}
    std::array<double, 5> penal_lookup;
    for (int p = 0; p <= 4; p++)
    {
        penal_lookup[p] = pow(Penal, (double)p);
    }

    // 主模拟循环（时间步循环）
    for (size_t ts = 0; ts < steps; ts++)
    {
        // 更新每个病毒当前的感染数等状态分布
        infe_cases = conv_to<vec>::from(sum(agent_status == 1, 0));
        vec sus_cases = conv_to<vec>::from(sum(agent_status == 0, 0));
        results.row(ts + 1).subvec(1, NUM_VIRUSES) = infe_cases.t();

        // 当前时间（单位与 dt 参数一致）
        double current_time = results(ts, 0);

        // 季节性效果
        double seasonal_force = (1.0 + beta_seasonal * cos(dt_seasonal * current_time - beta_phi));
        vec beta = beta_value * seasonal_force;
        if (NPI)
        {
            beta = beta * (1 - (NPI_value * exp(-decay_coef * current_time)));
        }
        // 如果病毒4为延时启动且当前时间未达到启动时间，则禁用其传播（将其 beta 设为 0）
        if (virus4_delayed && current_time < virus4_start_time)
        {
            beta(3) = 0.0;
        }

        // 根据 Lotka-Volterra 模型的思想计算易感性修正值 S_modify
        vec S_modify(NUM_VIRUSES, fill::zeros);
        for (size_t i = 0; i < NUM_VIRUSES; i++)
        {
            double competition_stress = 0.0;
            for (size_t j = 0; j < NUM_VIRUSES; j++)
            {
                if (i != j)
                {
                    competition_stress += (comp_value[j] - comp_value[i]) * infe_cases[j] * N_inv;
                }
            }
            S_modify(i) = sus_cases(i) * exp(-competition_stress);
        }

        // 计算每种病毒的“可被感染”概率
        vec prob_influence(NUM_VIRUSES, fill::zeros);
        for (size_t i = 0; i < NUM_VIRUSES; i++)
        {
            if (sus_cases(i) > 0)
                prob_influence(i) = S_modify(i) / sus_cases(i);
            else
                prob_influence(i) = 0.0;
        }

        // 生成随机数矩阵，判断哪些个体处于可感染状态
        mat rand_influence = randu<mat>(N, NUM_VIRUSES);
        umat Influence_flag(N, NUM_VIRUSES, fill::zeros);
#pragma omp parallel for num_threads(ncores) schedule(static)
        for (size_t n = 0; n < N; n++)
        {
            for (size_t i = 0; i < NUM_VIRUSES; i++)
            {
                if (agent_status(n, i) == 0)
                {
                    double r = rand_influence(n, i);
                    if (r < prob_influence(i))
                        Influence_flag(n, i) = 1;
                    else
                        Influence_flag(n, i) = 0;
                }
                else
                {
                    Influence_flag(n, i) = 1;
                }
            }
        }

// 并行计算状态转换概率（S->I, I->R, R->S）
#pragma omp parallel for num_threads(ncores) schedule(static)
        for (size_t n = 0; n < N; n++)
        {
            // 计算当前个体多重感染的惩罚因子
            double PenalNum = sum(agent_status.row(n) == 1);
            double penal_factor = penal_lookup[PenalNum];

            for (size_t i = 0; i < NUM_VIRUSES; i++)
            {
                double lambda = beta(i) * infe_cases(i) * N_inv * penal_factor;
                if (Influence_flag(n, i) == 0 && agent_status(n, i) == 0)
                {
                    State_probability(n, i) = 0.0;
                }
                else
                {
                    switch (static_cast<int>(agent_status(n, i)))
                    {
                    case 0: // S -> I
                        State_probability(n, i) = 1 - exp(-lambda * dt);
                        break;
                    case 1: // I -> R
                        State_probability(n, i) = gamma_prob(i);
                        break;
                    case -1: // R -> S
                        State_probability(n, i) = omega_prob(i);
                        break;
                    default:
                        State_probability(n, i) = 0.0;
                        break;
                    }
                }
            }
        }

        // 状态转换操作
        mat rand_draw = randu<mat>(N, NUM_VIRUSES);
        umat Comparing_matrix = (State_probability > rand_draw);
        umat condition_S = (agent_status == 0) % Comparing_matrix;
        umat condition_I = (agent_status == 1) % Comparing_matrix;
        umat condition_R = (agent_status == -1) % Comparing_matrix;
        agent_status.elem(find(condition_S)).fill(1);
        agent_status.elem(find(condition_I)).fill(-1);
        agent_status.elem(find(condition_R)).fill(0);

        // 增加新病例
        for (size_t col = 0; col < NUM_VIRUSES; col++)
        {
            // 如果是病毒4且为延时启动，且当前时间未到启动时间，则跳过添加新病例
            if (col == 3 && virus4_delayed && current_time < virus4_start_time)
                continue;
            uvec s_index = find(agent_status.col(col) == 0);
            size_t n_sus = s_index.n_elem;
            if (n_sus >= added_cases)
            {
                s_index = shuffle(s_index);
                for (size_t k = 0; k < added_cases; k++)
                {
                    agent_status(s_index(k), col) = 1;
                }
            }
        }
    }
    return results;
}