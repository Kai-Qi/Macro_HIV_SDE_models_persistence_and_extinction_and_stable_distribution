# __author__ = "Kai Qi (QiaoJianYe)"
# __email___ = "upcqikai@hotmail.com"
# __time____ = "2020.09.15"
'''
Using the Milstein’s higher order method in paper

D.J. Higham , An algorithmic introduction to numerical simulation of stochastic differential equations, SIAM Rev. 43 (3) (2001) 525–546

to provide examples and numerical simulations to verify  results of a HIV Stochastic Differnetial Equcations in paper

Kai Qi and Daqing Jiang, The impact of virus carrier screening and actively seeking treatment on dynamical behavior of a stochastic HIV/AIDS infection model, Applied Mathematical Modelling, 85 (2020) 378-404.
(see https://www.sciencedirect.com/science/article/pii/S0307904X2030161X?via%3Dihub)

'''
import random
import numpy as np       
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


def Marco_HIV_SDE_Stable_Distribution():
    Pi = 25000
    k = 1.2
    beta1 = 0.8*pow(10,-5) 
    beta2 = 0.6*pow(10,-5) 
    beta3 = 0.4*pow(10,-5) 
    beta4 = 0.7*pow(10,-5) 
    mu = 0.015; sigma = 0.2 ; rho1 = 0.001; gamma1 = 0.2; rho2 = 0.09
    psi = 0.2; 
    rho3 = 0.4; gamma2 = 0.15; rho4 = 0.2; gamma3 = 0.3
    delta1 = 0.33; delta2 = 0.2
    sigma_1 = 0.006** 0.5; sigma_2 = 0.006** 0.5; sigma_3 = 0.006** 0.5; sigma_4 = 0.006** 0.5 
    sigma_5 = 0.006** 0.5; sigma_6 = 0.006** 0.5; sigma_7 = 0.006** 0.5
    dt = np.array([0.01 ])

    random.seed( 10 )
    N=200000

    S = [350000]
    I = [6000]
    C = [7500]
    C_s = [8500]
    T = [12000]
    A = [8000]
    A_t = [7500]


    for t in range(0, 2*N - 1):

        Nor1 = np.random.randn(1); Nor2 = np.random.randn(1); Nor3 = np.random.randn(1)
        Nor4 = np.random.randn(1); Nor5 = np.random.randn(1)
        Nor6 = np.random.randn(1); Nor7 = np.random.randn(1)
                
        Tem_Var_1 = S[t] + [Pi - k * (beta1 * I[t] + beta2 * C[t] + beta3 * T[t] + beta4 * A[t]) * S[t] - mu * S[t]] * dt + sigma_1 * S[t] * pow(dt, 0.5) * Nor1 + 0.5 * (sigma_1 ** 2) * S[t] * pow(dt, 0.5) * (Nor1 ** 2 - 1)
        S.append(Tem_Var_1[0] )
        Tem_Var_2 = I[t] + [k * (beta1 * I[t] + beta2 * C[t] + beta3 * T[t] + beta4 * A[t]) * S[t] - (mu + sigma + rho1 + gamma1) * I[t]] * dt + sigma_2 * I[t] * pow(dt, 0.5) * Nor2 + 0.5 * (sigma_2 ** 2) * I[t] * pow(dt, 0.5) * (Nor2 ** 2 - 1)
        I.append(Tem_Var_2[0] )
        Tem_Var_3 = C[t] + [sigma * I[t] - (mu + rho2 + psi) * C[t]] * dt + sigma_3 * C[t] * pow(dt, 0.5) * Nor3 + 0.5 * (sigma_3 ** 2) * C[t] * pow(dt, 0.5) * (Nor3 ** 2 - 1)
        C.append(Tem_Var_3[0] )
        Tem_Var_4 = C_s[t] + [psi * C[t] - (mu + rho3 + gamma2) * C_s[t]] * dt + sigma_4 * C_s[t] * pow(dt, 0.5) * Nor4 + 0.5 * (sigma_4 ** 2) * C_s[t] * pow(dt, 0.5) * (Nor4 ** 2 - 1)
        C_s.append(Tem_Var_4[0]) 
        Tem_Var_5 = T[t] + [gamma1 * I[t] + gamma2 * C_s[t] - (mu + rho4) * T[t]] * dt + sigma_5 * T[t]* pow(dt, 0.5) * Nor5 + 0.5 * (sigma_5 ** 2) * T[t] * pow(dt, 0.5) * (Nor5 ** 2 - 1)
        T.append(Tem_Var_5[0])
        Tem_Var_6 = A[t] + [rho1 * I[t] + rho2 * C[t] + rho3 * C_s[t] + rho4 * T[t] - (mu + gamma3 + delta1) * A[t]] * dt + sigma_6 * A[t] * pow(dt, 0.5) * Nor6 + 0.5 * (sigma_6 ** 2) * A[t] * pow(dt, 0.5) * (Nor6 ** 2 - 1)
        A.append(Tem_Var_6[0])
        Tem_Var_7 = A_t[t] + [gamma3 * A[t] - (mu + delta2) * A_t[t]] * dt + sigma_7 * A_t[t] * pow(dt, 0.5) * Nor7 + 0.5 * (sigma_7 ** 2) * A_t[t] * pow(dt, 0.5) * (Nor7 ** 2 - 1)
        A_t.append(Tem_Var_7[0])

                
    for t in range(2*N-1, 3*N - 1):

        Nor1 = np.random.randn(1); Nor2 = np.random.randn(1); Nor3 = np.random.randn(1)
        Nor4 = np.random.randn(1); Nor5 = np.random.randn(1)
        Nor6 = np.random.randn(1); Nor7 = np.random.randn(1)
                
        Tem_Var_1 = S[t] + [Pi - k * (beta1 * I[t] + beta2 * C[t] + beta3 * T[t] + beta4 * A[t]) * S[t] - mu * S[t]] * dt + sigma_1 * S[t] * pow(dt, 0.5) * Nor1 + 0.5 * (sigma_1 ** 2) * S[t] * pow(dt, 0.5) * (Nor1 ** 2 - 1)
        S.append(Tem_Var_1[0] )
        Tem_Var_2 = I[t] + [k * (beta1 * I[t] + beta2 * C[t] + beta3 * T[t] + beta4 * A[t]) * S[t] - (mu + sigma + rho1 + gamma1) * I[t]] * dt + sigma_2 * I[t] * pow(dt, 0.5) * Nor2 + 0.5 * (sigma_2 ** 2) * I[t] * pow(dt, 0.5) * (Nor2 ** 2 - 1)
        I.append(Tem_Var_2[0] )
        Tem_Var_3 = C[t] + [sigma * I[t] - (mu + rho2 + psi) * C[t]] * dt + sigma_3 * C[t] * pow(dt, 0.5) * Nor3 + 0.5 * (sigma_3 ** 2) * C[t] * pow(dt, 0.5) * (Nor3 ** 2 - 1)
        C.append(Tem_Var_3[0] )
        Tem_Var_4 = C_s[t] + [psi * C[t] - (mu + rho3 + gamma2) * C_s[t]] * dt + sigma_4 * C_s[t] * pow(dt, 0.5) * Nor4 + 0.5 * (sigma_4 ** 2) * C_s[t] * pow(dt, 0.5) * (Nor4 ** 2 - 1)
        C_s.append(Tem_Var_4[0]) 
        Tem_Var_5 = T[t] + [gamma1 * I[t] + gamma2 * C_s[t] - (mu + rho4) * T[t]] * dt + sigma_5 * T[t]* pow(dt, 0.5) * Nor5 + 0.5 * (sigma_5 ** 2) * T[t] * pow(dt, 0.5) * (Nor5 ** 2 - 1)
        T.append(Tem_Var_5[0])
        Tem_Var_6 = A[t] + [rho1 * I[t] + rho2 * C[t] + rho3 * C_s[t] + rho4 * T[t] - (mu + gamma3 + delta1) * A[t]] * dt + sigma_6 * A[t] * pow(dt, 0.5) * Nor6 + 0.5 * (sigma_6 ** 2) * A[t] * pow(dt, 0.5) * (Nor6 ** 2 - 1)
        A.append(Tem_Var_6[0])
        Tem_Var_7 = A_t[t] + [gamma3 * A[t] - (mu + delta2) * A_t[t]] * dt + sigma_7 * A_t[t] * pow(dt, 0.5) * Nor7 + 0.5 * (sigma_7 ** 2) * A_t[t] * pow(dt, 0.5) * (Nor7 ** 2 - 1)
        A_t.append(Tem_Var_7[0])



    for t in range(3*N-1, 4*N - 1):

        Nor1 = np.random.randn(1); Nor2 = np.random.randn(1); Nor3 = np.random.randn(1)
        Nor4 = np.random.randn(1); Nor5 = np.random.randn(1)
        Nor6 = np.random.randn(1); Nor7 = np.random.randn(1)
                
        Tem_Var_1 = S[t] + [Pi - k * (beta1 * I[t] + beta2 * C[t] + beta3 * T[t] + beta4 * A[t]) * S[t] - mu * S[t]] * dt + sigma_1 * S[t] * pow(dt, 0.5) * Nor1 + 0.5 * (sigma_1 ** 2) * S[t] * pow(dt, 0.5) * (Nor1 ** 2 - 1)
        S.append(Tem_Var_1[0] )
        Tem_Var_2 = I[t] + [k * (beta1 * I[t] + beta2 * C[t] + beta3 * T[t] + beta4 * A[t]) * S[t] - (mu + sigma + rho1 + gamma1) * I[t]] * dt + sigma_2 * I[t] * pow(dt, 0.5) * Nor2 + 0.5 * (sigma_2 ** 2) * I[t] * pow(dt, 0.5) * (Nor2 ** 2 - 1)
        I.append(Tem_Var_2[0] )
        Tem_Var_3 = C[t] + [sigma * I[t] - (mu + rho2 + psi) * C[t]] * dt + sigma_3 * C[t] * pow(dt, 0.5) * Nor3 + 0.5 * (sigma_3 ** 2) * C[t] * pow(dt, 0.5) * (Nor3 ** 2 - 1)
        C.append(Tem_Var_3[0] )
        Tem_Var_4 = C_s[t] + [psi * C[t] - (mu + rho3 + gamma2) * C_s[t]] * dt + sigma_4 * C_s[t] * pow(dt, 0.5) * Nor4 + 0.5 * (sigma_4 ** 2) * C_s[t] * pow(dt, 0.5) * (Nor4 ** 2 - 1)
        C_s.append(Tem_Var_4[0]) 
        Tem_Var_5 = T[t] + [gamma1 * I[t] + gamma2 * C_s[t] - (mu + rho4) * T[t]] * dt + sigma_5 * T[t]* pow(dt, 0.5) * Nor5 + 0.5 * (sigma_5 ** 2) * T[t] * pow(dt, 0.5) * (Nor5 ** 2 - 1)
        T.append(Tem_Var_5[0])
        Tem_Var_6 = A[t] + [rho1 * I[t] + rho2 * C[t] + rho3 * C_s[t] + rho4 * T[t] - (mu + gamma3 + delta1) * A[t]] * dt + sigma_6 * A[t] * pow(dt, 0.5) * Nor6 + 0.5 * (sigma_6 ** 2) * A[t] * pow(dt, 0.5) * (Nor6 ** 2 - 1)
        A.append(Tem_Var_6[0])
        Tem_Var_7 = A_t[t] + [gamma3 * A[t] - (mu + delta2) * A_t[t]] * dt + sigma_7 * A_t[t] * pow(dt, 0.5) * Nor7 + 0.5 * (sigma_7 ** 2) * A_t[t] * pow(dt, 0.5) * (Nor7 ** 2 - 1)
        A_t.append(Tem_Var_7[0])

        
    plt.subplot(231)
    sns.kdeplot(S[N:2 * N - 1]); sns.kdeplot(S[2 * N - 1:3 * N - 1]); sns.kdeplot(S[3 * N - 1:4 * N - 1]); plt.xlabel("The pdf of S")

    plt.subplot(232)
    sns.kdeplot(I[N:2 * N - 1]); sns.kdeplot(I[2 * N - 1:3 * N - 1]); sns.kdeplot(I[3 * N - 1:4 * N - 1]); plt.xlabel("The pdf of I")

    plt.subplot(233)
    sns.kdeplot(C[N:2 * N - 1]); sns.kdeplot(C[2 * N - 1:3 * N - 1]); sns.kdeplot(C[3 * N - 1:4 * N - 1]); plt.xlabel("The pdf of C")

    plt.subplot(234)
    sns.kdeplot(C_s[N:2 * N - 1]); sns.kdeplot(C_s[2 * N - 1:3 * N - 1]); sns.kdeplot(C_s[3 * N - 1:4 * N - 1]); plt.xlabel("The pdf of C_s")

    plt.subplot(235)
    sns.kdeplot(T[N:2 * N - 1]); sns.kdeplot(T[2 * N - 1:3 * N - 1]); sns.kdeplot(T[3 * N - 1:4 * N - 1]); plt.xlabel("The pdf of T")

    plt.subplot(236)
    sns.kdeplot(A[N:2 * N - 1]); sns.kdeplot(A[2 * N - 1:3 * N - 1]); sns.kdeplot(A[3 * N - 1:4 * N - 1]); plt.xlabel("The pdf of A")

    plt.show()



if __name__ == '__main__':
    Marco_HIV_SDE_Stable_Distribution()