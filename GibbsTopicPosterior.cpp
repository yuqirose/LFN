/*
 * GibbsSamplerLFN.cpp
 *
 *  Created on: Apr 26, 2014
 *      Author: roseyu
 */


#include "mex.h"
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <math.h>

void GibbsSamplerTM(const double* THETA, const double* BETA, const double* Wp, double* tp, double* tpcount, double* F, const double* tau, const double* Gpfeature, const int K, const int C, const double* numToken, const int N, double* probs)
{
    int i,wid,wid0,k,k1,n,n1;
    double denominator, prior, likeli, fWeight, maxprob;
    double sampler;
    double *fWeight0;
    fWeight0 = (double*) mxCalloc(N, sizeof(double));

    int startID, endID;
//     for(n=0;n<N;n++)
//     {
//         startID = n*C;
//         endID = (int) n*C + (int)numToken[n];
//         for(i=startID;i<endID;i++)
//         {
//             printf("user %d, token %d has topic %d\n", n, i-startID+1, (int)*(tp+i));
//         }
//     }
    
    for(n=0;n<N;n++)
    {
//         printf("sampling user %d\n", n);
        startID = n*C;
        endID = (int) n*C + (int)numToken[n];

        for(n1=0;n1<N;n1++)
            fWeight0[n1] = 0;

        for(n1=0;n1<N;n1++)
        {
            if(n1!=n)
            {
                for(k1=0;k1<K;k1++)
                {
                    fWeight0[n1] = fWeight0[n1] + tpcount[n1*K+k1]*tpcount[n*K+k1];
                }
            }
        }

        for(i=startID;i<endID;i++)
        {
            tpcount[n*K+(int)tp[i]] = tpcount[n*K+(int)tp[i]]-1;
            wid = *(Wp+i)-1;

            if(i==startID || wid0!=wid )
            {
                wid0 = wid;
                for(k=0;k<K;k++)
                {
                    prior = THETA[n*K+k]; /* Theta(k,n) */
                    likeli = BETA[k+wid*K]; /* Beta(k, wid) */
                    probs[k] = prior*likeli;
                    probs[k] = log(probs[k]);
                    tpcount[n*K+k] = tpcount[n*K+k]+1;
                    for(n1=0;n1<N;n1++)
                    {
                        if(n!=n1)
                        {
                            fWeight = 0;
                            /* use BLAS to optimize the innder product */
                            fWeight = fWeight0[n1] - tpcount[n1*K+(int)tp[i]] + tpcount[n1*K+k];
                            fWeight = fWeight/(numToken[n1]*numToken[n]);
                            fWeight = 1/(1+exp(-tau[0]*fWeight - tau[1]*Gpfeature[n*N+n1] - tau[2]));
                            probs[k] = probs[k] + F[n*N+n1] * log(fWeight+1e-32) + (1-F[n*N+n1])*log(1-fWeight+1e-32);
                        }
                    }
                    tpcount[n*K+k] = tpcount[n*K+k]-1;
                }
                maxprob = probs[0];
                for(k=1;k<K;k++)
                {
                    if(maxprob<probs[k])
                        maxprob = probs[k];
                }

                denominator = 0;
                for(k=0;k<K;k++)
                {
                    probs[k] = exp(probs[k]-maxprob);
                    denominator = denominator+probs[k];
                }

                for(k=1;k<K;k++)
                {
                    probs[k] = probs[k]+probs[k-1];
                }

                for (k=0;k<K;k++)
                    probs[k] = probs[k]/denominator;

                probs[K-1] = 1;
            }

            sampler = ((double) rand() / (RAND_MAX));
            for (k=0;k<K;k++)
            {
                if(sampler<=probs[k])
                    break;
                // the topic sample is "k"
            }

            for(n1=0;n1<N;n1++)
            {
                if(n1!=n)
                {
                    fWeight0[n1] = fWeight0[n1] - tpcount[n1*K+(int)tp[i]] + tpcount[n1*K+k];
                }
            }
            *(tp+i) = (double) k;
            tpcount[n*K+k] = tpcount[n*K+k]+1;
        }
    }
}
    
    

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    /* attention */
    /* make sure Theta and Beta are processed to remove the 0 elements */
    /* Theta = bsxfun(@rdivide, Theta+1e-32, sum(Theta+1e-32)) */
    /* Beta = bsxfun(@rdivide, Beta+1e-32, sum(Beta+1e-32,2)) */
//     The input includes:
//     1. matrix of words in vector format: [C x N]
//     2. matrix of topics in vector format: [C x N]
//     3. number of tokens from each user: [1 x N] 
//     4. count of topics from each user: [K x N]
//     5. Grouping feature between users: [N x N]
//     6. F, observation of following relationship: [N x N]
//     7. tau, following distribution parameters: [1 X 3]
//     8. THETA, parameters of user-specific topics: [K x N]
//     9. BETA, topic conditioned distributions: [K X V]

    /* process the input variables & parameters */
    
    // input group 1: data
    /* 1-1: words (contents)*/
    double *Wp;
    Wp = mxGetPr(prhs[0]); 
    
    /* 1-2: topics */
    double *Tp, *numToken;
    Tp = mxGetPr(prhs[1]);
    
    /* 1-3: sizes ( number of tokens ) */
    numToken = mxGetPr(prhs[2]);
    
    /* 1-3: topics-statistics */
    double *Tpcount, *numWords;
    Tpcount = mxGetPr(prhs[3]);
    
    /* 1-4: grouping features */
    double *Gpfeature;
    Gpfeature = mxGetPr(prhs[4]);
    
    /* 1-5: connection data and parameters*/
    double *F, *tau;
    F = mxGetPr(prhs[5]);
    tau = mxGetPr(prhs[6]);
    
    // input group 2: parameters
    /* 2-1: topic weight parameters */
    double *THETA;
    THETA = mxGetPr(prhs[7]); 
    
    /* 2-2: word distribution parameters */
    double *BETA;
    BETA = mxGetPr(prhs[8]); // K*V x 1 vector
    
    /* constant values & parameters */
    int K,C,V,N;
    K = mxGetM(prhs[8]);
    V = mxGetN(prhs[8]);
    C = mxGetM(prhs[0]); /* the length of the longest user-doc */
    N = mxGetN(prhs[0]);
    
    
    int i,j, userid, n;
    
    /* set up the output varaibles */
    double *tp, *tpcount, *probs;
    plhs[0] =  mxCreateDoubleMatrix(C*N, 1, mxREAL);
    tp = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(K*N, 1, mxREAL);
    tpcount = mxGetPr(plhs[1]);
    
    probs = (double*) mxCalloc(K, sizeof(double));
    /*probs = new double[5];*/
    
    for (i=0;i<C*N;i++)
    {
        if(Wp[i]>0)
        {
            Wp[i] = Wp[i];
            Tp[i] = Tp[i];
            *(tp+i) = Tp[i]-1;
        }else{
            *(tp+i) = -1; /* trick needs attention if modifying code */
        }
    }
    
    for (i=0;i<K*N;i++)
    {
        *(tpcount+i) = Tpcount[i];
    }
    
    GibbsSamplerTM( THETA, BETA, Wp, tp, tpcount, F, tau, Gpfeature, K, C, numToken, N, probs);
    
//     
//     for (i=0;i<C;i=i+50)
//         printf("returned word %d topic is %f\n", i, *(tp+i));
    
    
}




