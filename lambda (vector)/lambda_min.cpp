#include <stdio.h>
#include <stdlib.h>
#include <math.h>





int main(void)
{
    
    const int n=1000;
    
    
    
    FILE *fp;
    FILE *fp_2;
    FILE *fp_t;
    fp=fopen("lambda_min.txt","w");
    //fp_2=fopen("sequence2.txt","w");
    //fp_t=fopen("sequencet.txt","w");
    int i,j,k;
    
    double delta=10.0/n;
    double lambda_1=1.0;
    double lambda_2=1.0;
    
    double lambda_max[2000], lambda_min[2000];
    
    double alfa=1.0;
    
    double m_1=1.0;
    double m_2=1.0;
    
    
    double bot_1[n],mid_1[n],top_1[n],diag_1[n];
    double bot_2[n],mid_2[n],top_2[n],diag_2[n];
    double bot_3[n],mid_3[n],top_3[n],diag_3[n];
    
    double norm_1=2.0,norm_p_1=2.0;
    double norm_2=2.0,norm_p_2=2.0;
    
    
    double col_1[n], psi_1[n],psi_p_1[n];
    double col_2[n], psi_2[n],psi_p_2[n];
    double col_3[n], theta[n], theta_p[n];
    
    double S1=0.0,S2=0.0;
    
    
    
    
    
    k=0;
    
    
    
    
    for(lambda_1=0.1; lambda_1<=5.0; lambda_1+=0.01)
    
    {
                      
    norm_2=2.0;
    
    lambda_2=lambda_1;
    while(norm_2>0.5)
    {
    
    lambda_2=lambda_2-0.01;
    
    for (i=0; i<n; i++)
        {    
        psi_1[i]=pow((i*delta),m_1)*pow(2.71, -i*delta*(i*delta));
        psi_2[i]=pow((i*delta),m_2)*pow(2.71, -i*delta*(i*delta));
        theta[i]=pow(2.71, -i*delta*i*delta);
        }
    
    S2=0.0;
    S1=0.0;
    
    for (i=1; i<n; i++)
        {
        bot_1[i]=-1.0/(delta*delta)+1.0/(2.0*i*delta*delta);
        mid_1[i]=lambda_1+2/(delta*delta)+m_1*m_1/(i*i*delta*delta);
        top_1[i]=-1.0/(delta*delta)-1.0/(2.0*i*delta*delta);
        
        bot_2[i]=-1.0/(delta*delta)+1.0/(2.0*i*delta*delta);
        mid_2[i]=lambda_2+2/(delta*delta)+m_2*m_2/(i*i*delta*delta);
        top_2[i]=-1.0/(delta*delta)-1.0/(2.0*i*delta*delta);
        
        bot_3[i]=-1.0/(delta*delta)+1.0/(2.0*i*delta*delta);
        top_3[i]=-1.0/(delta*delta)-1.0/(2.0*i*delta*delta);
        mid_3[i]=alfa+2/(delta*delta);
        
        }
    
    bot_1[0]=0.0;
    top_1[0]=-2/(delta*delta);
    
    bot_2[0]=0.0;
    top_2[0]=-2/(delta*delta);
    
    bot_3[0]=0.0;
    top_3[0]=-2/(delta*delta);
    mid_3[0]=alfa+2/(delta*delta); 
    
    
    
    
    
    
    while(((S1-1)*(S1-1))>0.0001)
    {
    norm_p_1=0.0;
    norm_1=0.0;
    norm_p_2=0.0;
    norm_2=0.0;
    
    
    
    for (i=0; i<n; i++)
        {    
         diag_1[i]=mid_1[i];
         
         col_1[i]=psi_1[i]*theta[i];
         
         psi_p_1[i]=psi_1[i];
         
         diag_2[i]=mid_2[i];
         
         col_2[i]=psi_2[i]*theta[i];
         
         psi_p_2[i]=psi_2[i];
         }
        
    for (i=2; i<n; i++)
        {
        diag_1[i]=diag_1[i]-top_1[i-1]*bot_1[i]/diag_1[i-1];
        col_1[i]=col_1[i]-col_1[i-1]*bot_1[i]/diag_1[i-1];
        
        diag_2[i]=diag_2[i]-top_2[i-1]*bot_2[i]/diag_2[i-1];
        col_2[i]=col_2[i]-col_2[i-1]*bot_2[i]/diag_2[i-1];
        }
    
    
    psi_1[n-1]=col_1[n-1]/diag_1[n-1];
    psi_2[n-1]=col_2[n-1]/diag_2[n-1];
    
    for (i=n-2; i>=1; i--)
        {
        col_1[i]=col_1[i]-col_1[i+1]*top_1[i]/diag_1[i+1];
        psi_1[i]=col_1[i]/diag_1[i];
        
        col_2[i]=col_2[i]-col_2[i+1]*top_2[i]/diag_2[i+1];
        psi_2[i]=col_2[i]/diag_2[i];
        }
    psi_1[0]=0.0;
    
    
    for (i=1; i<n; i++)
        {
        norm_p_1+=(psi_p_1[i]*psi_p_1[i]*i*delta*delta);
        norm_1+=(psi_1[i]*psi_1[i]*i*delta*delta);
        
        norm_p_2+=(psi_p_2[i]*psi_p_2[i]*i*delta*delta);
        norm_2+=(psi_2[i]*psi_2[i]*i*delta*delta);
        }
    
    
    S1=pow((norm_p_1/norm_1),0.75);
    S2=pow((norm_p_2/norm_2),0.75);
    
    
    for (i=0; i<n; i++)
        {
        psi_1[i]*=S1;
        psi_2[i]*=S2;
        }
    
    for (i=0; i<n; i++)
        {
        diag_3[i]=mid_3[i];
        col_3[i]=psi_1[i]*psi_1[i]+psi_2[i]*psi_2[i];
        theta_p[i]=theta[i];
        }
    
    for (i=1; i<n; i++)
        {
        diag_3[i]=diag_3[i]-top_3[i-1]*bot_3[i]/diag_3[i-1];
        col_3[i]=col_3[i]-col_3[i-1]*bot_3[i]/diag_3[i-1];
        }
    
    theta[n-1]=col_3[n-1]/diag_3[n-1];
    
    for (i=n-2; i>=0; i--)
        {
        col_3[i]=col_3[i]-col_3[i+1]*top_3[i]/diag_3[i+1];
        theta[i]=col_3[i]/diag_3[i];
        }
    
    
    
    
    //printf("\n%lf %lf\n",S1,S2);
    
    }
    //printf("%lf %lf\n",lambda_2, norm_2);
    //system("PAUSE");
    }
    printf("%lf %lf\n", lambda_1, lambda_2);
    fprintf(fp,"%lf %lf\n", lambda_1, lambda_2);
    }
    //printf("%lf %lf\n", norm_1,norm_2);
    
    /*for (i=0; i<n; i++)
        fprintf(fp_1, "%lf %lf\n", (i*delta), psi_1[i]);
        
    for (i=0; i<n; i++)
        fprintf(fp_2, "%lf %lf\n", (i*delta), psi_2[i]);
    for (i=0; i<n; i++)
        fprintf(fp_t, "%lf %lf\n", (i*delta), theta[i]);*/
        
        
        
    fclose(fp);
    
    system("PAUSE");
    return 0;
}
    
    
