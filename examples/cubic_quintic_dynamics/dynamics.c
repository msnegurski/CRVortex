#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <FFT_functions.h>


double double_abs(double x)
{
    if (x>=0.0)
        return x;
    else
        return -x;


}

double fract_double(double x)
{
    while (x>=1)
        x-=1;
    return x;
}

int int_double(double x)
{
    int k=0;
    while (x>=1)
    {
        x-=1;
        k+=1;
    }
    return k;


}

double psi_sq(double **Rdat, double **Idat, int i, int j)
{
    return (Rdat[i][j]*Rdat[i][j]+Idat[i][j]*Idat[i][j]);
}

double intence(double **Rdat,double **Idat,int i,int j)
{
    return sqrt(psi_sq(Rdat,Idat,i,j));
}


double gradx_sq(double **Rdat, double **Idat, int x, int y, double dx, int n, int logn)
{
    double **ReBuf,**ImBuf;
    int i,j;
    double p=2*3.141592653589/n;

    double t=p/dx;
    
    ReBuf = (double **) malloc (n*sizeof(double *));
    ImBuf = (double **) malloc (n*sizeof(double *));

    for (i=0; i<n; i++)
    {
        ReBuf[i] = (double *) malloc (n*sizeof(double));
        ImBuf[i] = (double *) malloc (n*sizeof(double));
    }
    
    
    for (i=0;i<n;i++)
    for (j=0;j<n;j++)
    {
        ReBuf[i][j]=Rdat[i][j];
        ImBuf[i][j]=Idat[i][j];
    }
    double buf;
    
    FFTX(ReBuf,ImBuf,n,logn,-1,y);
    
    for (i=0; i<n/2; i++)
    {
        buf=ReBuf[i][y];
        ReBuf[i][y]=t*i*ImBuf[i][y];
        ImBuf[i][y]=(-t*i)*buf;
    }
    
    for (i=n-1; i>=(n/2); i--)
    {
        buf=ReBuf[i][y];
        ReBuf[i][y]=-t*(n-i)*ImBuf[i][y];
        ImBuf[i][y]=(t*(n-i))*buf;
    }
  	
    FFTX(ReBuf,ImBuf,n,logn,1,y);
    
    buf=(ReBuf[x][y]*ReBuf[x][y]+ImBuf[x][y]*ImBuf[x][y]);
      
    for (i=0; i<n; i++)
    {
        free(ReBuf[i]);
        free(ImBuf[i]);
    }
    free(ReBuf);
    free(ImBuf);
    
    
    return buf;
    

}

double grady_sq(double **Rdat, double **Idat, int x, int y, double dx, int n, int logn)
{
    double **ReBuf,**ImBuf;
    int i,j;
    double p=2*3.141592653589/n;

    double t=p/dx;

    ReBuf = (double **) malloc (n*sizeof(double *));
    ImBuf = (double **) malloc (n*sizeof(double *));

    for (i=0; i<n; i++)
    {
        ReBuf[i] = (double *) malloc (n*sizeof(double));
        ImBuf[i] = (double *) malloc (n*sizeof(double));
    }

    for (i=0;i<n;i++)
    for (j=0;j<n;j++)
    {
        ReBuf[i][j]=Rdat[i][j];
        ImBuf[i][j]=Idat[i][j];
    }
    

    double buf;

    FFTY(ReBuf,ImBuf,n,logn,-1,x);

    for (i=0; i<n/2; i++)
    {
        buf=ReBuf[x][i];
        ReBuf[x][i]=t*i*ImBuf[x][i];
        ImBuf[x][i]=(-t*i)*buf;
    }
    for (i=n-1; i>=(n/2); i--)
    {
        buf=ReBuf[i][y];
        ReBuf[x][i]=(-t*(n-i))*ImBuf[x][i];
        ImBuf[x][i]=(t*(n-i))*buf;
    }

    FFTY(ReBuf,ImBuf,n,logn,1,x);
    
    buf=(ReBuf[x][y]*ReBuf[x][y]+ImBuf[x][y]*ImBuf[x][y]);
      
    for (i=0; i<n; i++)
    {
        free(ReBuf[i]);
        free(ImBuf[i]);
    }
    free(ReBuf);
    free(ImBuf);
    
    
    return buf;
}


double integrate_grad_xy_sq(double **Rdat, double **Idat, double dx, int n, int logn)
{
    double **ReBuf,**ImBuf;
    int i,j;
    double p=2*3.141592653589/n;
    double ky,kx;

    double t=p*p/dx/dx;

    ReBuf = (double **) malloc (n*sizeof(double *));
    ImBuf = (double **) malloc (n*sizeof(double *));

    for (i=0; i<n; i++)
    {
        ReBuf[i] = (double *) malloc (n*sizeof(double));
        ImBuf[i] = (double *) malloc (n*sizeof(double));
    }

    for (i=0;i<n;i++)
    for (j=0;j<n;j++)
    {
        ReBuf[i][j]=Rdat[i][j];
        ImBuf[i][j]=Idat[i][j];
    }
    
    double res=0.0;

    double buf; 
    
    FFT2D(ReBuf, ImBuf, n, logn, -1);


    for (i=0; i<n; i++)
        for (j=0; j<n; j++)
        {

            if (i<n/2)
                kx=i;
            else
                kx=(n-i);

            if (j<n/2)
                ky=j;
            else
                ky=(n-j);


            buf=-(kx*kx+ky*ky)*t;

            ReBuf[i][j]*=buf;
            ImBuf[i][j]*=buf;	
        }

    FFT2D(ReBuf, ImBuf, n, logn, 1);
    
    for (i=0; i<n; i++)
    for (j=0; j<n; j++)
    {
        res-=(Rdat[i][j]*ReBuf[i][j]+Idat[i][j]*ImBuf[i][j])*dx*dx;
    
    }
    
    return res;
       
       
}

double power_x(double **Rdat, double **Idat, int N, double dx, int j)
{
    double p=0.0;
    int i;
    for (i=1; i<N-1; i+=2)
        {
            //p+=psi_sq(Rdat,Idat,i,j);
            /*p+=psi_sq(Rdat,Idat,i,j);
            p+=(psi_sq(Rdat,Idat,i+1,j)+psi_sq(Rdat,Idat,i-1,j)+psi_sq(Rdat,Idat,i,j+1)+psi_sq(Rdat,Idat,i,j-1)-4*psi_sq(Rdat,Idat,i,j))*0.1666666666666666;
            p+=0.25*(psi_sq(Rdat,Idat,i+1,j)+psi_sq(Rdat,Idat,i,j+1)-psi_sq(Rdat,Idat,i-1,j)-psi_sq(Rdat,Idat,i,j-1));*/
            //p+=0.5*(psi_sq(Rdat,Idat,i+1,j)+psi_sq(Rdat,Idat,i,j+1));
            
            
            p+=0.3333333333333333333333333*(psi_sq(Rdat,Idat,i+1,j)+psi_sq(Rdat,Idat,i-1,j)+4*psi_sq(Rdat,Idat,i,j));
            
        }
    p*=dx;
    return p; 
       
}



double power(double **Rdat, double **Idat, int N, double dx)
{
    double p=0.0;
    int i;
    for (i=1; i<N-1; i+=2)
        {
            /*p+=0.16666666666666*(power_x(Rdat,Idat,N,dx,i+1)+power_x(Rdat,Idat,N,dx,i-1)-2*power_x(Rdat,Idat,N,dx,i));
            p+=0.25*(power_x(Rdat,Idat,N,dx,i+1)-power_x(Rdat,Idat,N,dx,i-1));
            p+=power_x(Rdat,Idat,N,dx,i);*/
            p+=0.33333333333333333333333333*(power_x(Rdat,Idat,N,dx,i+1)+power_x(Rdat,Idat,N,dx,i-1)+4*power_x(Rdat,Idat,N,dx,i));
        }
    p*=dx;
    return p;
}  



double energy(double **Rdat, double **Idat, int N, double dx, int LogN)
{
    double e=0.0;
    int i,j;
    for (i=1; i<N; i++)
        for (j=1; j<N; j++)
            e+=0.5*psi_sq(Rdat,Idat,i,j)*psi_sq(Rdat,Idat,i,j)-0.33333333333333333*psi_sq(Rdat,Idat,i,j)*psi_sq(Rdat,Idat,i,j)*psi_sq(Rdat,Idat,i,j);

    e*=dx*dx;
    
    //e+=integrate_grad_xy_sq(Rdat,Idat,dx,N,LogN);

    return e;
}

void find_stat(double **Rdat, double **Idat, int cnum, double dx, double size)
{
    const int n=10*cnum;

    int i,j,k;

    double delta=size/n/2;
    double lambda=0.85;

    FILE *fp;
    fp=fopen("sequence.txt","w");

    double a[n],b[n],c[n],diag[n];

    double v1[n], v2[n];

    double Nn=0.0,Np=0.0;

    double d[n], psi[n],psip[n];
    double S;

    int fl=0;

    for (i=0; i<n; i++)
    {    
        psi[i]=pow(2.71, -i*delta*(i*delta));
    }

    for (i=0; i<n; i++)
    {
        a[i]=-1.0/(delta*delta)+1.0/(2.0*i*delta*delta);
        c[i]=lambda+2/(delta*delta);//+m*m/(i*i*delta*delta);
        b[i]=-1.0/(delta*delta)-1.0/(2.0*i*delta*delta);
    }

    a[0]=0.0;
    b[0]=-2/(delta*delta);



    while(((S-1)*(S-1))>0.000000001)
    {
        Np=0.0;
        Nn=0.0;


        fl++;


        for (i=0; i<n; i++)
        {    
            diag[i]=c[i];

            d[i]=-2*psi[i]*psi[i]*psi[i]+4*psi[i]*psi[i]*psi[i]*psi[i]*psi[i];

            psip[i]=psi[i];
        }



        for (i=1; i<n; i++)
        {
            diag[i]=diag[i]-b[i-1]*a[i]/diag[i-1];
            d[i]=d[i]-d[i-1]*a[i]/diag[i-1];
        }

        psi[n-1]=d[n-1]/diag[n-1];

        for (i=n-2; i>=0; i--)
        {
            d[i]=d[i]-d[i+1]*b[i]/diag[i+1];
            psi[i]=d[i]/diag[i];
        }



        for (i=0; i<n; i++)
        {
            Np+=(psip[i]*psip[i]*i*delta*delta);
            Nn+=(psi[i]*psi[i]*i*delta*delta);
        }


        S=pow((Np/Nn),0.6);


        for (i=0; i<n; i++)
        {

            psi[i]*=S;

        }

        //printf("%lf\n",S);

    }
    printf("%lf\n",psi[0]);


    for (i=0; i<n; i++)
        fprintf(fp,"%lf %lf\n", i*delta, psi[i]);


    int r;
    double buf;

    for (i=0;i<cnum;i++)
        for (j=0;j<cnum;j++)
        {
            
            buf=sqrt((i-cnum/2)*(i-cnum/2)+(j-cnum/2)*(j-cnum/2)*1.0)*20;

            r=int_double(buf);

            if ((i==cnum/2) && (j==cnum/2))
                printf("%d\n",r);
            if ((r*delta*2)>size)
                Rdat[i][j]=0.0;
            else
                Rdat[i][j]=psi[r]; 
            Idat[i][j]=Rdat[i][j];
        }


}




void cubic_evol(double **Rdat, double **Idat, int N, double dt)
{
    int i,j;
    double buf1,buf2,buf3;

    for (i=0; i<N; i++)
        for (j=0; j<N; j++)
        {
            buf1=Rdat[i][j];
            buf2=Idat[i][j];
            buf3=(psi_sq(Rdat,Idat,i,j)-psi_sq(Rdat,Idat,i,j)*psi_sq(Rdat,Idat,i,j))*dt;
            //printf("%lf\n",buf3);
            //system("PAUSE");
            Rdat[i][j]=buf1*cos(buf3)-buf2*sin(buf3);
            Idat[i][j]=buf1*sin(buf3)+buf2*cos(buf3);
        }




}


void linear_evol(double **Rdat, double **Idat, int N, int LogN, double dt, double dx)
{

    int i,j;
    double kx,ky;
    double p=2*3.141592653589/N;

    double t=p*p/dx/dx;

    double buf1,buf2,buf3;


    FFT2D(Rdat, Idat, N, LogN, -1);


    for (i=0; i<N; i++)
        for (j=0; j<N; j++)
        {

            if (i<N/2)
                kx=i;
            else
                kx=-(N-i);

            if (j<N/2)
                ky=j;
            else
                ky=-(N-j);


            buf1=-(kx*kx+ky*ky)*t*dt;
            buf2=Rdat[i][j];
            buf3=Idat[i][j];

            Rdat[i][j]=buf2*cos(buf1)-buf3*sin(buf1);
            Idat[i][j]=buf2*sin(buf1)+buf3*cos(buf1);	
        }

    FFT2D(Rdat, Idat, N, LogN, 1);

}


void total_evol(double **Rdat, double **Idat, int N, int LogN, double dt, double dx,int num_step)
{
    int i;

    //linear_evol(Rdat,Idat,N,LogN,dt/2,dx);
    cubic_evol(Rdat,Idat,N,dt/2);

    for (i=0; i<num_step-1; i++)
    {

        linear_evol(Rdat,Idat,N,LogN,dt,dx);
        cubic_evol(Rdat,Idat,N,dt);
        
    }

    linear_evol(Rdat,Idat,N,LogN,dt,dx);
    cubic_evol(Rdat,Idat,N,dt/2);
    



}



double average_width(double **Rdat, double **Idat, int N, double dx)
{
    double s=0.0;
    double p=0.0;
    int i,j;
    for (i=0;i<N;i++)
        for (j=0;j<N;j++)
        {
            s+=((i-N/2)*(i-N/2)+(j-N/2)*(j-N/2))*dx*dx*psi_sq(Rdat,Idat,i,j)*dx*dx;
            p+=psi_sq(Rdat,Idat,i,j)*dx*dx;
        }

    s=s/p;

    return s;
}




int main()
{

    int i,j;
    int deg=8;
    int cnum=1;
    double **RePsi,**ImPsi;
    double size=60;
    double ampl=10.0;
    FILE *final, *start,*ham;
    double dt=0.00000000000001;
    int n=10;
    FILE *x,*y,*z;

    final=fopen("final.txt","w");
    start=fopen("start.txt","w");
    ham=fopen("hamilton.txt","w");


    x=fopen("x.txt","w");
    y=fopen("y.txt","w");
    z=fopen("z.txt","w");

    for (i=0; i<deg; i++)
        cnum*=2;

    double dx=size/cnum;

    RePsi = (double **) malloc (cnum*sizeof(double *));
    ImPsi = (double **) malloc (cnum*sizeof(double *));

    for (i=0; i<cnum; i++)
    {
        RePsi[i] = (double *) malloc (cnum*sizeof(double));
        ImPsi[i] = (double *) malloc (cnum*sizeof(double));
    }

    //----------------------------------------------------------------------------------------



    find_stat(RePsi,ImPsi,cnum,dx,size);
    /*
       for (i=0;i<cnum;i++)
       fprintf(start,"%f %f\n",(i-cnum/2)*dx,sqrt(psi_sq(RePsi,ImPsi,i,cnum/2)));
       i=0;

       printf("%f %f %f\n",i*n*dt,average_width(RePsi,ImPsi,cnum,dx),energy(RePsi,ImPsi,cnum,dx));


       for (i=0;i<1;i++)
       {
       total_evol(RePsi,ImPsi,cnum,deg,dt,dx,n);
       printf("%lf %lf %lf\n", (i+1)*n*dt, average_width(RePsi,ImPsi,cnum,dx),energy(RePsi,ImPsi,cnum,dx));
       }

       for (i=0;i<cnum;i++)
       fprintf(final,"%d %f\n",i,sqrt(psi_sq(RePsi,ImPsi,cnum/2,i)));

*/
    /*
       for (i=0;i<cnum;i++)
       {
       for (j=0;j<cnum;j++)
       fprintf(start,"%f ",sqrt(psi_sq(RePsi,ImPsi,i,j)));
       fprintf(start,"\n");
       }

*/

    /*
       for (i=0;i<cnum;i++)
       fprintf(start,"%f %f\n",(i-cnum/2)*dx,sqrt(psi_sq(RePsi,ImPsi,i,cnum/2)));

       printf("%lf\n",RePsi[cnum-1][cnum-1]);*/

    for (i=0; i<cnum; i++)
    {
        for (j=0; j<cnum; j++)
        {
            fprintf(x,"%lf ",(i-cnum/2)*dx);
            fprintf(x,"%lf ",(j-cnum/2)*dx);
            fprintf(x,"%lf ",psi_sq(RePsi,ImPsi,i,j));
            fprintf(x,"\n");
        }

        //fprintf(y,"\n");
        //fprintf(z,"\n");
    }

    for (i=0;i<cnum;i++)
        fprintf(start,"%f %f\n",(i-cnum/2)*dx,sqrt(psi_sq(RePsi,ImPsi,i,cnum/2)));
    i=0;

    //printf("%f %f %f\n",i*n*dt,average_width(RePsi,ImPsi,cnum,dx),energy(RePsi,ImPsi,cnum,dx,deg));
    printf("%lf %lf %lf\n",i*n*dt,average_width(RePsi,ImPsi,cnum,dx), power(RePsi,ImPsi,cnum,dx));
    
    for (i=0;i<10;i++)
    {
        total_evol(RePsi,ImPsi,cnum,deg,dt,dx,n);
        //printf("%lf %lf %lf\n", (i+1)*n*dt, average_width(RePsi,ImPsi,cnum,dx),energy(RePsi,ImPsi,cnum,dx,deg));
        printf("%lf %lf %lf\n", (i+1)*n*dt, average_width(RePsi,ImPsi,cnum,dx), power(RePsi,ImPsi,cnum,dx));
    }
    for (i=0;i<cnum;i++)
        fprintf(final,"%lf %lf\n",(i-cnum/2)*dx,sqrt(psi_sq(RePsi,ImPsi,i,cnum/2)));
    //fprintf(final,"%f %f\n",(i-cnum/2)*dx,sqrt(RePsi[i][cnum/2]*RePsi[i][cnum/2]));
    //----------------------------------------------------------------------------------------
    fclose(start);
    fclose(final);

    for (i=0; i<cnum; i++)
    {
        free(RePsi[i]);
        free(ImPsi[i]);
    }

    free(RePsi);
    free(ImPsi);
    
    return 0;

}
