#include<stdio.h>
#include<stdlib.h>
#include<float.h>
#include<math.h>
#include<mpi.h>
#define llu unsigned long long

int main(int argc, char*argv[]){
    MPI_Init(&argc,&argv);
    MPI_Status status;
    llu M=atoi(argv[1]);
    int D1=atoi(argv[2]);
    int D2=atoi(argv[3]);
    int T=atoi(argv[4]);
    int seed=atoi(argv[5]);
    
    int myrank, P;
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
    //Allocate buffers:
//data_received      - holds own/initial data
//data_from_send     - temporary receive buffer for incoming message
//data_at_D1/D2      - processed data to send to neighbours (squared / log)
//buffer_updated_*   - transformed buffers used in subsequent iterations
//data_from_D1/D2    - final results received back from neighbours

    double *data_received = (double*)malloc(M * sizeof(double));
    double *data_from_send = (double*)malloc(M * sizeof(double));
    double *data_at_D1 = (double*)malloc(M * sizeof(double));
    double *data_at_D2 = (double*)malloc(M * sizeof(double));
    double *buffer_updated_for_D1 = (double*)malloc(M * sizeof(double));
    double *buffer_updated_for_D2 = (double*)malloc(M * sizeof(double));
    double *data_from_D1 = (double*)malloc(M * sizeof(double));
    double *data_from_D2 = (double*)malloc(M * sizeof(double));
    double maxD1=0;
    double maxD2=0;
    double stime,etime,time;
    int t=1;  
    
    // Initialize data_from_send to avoid garbage math if a Recv is skipped
    for(llu i=0; i<M; i++) data_from_send[i] = 0.0;

    stime=MPI_Wtime(); // Start clock timer 
    while(t<=T){
        if(t==1){
            srand(seed);
            for(llu i=0;i<M;i++){
                data_received[i]=(double)rand()*(myrank+1)/10000.0;
            }
            
            int flag1=myrank/D1;
            // Communication for D1 - Corrected to prevent deadlock
            if(flag1 % 2 == 0){
                if((myrank+D1) < P) MPI_Send(data_received, M, MPI_DOUBLE, myrank+D1, 10, MPI_COMM_WORLD);
                if((myrank-D1) >= 0) MPI_Recv(data_from_send, M, MPI_DOUBLE, myrank-D1, 10, MPI_COMM_WORLD, &status);
            } else {
                if((myrank-D1) >= 0) MPI_Recv(data_from_send, M, MPI_DOUBLE, myrank-D1, 10, MPI_COMM_WORLD, &status);
                if((myrank+D1) < P) MPI_Send(data_received, M, MPI_DOUBLE, myrank+D1, 10, MPI_COMM_WORLD);
            }

            for(llu i=0;i<M;i++) data_at_D1[i]=data_from_send[i]*data_from_send[i];

            if(flag1 % 2 == 1){
                if((myrank-D1) >= 0) MPI_Send(data_at_D1, M, MPI_DOUBLE, myrank-D1, 20, MPI_COMM_WORLD);
                if((myrank+D1) < P) MPI_Recv(data_from_D1, M, MPI_DOUBLE, myrank+D1, 20, MPI_COMM_WORLD, &status);
            } else {
                if((myrank+D1) < P) MPI_Recv(data_from_D1, M, MPI_DOUBLE, myrank+D1, 20, MPI_COMM_WORLD, &status);
                if((myrank-D1) >= 0) MPI_Send(data_at_D1, M, MPI_DOUBLE, myrank-D1, 20, MPI_COMM_WORLD);
            }

            int flag2=myrank/D2;
            // Communication for D2
            if(flag2 % 2 == 0){
                if((myrank+D2) < P) MPI_Send(data_received, M, MPI_DOUBLE, myrank+D2, 30, MPI_COMM_WORLD);
                if((myrank-D2) >= 0) MPI_Recv(data_from_send, M, MPI_DOUBLE, myrank-D2, 30, MPI_COMM_WORLD, &status);
            } else {
                if((myrank-D2) >= 0) MPI_Recv(data_from_send, M, MPI_DOUBLE, myrank-D2, 30, MPI_COMM_WORLD, &status);
                if((myrank+D2) < P) MPI_Send(data_received, M, MPI_DOUBLE, myrank+D2, 30, MPI_COMM_WORLD);
            }

            for(llu i=0;i<M;i++) data_at_D2[i]=log(data_from_send[i]+DBL_MIN); // +1 to avoid log(0)

            if(flag2 % 2 == 1){
                if((myrank-D2) >= 0) MPI_Send(data_at_D2, M, MPI_DOUBLE, myrank-D2, 40, MPI_COMM_WORLD);
                if((myrank+D2) < P) MPI_Recv(data_from_D2, M, MPI_DOUBLE, myrank+D2, 40, MPI_COMM_WORLD, &status);
            } else {
                if((myrank+D2) < P) MPI_Recv(data_from_D2, M, MPI_DOUBLE, myrank+D2, 40, MPI_COMM_WORLD, &status);
                if((myrank-D2) >= 0) MPI_Send(data_at_D2, M, MPI_DOUBLE, myrank-D2, 40, MPI_COMM_WORLD);
            }

            for(llu i=0;i<M;i++){
                buffer_updated_for_D1[i]=(unsigned long long)data_from_D1[i]%100000;
                buffer_updated_for_D2[i]=data_from_D2[i]*100000;
            }
        }
        else {
            // Logic for t>1 (Using buffer_updated)
            int flag1=myrank/D1;
            if(flag1 % 2 == 0){
                if((myrank+D1) < P) MPI_Send(buffer_updated_for_D1, M, MPI_DOUBLE, myrank+D1, 50, MPI_COMM_WORLD);
                if((myrank-D1) >= 0) MPI_Recv(data_from_send, M, MPI_DOUBLE, myrank-D1, 50, MPI_COMM_WORLD, &status);
            } else {
                if((myrank-D1) >= 0) MPI_Recv(data_from_send, M, MPI_DOUBLE, myrank-D1, 50, MPI_COMM_WORLD, &status);
                if((myrank+D1) < P) MPI_Send(buffer_updated_for_D1, M, MPI_DOUBLE, myrank+D1, 50, MPI_COMM_WORLD);
            }

            for(llu i=0;i<M;i++) data_at_D1[i]=data_from_send[i]*data_from_send[i];

            if(flag1 % 2 == 1){
                if((myrank-D1) >= 0) MPI_Send(data_at_D1, M, MPI_DOUBLE, myrank-D1, 60, MPI_COMM_WORLD);
                if((myrank+D1) < P) MPI_Recv(data_from_D1, M, MPI_DOUBLE, myrank+D1, 60, MPI_COMM_WORLD, &status);
            } else {
                if((myrank+D1) < P) MPI_Recv(data_from_D1, M, MPI_DOUBLE, myrank+D1, 60, MPI_COMM_WORLD, &status);
                if((myrank-D1) >= 0) MPI_Send(data_at_D1, M, MPI_DOUBLE, myrank-D1, 60, MPI_COMM_WORLD);
            }

            int flag2=myrank/D2;
            if(flag2 % 2 == 0){
                if((myrank+D2) < P) MPI_Send(buffer_updated_for_D2, M, MPI_DOUBLE, myrank+D2, 70, MPI_COMM_WORLD);
                if((myrank-D2) >= 0) MPI_Recv(data_from_send, M, MPI_DOUBLE, myrank-D2, 70, MPI_COMM_WORLD, &status);
            } else {
                if((myrank-D2) >= 0) MPI_Recv(data_from_send, M, MPI_DOUBLE, myrank-D2, 70, MPI_COMM_WORLD, &status);
                if((myrank+D2) < P) MPI_Send(buffer_updated_for_D2, M, MPI_DOUBLE, myrank+D2, 70, MPI_COMM_WORLD);
            }

            for(llu i=0;i<M;i++) data_at_D2[i]=log(data_from_send[i]+DBL_MIN);

            if(flag2 % 2 == 1){
                if((myrank-D2) >= 0) MPI_Send(data_at_D2, M, MPI_DOUBLE, myrank-D2, 80, MPI_COMM_WORLD);
                if((myrank+D2) < P) MPI_Recv(data_from_D2, M, MPI_DOUBLE, myrank+D2, 80, MPI_COMM_WORLD, &status);
            } else {
                if((myrank+D2) < P) MPI_Recv(data_from_D2, M, MPI_DOUBLE, myrank+D2, 80, MPI_COMM_WORLD, &status);
                if((myrank-D2) >= 0) MPI_Send(data_at_D2, M, MPI_DOUBLE, myrank-D2, 80, MPI_COMM_WORLD);
            }

            for(llu i=0;i<M;i++){
                buffer_updated_for_D1[i]=(unsigned long long)data_from_D1[i]%100000;
                buffer_updated_for_D2[i]=data_from_D2[i]*100000;
            }
        }
        t++;
    }

    // Final reduction logic moved OUTSIDE the while loop
    double maxD1_local=0;
    double maxD2_local=0;
    for(llu i=0;i<M;i++) {
        if(data_from_D1[i]>maxD1_local){
            maxD1_local=data_from_D1[i];
        }
        if(data_from_D2[i]>maxD2_local){
            maxD2_local=data_from_D2[i];
        }
    }

    if(myrank==0){
        maxD1=maxD1_local;
        maxD2=maxD2_local;
        double remote_max1, remote_max2;
        for(int i=1;i<P;i++) {
            MPI_Recv(&remote_max1,1,MPI_DOUBLE,i,100,MPI_COMM_WORLD,&status);
            MPI_Recv(&remote_max2,1,MPI_DOUBLE,i,101,MPI_COMM_WORLD,&status);
            if(remote_max1>maxD1){
                maxD1=remote_max1;
            }
            if(remote_max2>maxD2){
                maxD2=remote_max2;
            }
        }
    } 
    else{
        MPI_Send(&maxD1_local,1,MPI_DOUBLE,0,100,MPI_COMM_WORLD);
        MPI_Send(&maxD2_local,1,MPI_DOUBLE,0,101,MPI_COMM_WORLD);
    }

    etime=MPI_Wtime();
    
    double local_time,max_time;
    local_time=etime-stime;
    
    MPI_Reduce(&local_time,&max_time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

    if(myrank==0){
        printf("%f %f %f\n",maxD1,maxD2,max_time);
    }
    
    MPI_Finalize();
    return 0;

}