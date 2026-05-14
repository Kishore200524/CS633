#include <mpi.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Global variables declaration:
int nx, ny, nz;
int px, py, pz; 
int F, T, h, halo_w; 
int myrank, nprocs;

int left_rank, right_rank, down_rank, up_rank, back_rank, front_rank; 
int has_left, has_right, has_down, has_up, has_back, has_front; // flags for neighbor existence
int gx, gy, gz; // local padded dimensions for one rank (including ghost layers)
size_t sx, sy, field_size; // strides and total size for padded grid -> we did precomputation to avoid doing it repeatedly in the inner loops
size_t face_x, face_y, face_z; // number of elements in each face for halo exchange
double *current_data, *next_data; // pointers to the current and next time step data arrays (each holds F fields in a contiguous block)

// Buffers for sending and receiving halo data along each axis(allocated only if face size is non zero)
double *send_left, *send_right, *recv_left, *recv_right; // X axis
double *send_down, *send_up, *recv_down, *recv_up; // Y axis
double *send_back, *send_front, *recv_back, *recv_front; // Z axis

int *valid_minus[3], *valid_plus[3]; // stores count of valid stencil neighbors available, minus is for left/down/back and plus is for right/up/front

// Getting 1D index from 3D coordinates (including halo)
size_t idx(int x, int y, int z) {
    return (size_t)x*sx + (size_t)y*sy + (size_t)z;
}

// Checks for crossing of isovalue
unsigned long long crosses(double a, double b, double iso) {
    int k1 = (a<=iso && iso<=b);
    int k2 = (b<=iso && iso<=a);
    return (k1||k2) ? 1ULL : 0ULL;
}

// Packs the face values of a field into a buffer for sending to a neighbor along a given axis and side
void pack_face(double *field, double *buf, int axis, int side) {
    size_t p=0;
    int layer, x, y, z;
    // getting exact layer of local grid that neighbor needs to fill halo
    // 0 gets the first local layers, 1 gets the last local layers
    if(axis == 0){
        for (layer=0; layer<halo_w; layer++) {
            int x_side;
            if (side==0){
                x_side = halo_w+layer;
            }
            else {
                x_side = halo_w+nx-halo_w+layer;
            }

            for (y=0; y<ny; y++) {
                for (z=0; z<nz; z++) {
                    buf[p++]=field[idx(x_side,halo_w+y,halo_w+z)];
                }
            }
        }
    } else if(axis==1){
        for (layer=0; layer<halo_w; layer++) {
            int y_side;
            if (side==0){
                y_side = halo_w+layer;
            }
            else {
                y_side = halo_w+ny-halo_w+layer;
            }

            for (x=0; x<nx; x++) {
                for (z=0; z<nz; z++) {
                    buf[p++]=field[idx(halo_w+x,y_side,halo_w+z)];
                }
            }
        }
    } else {
        for (layer=0; layer<halo_w; layer++) {
            int z_side;
            if (side==0){
                z_side = halo_w+layer;
            }
            else {
                z_side = halo_w+nz-halo_w+layer;
            }

            for (x=0; x<nx; x++) {
                for (y=0; y<ny; y++) {
                    buf[p++]=field[idx(halo_w+x,halo_w+y,z_side)];
                }
            }
        }
    }
}

// Writes received halo values from a neighbor into the ghost-cell layers of the local grid along a given axis
void unpack_face(double *field, double *buf, int axis, int side) {
    size_t p=0;
    int layer, x, y, z;

    if(axis==0){
        for (layer=0; layer<halo_w; layer++) {
            int x_side;
            if (side==0){
                x_side = halo_w-1-layer;
            } else {
                x_side = halo_w+nx+layer;
            }

            for (y=0; y<ny; y++) {
                for (z=0; z<nz; z++) {
                    field[idx(x_side, halo_w+y, halo_w+z)]=buf[p++];
                }
            }
        }
    } else if(axis==1){
        for (layer=0; layer<halo_w; layer++) {
            int y_side;
            if (side==0){
                y_side = halo_w-1-layer;
            } else {
                y_side = halo_w+ny+layer;
            }
            for (x=0; x<nx; x++) {
                for (z=0; z<nz; z++) {
                    field[idx(halo_w+x, y_side, halo_w+z)]=buf[p++];
                }
            }
        }
    } else {
        for (layer=0; layer<halo_w; layer++) {
            int z_side;
            if (side==0){
                z_side = halo_w-1-layer;
            } else {
                z_side = halo_w + nz + layer;
            }
            for (x=0; x<nx; x++) {
                for (y=0; y<ny; y++) {
                    field[idx(halo_w+x, halo_w+y, z_side)]=buf[p++];
                }
            }
        }
    }
}

// This is non blocking halo communication (so that overlap can be done)
int start_halo(double *data, MPI_Request *reqs) {
    int nreq=0;
    int f;

    if (halo_w==0) return 0;

    // Doing all the receives first so that communication can be overlapped with computation of interior points, which does not depend on halo data
    if (has_left)  MPI_Irecv(recv_left, (int)((size_t)F*face_x), MPI_DOUBLE, left_rank, 0, MPI_COMM_WORLD, &reqs[nreq++]);
    if (has_right) MPI_Irecv(recv_right, (int)((size_t)F*face_x), MPI_DOUBLE, right_rank, 1, MPI_COMM_WORLD, &reqs[nreq++]);
    if (has_down)  MPI_Irecv(recv_down, (int)((size_t)F*face_y), MPI_DOUBLE, down_rank, 2, MPI_COMM_WORLD, &reqs[nreq++]);
    if (has_up)    MPI_Irecv(recv_up, (int)((size_t)F*face_y), MPI_DOUBLE, up_rank, 3, MPI_COMM_WORLD, &reqs[nreq++]);
    if (has_back)  MPI_Irecv(recv_back, (int)((size_t)F*face_z), MPI_DOUBLE, back_rank, 4, MPI_COMM_WORLD, &reqs[nreq++]);
    if (has_front) MPI_Irecv(recv_front, (int)((size_t)F*face_z), MPI_DOUBLE, front_rank, 5, MPI_COMM_WORLD, &reqs[nreq++]);

    for (f=0; f<F; f++) {
        double *field=data+(size_t)f*field_size;
        if (has_left)  pack_face(field,send_left+(size_t)f*face_x,0,0);
        if (has_right) pack_face(field,send_right+(size_t)f*face_x,0,1);
        if (has_down)  pack_face(field,send_down+(size_t)f*face_y,1,0);
        if (has_up)    pack_face(field,send_up+(size_t)f*face_y,1,1);
        if (has_back)  pack_face(field,send_back+(size_t)f*face_z,2,0);
        if (has_front) pack_face(field,send_front+(size_t)f*face_z,2,1);
    }

    if (has_left)  MPI_Isend(send_left, (int)((size_t)F*face_x), MPI_DOUBLE, left_rank, 1, MPI_COMM_WORLD, &reqs[nreq++]);
    if (has_right) MPI_Isend(send_right, (int)((size_t)F*face_x), MPI_DOUBLE, right_rank, 0, MPI_COMM_WORLD, &reqs[nreq++]);
    if (has_down)  MPI_Isend(send_down, (int)((size_t)F*face_y), MPI_DOUBLE, down_rank, 3, MPI_COMM_WORLD, &reqs[nreq++]);
    if (has_up)    MPI_Isend(send_up, (int)((size_t)F*face_y), MPI_DOUBLE, up_rank, 2, MPI_COMM_WORLD, &reqs[nreq++]);
    if (has_back)  MPI_Isend(send_back, (int)((size_t)F*face_z), MPI_DOUBLE, back_rank, 5, MPI_COMM_WORLD, &reqs[nreq++]);
    if (has_front) MPI_Isend(send_front, (int)((size_t)F*face_z), MPI_DOUBLE, front_rank, 4, MPI_COMM_WORLD, &reqs[nreq++]);

    return nreq;
}

// Waits for all halo communication to complete and stores the received data into ghost-cell regions
void finish_halo(double *data, MPI_Request *reqs, int nreq) {
    int i;
    if (nreq>0) MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);
    for (i=0; i<F; i++) {
        double *field=data+(size_t)i*field_size;
        if (has_left)  unpack_face(field,recv_left+(size_t)i*face_x,0,0);
        if (has_right) unpack_face(field,recv_right+(size_t)i*face_x,0,1);
        if (has_down)  unpack_face(field,recv_down+(size_t)i*face_y,1,0);
        if (has_up)    unpack_face(field,recv_up+(size_t)i*face_y,1,1);
        if (has_back)  unpack_face(field,recv_back+(size_t)i*face_z,2,0);
        if (has_front) unpack_face(field,recv_front+(size_t)i*face_z,2,1);
    }
}

// Initializing the data array as mentioned in the assignment
void init_data(int seed) {
    int i, x, y, z;
    srand((unsigned int)seed);
    for (i=0; i<F; i++) {
        size_t t=0;
        for (x=0; x<nx; x++) {
            for (y=0; y<ny; y++) {
                for (z=0; z<nz; z++,t++) {
                    current_data[(size_t)i*field_size + idx(halo_w+x,halo_w+y,halo_w+z)]=
                        (double)rand()*(myrank+1)/(110426.0+i+t);
                }
            }
        }
    }
}

// Precomputes how many stencil points are valid on each side of each local index
void build_valid_counts() {
    int i;
    // 0 is for local layers adjacent to 0 neighbor (left/down/back side), 
    // 1 is for local layers adjacent to 1 neighbor (right/up/front side)
    valid_minus[0]=(int*)malloc((size_t)nx*sizeof(int));
    valid_plus[0]=(int*)malloc((size_t)nx*sizeof(int));
    valid_minus[1]=(int*)malloc((size_t)ny*sizeof(int));
    valid_plus[1]=(int*)malloc((size_t)ny*sizeof(int));
    valid_minus[2]=(int*)malloc((size_t)nz*sizeof(int));
    valid_plus[2]=(int*)malloc((size_t)nz*sizeof(int));

    if (!valid_minus[0] || !valid_plus[0] || !valid_minus[1] || !valid_plus[1] || !valid_minus[2] || !valid_plus[2]) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (i=0; i<nx; i++) {
        valid_minus[0][i]=has_left? h : (i<h?i:h);
        valid_plus[0][i]=has_right? h : ((nx-1-i<h)?(nx-1-i):h);
    }
    for (i=0; i<ny; i++) {
        valid_minus[1][i]=has_down? h : (i<h?i:h);
        valid_plus[1][i]=has_up? h : ((ny-1-i<h)?(ny-1-i):h);
    }
    for (i=0; i<nz; i++) {
        valid_minus[2][i]=has_back? h : (i<h?i:h);
        valid_plus[2][i]=has_front? h : ((nz-1-i<h)?(nz-1-i):h);
    }
}

// Updates grid points that do not require halo data
void compute_interior(double *src, double *dst) {
    int f, x, y, z, i;
    long long step_x=(long long)sx;
    long long step_y=(long long)sy;

    if (h==0) {
        for (f=0; f<F; f++) {
            double *src_field= src + (size_t)f*field_size;
            double *dst_field= dst + (size_t)f*field_size;
            for (x=0; x<nx; x++) {
                for (y=0; y<ny; y++) {
                    memcpy(dst_field+idx(halo_w+x,halo_w+y,halo_w), src_field+idx(halo_w+x,halo_w+y,halo_w), (size_t)nz*sizeof(double));
                }
            }
        }
        return;
    }

    if (nx<=2*h || ny<=2*h || nz<=2*h) return;
    for (f=0; f<F; f++){
        double *src_field=src+(size_t)f*field_size;
        double *dst_field=dst+(size_t)f*field_size;
        int cnt = 1 + 6*h;
        for (x=h; x<nx-h; x++) {
            for (y=h; y<ny-h; y++) {
                // used the pointers to avoid calculating the index repeatedly in the innermost loop, which is the bottleneck
                double *p = src_field + idx(halo_w+x,halo_w+y,halo_w+h);
                double *q = dst_field + idx(halo_w+x,halo_w+y,halo_w+h);
                for (z=h; z<nz-h; z++,p++,q++) {
                    double sum=*p;
                    for (i=1; i<=h; i++) {
                        sum+=p[-(long long)i*step_x];
                        sum+=p[(long long)i*step_x];
                        sum+=p[-(long long)i*step_y];
                        sum+=p[(long long)i*step_y];
                        sum+=p[-i];
                        sum+=p[i];
                    }
                    *q = sum/(double)cnt;
                }
            }
        }
    }
}

// FOr updating grid points near edges (stencil might go out of local domain)
void update_box(double *src_field, double *dst_field,int x0, int x1, int y0, int y1, int z0, int z1) {
    int x, y, z, i;
    // checking validity of box
    for (x=x0; x<x1; x++) {

        int x_minus=valid_minus[0][x];
        int x_plus=valid_plus[0][x];

        for (y=y0; y<y1; y++) {
            int y_minus = valid_minus[1][y];
            int y_plus = valid_plus[1][y];
            for (z=z0; z<z1; z++) {
                int z_minus = valid_minus[2][z];
                int z_plus = valid_plus[2][z];
                double sum = src_field[idx(halo_w+x,halo_w+y,halo_w+z)];
                int cnt = 1 + x_minus+x_plus+y_minus+y_plus+z_minus+z_plus;

                for(i=1; i<=x_minus; i++) sum+=src_field[idx(halo_w+x-i,halo_w+y,halo_w+z)];
                for(i=1; i<=x_plus; i++) sum+=src_field[idx(halo_w+x+i,halo_w+y,halo_w+z)];
                for(i=1; i<=y_minus; i++) sum+=src_field[idx(halo_w+x,halo_w+y-i,halo_w+z)];
                for(i=1; i<=y_plus; i++) sum+=src_field[idx(halo_w+x,halo_w+y+i,halo_w+z)];
                for(i=1; i<=z_minus; i++) sum+=src_field[idx(halo_w+x,halo_w+y,halo_w+z-i)];
                for(i=1; i<=z_plus; i++) sum+=src_field[idx(halo_w+x,halo_w+y,halo_w+z+i)];

                dst_field[idx(halo_w+x,halo_w+y,halo_w+z)]=sum/(double)cnt;
            }
        }
    }
}

// Updates grid points near subdomain edges where the stencil go out of local domain
void compute_boundary (double *src, double *dst) {
    int f;
    if (h==0) return;
    for (f=0; f<F; f++) {
        double *src_field = src+ (size_t)f*field_size;
        double *dst_field = dst+ (size_t)f*field_size;
        update_box(src_field,dst_field,0,h,0,ny,0,nz);
        update_box(src_field,dst_field,nx-h,nx,0,ny,0,nz);
        update_box(src_field,dst_field,h,nx-h,0,h,0,nz);
        update_box(src_field,dst_field,h,nx-h,ny-h,ny,0,nz);
        update_box(src_field,dst_field,h,nx-h,h,ny-h,0,h);
        update_box(src_field,dst_field,h,nx-h,h,ny-h,nz-h,nz);
    }
}

// Counts how many edges between neighboring grid points cross the given isovalue, including edges that connect to neighboring subdomains
// Double counting is done here
void count_iso(double *data, double iso, unsigned long long *my_count) {
    int i, x, y, z;
    long long step_x = (long long)sx;
    long long step_y = (long long)sy;

    for (i=0; i<F; i++) {
        double *src_field = data + (size_t)i*field_size;
        unsigned long long count=0ULL;
        for (x=0; x<nx; x++) {
            int check_x_plus= (x+1<nx || has_right);
            int check_x_minus= (x==0 && has_left);
            for (y=0; y<ny; y++) {

                double *row = src_field + idx(halo_w+x,halo_w+y,halo_w);
                int check_y_plus = (y+1 < ny || has_up);
                int check_y_minus = (y==0 && has_down);
                for (z=0; z<nz; z++) {
                    double *p = row+z;
                    double v=*p;
                    if (check_x_plus)  count+=crosses(v,*(p+step_x),iso);
                    if (check_x_minus) count+=crosses(v,*(p-step_x),iso);
                    if (check_y_plus)  count+=crosses(v,*(p+step_y),iso);
                    if (check_y_minus) count+=crosses(v,*(p-step_y),iso);
                    if (z+1<nz||has_front) count+=crosses(v,*(p+1),iso);
                    if (z==0&&has_back) count+=crosses(v,*(p-1),iso);
                }
            }
        }

        my_count[i]=count;
    }
}

// Now using the above functions to implement the main logic
// Flow: start halo -> interior -> wait -> count previous -> boundary -> swap 
int main(int argc, char **argv) {
    int d, ppn, seed, step, f, nreq;
    double isovalue;
    unsigned long long *my_count, *total_count;
    double start_time, end_time, local_time, max_time;
    MPI_Request reqs[12];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc!=13) {
            MPI_Abort(MPI_COMM_WORLD, 1);
    }

    d=atoi(argv[1]);
    ppn=atoi(argv[2]);
    px=atoi(argv[3]);
    py=atoi(argv[4]);
    pz=atoi(argv[5]);
    nx=atoi(argv[6]);
    ny=atoi(argv[7]);
    nz=atoi(argv[8]);
    T=atoi(argv[9]);
    seed=atoi(argv[10]);
    F=atoi(argv[11]);
    isovalue=atof(argv[12]);
    (void)ppn;
    
    // Error checking for input
    if (px<=0||py<=0||pz<=0||nx<=0||ny<=0||nz<=0||F<=0||T<0) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (nprocs!=px*py*pz) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (d<1||(d-1)%6!=0) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (d>=nx||d>=ny||d>=nz) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Done for stencil of form 1+6r only
    h=(d-1)/6;
    halo_w=(h>0)?h:1;

    // Determine neighbors in the 3D process grid
    int pi=myrank%px;
    int pj=(myrank/px)%py;
    int pk=myrank/(px*py);
    left_rank=(pi>0)?myrank-1:-1;
    right_rank=(pi<px-1)?myrank+1:-1;
    down_rank=(pj>0)?myrank-px:-1;
    up_rank=(pj<py-1)?myrank+px:-1;
    back_rank=(pk>0)?myrank-px*py:-1;
    front_rank=(pk<pz-1)?myrank+px*py:-1;
    has_left=(left_rank>=0);
    has_right=(right_rank>=0);
    has_down=(down_rank>=0);
    has_up=(up_rank>=0);
    has_back=(back_rank>=0);
    has_front=(front_rank>=0);
    gx=nx+2*halo_w;
    gy=ny+2*halo_w;
    gz=nz+2*halo_w;
    sy=(size_t)gz;
    sx=(size_t)gy*(size_t)gz;
    field_size=(size_t)gx*(size_t)gy*(size_t)gz;

    // Allocating memory for current and next data arrays, count arrays, and halo buffers
    // Here we kept as one contiguous block instead of double **.
    current_data=(double*)calloc((size_t)F*field_size, sizeof(double));
    next_data=(double*)calloc((size_t)F*field_size, sizeof(double));
    my_count=(unsigned long long*)malloc((size_t)F*sizeof(unsigned long long));
    total_count=(myrank==0)?(unsigned long long*)malloc((size_t)F*sizeof(unsigned long long)):NULL;

    if (!current_data || !next_data ||!my_count||(myrank==0&&!total_count)) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    face_x = (size_t)halo_w*(size_t)ny*(size_t)nz;
    face_y = (size_t)nx*(size_t)halo_w*(size_t)nz;
    face_z = (size_t)nx*(size_t)ny*(size_t)halo_w;

    // Halo buffers for each face and field
    send_left=face_x?(double*)malloc((size_t)F*face_x*sizeof(double)):NULL;
    send_right=face_x?(double*)malloc((size_t)F*face_x*sizeof(double)):NULL;
    recv_left=face_x?(double*)malloc((size_t)F*face_x*sizeof(double)):NULL;
    recv_right=face_x?(double*)malloc((size_t)F*face_x*sizeof(double)):NULL;

    send_down=face_y?(double*)malloc((size_t)F*face_y*sizeof(double)):NULL;
    send_up=face_y?(double*)malloc((size_t)F*face_y*sizeof(double)):NULL;
    recv_down=face_y?(double*)malloc((size_t)F*face_y*sizeof(double)):NULL;
    recv_up=face_y?(double*)malloc((size_t)F*face_y*sizeof(double)):NULL;

    send_back=face_z?(double*)malloc((size_t)F*face_z*sizeof(double)):NULL;
    send_front=face_z?(double*)malloc((size_t)F*face_z*sizeof(double)):NULL;
    recv_back=face_z?(double*)malloc((size_t)F*face_z*sizeof(double)):NULL;
    recv_front=face_z?(double*)malloc((size_t)F*face_z*sizeof(double)):NULL;

    if ((face_x&&(!send_left||!send_right||!recv_left||!recv_right))||(face_y&&(!send_down||!send_up||!recv_down||!recv_up))||(face_z&&(!send_back||!send_front||!recv_back||!recv_front))) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    init_data(seed);
    build_valid_counts();

    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    for (step=1; step<=T; step++) { // each step corresponds to computing next time step from current time step
        double *temp;

        // Here we start the halo exchange and in the mean time we computed the interior points so that communication and computation are overlapped.
        nreq=start_halo(current_data,reqs);
        compute_interior(current_data,next_data);
        finish_halo(current_data,reqs,nreq); // waiting for halo to complete

        if (step>1) {
            count_iso(current_data,isovalue,my_count);
            MPI_Reduce(my_count,total_count,F,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);

            if (myrank==0) {
                for (f=0; f<F; f++) {
                    if (f) printf(" ");
                    printf("%llu", total_count[f]);
                }
                printf("\n");
            }
        }

        // Boundary points need the halo values, so doing them after the wait.
        compute_boundary(current_data,next_data);
        temp=current_data;
        current_data=next_data;
        next_data=temp;
    }

    if (T>0) { // counting for the last time step, which is not counted inside the loop
        nreq=start_halo(current_data,reqs);
        finish_halo(current_data,reqs,nreq);
        count_iso(current_data,isovalue,my_count);
        MPI_Reduce(my_count,total_count,F,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);

        if (myrank==0) {
            for (f=0; f<F; f++) {
                if (f) printf(" ");
                printf("%llu", total_count[f]);
            }
            printf("\n");
        }
    }

    end_time=MPI_Wtime();
    local_time=end_time-start_time;
    max_time=0.0;
    MPI_Reduce(&local_time,&max_time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

    if (myrank == 0) {
        printf("%.12f\n", max_time);
    }

    free(current_data);
    free(next_data);
    free(my_count);
    free(total_count);

    free(send_left);  free(send_right); free(recv_left);  free(recv_right);
    free(send_down);  free(send_up);    free(recv_down);  free(recv_up);
    free(send_back);  free(send_front); free(recv_back);  free(recv_front);
    free(valid_minus[0]); free(valid_plus[0]);
    free(valid_minus[1]); free(valid_plus[1]);
    free(valid_minus[2]); free(valid_plus[2]);

    MPI_Finalize();
    return 0;
}

