#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<string.h>
#include<math.h>
#include<time.h>
#define nt 10
#define mt 3
/*#define BLOCK_SIZE_xf 16
#define BLOCK_SIZE_x  32
#define BLOCK_SIZE_yf 4
#define BLOCK_SIZE_yg 4
#define BLOCK_SIZE_yh 2*/
#define BLOCK_SIZE_xf 16  // can not change
#define BLOCK_SIZE_x  32
#define BLOCK_SIZE_yf 4
#define BLOCK_SIZE_yg 4
#define BLOCK_SIZE_yh 2

#define ama 1.1
#define epweno 1.0e-6

#define checkCudaErrors(err)  __checkCudaErrors (err, __FILE__, __LINE__)

inline void __checkCudaErrors(cudaError err, const char *file, const int line )
{
    if(cudaSuccess != err)
    {
        fprintf(stderr, "%s(%i) : CUDA Runtime API error %d: %s.\n",file, line, (int)err, cudaGetErrorString( err ) );
        exit(-1);        
    }
}

__global__ void fx_gpu(int,int,int,int,double,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,cudaTextureObject_t,cudaTextureObject_t,cudaTextureObject_t,cudaTextureObject_t,cudaTextureObject_t,cudaTextureObject_t);
__global__ void fxleft_gpu(double*,int,int,int,double*);
__global__ void gy_gpu(double*,int,int,int,int,double,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,cudaTextureObject_t,cudaTextureObject_t,cudaTextureObject_t,cudaTextureObject_t,cudaTextureObject_t,cudaTextureObject_t);
__global__ void hz_gpu(int,int,int,int,double,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,cudaTextureObject_t,cudaTextureObject_t,cudaTextureObject_t,cudaTextureObject_t,cudaTextureObject_t,cudaTextureObject_t);
__global__ void hzleft_gpu(double*,int,int,int,double*);
__global__ void ucIn(double*,double*,int,int,int,int,int);
__global__ void rhsIn(double*,double*,int,int,int,int,double);
__global__ void rhsOut(double*,double*,int,int,int,int);

__inline__ __device__ double fetch(cudaTextureObject_t tex,float i,float j){
	int2 data=tex2D<int2>(tex,i,j);
	return __hiloint2double(data.y,data.x);	
}

extern "C"
void fghx_gpuweno_(double* f_uc,double* f_rhs,double* f_ark,long* f_io,double* f_gamma,long* f_nLx,long* f_nLy,long* f_nLz,long* f_mnp,long* f_myid)
{
	clock_t tstart,tdec,tmalloc,tkernel,tdata,tres,tfree;
	tstart=clock();
	static int step=1;
	//cpu
	int nLx=*f_nLx;//!
	int nLy=*f_nLy;//!
	int nLz=*f_nLz;//!
    int io=*f_io;
	int mnp=*f_mnp;//!
    int myid=*f_myid;
	double* uc=f_uc;//!
	double* rhs=f_rhs;//!

    double* ark=f_ark;
	double gamma=*f_gamma;//!
    double qk=*(ark+io);//!	
	
	static int LLx,LLy,LLz;
	size_t shmem;
    //gpu
	static size_t ucsize,rhssize; 
	static double* d_uc;
    static double* d_rhs;
    static double* d_fh;
    static double* d_fd,*d_fdr;
    static double *d_gg1,* d_gg2;
    static double *d_cs,*d_vx,*d_vy,*d_vz,*d_w,*d_h;
	static cudaEvent_t start,stop;
	float kxtime,kytime,kztime,kuit,krit,krot;
	static dim3 dimGrid,dimBlock;
	size_t offset;
	static size_t pitch;
	static const cudaChannelFormatDesc channelDesc=cudaCreateChannelDesc<int2>();
//	printf("myid:%d\n",myid);
    cudaSetDevice(myid%2);
//    cudaSetDevice(0);
	static double* h_uc;
	static double* h_rhs;
	static cudaStream_t stream[2];
    tdec=clock();
if(step==1){
	LLx=nLx+10;
	LLy=nLy+10;
	LLz=nLz+10;
	int ls=min(nLx,nLy);
	ls=min(ls,nLz);
    int msize=nLx*nLy*nLz/ls;//small
	int LL=ls+10;
	int LL1=ls+5;
	int LL2=ls+1;
	ucsize=mnp*LLx*LLy*LLz*sizeof(double); // 10=-4~0, +1~+5 5=-2~0, +1~+2//////////////////////////////
	rhssize=mnp*nLx*nLy*nLz*sizeof(double);
	size_t fhsize=mnp*LL2*msize*sizeof(double); 
	size_t ggsize=mnp*LL1*msize*sizeof(double);
	size_t fdsize=mnp*LLx*LLy*LLz*sizeof(double);//used to change the data
    size_t othersize=LL*msize*sizeof(double);
	//memory copying
	assert(cudaMallocPitch((void**)&d_uc,&pitch,LLx*sizeof(double),LLy*LLz*mnp)==cudaSuccess);
	pitch=(int)(pitch/sizeof(double));
	assert(cudaMalloc((void**)&d_rhs,rhssize)==cudaSuccess);
	assert(cudaMalloc((void**)&d_fh,fhsize)==cudaSuccess);
	assert(cudaMalloc((void**)&d_fd,fdsize)==cudaSuccess);
	assert(cudaMalloc((void**)&d_fdr,fdsize)==cudaSuccess);
	assert(cudaMalloc((void**)&d_gg1,ggsize)==cudaSuccess);
	assert(cudaMalloc((void**)&d_gg2,ggsize)==cudaSuccess);
	assert(cudaMalloc((void**)&d_cs,othersize)==cudaSuccess);
	assert(cudaMalloc((void**)&d_vx,othersize)==cudaSuccess);
	assert(cudaMalloc((void**)&d_vy,othersize)==cudaSuccess);
	assert(cudaMalloc((void**)&d_vz,othersize)==cudaSuccess);
	assert(cudaMalloc((void**)&d_w,othersize)==cudaSuccess);
	assert(cudaMalloc((void**)&d_h,othersize)==cudaSuccess);
    cudaEventCreate(&start);
	cudaEventCreate(&stop);

	assert(cudaMallocHost((void**)&h_uc,ucsize)==cudaSuccess);
	assert(cudaMallocHost((void**)&h_rhs,rhssize)==cudaSuccess);
	cudaStreamCreate(&stream[0]);
	cudaStreamCreate(&stream[1]);
}
   tmalloc=clock();
	memcpy(h_uc,uc,ucsize);
	memcpy(h_rhs,rhs,rhssize);
   tdata=clock();
	assert(cudaMemcpyAsync(d_fd,h_uc,ucsize,cudaMemcpyHostToDevice,stream[0])==cudaSuccess);
	assert(cudaMemcpyAsync(d_fdr,h_rhs,rhssize,cudaMemcpyHostToDevice,stream[1])==cudaSuccess);
	
		cudaEventRecord(start,0);
	    dimGrid.x=LLz;
		dimGrid.y=1;
        dimBlock.x=LLx;
		dimBlock.y=1;
		shmem=LLx*mnp*sizeof(double);
   	    ucIn<<<dimGrid,dimBlock,shmem,stream[0]>>>(d_fd,d_uc,LLx,LLy,LLz,mnp,pitch);
	size_t width=LLx;
	size_t height=LLy*LLz;

	// create texture object
	cudaResourceDesc resDescr;
	cudaResourceDesc resDescx;
	cudaResourceDesc resDescy;
	cudaResourceDesc resDescz;
	cudaResourceDesc resDesce;
	cudaResourceDesc resDesct;
	memset(&resDescr, 0, sizeof(resDescr));
	memset(&resDescx, 0, sizeof(resDescx));
	memset(&resDescy, 0, sizeof(resDescy));
	memset(&resDescz, 0, sizeof(resDescz));
	memset(&resDesce, 0, sizeof(resDesce));
	memset(&resDesct, 0, sizeof(resDesct));
	resDescr.resType = cudaResourceTypePitch2D;
	resDescx.resType = cudaResourceTypePitch2D;
	resDescy.resType = cudaResourceTypePitch2D;
	resDescz.resType = cudaResourceTypePitch2D;
	resDesce.resType = cudaResourceTypePitch2D;
	resDesct.resType = cudaResourceTypePitch2D;
	resDescr.res.pitch2D.pitchInBytes = pitch*sizeof(double);
	resDescx.res.pitch2D.pitchInBytes = pitch*sizeof(double);
	resDescy.res.pitch2D.pitchInBytes = pitch*sizeof(double);
	resDescz.res.pitch2D.pitchInBytes = pitch*sizeof(double);
	resDesce.res.pitch2D.pitchInBytes = pitch*sizeof(double);
	resDesct.res.pitch2D.pitchInBytes = pitch*sizeof(double);
	resDescr.res.pitch2D.width = width;
	resDescx.res.pitch2D.width = width;
	resDescy.res.pitch2D.width = width;
	resDescz.res.pitch2D.width = width;
	resDesce.res.pitch2D.width = width;
	resDesct.res.pitch2D.width = width;
	resDescr.res.pitch2D.height = height;
	resDescx.res.pitch2D.height = height;
	resDescy.res.pitch2D.height = height;
	resDescz.res.pitch2D.height = height;
	resDesce.res.pitch2D.height = height;
	resDesct.res.pitch2D.height = height;
	resDescr.res.pitch2D.devPtr = d_uc;
	resDescx.res.pitch2D.devPtr = d_uc+pitch*height;
	resDescy.res.pitch2D.devPtr = d_uc+2*pitch*height;
	resDescz.res.pitch2D.devPtr = d_uc+3*pitch*height;
	resDesce.res.pitch2D.devPtr = d_uc+4*pitch*height;
	resDesct.res.pitch2D.devPtr = d_uc+5*pitch*height;
	resDescr.res.pitch2D.desc = channelDesc;
	resDescx.res.pitch2D.desc = channelDesc;
	resDescy.res.pitch2D.desc = channelDesc;
	resDescz.res.pitch2D.desc = channelDesc;
	resDesce.res.pitch2D.desc = channelDesc;
	resDesct.res.pitch2D.desc = channelDesc;

	cudaTextureDesc texDescr;
	cudaTextureDesc texDescx;
	cudaTextureDesc texDescy;
	cudaTextureDesc texDescz;
	cudaTextureDesc texDesce;
	cudaTextureDesc texDesct;
	memset(&texDescr, 0, sizeof(texDescr));
	memset(&texDescx, 0, sizeof(texDescx));
	memset(&texDescy, 0, sizeof(texDescy));
	memset(&texDescz, 0, sizeof(texDescz));
	memset(&texDesce, 0, sizeof(texDesce));
	memset(&texDesct, 0, sizeof(texDesct));
	texDescr.readMode = cudaReadModeElementType;
	texDescx.readMode = cudaReadModeElementType;
	texDescy.readMode = cudaReadModeElementType;
	texDescz.readMode = cudaReadModeElementType;
	texDesce.readMode = cudaReadModeElementType;
	texDesct.readMode = cudaReadModeElementType;
	texDescr.addressMode[0] = cudaAddressModeWrap;
	texDescx.addressMode[0] = cudaAddressModeWrap;
	texDescy.addressMode[0] = cudaAddressModeWrap;
	texDescz.addressMode[0] = cudaAddressModeWrap;
	texDesce.addressMode[0] = cudaAddressModeWrap;
	texDesct.addressMode[0] = cudaAddressModeWrap;
	texDescr.addressMode[1] = cudaAddressModeWrap;
	texDescx.addressMode[1] = cudaAddressModeWrap;
	texDescy.addressMode[1] = cudaAddressModeWrap;
	texDescz.addressMode[1] = cudaAddressModeWrap;
	texDesce.addressMode[1] = cudaAddressModeWrap;
	texDesct.addressMode[1] = cudaAddressModeWrap;
	texDescr.normalizedCoords = false;
	texDescx.normalizedCoords = false;
	texDescy.normalizedCoords = false;
	texDescz.normalizedCoords = false;
	texDesce.normalizedCoords = false;
	texDesct.normalizedCoords = false;
	// create texture object: we only have to do this once!
	cudaTextureObject_t texro=0;
	cudaTextureObject_t texvx=0;
	cudaTextureObject_t texvy=0;
	cudaTextureObject_t texvz=0;
	cudaTextureObject_t texeng=0;
	cudaTextureObject_t textracer=0;

	cudaCreateTextureObject(&texro, &resDescr, &texDescr, NULL);
	cudaCreateTextureObject(&texvx, &resDescx, &texDescx, NULL);
	cudaCreateTextureObject(&texvy, &resDescy, &texDescy, NULL);
	cudaCreateTextureObject(&texvz, &resDescz, &texDescz, NULL);
	cudaCreateTextureObject(&texeng, &resDesce, &texDesce, NULL);
	cudaCreateTextureObject(&textracer, &resDesct, &texDesct, NULL);


        dimGrid.x=nLy/BLOCK_SIZE_xf;//////////////////////////////////////
		dimGrid.y=nLz/BLOCK_SIZE_yf;//////////////////////////////////////
        dimBlock.x=BLOCK_SIZE_xf;
		dimBlock.y=BLOCK_SIZE_yf;
		fx_gpu<<<dimGrid,dimBlock,0,stream[0]>>>(nLx,nLy,nLz,mnp,gamma,d_fd,d_cs,d_vx,d_vy,d_vz,d_w,d_h,d_gg1,d_gg2,d_fh,texro,texvx,texvy,texvz,texeng,textracer);
        dimGrid.x=nLz;
		dimGrid.y=1;
        dimBlock.x=nLx;
		dimBlock.y=1;
	    shmem=nLx*mnp*sizeof(double);
	    rhsIn<<<dimGrid,dimBlock,shmem,stream[1]>>>(d_fdr,d_rhs,nLx,nLy,nLz,mnp,qk);
        
		dimGrid.x=nLy/BLOCK_SIZE_xf;//////////////////////////////////////
		dimGrid.y=nLz/BLOCK_SIZE_yf;//////////////////////////////////////
        dimBlock.x=BLOCK_SIZE_xf;
		dimBlock.y=BLOCK_SIZE_yf;
		fxleft_gpu<<<dimGrid,dimBlock>>>(d_rhs,nLx,nLy,nLz,d_fh);
	    cudaEventRecord(stop,0);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&kxtime,start,stop);
	
		cudaEventRecord(start,0);
        dimGrid.x=nLx/BLOCK_SIZE_x;
		dimGrid.y=nLz/BLOCK_SIZE_yg;//////////////////////////////////////
        dimBlock.x=BLOCK_SIZE_x;
		dimBlock.y=BLOCK_SIZE_yg;//////////////////////////////////////
		gy_gpu<<<dimGrid,dimBlock>>>(d_rhs,nLx,nLy,nLz,mnp,gamma,d_fd,d_cs,d_vx,d_vy,d_vz,d_w,d_h,d_gg1,d_gg2,d_fh,texro,texvx,texvy,texvz,texeng,textracer);
	    cudaEventRecord(stop,0);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&kytime,start,stop);

		cudaEventRecord(start,0);
        dimGrid.x=nLx/BLOCK_SIZE_x;
		dimGrid.y=nLy/BLOCK_SIZE_yh;//////////////////////////////////////
        dimBlock.x=BLOCK_SIZE_x;
		dimBlock.y=BLOCK_SIZE_yh;//////////////////////////////////////
		hz_gpu<<<dimGrid,dimBlock>>>(nLx,nLy,nLz,mnp,gamma,d_fd,d_cs,d_vx,d_vy,d_vz,d_w,d_h,d_gg1,d_gg2,d_fh,texro,texvx,texvy,texvz,texeng,textracer);
		hzleft_gpu<<<dimGrid,dimBlock>>>(d_rhs,nLx,nLy,nLz,d_fh);
	    cudaEventRecord(stop,0);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&kztime,start,stop);
	// get the result
	tkernel=clock();

    dimGrid.x=nLz;
	dimGrid.y=1;
    dimBlock.x=nLx;
	dimBlock.y=1;
	shmem=nLx*mnp*sizeof(double);
		cudaEventRecord(start,0);
    	rhsOut<<<dimGrid,dimBlock,shmem>>>(d_rhs,d_fdr,nLx,nLy,nLz,mnp);

	// destroy texture object
	cudaDestroyTextureObject(texro);
	cudaDestroyTextureObject(texvx);
	cudaDestroyTextureObject(texvy);
	cudaDestroyTextureObject(texvz);
	cudaDestroyTextureObject(texeng);
	cudaDestroyTextureObject(textracer);

	    cudaEventRecord(stop,0);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&krot,start,stop);
	assert(cudaMemcpy(h_rhs,d_fdr,rhssize,cudaMemcpyDeviceToHost)==cudaSuccess);
	memcpy(rhs,h_rhs,rhssize);
    tres=clock();
if(step==mt*nt+1){
    assert(cudaFree(d_uc)==cudaSuccess);
    assert(cudaFree(d_rhs)==cudaSuccess);
    assert(cudaFree(d_fh)==cudaSuccess);
    assert(cudaFree(d_fd)==cudaSuccess);
    assert(cudaFree(d_fdr)==cudaSuccess);
    assert(cudaFree(d_gg1)==cudaSuccess);
    assert(cudaFree(d_gg2)==cudaSuccess);
    assert(cudaFree(d_cs)==cudaSuccess);
    assert(cudaFree(d_vx)==cudaSuccess);
    assert(cudaFree(d_vy)==cudaSuccess);
    assert(cudaFree(d_vz)==cudaSuccess);
    assert(cudaFree(d_w)==cudaSuccess);
    assert(cudaFree(d_h)==cudaSuccess);
        cudaEventDestroy(start);
        cudaEventDestroy(stop);
 
    assert(cudaFreeHost(h_uc)==cudaSuccess);
    assert(cudaFreeHost(h_rhs)==cudaSuccess);
    cudaStreamDestroy(stream[0]);
    cudaStreamDestroy(stream[1]);
}
	
    tfree=clock();
	printf("dec:%f\nmalloc:%f\ndataready:%f\nucin:%f\nrhsin:%f\nkernelfx:%f\nkernelgy:%f\nkernelhz:%f\nkernel:%f\nres:%f\nrhsout:%f\nfree:%f\n",(float)(tdec-tstart)/CLOCKS_PER_SEC*1000,(float)(tmalloc-tdec)/CLOCKS_PER_SEC*1000,(float)(tdata-tmalloc)/CLOCKS_PER_SEC*1000,kuit,krit,kxtime,kytime,kztime,(float)(tkernel-tdata)/CLOCKS_PER_SEC*1000,(float)(tres-tkernel)/CLOCKS_PER_SEC*1000,krot,(float)(tfree-tres)/CLOCKS_PER_SEC*1000);
	step++;
} 

__global__ void ucIn(double* org,double* cur,int Lx,int Ly,int Lz,int mnp,int pitch)
{
	extern __shared__ double middle[];//blsx*mnp
	unsigned int tx=threadIdx.x;
	unsigned int bx=blockIdx.x;
	unsigned int block=blockDim.x;
	unsigned int index;
	for(unsigned int i=0;i<Ly;i++){
		index=tx;
	    __syncthreads();	
		middle[(index%mnp)*block+index/mnp]=org[i*mnp*Lx+bx*mnp*Lx*Ly+index];
        index+=block;		
		middle[(index%mnp)*block+index/mnp]=org[i*mnp*Lx+bx*mnp*Lx*Ly+index];
        index+=block;		
		middle[(index%mnp)*block+index/mnp]=org[i*mnp*Lx+bx*mnp*Lx*Ly+index];
        index+=block;		
		middle[(index%mnp)*block+index/mnp]=org[i*mnp*Lx+bx*mnp*Lx*Ly+index];
        index+=block;		
		middle[(index%mnp)*block+index/mnp]=org[i*mnp*Lx+bx*mnp*Lx*Ly+index];
        index+=block;		
		middle[(index%mnp)*block+index/mnp]=org[i*mnp*Lx+bx*mnp*Lx*Ly+index];
	    __syncthreads();	
		index=tx;
		cur[0*pitch*Ly*Lz+bx*pitch*Ly+i*pitch+index]=middle[0*block+index];
		cur[pitch*Ly*Lz+bx*pitch*Ly+i*pitch+index]=middle[block+index];
		cur[2*pitch*Ly*Lz+bx*pitch*Ly+i*pitch+index]=middle[2*block+index];
		cur[3*pitch*Ly*Lz+bx*pitch*Ly+i*pitch+index]=middle[3*block+index];
		cur[4*pitch*Ly*Lz+bx*pitch*Ly+i*pitch+index]=middle[4*block+index];
		cur[5*pitch*Ly*Lz+bx*pitch*Ly+i*pitch+index]=middle[5*block+index];
	}
}
__global__ void rhsIn(double* org,double* cur,int Lx,int Ly,int Lz,int mnp,double qk)
{
	extern __shared__ double middle[];//blsx*mnp
	unsigned int tx=threadIdx.x;
	unsigned int bx=blockIdx.x;
	unsigned int block=blockDim.x;
	unsigned int i;
	unsigned int index;
	for(i=0;i<Ly;i++){
		index=tx;
	    __syncthreads();	
		middle[(index%mnp)*block+index/mnp]=org[i*mnp*Lx+bx*mnp*Lx*Ly+index]*qk;
        index+=block;		
		middle[(index%mnp)*block+index/mnp]=org[i*mnp*Lx+bx*mnp*Lx*Ly+index]*qk;
        index+=block;		
		middle[(index%mnp)*block+index/mnp]=org[i*mnp*Lx+bx*mnp*Lx*Ly+index]*qk;
        index+=block;		
		middle[(index%mnp)*block+index/mnp]=org[i*mnp*Lx+bx*mnp*Lx*Ly+index]*qk;
        index+=block;		
		middle[(index%mnp)*block+index/mnp]=org[i*mnp*Lx+bx*mnp*Lx*Ly+index]*qk;
        index+=block;		
		middle[(index%mnp)*block+index/mnp]=org[i*mnp*Lx+bx*mnp*Lx*Ly+index]*qk;
	    __syncthreads();	
		index=tx;
		cur[0*Lx*Ly*Lz+bx*Lx*Ly+i*Lx+index]=middle[0*block+index];
		cur[Lx*Ly*Lz+bx*Lx*Ly+i*Lx+index]=middle[block+index];
		cur[2*Lx*Ly*Lz+bx*Lx*Ly+i*Lx+index]=middle[2*block+index];
		cur[3*Lx*Ly*Lz+bx*Lx*Ly+i*Lx+index]=middle[3*block+index];
		cur[4*Lx*Ly*Lz+bx*Lx*Ly+i*Lx+index]=middle[4*block+index];
		cur[5*Lx*Ly*Lz+bx*Lx*Ly+i*Lx+index]=middle[5*block+index];
	}
}

__global__ void rhsOut(double* org,double* cur,int Lx,int Ly,int Lz,int mnp)
{
	extern __shared__ double middle[];//blsx*mnp
	unsigned int tx=threadIdx.x;
	unsigned int bx=blockIdx.x;
	unsigned int block=blockDim.x;
	unsigned int i;
	unsigned int index;
	for(i=0;i<Ly;i++){
		index=tx;
	    __syncthreads();	
		middle[index*6]=org[0*Lx*Ly*Lz+bx*Lx*Ly+i*Lx+index];
		middle[index*6+1]=org[Lx*Ly*Lz+bx*Lx*Ly+i*Lx+index];
		middle[index*6+2]=org[2*Lx*Ly*Lz+bx*Lx*Ly+i*Lx+index];
		middle[index*6+3]=org[3*Lx*Ly*Lz+bx*Lx*Ly+i*Lx+index];
		middle[index*6+4]=org[4*Lx*Ly*Lz+bx*Lx*Ly+i*Lx+index];
		middle[index*6+5]=org[5*Lx*Ly*Lz+bx*Lx*Ly+i*Lx+index];
	    __syncthreads();	
		index=tx;
		cur[i*mnp*Lx+bx*mnp*Lx*Ly+index]=middle[index];
        index+=block;		
		cur[i*mnp*Lx+bx*mnp*Lx*Ly+index]=middle[index];
        index+=block;		
		cur[i*mnp*Lx+bx*mnp*Lx*Ly+index]=middle[index];
        index+=block;		
		cur[i*mnp*Lx+bx*mnp*Lx*Ly+index]=middle[index];
        index+=block;		
		cur[i*mnp*Lx+bx*mnp*Lx*Ly+index]=middle[index];
        index+=block;		
		cur[i*mnp*Lx+bx*mnp*Lx*Ly+index]=middle[index];
	}
}


__global__ void fx_gpu(int nLx,int nLy,int nLz,int mnp,double gamma,double* fd,double* cs,double* vx,double* vy,double* vz,double* w,double* h,double* gg1,double* gg2,double* fh,cudaTextureObject_t texro,cudaTextureObject_t texvx,cudaTextureObject_t texvy,cudaTextureObject_t texvz,cudaTextureObject_t texeng,cudaTextureObject_t textracer)
{
	int i;
	double gm1=gamma-1;
	double den,xmt,ymt,zmt,eng,tracer;
	double am1,am2,am5;
	double aden,pre,vtemp,vxm,vym,vzm;
	int ip4,ip3,ip2,ip1,i0,in1,in2;
	double ta1,ta2,ta3,ta4,ta5,ta6;
	double tg1,tg2,tg3,tg4,tg5,tg6;
	double ttg1,ttg2,ttg3,ttg4,ttg5,ttg6;
	double te0,te1,te2,te3,te4,te5,te6,te7,te8,te9,tea,teb,tec;
	double evr15,evr55,evl11,evl21,evl12,evl13,evl14,evl15,evl25;
	double hm,qm,cm,cm2,cmm;
/*	double h11,h21,h31,h41,h51,h61;
	double h12,h22,h32,h42,h52,h62;
	double h13,h23,h33,h43,h53,h63;
	double h14,h24,h34,h44,h54,h64;
	double h15,h25,h35,h45,h55,h65;
	double h16,h26,h36,h46,h56,h66;*/
	double h17,h27,h37,h47,h57,h67;
	double h18,h28,h38,h48,h58,h68;

	double tt11,tt12,tt13,tt14,tt15,tt16;
	double ts11,ts12,ts13,ts14,ts15,ts16;

	double tfh1;

	int LL=nLx+10,LL1=nLx+5,LL2=nLx+1;
	int LLy=nLy+10;
//	int LLz=nLz+10;
	int dist3dfd=LL*nLy*nLz,dist3dgg=LL1*nLy*nLz,dist3dfh=LL2*nLy*nLz;
	unsigned int tx=threadIdx.x;
	unsigned int ty=threadIdx.y;
	unsigned int xid=blockIdx.x*BLOCK_SIZE_xf+tx;
	unsigned int yid=blockIdx.y*BLOCK_SIZE_yf+ty;
//	uc+=((yid+5)*LLy+xid+5)*LL;
	fd+=yid*nLy*LL+xid;
    cs+=yid*nLy*LL+xid;
    vx+=yid*nLy*LL+xid;
    vy+=yid*nLy*LL+xid;
    vz+=yid*nLy*LL+xid;
    w+=yid*nLy*LL+xid;
    h+=yid*nLy*LL+xid;
	gg1+=yid*nLy*LL1+xid;
	gg2+=yid*nLy*LL1+xid;
    fh+=yid*nLy*LL2+xid;
//	rhs+=(yid*nLy+xid-tx)*nLx+tx;
	float texj=(yid+5)*LLy+xid+5+0.5f;
    __shared__ double shm[BLOCK_SIZE_xf*BLOCK_SIZE_yf][17];
	if(tx==0) shm[0][0]=1;

    for(i=0;i<5;i++)
	{
		den=fetch(texro,i,texj);
		xmt=fetch(texvx,i,texj);
		ymt=fetch(texvy,i,texj);
		zmt=fetch(texvz,i,texj);
		eng=fetch(texeng,i,texj);
		tracer=fetch(textracer,i,texj);

		aden=1.0/den;   ///////////////////////////cuda dividing is slow?????????????????????
		vxm=xmt*aden;
		vym=ymt*aden;
		vzm=zmt*aden;
        vtemp=0.5*(vxm*vxm+vym*vym+vzm*vzm);
		pre=gm1*(eng-vtemp*den);
         
		*(fd+i*nLy)=xmt;
		*(fd+dist3dfd+i*nLy)=xmt*vxm+pre;
		*(fd+2*dist3dfd+i*nLy)=ymt*vxm;
		*(fd+3*dist3dfd+i*nLy)=zmt*vxm;
		*(fd+4*dist3dfd+i*nLy)=vxm*(pre+eng);
		*(fd+5*dist3dfd+i*nLy)=tracer*vxm;

        *(cs+i*nLy)=sqrt(fabs(gamma*pre*aden)); 
        *(h+i*nLy)=fabs(gamma*pre*aden)/gm1+vtemp;
        *(vx+i*nLy)=vxm;
        *(vy+i*nLy)=vym;
        *(vz+i*nLy)=vzm;
        *(w+i*nLy)=sqrt(fabs(den));

    }//0~4

    for(i=5;i<9;i++)
	{
        ip2=i-5;
		ip1=i-4;
		i0=i-3;
		in1=i-2;
		in2=i-1;

		den=fetch(texro,i,texj);
		xmt=fetch(texvx,i,texj);
		ymt=fetch(texvy,i,texj);
		zmt=fetch(texvz,i,texj);
		eng=fetch(texeng,i,texj);
		tracer=fetch(textracer,i,texj);

		aden=1.0/den;   ///////////////////////////cuda dividing is slow?????????????????????
		vxm=xmt*aden;
		vym=ymt*aden;
		vzm=zmt*aden;
        vtemp=0.5*(vxm*vxm+vym*vym+vzm*vzm);
		pre=gm1*(eng-vtemp*den);
         
		*(fd+i*nLy)=xmt;
		*(fd+dist3dfd+i*nLy)=xmt*vxm+pre;
		*(fd+2*dist3dfd+i*nLy)=ymt*vxm;
		*(fd+3*dist3dfd+i*nLy)=zmt*vxm;
		*(fd+4*dist3dfd+i*nLy)=vxm*(pre+eng);
		*(fd+5*dist3dfd+i*nLy)=tracer*vxm;

        *(cs+i*nLy)=sqrt(fabs(gamma*pre*aden)); 
        *(h+i*nLy)=fabs(gamma*pre*aden)/gm1+vtemp;
        *(vx+i*nLy)=vxm;
        *(vy+i*nLy)=vym;
        *(vz+i*nLy)=vzm;
        *(w+i*nLy)=sqrt(fabs(den));
        ta1=fmax(fabs(*(vx+ip2*nLy)-*(cs+ip2*nLy)),fabs(*(vx+ip1*nLy)-*(cs+ip1*nLy)));
        ta2=fmax(fabs(*(vx+ip2*nLy)),fabs(*(vx+ip1*nLy)));
        ta3=fmax(fabs(*(vx+ip2*nLy)+*(cs+ip2*nLy)),fabs(*(vx+ip1*nLy)+*(cs+ip1*nLy)));
        ta4=fmax(fabs(*(vx+i0*nLy)-*(cs+i0*nLy)),fabs(*(vx+in1*nLy)-*(cs+in1*nLy)));
        ta5=fmax(fabs(*(vx+i0*nLy)),fabs(*(vx+in1*nLy)));
        ta6=fmax(fabs(*(vx+i0*nLy)+*(cs+i0*nLy)),fabs(*(vx+in1*nLy)+*(cs+in1*nLy)));
        ta1=fmax(fabs(*(vx+in2*nLy)-*(cs+in2*nLy)),ta1);
        ta2=fmax(fabs(*(vx+in2*nLy)),ta2);
        ta3=fmax(fabs(*(vx+in2*nLy)+*(cs+in2*nLy)),ta3);
        ta4=fmax(fabs(*(vx+i*nLy)-*(cs+i*nLy)),ta4);
        ta5=fmax(fabs(*(vx+i*nLy)),ta5);
        ta6=fmax(fabs(*(vx+i*nLy)+*(cs+i*nLy)),ta6);
		am1=fmax(ta1,ta4)*ama;
		am2=fmax(ta2,ta5)*ama;
		am5=fmax(ta3,ta6)*ama;///////////////////////////////////////////////////////////!!!!!!!!!!!!!!!!!!!!!!
        
		tg1=*(fd+in1*nLy)-*(fd+i0*nLy);
		tg2=*(fd+dist3dfd+in1*nLy)-*(fd+dist3dfd+i0*nLy);
		tg3=*(fd+2*dist3dfd+in1*nLy)-*(fd+2*dist3dfd+i0*nLy);
		tg4=*(fd+3*dist3dfd+in1*nLy)-*(fd+3*dist3dfd+i0*nLy);
		tg5=*(fd+4*dist3dfd+in1*nLy)-*(fd+4*dist3dfd+i0*nLy);
		tg6=*(fd+5*dist3dfd+in1*nLy)-*(fd+5*dist3dfd+i0*nLy);

		ttg1=0.5*(tg1+am1*(fetch(texro,in1,texj)-fetch(texro,i0,texj)));
		ttg2=0.5*(tg2+am2*(fetch(texvx,in1,texj)-fetch(texvx,i0,texj)));
		ttg3=0.5*(tg3+am2*(fetch(texvy,in1,texj)-fetch(texvy,i0,texj)));
		ttg4=0.5*(tg4+am2*(fetch(texvz,in1,texj)-fetch(texvz,i0,texj)));
		ttg5=0.5*(tg5+am5*(fetch(texeng,in1,texj)-fetch(texeng,i0,texj)));
		ttg6=0.5*(tg6+am2*(fetch(textracer,in1,texj)-fetch(textracer,i0,texj)));
	
		*(gg1+(i0-2)*nLy)=ttg1;
		*(gg1+dist3dgg+(i0-2)*nLy)=ttg2;
		*(gg1+2*dist3dgg+(i0-2)*nLy)=ttg3;
		*(gg1+3*dist3dgg+(i0-2)*nLy)=ttg4;
		*(gg1+4*dist3dgg+(i0-2)*nLy)=ttg5;
		*(gg1+5*dist3dgg+(i0-2)*nLy)=ttg6;

		*(gg2+(i0-2)*nLy)=ttg1-tg1;
		*(gg2+dist3dgg+(i0-2)*nLy)=ttg2-tg2;
		*(gg2+2*dist3dgg+(i0-2)*nLy)=ttg3-tg3;
		*(gg2+3*dist3dgg+(i0-2)*nLy)=ttg4-tg4;
		*(gg2+4*dist3dgg+(i0-2)*nLy)=ttg5-tg5;
		*(gg2+5*dist3dgg+(i0-2)*nLy)=ttg6-tg6;

    }//5~8

    for(i=9;i<LL;i++)
	{
		ip4=i-7;
		ip3=i-6;
        ip2=i-5;
		ip1=i-4;
		i0=i-3;
		in1=i-2;
		in2=i-1;

		den=fetch(texro,i,texj);
		xmt=fetch(texvx,i,texj);
		ymt=fetch(texvy,i,texj);
		zmt=fetch(texvz,i,texj);
		eng=fetch(texeng,i,texj);
		tracer=fetch(textracer,i,texj);

		aden=1.0/den;   ///////////////////////////cuda dividing is slow?????????????????????
		vxm=xmt*aden;
		vym=ymt*aden;
		vzm=zmt*aden;
        vtemp=0.5*(vxm*vxm+vym*vym+vzm*vzm);
		pre=gm1*(eng-vtemp*den);
         
		*(fd+i*nLy)=xmt;
		*(fd+dist3dfd+i*nLy)=xmt*vxm+pre;
		*(fd+2*dist3dfd+i*nLy)=ymt*vxm;
		*(fd+3*dist3dfd+i*nLy)=zmt*vxm;
		*(fd+4*dist3dfd+i*nLy)=vxm*(pre+eng);
		*(fd+5*dist3dfd+i*nLy)=tracer*vxm;

        *(cs+i*nLy)=sqrt(fabs(gamma*pre*aden)); 
        *(h+i*nLy)=fabs(gamma*pre*aden)/gm1+vtemp;
        *(vx+i*nLy)=vxm;
        *(vy+i*nLy)=vym;
        *(vz+i*nLy)=vzm;
        *(w+i*nLy)=sqrt(fabs(den));

        ta1=fmax(fabs(*(vx+ip2*nLy)-*(cs+ip2*nLy)),fabs(*(vx+ip1*nLy)-*(cs+ip1*nLy)));
        ta2=fmax(fabs(*(vx+ip2*nLy)),fabs(*(vx+ip1*nLy)));
        ta3=fmax(fabs(*(vx+ip2*nLy)+*(cs+ip2*nLy)),fabs(*(vx+ip1*nLy)+*(cs+ip1*nLy)));
        ta4=fmax(fabs(*(vx+i0*nLy)-*(cs+i0*nLy)),fabs(*(vx+in1*nLy)-*(cs+in1*nLy)));
        ta5=fmax(fabs(*(vx+i0*nLy)),fabs(*(vx+in1*nLy)));
        ta6=fmax(fabs(*(vx+i0*nLy)+*(cs+i0*nLy)),fabs(*(vx+in1*nLy)+*(cs+in1*nLy)));
        ta1=fmax(fabs(*(vx+in2*nLy)-*(cs+in2*nLy)),ta1);
        ta2=fmax(fabs(*(vx+in2*nLy)),ta2);
        ta3=fmax(fabs(*(vx+in2*nLy)+*(cs+in2*nLy)),ta3);
        ta4=fmax(fabs(*(vx+i*nLy)-*(cs+i*nLy)),ta4);
        ta5=fmax(fabs(*(vx+i*nLy)),ta5);
        ta6=fmax(fabs(*(vx+i*nLy)+*(cs+i*nLy)),ta6);
		am1=fmax(ta1,ta4)*ama;
		am2=fmax(ta2,ta5)*ama;
		am5=fmax(ta3,ta6)*ama;///////////////////////////////////////////////////////////!!!!!!!!!!!!!!!!!!!!!!
        
		tg1=*(fd+in1*nLy)-*(fd+i0*nLy);
		tg2=*(fd+dist3dfd+in1*nLy)-*(fd+dist3dfd+i0*nLy);
		tg3=*(fd+2*dist3dfd+in1*nLy)-*(fd+2*dist3dfd+i0*nLy);
		tg4=*(fd+3*dist3dfd+in1*nLy)-*(fd+3*dist3dfd+i0*nLy);
		tg5=*(fd+4*dist3dfd+in1*nLy)-*(fd+4*dist3dfd+i0*nLy);
		tg6=*(fd+5*dist3dfd+in1*nLy)-*(fd+5*dist3dfd+i0*nLy);

		ttg1=0.5*(tg1+am1*(fetch(texro,in1,texj)-fetch(texro,i0,texj)));
		ttg2=0.5*(tg2+am2*(fetch(texvx,in1,texj)-fetch(texvx,i0,texj)));
		ttg3=0.5*(tg3+am2*(fetch(texvy,in1,texj)-fetch(texvy,i0,texj)));
		ttg4=0.5*(tg4+am2*(fetch(texvz,in1,texj)-fetch(texvz,i0,texj)));
		ttg5=0.5*(tg5+am5*(fetch(texeng,in1,texj)-fetch(texeng,i0,texj)));
		ttg6=0.5*(tg6+am2*(fetch(textracer,in1,texj)-fetch(textracer,i0,texj)));
		*(gg1+(i0-2)*nLy)=ttg1;
		*(gg1+dist3dgg+(i0-2)*nLy)=ttg2;
		*(gg1+2*dist3dgg+(i0-2)*nLy)=ttg3;
		*(gg1+3*dist3dgg+(i0-2)*nLy)=ttg4;
		*(gg1+4*dist3dgg+(i0-2)*nLy)=ttg5;
		*(gg1+5*dist3dgg+(i0-2)*nLy)=ttg6;

		*(gg2+(i0-2)*nLy)=ttg1-tg1;
		*(gg2+dist3dgg+(i0-2)*nLy)=ttg2-tg2;
		*(gg2+2*dist3dgg+(i0-2)*nLy)=ttg3-tg3;
		*(gg2+3*dist3dgg+(i0-2)*nLy)=ttg4-tg4;
		*(gg2+4*dist3dgg+(i0-2)*nLy)=ttg5-tg5;
		*(gg2+5*dist3dgg+(i0-2)*nLy)=ttg6-tg6;
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
     //   te0=sqrt(fabs(*(w+ip2*nLy)))/(sqrt(fabs(*(w+ip2*nLy)))+sqrt(fabs(*(w+ip1*nLy))));
        te0=*(w+ip2*nLy)/(*(w+ip2*nLy)+*(w+ip1*nLy));

		te1=1-te0;
		te2=*(vx+ip2*nLy);
		te3=*(vy+ip2*nLy);
		te4=*(vz+ip2*nLy);
		te5=*(vx+ip1*nLy);
		te6=*(vy+ip1*nLy);
		te7=*(vz+ip1*nLy);
		tea=te0*te2+te1*te5;
		teb=te0*te3+te1*te6;
		tec=te0*te4+te1*te7;
        
		hm=te0*(*(h+ip2*nLy))+te1*(*(h+ip1*nLy));
		qm=0.5*(tea*tea+teb*teb+tec*tec);
		te8=gm1*qm;
        cm2=te0*pow(*(cs+ip2*nLy),2)+te1*(pow(*(cs+ip1*nLy),2))+0.5*te0*te1*gm1*((te2-te5)*(te2-te5)+(te3-te6)*(te3-te6)+(te4-te7)*(te4-te7));
		cm=sqrt(cm2);
		cmm=1.0/(3.0*cm2);
		te9=tea*cm;
		tea*=gm1;
		teb*=gm1;
		tec*=gm1;

		evr15=hm-te9;
		evr55=hm+te9;

		evl11=te8+te9;
		evl21=-(cm+tea);
		evl12=-teb/gm1*cm2;
		evl13=-tec/gm1*cm2;
		evl14=cm2-te8;
		evl15=te8-te9;
		evl25=cm-tea;

		tg1=*(gg1+(ip4-2)*nLy);
		tg2=*(gg1+dist3dgg+(ip4-2)*nLy);
		tg3=*(gg1+2*dist3dgg+(ip4-2)*nLy);
		tg4=*(gg1+3*dist3dgg+(ip4-2)*nLy);
		tg5=*(gg1+4*dist3dgg+(ip4-2)*nLy);
		tg6=*(gg1+5*dist3dgg+(ip4-2)*nLy);

	    den=(evl11*tg1+evl21*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		aden=evl12*tg1+cm2*tg3;
		ta1=evl13*tg1+cm2*tg4;
		ttg1=evl14*tg1+tea*tg2+teb*tg3+tec*tg4-gm1*tg5;
	    te0=(evl15*tg1+evl25*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		te6=tg6;

		tg1=*(gg1+(ip3-2)*nLy);
		tg2=*(gg1+dist3dgg+(ip3-2)*nLy);
		tg3=*(gg1+2*dist3dgg+(ip3-2)*nLy);
		tg4=*(gg1+3*dist3dgg+(ip3-2)*nLy);
		tg5=*(gg1+4*dist3dgg+(ip3-2)*nLy);
		tg6=*(gg1+5*dist3dgg+(ip3-2)*nLy);
	  
	    xmt=(evl11*tg1+evl21*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		pre=evl12*tg1+cm2*tg3;
		ta2=evl13*tg1+cm2*tg4;
		ttg2=evl14*tg1+tea*tg2+teb*tg3+tec*tg4-gm1*tg5;
	    te1=(evl15*tg1+evl25*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		te7=tg6;
		
		tg1=*(gg1+(ip2-2)*nLy);
		tg2=*(gg1+dist3dgg+(ip2-2)*nLy);
		tg3=*(gg1+2*dist3dgg+(ip2-2)*nLy);
		tg4=*(gg1+3*dist3dgg+(ip2-2)*nLy);
		tg5=*(gg1+4*dist3dgg+(ip2-2)*nLy);
		tg6=*(gg1+5*dist3dgg+(ip2-2)*nLy);

	    ymt=(evl11*tg1+evl21*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		vtemp=evl12*tg1+cm2*tg3;
		ta3=evl13*tg1+cm2*tg4;
		ttg3=evl14*tg1+tea*tg2+teb*tg3+tec*tg4-gm1*tg5;
	    te2=(evl15*tg1+evl25*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		te8=tg6;

		tg1=*(gg1+(ip1-2)*nLy);
		tg2=*(gg1+dist3dgg+(ip1-2)*nLy);
		tg3=*(gg1+2*dist3dgg+(ip1-2)*nLy);
		tg4=*(gg1+3*dist3dgg+(ip1-2)*nLy);
		tg5=*(gg1+4*dist3dgg+(ip1-2)*nLy);
		tg6=*(gg1+5*dist3dgg+(ip1-2)*nLy);

	    zmt=(evl11*tg1+evl21*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		vxm=evl12*tg1+cm2*tg3;
		ta4=evl13*tg1+cm2*tg4;
		ttg4=evl14*tg1+tea*tg2+teb*tg3+tec*tg4-gm1*tg5;
	    te3=(evl15*tg1+evl25*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		te9=tg6;

		tg1=*(gg2+(i0-2)*nLy);
		tg2=*(gg2+dist3dgg+(i0-2)*nLy);
		tg3=*(gg2+2*dist3dgg+(i0-2)*nLy);
		tg4=*(gg2+3*dist3dgg+(i0-2)*nLy);
		tg5=*(gg2+4*dist3dgg+(i0-2)*nLy);
		tg6=*(gg2+5*dist3dgg+(i0-2)*nLy);

	    eng=(evl11*tg1+evl21*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		vym=evl12*tg1+cm2*tg3;
		ta5=evl13*tg1+cm2*tg4;
		ttg5=evl14*tg1+tea*tg2+teb*tg3+tec*tg4-gm1*tg5;
	    te4=(evl15*tg1+evl25*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		am1=tg6;

		tg1=*(gg2+(ip1-2)*nLy);
		tg2=*(gg2+dist3dgg+(ip1-2)*nLy);
		tg3=*(gg2+2*dist3dgg+(ip1-2)*nLy);
		tg4=*(gg2+3*dist3dgg+(ip1-2)*nLy);
		tg5=*(gg2+4*dist3dgg+(ip1-2)*nLy);
		tg6=*(gg2+5*dist3dgg+(ip1-2)*nLy);

	    tracer=(evl11*tg1+evl21*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		vzm=evl12*tg1+cm2*tg3;
		ta6=evl13*tg1+cm2*tg4;
		ttg6=evl14*tg1+tea*tg2+teb*tg3+tec*tg4-gm1*tg5;
	    te5=(evl15*tg1+evl25*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		am2=tg6;

		tg1=*(gg2+(ip2-2)*nLy);
		tg2=*(gg2+dist3dgg+(ip2-2)*nLy);
		tg3=*(gg2+2*dist3dgg+(ip2-2)*nLy);
		tg4=*(gg2+3*dist3dgg+(ip2-2)*nLy);
		tg5=*(gg2+4*dist3dgg+(ip2-2)*nLy);
		tg6=*(gg2+5*dist3dgg+(ip2-2)*nLy);

	    h17=(evl11*tg1+evl21*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h27=evl12*tg1+cm2*tg3;
		h37=evl13*tg1+cm2*tg4;
		h47=evl14*tg1+tea*tg2+teb*tg3+tec*tg4-gm1*tg5;
	    h57=(evl15*tg1+evl25*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h67=tg6;

		tg1=*(gg2+(ip3-2)*nLy);
		tg2=*(gg2+dist3dgg+(ip3-2)*nLy);
		tg3=*(gg2+2*dist3dgg+(ip3-2)*nLy);
		tg4=*(gg2+3*dist3dgg+(ip3-2)*nLy);
		tg5=*(gg2+4*dist3dgg+(ip3-2)*nLy);
		tg6=*(gg2+5*dist3dgg+(ip3-2)*nLy);

	    h18=(evl11*tg1+evl21*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h28=evl12*tg1+cm2*tg3;
		h38=evl13*tg1+cm2*tg4;
		h48=evl14*tg1+tea*tg2+teb*tg3+tec*tg4-gm1*tg5;
	    h58=(evl15*tg1+evl25*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h68=tg6;

		ts11=den-xmt;
		ts12=xmt-ymt;
		ts13=ymt-zmt;
		ts14=eng-tracer;
		ts15=tracer-h17;
		ts16=h17-h18;

		tt11=13.0*ts11*ts11+3.0*(den-3.0*xmt)*(den-3.0*xmt)+epweno;
		tt12=13.0*ts12*ts12+3.0*(xmt+ymt)*(xmt+ymt)+epweno;
		tt13=13.0*ts13*ts13+3.0*(3.0*ymt-zmt)*(3.0*ymt-zmt)+epweno;
		tt14=13.0*ts14*ts14+3.0*(eng-3.0*tracer)*(eng-3.0*tracer)+epweno;
		tt15=13.0*ts15*ts15+3.0*(tracer+h17)*(tracer+h17)+epweno;
		tt16=13.0*ts16*ts16+3.0*(3.0*h17-h18)*(3.0*h17-h18)+epweno;

		den=(ts12-ts11)/(1.0+6.0*(tt11/tt12)*(tt11/tt12)+3.0*(tt11/tt13)*(tt11/tt13));
		xmt=(ts13-ts12)/(2.0+(tt13/tt11)*(tt13/tt11)/1.5+4.0*(tt13/tt12)*(tt13/tt12))-0.25*(ts13-ts12);
		ymt=(ts15-ts14)/(1.0+6.0*(tt14/tt15)*(tt14/tt15)+3.0*(tt14/tt16)*(tt14/tt16));
		zmt=(ts16-ts15)/(2.0+(tt16/tt14)*(tt16/tt14)/1.5+4.0*(tt16/tt15)*(tt16/tt15))-0.25*(ts16-ts15);
		eng=(den+xmt+ymt+zmt)*cmm;

		ts11=aden-pre;
		ts12=pre-vtemp;
		ts13=vtemp-vxm;
		ts14=vym-vzm;
		ts15=vzm-h27;
		ts16=h27-h28;

		tt11=13.0*ts11*ts11+3.0*(aden-3.0*pre)*(aden-3.0*pre)+epweno;
		tt12=13.0*ts12*ts12+3.0*(pre+vtemp)*(pre+vtemp)+epweno;
		tt13=13.0*ts13*ts13+3.0*(3.0*vtemp-vxm)*(3.0*vtemp-vxm)+epweno;
		tt14=13.0*ts14*ts14+3.0*(vym-3.0*vzm)*(vym-3.0*vzm)+epweno;
		tt15=13.0*ts15*ts15+3.0*(vzm+h27)*(vzm+h27)+epweno;
		tt16=13.0*ts16*ts16+3.0*(3.0*h27-h28)*(3.0*h27-h28)+epweno;

		aden=(ts12-ts11)/(1.0+6.0*(tt11/tt12)*(tt11/tt12)+3.0*(tt11/tt13)*(tt11/tt13));
		pre=(ts13-ts12)/(2.0+(tt13/tt11)*(tt13/tt11)/1.5+4.0*(tt13/tt12)*(tt13/tt12))-0.25*(ts13-ts12);
		vtemp=(ts15-ts14)/(1.0+6.0*(tt14/tt15)*(tt14/tt15)+3.0*(tt14/tt16)*(tt14/tt16));
		vxm=(ts16-ts15)/(2.0+(tt16/tt14)*(tt16/tt14)/1.5+4.0*(tt16/tt15)*(tt16/tt15))-0.25*(ts16-ts15);
		vym=(aden+pre+vtemp+vxm)*cmm;

		ts11=ta1-ta2;
		ts12=ta2-ta3;
		ts13=ta3-ta4;
		ts14=ta5-ta6;
		ts15=ta6-h37;
		ts16=h37-h38;

		tt11=13.0*ts11*ts11+3.0*(ta1-3.0*ta2)*(ta1-3.0*ta2)+epweno;
		tt12=13.0*ts12*ts12+3.0*(ta2+ta3)*(ta2+ta3)+epweno;
		tt13=13.0*ts13*ts13+3.0*(3.0*ta3-ta4)*(3.0*ta3-ta4)+epweno;
		tt14=13.0*ts14*ts14+3.0*(ta5-3.0*ta6)*(ta5-3.0*ta6)+epweno;
		tt15=13.0*ts15*ts15+3.0*(ta6+h37)*(ta6+h37)+epweno;
		tt16=13.0*ts16*ts16+3.0*(3.0*h37-h38)*(3.0*h37-h38)+epweno;

		ta1=(ts12-ts11)/(1.0+6.0*(tt11/tt12)*(tt11/tt12)+3.0*(tt11/tt13)*(tt11/tt13));
		ta2=(ts13-ts12)/(2.0+(tt13/tt11)*(tt13/tt11)/1.5+4.0*(tt13/tt12)*(tt13/tt12))-0.25*(ts13-ts12);
		ta3=(ts15-ts14)/(1.0+6.0*(tt14/tt15)*(tt14/tt15)+3.0*(tt14/tt16)*(tt14/tt16));
		ta4=(ts16-ts15)/(2.0+(tt16/tt14)*(tt16/tt14)/1.5+4.0*(tt16/tt15)*(tt16/tt15))-0.25*(ts16-ts15);
		ta5=(ta1+ta2+ta3+ta4)*cmm;

		ts11=ttg1-ttg2;
		ts12=ttg2-ttg3;
		ts13=ttg3-ttg4;
		ts14=ttg5-ttg6;
		ts15=ttg6-h47;
		ts16=h47-h48;

		tt11=13.0*ts11*ts11+3.0*(ttg1-3.0*ttg2)*(ttg1-3.0*ttg2)+epweno;
		tt12=13.0*ts12*ts12+3.0*(ttg2+ttg3)*(ttg2+ttg3)+epweno;
		tt13=13.0*ts13*ts13+3.0*(3.0*ttg3-ttg4)*(3.0*ttg3-ttg4)+epweno;
		tt14=13.0*ts14*ts14+3.0*(ttg5-3.0*ttg6)*(ttg5-3.0*ttg6)+epweno;
		tt15=13.0*ts15*ts15+3.0*(ttg6+h47)*(ttg6+h47)+epweno;
		tt16=13.0*ts16*ts16+3.0*(3.0*h47-h48)*(3.0*h47-h48)+epweno;

		ttg1=(ts12-ts11)/(1.0+6.0*(tt11/tt12)*(tt11/tt12)+3.0*(tt11/tt13)*(tt11/tt13));
		ttg2=(ts13-ts12)/(2.0+(tt13/tt11)*(tt13/tt11)/1.5+4.0*(tt13/tt12)*(tt13/tt12))-0.25*(ts13-ts12);
		ttg3=(ts15-ts14)/(1.0+6.0*(tt14/tt15)*(tt14/tt15)+3.0*(tt14/tt16)*(tt14/tt16));
		ttg4=(ts16-ts15)/(2.0+(tt16/tt14)*(tt16/tt14)/1.5+4.0*(tt16/tt15)*(tt16/tt15))-0.25*(ts16-ts15);
		ttg5=(ttg1+ttg2+ttg3+ttg4)*cmm;

		ts11=te0-te1;
		ts12=te1-te2;
		ts13=te2-te3;
		ts14=te4-te5;
		ts15=te5-h57;
		ts16=h57-h58;

		tt11=13.0*ts11*ts11+3.0*(te0-3.0*te1)*(te0-3.0*te1)+epweno;
		tt12=13.0*ts12*ts12+3.0*(te1+te2)*(te1+te2)+epweno;
		tt13=13.0*ts13*ts13+3.0*(3.0*te2-te3)*(3.0*te2-te3)+epweno;
		tt14=13.0*ts14*ts14+3.0*(te4-3.0*te5)*(te4-3.0*te5)+epweno;
		tt15=13.0*ts15*ts15+3.0*(te5+h57)*(te5+h57)+epweno;
		tt16=13.0*ts16*ts16+3.0*(3.0*h57-h58)*(3.0*h57-h58)+epweno;

		te0=(ts12-ts11)/(1.0+6.0*(tt11/tt12)*(tt11/tt12)+3.0*(tt11/tt13)*(tt11/tt13));
		te1=(ts13-ts12)/(2.0+(tt13/tt11)*(tt13/tt11)/1.5+4.0*(tt13/tt12)*(tt13/tt12))-0.25*(ts13-ts12);
		te2=(ts15-ts14)/(1.0+6.0*(tt14/tt15)*(tt14/tt15)+3.0*(tt14/tt16)*(tt14/tt16));
		te3=(ts16-ts15)/(2.0+(tt16/tt14)*(tt16/tt14)/1.5+4.0*(tt16/tt15)*(tt16/tt15))-0.25*(ts16-ts15);
		te4=(te0+te1+te2+te3)*cmm;

		ts11=te6-te7;
		ts12=te7-te8;
		ts13=te8-te9;
		ts14=am1-am2;
		ts15=am2-h67;
		ts16=h67-h68;

		tt11=13.0*ts11*ts11+3.0*(te6-3.0*te7)*(te6-3.0*te7)+epweno;
		tt12=13.0*ts12*ts12+3.0*(te7+te8)*(te7+te8)+epweno;
		tt13=13.0*ts13*ts13+3.0*(3.0*te8-te9)*(3.0*te8-te9)+epweno;
		tt14=13.0*ts14*ts14+3.0*(am1-3.0*am2)*(am1-3.0*am2)+epweno;
		tt15=13.0*ts15*ts15+3.0*(am2+h67)*(am2+h67)+epweno;
		tt16=13.0*ts16*ts16+3.0*(3.0*h67-h68)*(3.0*h67-h68)+epweno;

		te6=(ts12-ts11)/(3.0+18.0*(tt11/tt12)*(tt11/tt12)+9.0*(tt11/tt13)*(tt11/tt13));
		te7=(ts13-ts12)/(6.0+(tt13/tt11)*(tt13/tt11)*2.0+12.0*(tt13/tt12)*(tt13/tt12))-0.25*(ts13-ts12)/3.0;
		te8=(ts15-ts14)/(3.0+18.0*(tt14/tt15)*(tt14/tt15)+9.0*(tt14/tt16)*(tt14/tt16));
		te9=(ts16-ts15)/(6.0+(tt16/tt14)*(tt16/tt14)*2.0+12.0*(tt16/tt15)*(tt16/tt15))-0.25*(ts16-ts15)/3.0;
		am1=te6+te7+te8+te9;

		tfh1=eng+ttg5+te4;

		*(fh+(ip2-4)*nLy)=(-*(fd+ip3*nLy)+7.0*(*(fd+ip2*nLy)+*(fd+ip1*nLy))-*(fd+i0*nLy))/12.0+tfh1;
		*(fh+dist3dfh+(ip2-4)*nLy)=(-*(fd+dist3dfd+ip3*nLy)+7.0*(*(fd+dist3dfd+ip2*nLy)+*(fd+dist3dfd+ip1*nLy))-*(fd+dist3dfd+i0*nLy))/12.0+tfh1*tea/gm1-(eng-te4)*cm;
		*(fh+2*dist3dfh+(ip2-4)*nLy)=(-*(fd+2*dist3dfd+ip3*nLy)+7.0*(*(fd+2*dist3dfd+ip2*nLy)+*(fd+2*dist3dfd+ip1*nLy))-*(fd+2*dist3dfd+i0*nLy))/12.0+tfh1*teb/gm1+vym;
		*(fh+3*dist3dfh+(ip2-4)*nLy)=(-*(fd+3*dist3dfd+ip3*nLy)+7.0*(*(fd+3*dist3dfd+ip2*nLy)+*(fd+3*dist3dfd+ip1*nLy))-*(fd+3*dist3dfd+i0*nLy))/12.0+tfh1*tec/gm1+ta5;
		*(fh+4*dist3dfh+(ip2-4)*nLy)=(-*(fd+4*dist3dfd+ip3*nLy)+7.0*(*(fd+4*dist3dfd+ip2*nLy)+*(fd+4*dist3dfd+ip1*nLy))-*(fd+4*dist3dfd+i0*nLy))/12.0+evr15*eng+teb/gm1*vym+tec/gm1*ta5+qm*ttg5+evr55*te4;
		*(fh+5*dist3dfh+(ip2-4)*nLy)=(-*(fd+5*dist3dfd+ip3*nLy)+7.0*(*(fd+5*dist3dfd+ip2*nLy)+*(fd+5*dist3dfd+ip1*nLy))-*(fd+5*dist3dfd+i0*nLy))/12.0+am1;
    }//9~nLx+9
}

__global__ void gy_gpu(double* rhs,int nLx,int nLy,int nLz,int mnp,double gamma,double* fd,double* cs,double* vx,double* vy,double* vz,double* w,double* h,double* gg1,double* gg2,double* fh,cudaTextureObject_t texro,cudaTextureObject_t texvx,cudaTextureObject_t texvy,cudaTextureObject_t texvz,cudaTextureObject_t texeng,cudaTextureObject_t textracer)
{
	int i;
	double gm1=gamma-1;
	double den,xmt,ymt,zmt,eng,tracer;
	double am1,am2,am5;
	double aden,pre,vtemp,vxm,vym,vzm;
	int ip4,ip3,ip2,ip1,i0,in1,in2;
	double ta1,ta2,ta3,ta4,ta5,ta6;
	double tg1,tg2,tg3,tg4,tg5,tg6;
	double ttg1,ttg2,ttg3,ttg4,ttg5,ttg6;
	double te0,te1,te2,te3,te4,te5,te6,te7,te8,te9,tea,teb,tec;
	double evr15,evr55,evl11,evl21,evl12,evl13,evl14,evl15,evl25;
	double hm,qm,cm,cm2,cmm;
	double h11,h21,h31,h41,h51,h61;
	double h12,h22,h32,h42,h52,h62;
	double h13,h23,h33,h43,h53,h63;
	double h14,h24,h34,h44,h54,h64;
	double h15,h25,h35,h45,h55,h65;
	double h16,h26,h36,h46,h56,h66;
	double h17,h27,h37,h47,h57,h67;
	double h18,h28,h38,h48,h58,h68;

	double tt11,tt12,tt13,tt14,tt15,tt16;
	double ts11,ts12,ts13,ts14,ts15,ts16;

	double tfh1;

	int LL=nLy+10,LL1=nLy+5,LL2=nLy+1;
//	int LLx=nLx+10;
//	int LLz=nLz+10;
	int dist3dfd=LL*nLx*nLz,dist3dgg=LL1*nLx*nLz,dist3dfh=LL2*nLx*nLz,dist3drhs=nLx*nLy*nLz;
	int tx=threadIdx.x;
	int ty=threadIdx.y;
	int bx=blockIdx.x;
	int by=blockIdx.y;
	int xid=bx*BLOCK_SIZE_x+tx;
	int yid=by*BLOCK_SIZE_yg+ty;
//	int id=yid*nLx+xid;
//	uc+=(yid+5)*LL*LLx+xid+5;
	fd+=yid*nLx*LL+xid;
    cs+=yid*nLx*LL+xid;
    vx+=yid*nLx*LL+xid;
    vy+=yid*nLx*LL+xid;
    vz+=yid*nLx*LL+xid;
    w+=yid*nLx*LL+xid;
    h+=yid*nLx*LL+xid;
	gg1+=yid*nLx*LL1+xid;
	gg2+=yid*nLx*LL1+xid;
    fh+=yid*nLx*LL2+xid;
	rhs+=yid*nLx*nLy+xid;
	float texj=(yid+5)*LL+0.5f;
	float texi=xid+5+0.5f;
    for(i=0;i<5;i++)
	{
		den=fetch(texro,texi,texj+i);
		xmt=fetch(texvy,texi,texj+i);
		ymt=fetch(texvx,texi,texj+i);
		zmt=fetch(texvz,texi,texj+i);
		eng=fetch(texeng,texi,texj+i);
		tracer=fetch(textracer,texi,texj+i);

		aden=1.0/den;   ///////////////////////////cuda dividing is slow?????????????????????
		vxm=xmt*aden;
		vym=ymt*aden;
		vzm=zmt*aden;
        vtemp=0.5*(vxm*vxm+vym*vym+vzm*vzm);
		pre=gm1*(eng-vtemp*den);
         
		*(fd+i*nLx)=xmt;
		*(fd+dist3dfd+i*nLx)=xmt*vxm+pre;
		*(fd+2*dist3dfd+i*nLx)=ymt*vxm;
		*(fd+3*dist3dfd+i*nLx)=zmt*vxm;
		*(fd+4*dist3dfd+i*nLx)=vxm*(pre+eng);
		*(fd+5*dist3dfd+i*nLx)=tracer*vxm;

        *(cs+i*nLx)=sqrt(fabs(gamma*pre*aden)); 
        *(vx+i*nLx)=vxm;
        *(vy+i*nLx)=vym;
        *(vz+i*nLx)=vzm;
        *(w+i*nLx)=sqrt(fabs(den));
        *(h+i*nLx)=fabs(gamma*pre*aden)/gm1+vtemp;

    }//0~4

    for(i=5;i<9;i++)
	{
        ip2=i-5;
		ip1=i-4;
		i0=i-3;
		in1=i-2;
		in2=i-1;

		den=fetch(texro,texi,texj+i);///ro
		xmt=fetch(texvy,texi,texj+i);///x
		ymt=fetch(texvx,texi,texj+i);///y
		zmt=fetch(texvz,texi,texj+i);///z
		eng=fetch(texeng,texi,texj+i);///e
		tracer=fetch(textracer,texi,texj+i);///t

		aden=1.0/den;   ///ro////////////////////////cuda dividing is slow?????????????????????
		vxm=xmt*aden; ///x
		vym=ymt*aden; ///y
		vzm=zmt*aden; ///
        vtemp=0.5*(vxm*vxm+vym*vym+vzm*vzm);
		pre=gm1*(eng-vtemp*den);
         
		*(fd+i*nLx)=xmt;
		*(fd+dist3dfd+i*nLx)=xmt*vxm+pre;
		*(fd+2*dist3dfd+i*nLx)=ymt*vxm;
		*(fd+3*dist3dfd+i*nLx)=zmt*vxm;
		*(fd+4*dist3dfd+i*nLx)=vxm*(pre+eng);
		*(fd+5*dist3dfd+i*nLx)=tracer*vxm;

        *(cs+i*nLx)=sqrt(fabs(gamma*pre*aden)); 
        *(w+i*nLx)=sqrt(fabs(den));
        *(vx+i*nLx)=vxm;
        *(vy+i*nLx)=vym;
        *(vz+i*nLx)=vzm;
        *(h+i*nLx)=fabs(gamma*pre*aden)/gm1+vtemp;

        ta1=fmax(fabs(*(vx+ip2*nLx)-*(cs+ip2*nLx)),fabs(*(vx+ip1*nLx)-*(cs+ip1*nLx)));
        ta2=fmax(fabs(*(vx+ip2*nLx)),fabs(*(vx+ip1*nLx)));
        ta3=fmax(fabs(*(vx+ip2*nLx)+*(cs+ip2*nLx)),fabs(*(vx+ip1*nLx)+*(cs+ip1*nLx)));
        ta4=fmax(fabs(*(vx+i0*nLx)-*(cs+i0*nLx)),fabs(*(vx+in1*nLx)-*(cs+in1*nLx)));
        ta5=fmax(fabs(*(vx+i0*nLx)),fabs(*(vx+in1*nLx)));
        ta6=fmax(fabs(*(vx+i0*nLx)+*(cs+i0*nLx)),fabs(*(vx+in1*nLx)+*(cs+in1*nLx)));
        ta1=fmax(fabs(*(vx+in2*nLx)-*(cs+in2*nLx)),ta1);
        ta2=fmax(fabs(*(vx+in2*nLx)),ta2);
        ta3=fmax(fabs(*(vx+in2*nLx)+*(cs+in2*nLx)),ta3);
        ta4=fmax(fabs(*(vx+i*nLx)-*(cs+i*nLx)),ta4);
        ta5=fmax(fabs(*(vx+i*nLx)),ta5);
        ta6=fmax(fabs(*(vx+i*nLx)+*(cs+i*nLx)),ta6);
		am1=fmax(ta1,ta4)*ama;
		am2=fmax(ta2,ta5)*ama;
		am5=fmax(ta3,ta6)*ama;///////////////////////////////////////////////////////////!!!!!!!!!!!!!!!!!!!!!!
        
		tg1=*(fd+in1*nLx)-*(fd+i0*nLx);
		tg2=*(fd+dist3dfd+in1*nLx)-*(fd+dist3dfd+i0*nLx);
		tg3=*(fd+2*dist3dfd+in1*nLx)-*(fd+2*dist3dfd+i0*nLx);
		tg4=*(fd+3*dist3dfd+in1*nLx)-*(fd+3*dist3dfd+i0*nLx);
		tg5=*(fd+4*dist3dfd+in1*nLx)-*(fd+4*dist3dfd+i0*nLx);
		tg6=*(fd+5*dist3dfd+in1*nLx)-*(fd+5*dist3dfd+i0*nLx);

		ttg1=0.5*(tg1+am1*(fetch(texro,texi,texj+in1)-fetch(texro,texi,texj+i0)));
		ttg2=0.5*(tg2+am2*(fetch(texvy,texi,texj+in1)-fetch(texvy,texi,texj+i0)));
		ttg3=0.5*(tg3+am2*(fetch(texvx,texi,texj+in1)-fetch(texvx,texi,texj+i0)));
		ttg4=0.5*(tg4+am2*(fetch(texvz,texi,texj+in1)-fetch(texvz,texi,texj+i0)));
		ttg5=0.5*(tg5+am5*(fetch(texeng,texi,texj+in1)-fetch(texeng,texi,texj+i0)));
		ttg6=0.5*(tg6+am2*(fetch(textracer,texi,texj+in1)-fetch(textracer,texi,texj+i0)));
		*(gg1+(i0-2)*nLx)=ttg1;
		*(gg1+dist3dgg+(i0-2)*nLx)=ttg2;
		*(gg1+2*dist3dgg+(i0-2)*nLx)=ttg3;
		*(gg1+3*dist3dgg+(i0-2)*nLx)=ttg4;
		*(gg1+4*dist3dgg+(i0-2)*nLx)=ttg5;
		*(gg1+5*dist3dgg+(i0-2)*nLx)=ttg6;

		*(gg2+(i0-2)*nLx)=ttg1-tg1;
		*(gg2+dist3dgg+(i0-2)*nLx)=ttg2-tg2;
		*(gg2+2*dist3dgg+(i0-2)*nLx)=ttg3-tg3;
		*(gg2+3*dist3dgg+(i0-2)*nLx)=ttg4-tg4;
		*(gg2+4*dist3dgg+(i0-2)*nLx)=ttg5-tg5;
		*(gg2+5*dist3dgg+(i0-2)*nLx)=ttg6-tg6;

    }//5~8

    for(i=9;i<LL;i++)
	{
		ip4=i-7;
		ip3=i-6;
        ip2=i-5;
		ip1=i-4;
		i0=i-3;
		in1=i-2;
		in2=i-1;

		den=fetch(texro,texi,texj+i);///ro
		xmt=fetch(texvy,texi,texj+i);///x
		ymt=fetch(texvx,texi,texj+i);///y
		zmt=fetch(texvz,texi,texj+i);///z
		eng=fetch(texeng,texi,texj+i);///e
		tracer=fetch(textracer,texi,texj+i);///t

		aden=1.0/den;   ///ro////////////////////////cuda dividing is slow?????????????????????
		vxm=xmt*aden; ///x
		vym=ymt*aden; ///y
		vzm=zmt*aden; ///
        vtemp=0.5*(vxm*vxm+vym*vym+vzm*vzm);
		pre=gm1*(eng-vtemp*den);
         
		*(fd+i*nLx)=xmt;
		*(fd+dist3dfd+i*nLx)=xmt*vxm+pre;
		*(fd+2*dist3dfd+i*nLx)=ymt*vxm;
		*(fd+3*dist3dfd+i*nLx)=zmt*vxm;
		*(fd+4*dist3dfd+i*nLx)=vxm*(pre+eng);
		*(fd+5*dist3dfd+i*nLx)=tracer*vxm;

        *(cs+i*nLx)=sqrt(fabs(gamma*pre*aden)); 
        *(w+i*nLx)=sqrt(fabs(den));
        *(vx+i*nLx)=vxm;
        *(vy+i*nLx)=vym;
        *(vz+i*nLx)=vzm;
        *(h+i*nLx)=fabs(gamma*pre*aden)/gm1+vtemp;

        ta1=fmax(fabs(*(vx+ip2*nLx)-*(cs+ip2*nLx)),fabs(*(vx+ip1*nLx)-*(cs+ip1*nLx)));
        ta2=fmax(fabs(*(vx+ip2*nLx)),fabs(*(vx+ip1*nLx)));
        ta3=fmax(fabs(*(vx+ip2*nLx)+*(cs+ip2*nLx)),fabs(*(vx+ip1*nLx)+*(cs+ip1*nLx)));
        ta4=fmax(fabs(*(vx+i0*nLx)-*(cs+i0*nLx)),fabs(*(vx+in1*nLx)-*(cs+in1*nLx)));
        ta5=fmax(fabs(*(vx+i0*nLx)),fabs(*(vx+in1*nLx)));
        ta6=fmax(fabs(*(vx+i0*nLx)+*(cs+i0*nLx)),fabs(*(vx+in1*nLx)+*(cs+in1*nLx)));
        ta1=fmax(fabs(*(vx+in2*nLx)-*(cs+in2*nLx)),ta1);
        ta2=fmax(fabs(*(vx+in2*nLx)),ta2);
        ta3=fmax(fabs(*(vx+in2*nLx)+*(cs+in2*nLx)),ta3);
        ta4=fmax(fabs(*(vx+i*nLx)-*(cs+i*nLx)),ta4);
        ta5=fmax(fabs(*(vx+i*nLx)),ta5);
        ta6=fmax(fabs(*(vx+i*nLx)+*(cs+i*nLx)),ta6);
		am1=fmax(ta1,ta4)*ama;
		am2=fmax(ta2,ta5)*ama;
		am5=fmax(ta3,ta6)*ama;///////////////////////////////////////////////////////////!!!!!!!!!!!!!!!!!!!!!!
        
		tg1=*(fd+in1*nLx)-*(fd+i0*nLx);
		tg2=*(fd+dist3dfd+in1*nLx)-*(fd+dist3dfd+i0*nLx);
		tg3=*(fd+2*dist3dfd+in1*nLx)-*(fd+2*dist3dfd+i0*nLx);
		tg4=*(fd+3*dist3dfd+in1*nLx)-*(fd+3*dist3dfd+i0*nLx);
		tg5=*(fd+4*dist3dfd+in1*nLx)-*(fd+4*dist3dfd+i0*nLx);
		tg6=*(fd+5*dist3dfd+in1*nLx)-*(fd+5*dist3dfd+i0*nLx);

		ttg1=0.5*(tg1+am1*(fetch(texro,texi,texj+in1)-fetch(texro,texi,texj+i0)));
		ttg2=0.5*(tg2+am2*(fetch(texvy,texi,texj+in1)-fetch(texvy,texi,texj+i0)));
		ttg3=0.5*(tg3+am2*(fetch(texvx,texi,texj+in1)-fetch(texvx,texi,texj+i0)));
		ttg4=0.5*(tg4+am2*(fetch(texvz,texi,texj+in1)-fetch(texvz,texi,texj+i0)));
		ttg5=0.5*(tg5+am5*(fetch(texeng,texi,texj+in1)-fetch(texeng,texi,texj+i0)));
		ttg6=0.5*(tg6+am2*(fetch(textracer,texi,texj+in1)-fetch(textracer,texi,texj+i0)));
		*(gg1+(i0-2)*nLx)=ttg1;
		*(gg1+dist3dgg+(i0-2)*nLx)=ttg2;
		*(gg1+2*dist3dgg+(i0-2)*nLx)=ttg3;
		*(gg1+3*dist3dgg+(i0-2)*nLx)=ttg4;
		*(gg1+4*dist3dgg+(i0-2)*nLx)=ttg5;
		*(gg1+5*dist3dgg+(i0-2)*nLx)=ttg6;

		*(gg2+(i0-2)*nLx)=ttg1-tg1;
		*(gg2+dist3dgg+(i0-2)*nLx)=ttg2-tg2;
		*(gg2+2*dist3dgg+(i0-2)*nLx)=ttg3-tg3;
		*(gg2+3*dist3dgg+(i0-2)*nLx)=ttg4-tg4;
		*(gg2+4*dist3dgg+(i0-2)*nLx)=ttg5-tg5;
		*(gg2+5*dist3dgg+(i0-2)*nLx)=ttg6-tg6;
        te0=*(w+ip2*nLx)/(*(w+ip2*nLx)+*(w+ip1*nLx)); //////////////////////////////weight 0

		te1=1-te0;                                    //////////////////////////////weight 1
		te2=*(vx+ip2*nLx);
		te3=*(vy+ip2*nLx);
		te4=*(vz+ip2*nLx);
		te5=*(vx+ip1*nLx);
		te6=*(vy+ip1*nLx);
		te7=*(vz+ip1*nLx);                             /////////////////////////////vx,vy,vz
		tea=te0*te2+te1*te5;                
		teb=te0*te3+te1*te6;
		tec=te0*te4+te1*te7;                           /////////////middle vx,middle vy,middle vz
        
		hm=te0*(*(h+ip2*nLx))+te1*(*(h+ip1*nLx));        //////middle h
		qm=0.5*(tea*tea+teb*teb+tec*tec);                
		te8=gm1*qm;
        cm2=te0*pow(*(cs+ip2*nLx),2)+te1*(pow(*(cs+ip1*nLx),2))+0.5*te0*te1*gm1*((te2-te5)*(te2-te5)+(te3-te6)*(te3-te6)+(te4-te7)*(te4-te7));
		cm=sqrt(cm2);
		cmm=1.0/(3.0*cm2);
		te9=tea*cm;
		tea*=gm1;
		teb*=gm1;
		tec*=gm1;

		evr15=hm-te9;
		evr55=hm+te9;

		evl11=te8+te9;
		evl21=-(cm+tea);
		evl12=-teb/gm1*cm2;
		evl13=-tec/gm1*cm2;
		evl14=cm2-te8;
		evl15=te8-te9;
		evl25=cm-tea;

		tg1=*(gg1+(ip4-2)*nLx);
		tg2=*(gg1+dist3dgg+(ip4-2)*nLx);
		tg3=*(gg1+2*dist3dgg+(ip4-2)*nLx);
		tg4=*(gg1+3*dist3dgg+(ip4-2)*nLx);
		tg5=*(gg1+4*dist3dgg+(ip4-2)*nLx);
		tg6=*(gg1+5*dist3dgg+(ip4-2)*nLx);

	    h11=(evl11*tg1+evl21*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h21=evl12*tg1+cm2*tg3;
		h31=evl13*tg1+cm2*tg4;
		h41=evl14*tg1+tea*tg2+teb*tg3+tec*tg4-gm1*tg5;
	    h51=(evl15*tg1+evl25*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h61=tg6;

		tg1=*(gg1+(ip3-2)*nLx);
		tg2=*(gg1+dist3dgg+(ip3-2)*nLx);
		tg3=*(gg1+2*dist3dgg+(ip3-2)*nLx);
		tg4=*(gg1+3*dist3dgg+(ip3-2)*nLx);
		tg5=*(gg1+4*dist3dgg+(ip3-2)*nLx);
		tg6=*(gg1+5*dist3dgg+(ip3-2)*nLx);
	  
	    h12=(evl11*tg1+evl21*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h22=evl12*tg1+cm2*tg3;
		h32=evl13*tg1+cm2*tg4;
		h42=evl14*tg1+tea*tg2+teb*tg3+tec*tg4-gm1*tg5;
	    h52=(evl15*tg1+evl25*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h62=tg6;
		
		tg1=*(gg1+(ip2-2)*nLx);
		tg2=*(gg1+dist3dgg+(ip2-2)*nLx);
		tg3=*(gg1+2*dist3dgg+(ip2-2)*nLx);
		tg4=*(gg1+3*dist3dgg+(ip2-2)*nLx);
		tg5=*(gg1+4*dist3dgg+(ip2-2)*nLx);
		tg6=*(gg1+5*dist3dgg+(ip2-2)*nLx);

	    h13=(evl11*tg1+evl21*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h23=evl12*tg1+cm2*tg3;
		h33=evl13*tg1+cm2*tg4;
		h43=evl14*tg1+tea*tg2+teb*tg3+tec*tg4-gm1*tg5;
	    h53=(evl15*tg1+evl25*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h63=tg6;

		tg1=*(gg1+(ip1-2)*nLx);
		tg2=*(gg1+dist3dgg+(ip1-2)*nLx);
		tg3=*(gg1+2*dist3dgg+(ip1-2)*nLx);
		tg4=*(gg1+3*dist3dgg+(ip1-2)*nLx);
		tg5=*(gg1+4*dist3dgg+(ip1-2)*nLx);
		tg6=*(gg1+5*dist3dgg+(ip1-2)*nLx);

	    h14=(evl11*tg1+evl21*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h24=evl12*tg1+cm2*tg3;
		h34=evl13*tg1+cm2*tg4;
		h44=evl14*tg1+tea*tg2+teb*tg3+tec*tg4-gm1*tg5;
	    h54=(evl15*tg1+evl25*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h64=tg6;

		tg1=*(gg2+(i0-2)*nLx);
		tg2=*(gg2+dist3dgg+(i0-2)*nLx);
		tg3=*(gg2+2*dist3dgg+(i0-2)*nLx);
		tg4=*(gg2+3*dist3dgg+(i0-2)*nLx);
		tg5=*(gg2+4*dist3dgg+(i0-2)*nLx);
		tg6=*(gg2+5*dist3dgg+(i0-2)*nLx);

	    h15=(evl11*tg1+evl21*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h25=evl12*tg1+cm2*tg3;
		h35=evl13*tg1+cm2*tg4;
		h45=evl14*tg1+tea*tg2+teb*tg3+tec*tg4-gm1*tg5;
	    h55=(evl15*tg1+evl25*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h65=tg6;

		tg1=*(gg2+(ip1-2)*nLx);
		tg2=*(gg2+dist3dgg+(ip1-2)*nLx);
		tg3=*(gg2+2*dist3dgg+(ip1-2)*nLx);
		tg4=*(gg2+3*dist3dgg+(ip1-2)*nLx);
		tg5=*(gg2+4*dist3dgg+(ip1-2)*nLx);
		tg6=*(gg2+5*dist3dgg+(ip1-2)*nLx);

	    h16=(evl11*tg1+evl21*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h26=evl12*tg1+cm2*tg3;
		h36=evl13*tg1+cm2*tg4;
		h46=evl14*tg1+tea*tg2+teb*tg3+tec*tg4-gm1*tg5;
	    h56=(evl15*tg1+evl25*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h66=tg6;

		tg1=*(gg2+(ip2-2)*nLx);
		tg2=*(gg2+dist3dgg+(ip2-2)*nLx);
		tg3=*(gg2+2*dist3dgg+(ip2-2)*nLx);
		tg4=*(gg2+3*dist3dgg+(ip2-2)*nLx);
		tg5=*(gg2+4*dist3dgg+(ip2-2)*nLx);
		tg6=*(gg2+5*dist3dgg+(ip2-2)*nLx);

	    h17=(evl11*tg1+evl21*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h27=evl12*tg1+cm2*tg3;
		h37=evl13*tg1+cm2*tg4;
		h47=evl14*tg1+tea*tg2+teb*tg3+tec*tg4-gm1*tg5;
	    h57=(evl15*tg1+evl25*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h67=tg6;

		tg1=*(gg2+(ip3-2)*nLx);
		tg2=*(gg2+dist3dgg+(ip3-2)*nLx);
		tg3=*(gg2+2*dist3dgg+(ip3-2)*nLx);
		tg4=*(gg2+3*dist3dgg+(ip3-2)*nLx);
		tg5=*(gg2+4*dist3dgg+(ip3-2)*nLx);
		tg6=*(gg2+5*dist3dgg+(ip3-2)*nLx);

	    h18=(evl11*tg1+evl21*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h28=evl12*tg1+cm2*tg3;
		h38=evl13*tg1+cm2*tg4;
		h48=evl14*tg1+tea*tg2+teb*tg3+tec*tg4-gm1*tg5;
	    h58=(evl15*tg1+evl25*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h68=tg6;

		ts11=h11-h12;
		ts12=h12-h13;
		ts13=h13-h14;
		ts14=h15-h16;
		ts15=h16-h17;
		ts16=h17-h18;

		tt11=13.0*ts11*ts11+3.0*(h11-3.0*h12)*(h11-3.0*h12)+epweno;
		tt12=13.0*ts12*ts12+3.0*(h12+h13)*(h12+h13)+epweno;
		tt13=13.0*ts13*ts13+3.0*(3.0*h13-h14)*(3.0*h13-h14)+epweno;
		tt14=13.0*ts14*ts14+3.0*(h15-3.0*h16)*(h15-3.0*h16)+epweno;
		tt15=13.0*ts15*ts15+3.0*(h16+h17)*(h16+h17)+epweno;
		tt16=13.0*ts16*ts16+3.0*(3.0*h17-h18)*(3.0*h17-h18)+epweno;

		h11=(ts12-ts11)/(1.0+6.0*(tt11/tt12)*(tt11/tt12)+3.0*(tt11/tt13)*(tt11/tt13));
		h12=(ts13-ts12)/(2.0+(tt13/tt11)*(tt13/tt11)/1.5+4.0*(tt13/tt12)*(tt13/tt12))-0.25*(ts13-ts12);
		h13=(ts15-ts14)/(1.0+6.0*(tt14/tt15)*(tt14/tt15)+3.0*(tt14/tt16)*(tt14/tt16));
		h14=(ts16-ts15)/(2.0+(tt16/tt14)*(tt16/tt14)/1.5+4.0*(tt16/tt15)*(tt16/tt15))-0.25*(ts16-ts15);
		h15=(h11+h12+h13+h14)*cmm;

		ts11=h21-h22;
		ts12=h22-h23;
		ts13=h23-h24;
		ts14=h25-h26;
		ts15=h26-h27;
		ts16=h27-h28;

		tt11=13.0*ts11*ts11+3.0*(h21-3.0*h22)*(h21-3.0*h22)+epweno;
		tt12=13.0*ts12*ts12+3.0*(h22+h23)*(h22+h23)+epweno;
		tt13=13.0*ts13*ts13+3.0*(3.0*h23-h24)*(3.0*h23-h24)+epweno;
		tt14=13.0*ts14*ts14+3.0*(h25-3.0*h26)*(h25-3.0*h26)+epweno;
		tt15=13.0*ts15*ts15+3.0*(h26+h27)*(h26+h27)+epweno;
		tt16=13.0*ts16*ts16+3.0*(3.0*h27-h28)*(3.0*h27-h28)+epweno;

		h21=(ts12-ts11)/(1.0+6.0*(tt11/tt12)*(tt11/tt12)+3.0*(tt11/tt13)*(tt11/tt13));
		h22=(ts13-ts12)/(2.0+(tt13/tt11)*(tt13/tt11)/1.5+4.0*(tt13/tt12)*(tt13/tt12))-0.25*(ts13-ts12);
		h23=(ts15-ts14)/(1.0+6.0*(tt14/tt15)*(tt14/tt15)+3.0*(tt14/tt16)*(tt14/tt16));
		h24=(ts16-ts15)/(2.0+(tt16/tt14)*(tt16/tt14)/1.5+4.0*(tt16/tt15)*(tt16/tt15))-0.25*(ts16-ts15);
		h25=(h21+h22+h23+h24)*cmm;

		ts11=h31-h32;
		ts12=h32-h33;
		ts13=h33-h34;
		ts14=h35-h36;
		ts15=h36-h37;
		ts16=h37-h38;

		tt11=13.0*ts11*ts11+3.0*(h31-3.0*h32)*(h31-3.0*h32)+epweno;
		tt12=13.0*ts12*ts12+3.0*(h32+h33)*(h32+h33)+epweno;
		tt13=13.0*ts13*ts13+3.0*(3.0*h33-h34)*(3.0*h33-h34)+epweno;
		tt14=13.0*ts14*ts14+3.0*(h35-3.0*h36)*(h35-3.0*h36)+epweno;
		tt15=13.0*ts15*ts15+3.0*(h36+h37)*(h36+h37)+epweno;
		tt16=13.0*ts16*ts16+3.0*(3.0*h37-h38)*(3.0*h37-h38)+epweno;

		h31=(ts12-ts11)/(1.0+6.0*(tt11/tt12)*(tt11/tt12)+3.0*(tt11/tt13)*(tt11/tt13));
		h32=(ts13-ts12)/(2.0+(tt13/tt11)*(tt13/tt11)/1.5+4.0*(tt13/tt12)*(tt13/tt12))-0.25*(ts13-ts12);
		h33=(ts15-ts14)/(1.0+6.0*(tt14/tt15)*(tt14/tt15)+3.0*(tt14/tt16)*(tt14/tt16));
		h34=(ts16-ts15)/(2.0+(tt16/tt14)*(tt16/tt14)/1.5+4.0*(tt16/tt15)*(tt16/tt15))-0.25*(ts16-ts15);
		h35=(h31+h32+h33+h34)*cmm;

		ts11=h41-h42;
		ts12=h42-h43;
		ts13=h43-h44;
		ts14=h45-h46;
		ts15=h46-h47;
		ts16=h47-h48;

		tt11=13.0*ts11*ts11+3.0*(h41-3.0*h42)*(h41-3.0*h42)+epweno;
		tt12=13.0*ts12*ts12+3.0*(h42+h43)*(h42+h43)+epweno;
		tt13=13.0*ts13*ts13+3.0*(3.0*h43-h44)*(3.0*h43-h44)+epweno;
		tt14=13.0*ts14*ts14+3.0*(h45-3.0*h46)*(h45-3.0*h46)+epweno;
		tt15=13.0*ts15*ts15+3.0*(h46+h47)*(h46+h47)+epweno;
		tt16=13.0*ts16*ts16+3.0*(3.0*h47-h48)*(3.0*h47-h48)+epweno;

		h41=(ts12-ts11)/(1.0+6.0*(tt11/tt12)*(tt11/tt12)+3.0*(tt11/tt13)*(tt11/tt13));
		h42=(ts13-ts12)/(2.0+(tt13/tt11)*(tt13/tt11)/1.5+4.0*(tt13/tt12)*(tt13/tt12))-0.25*(ts13-ts12);
		h43=(ts15-ts14)/(1.0+6.0*(tt14/tt15)*(tt14/tt15)+3.0*(tt14/tt16)*(tt14/tt16));
		h44=(ts16-ts15)/(2.0+(tt16/tt14)*(tt16/tt14)/1.5+4.0*(tt16/tt15)*(tt16/tt15))-0.25*(ts16-ts15);
		h45=(h41+h42+h43+h44)*cmm;

		ts11=h51-h52;
		ts12=h52-h53;
		ts13=h53-h54;
		ts14=h55-h56;
		ts15=h56-h57;
		ts16=h57-h58;

		tt11=13.0*ts11*ts11+3.0*(h51-3.0*h52)*(h51-3.0*h52)+epweno;
		tt12=13.0*ts12*ts12+3.0*(h52+h53)*(h52+h53)+epweno;
		tt13=13.0*ts13*ts13+3.0*(3.0*h53-h54)*(3.0*h53-h54)+epweno;
		tt14=13.0*ts14*ts14+3.0*(h55-3.0*h56)*(h55-3.0*h56)+epweno;
		tt15=13.0*ts15*ts15+3.0*(h56+h57)*(h56+h57)+epweno;
		tt16=13.0*ts16*ts16+3.0*(3.0*h57-h58)*(3.0*h57-h58)+epweno;

		h51=(ts12-ts11)/(1.0+6.0*(tt11/tt12)*(tt11/tt12)+3.0*(tt11/tt13)*(tt11/tt13));
		h52=(ts13-ts12)/(2.0+(tt13/tt11)*(tt13/tt11)/1.5+4.0*(tt13/tt12)*(tt13/tt12))-0.25*(ts13-ts12);
		h53=(ts15-ts14)/(1.0+6.0*(tt14/tt15)*(tt14/tt15)+3.0*(tt14/tt16)*(tt14/tt16));
		h54=(ts16-ts15)/(2.0+(tt16/tt14)*(tt16/tt14)/1.5+4.0*(tt16/tt15)*(tt16/tt15))-0.25*(ts16-ts15);
		h55=(h51+h52+h53+h54)*cmm;

		ts11=h61-h62;
		ts12=h62-h63;
		ts13=h63-h64;
		ts14=h65-h66;
		ts15=h66-h67;
		ts16=h67-h68;

		tt11=13.0*ts11*ts11+3.0*(h61-3.0*h62)*(h61-3.0*h62)+epweno;
		tt12=13.0*ts12*ts12+3.0*(h62+h63)*(h62+h63)+epweno;
		tt13=13.0*ts13*ts13+3.0*(3.0*h63-h64)*(3.0*h63-h64)+epweno;
		tt14=13.0*ts14*ts14+3.0*(h65-3.0*h66)*(h65-3.0*h66)+epweno;
		tt15=13.0*ts15*ts15+3.0*(h66+h67)*(h66+h67)+epweno;
		tt16=13.0*ts16*ts16+3.0*(3.0*h67-h68)*(3.0*h67-h68)+epweno;

		h61=(ts12-ts11)/(3.0+18.0*(tt11/tt12)*(tt11/tt12)+9.0*(tt11/tt13)*(tt11/tt13));
		h62=(ts13-ts12)/(6.0+(tt13/tt11)*(tt13/tt11)*2.0+12.0*(tt13/tt12)*(tt13/tt12))-0.25*(ts13-ts12)/3.0;
		h63=(ts15-ts14)/(3.0+18.0*(tt14/tt15)*(tt14/tt15)+9.0*(tt14/tt16)*(tt14/tt16));
		h64=(ts16-ts15)/(6.0+(tt16/tt14)*(tt16/tt14)*2.0+12.0*(tt16/tt15)*(tt16/tt15))-0.25*(ts16-ts15)/3.0;
		h65=h61+h62+h63+h64;

		tfh1=h15+h45+h55;

		*(fh+(ip2-4)*nLx)=(-*(fd+ip3*nLx)+7.0*(*(fd+ip2*nLx)+*(fd+ip1*nLx))-*(fd+i0*nLx))/12.0+tfh1;
		*(fh+dist3dfh+(ip2-4)*nLx)=(-*(fd+dist3dfd+ip3*nLx)+7.0*(*(fd+dist3dfd+ip2*nLx)+*(fd+dist3dfd+ip1*nLx))-*(fd+dist3dfd+i0*nLx))/12.0+tfh1*tea/gm1-(h15-h55)*cm;
		*(fh+2*dist3dfh+(ip2-4)*nLx)=(-*(fd+2*dist3dfd+ip3*nLx)+7.0*(*(fd+2*dist3dfd+ip2*nLx)+*(fd+2*dist3dfd+ip1*nLx))-*(fd+2*dist3dfd+i0*nLx))/12.0+tfh1*teb/gm1+h25;
		*(fh+3*dist3dfh+(ip2-4)*nLx)=(-*(fd+3*dist3dfd+ip3*nLx)+7.0*(*(fd+3*dist3dfd+ip2*nLx)+*(fd+3*dist3dfd+ip1*nLx))-*(fd+3*dist3dfd+i0*nLx))/12.0+tfh1*tec/gm1+h35;
		*(fh+4*dist3dfh+(ip2-4)*nLx)=(-*(fd+4*dist3dfd+ip3*nLx)+7.0*(*(fd+4*dist3dfd+ip2*nLx)+*(fd+4*dist3dfd+ip1*nLx))-*(fd+4*dist3dfd+i0*nLx))/12.0+evr15*h15+teb/gm1*h25+tec/gm1*h35+qm*h45+evr55*h55;
		*(fh+5*dist3dfh+(ip2-4)*nLx)=(-*(fd+5*dist3dfd+ip3*nLx)+7.0*(*(fd+5*dist3dfd+ip2*nLx)+*(fd+5*dist3dfd+ip1*nLx))-*(fd+5*dist3dfd+i0*nLx))/12.0+h65;

    }//9~nLx+9
    		for(i=0;i<nLy;i++){
					*(rhs+i*nLx)+=*(fh+i*nLx)-*(fh+(i+1)*nLx);
					*(rhs+dist3drhs+i*nLx)+=*(fh+2*dist3dfh+i*nLx)-*(fh+2*dist3dfh+(i+1)*nLx);
					*(rhs+2*dist3drhs+i*nLx)+=*(fh+dist3dfh+i*nLx)-*(fh+dist3dfh+(i+1)*nLx);
					*(rhs+3*dist3drhs+i*nLx)+=*(fh+3*dist3dfh+i*nLx)-*(fh+3*dist3dfh+(i+1)*nLx);
					*(rhs+4*dist3drhs+i*nLx)+=*(fh+4*dist3dfh+i*nLx)-*(fh+4*dist3dfh+(i+1)*nLx);
					*(rhs+5*dist3drhs+i*nLx)+=*(fh+5*dist3dfh+i*nLx)-*(fh+5*dist3dfh+(i+1)*nLx);
			}
}


__global__ void hz_gpu(int nLx,int nLy,int nLz,int mnp,double gamma,double* fd,double* cs,double* vx,double* vy,double* vz,double* w,double* h,double* gg1,double* gg2,double* fh,cudaTextureObject_t texro,cudaTextureObject_t texvx,cudaTextureObject_t texvy,cudaTextureObject_t texvz,cudaTextureObject_t texeng,cudaTextureObject_t textracer)
{
	int i;
	double gm1=gamma-1;
	double den,xmt,ymt,zmt,eng,tracer;
	double am1,am2,am5;
	double aden,pre,vtemp,vxm,vym,vzm;
	int ip4,ip3,ip2,ip1,i0,in1,in2;
	double ta1,ta2,ta3,ta4,ta5,ta6;
	double tg1,tg2,tg3,tg4,tg5,tg6;
	double ttg1,ttg2,ttg3,ttg4,ttg5,ttg6;
	double te0,te1,te2,te3,te4,te5,te6,te7,te8,te9,tea,teb,tec;
	double evr15,evr55,evl11,evl21,evl12,evl13,evl14,evl15,evl25;
	double hm,qm,cm,cm2,cmm;
	double h11,h21,h31,h41,h51,h61;
	double h12,h22,h32,h42,h52,h62;
	double h13,h23,h33,h43,h53,h63;
	double h14,h24,h34,h44,h54,h64;
	double h15,h25,h35,h45,h55,h65;
	double h16,h26,h36,h46,h56,h66;
	double h17,h27,h37,h47,h57,h67;
	double h18,h28,h38,h48,h58,h68;

	double tt11,tt12,tt13,tt14,tt15,tt16;
	double ts11,ts12,ts13,ts14,ts15,ts16;
	double tfh1;

	int LL=nLz+10,LL1=nLz+5,LL2=nLz+1;
//	int LLx=nLx+10;
	int LLy=nLy+10;
	int offset=nLy*nLx;
	int dist3dfd=LL*nLx*nLy,dist3dgg=LL1*nLx*nLy,dist3dfh=LL2*nLx*nLy;
	unsigned int tx=threadIdx.x;
	unsigned int ty=threadIdx.y;
	unsigned int bx=blockIdx.x;
	unsigned int by=blockIdx.y;
	unsigned int xid=bx*BLOCK_SIZE_x+tx;
	unsigned int yid=by*BLOCK_SIZE_yh+ty;
//	int id=yid*nLx+xid;
//	uc+=(yid+5)*LLx+xid+5;
	fd+=yid*nLx+xid;
    cs+=yid*nLx+xid;
    vx+=yid*nLx+xid;
    vy+=yid*nLx+xid;
    vz+=yid*nLx+xid;
    w+=yid*nLx+xid;
    h+=yid*nLx+xid;
	gg1+=yid*nLx+xid;
	gg2+=yid*nLx+xid;
    fh+=yid*nLx+xid;
//	rhs+=yid*nLx+xid;
	float texi=xid+5+0.5f;
	float texj=yid+5+0.5f;	
    for(i=0;i<5;i++)
	{
	/*	den=*(uc+i*LLx*LLy);
		xmt=*(uc+3*dist3d+i*LLx*LLy);
		ymt=*(uc+2*dist3d+i*LLx*LLy);
		zmt=*(uc+dist3d+i*LLx*LLy);
		eng=*(uc+4*dist3d+i*LLx*LLy);
		tracer=*(uc+5*dist3d+i*LLx*LLy);*/
		den=fetch(texro,texi,texj+i*LLy);
		xmt=fetch(texvz,texi,texj+i*LLy);
		ymt=fetch(texvy,texi,texj+i*LLy);
		zmt=fetch(texvx,texi,texj+i*LLy);
		eng=fetch(texeng,texi,texj+i*LLy);
		tracer=fetch(textracer,texi,texj+i*LLy);

		aden=1.0/den;   ///////////////////////////cuda dividing is slow?????????????????????
		vxm=xmt*aden;
		vym=ymt*aden;
		vzm=zmt*aden;
        vtemp=0.5*(vxm*vxm+vym*vym+vzm*vzm);
		pre=gm1*(eng-vtemp*den);
         
		*(fd+i*offset)=xmt;
		*(fd+dist3dfd+i*offset)=xmt*vxm+pre;
		*(fd+2*dist3dfd+i*offset)=ymt*vxm;
		*(fd+3*dist3dfd+i*offset)=zmt*vxm;
		*(fd+4*dist3dfd+i*offset)=vxm*(pre+eng);
		*(fd+5*dist3dfd+i*offset)=tracer*vxm;

        *(cs+i*offset)=sqrt(fabs(gamma*pre*aden)); 
        *(vx+i*offset)=vxm;
        *(vy+i*offset)=vym;
        *(vz+i*offset)=vzm;
        *(w+i*offset)=sqrt(fabs(den));
        *(h+i*offset)=fabs(gamma*pre*aden)/gm1+vtemp;

    }//0~4

    for(i=5;i<9;i++)
	{
        ip2=i-5;
		ip1=i-4;
		i0=i-3;
		in1=i-2;
		in2=i-1;

	/*	den=*(uc+i*LLx*LLy);
		xmt=*(uc+3*dist3d+i*LLx*LLy);
		ymt=*(uc+2*dist3d+i*LLx*LLy);
		zmt=*(uc+dist3d+i*LLx*LLy);
		eng=*(uc+4*dist3d+i*LLx*LLy);
		tracer=*(uc+5*dist3d+i*LLx*LLy);*/
		den=fetch(texro,texi,texj+i*LLy);
		xmt=fetch(texvz,texi,texj+i*LLy);
		ymt=fetch(texvy,texi,texj+i*LLy);
		zmt=fetch(texvx,texi,texj+i*LLy);
		eng=fetch(texeng,texi,texj+i*LLy);
		tracer=fetch(textracer,texi,texj+i*LLy);

		aden=1.0/den;   ///////////////////////////cuda dividing is slow?????????????????????
		vxm=xmt*aden;
		vym=ymt*aden;
		vzm=zmt*aden;
        vtemp=0.5*(vxm*vxm+vym*vym+vzm*vzm);
		pre=gm1*(eng-vtemp*den);
         
		*(fd+i*offset)=xmt;
		*(fd+dist3dfd+i*offset)=xmt*vxm+pre;
		*(fd+2*dist3dfd+i*offset)=ymt*vxm;
		*(fd+3*dist3dfd+i*offset)=zmt*vxm;
		*(fd+4*dist3dfd+i*offset)=vxm*(pre+eng);
		*(fd+5*dist3dfd+i*offset)=tracer*vxm;

        *(cs+i*offset)=sqrt(fabs(gamma*pre*aden)); 
        *(vx+i*offset)=vxm;
        *(vy+i*offset)=vym;
        *(vz+i*offset)=vzm;
        *(w+i*offset)=sqrt(fabs(den));
        *(h+i*offset)=fabs(gamma*pre*aden)/gm1+vtemp;

        ta1=fmax(fabs(*(vx+ip2*offset)-*(cs+ip2*offset)),fabs(*(vx+ip1*offset)-*(cs+ip1*offset)));
        ta2=fmax(fabs(*(vx+ip2*offset)),fabs(*(vx+ip1*offset)));
        ta3=fmax(fabs(*(vx+ip2*offset)+*(cs+ip2*offset)),fabs(*(vx+ip1*offset)+*(cs+ip1*offset)));
        ta4=fmax(fabs(*(vx+i0*offset)-*(cs+i0*offset)),fabs(*(vx+in1*offset)-*(cs+in1*offset)));
        ta5=fmax(fabs(*(vx+i0*offset)),fabs(*(vx+in1*offset)));
        ta6=fmax(fabs(*(vx+i0*offset)+*(cs+i0*offset)),fabs(*(vx+in1*offset)+*(cs+in1*offset)));
        ta1=fmax(fabs(*(vx+in2*offset)-*(cs+in2*offset)),ta1);
        ta2=fmax(fabs(*(vx+in2*offset)),ta2);
        ta3=fmax(fabs(*(vx+in2*offset)+*(cs+in2*offset)),ta3);
        ta4=fmax(fabs(*(vx+i*offset)-*(cs+i*offset)),ta4);
        ta5=fmax(fabs(*(vx+i*offset)),ta5);
        ta6=fmax(fabs(*(vx+i*offset)+*(cs+i*offset)),ta6);
		am1=fmax(ta1,ta4)*ama;
		am2=fmax(ta2,ta5)*ama;
		am5=fmax(ta3,ta6)*ama;///////////////////////////////////////////////////////////!!!!!!!!!!!!!!!!!!!!!!
        
		tg1=*(fd+in1*offset)-*(fd+i0*offset);
		tg2=*(fd+dist3dfd+in1*offset)-*(fd+dist3dfd+i0*offset);
		tg3=*(fd+2*dist3dfd+in1*offset)-*(fd+2*dist3dfd+i0*offset);
		tg4=*(fd+3*dist3dfd+in1*offset)-*(fd+3*dist3dfd+i0*offset);
		tg5=*(fd+4*dist3dfd+in1*offset)-*(fd+4*dist3dfd+i0*offset);
		tg6=*(fd+5*dist3dfd+in1*offset)-*(fd+5*dist3dfd+i0*offset);

/*		ttg1=0.5*(tg1+am1*(*(uc+in1*LLx*LLy)-*(uc+i0*LLx*LLy)));
        ttg2=0.5*(tg2+am2*(*(uc+3*dist3d+in1*LLx*LLy)-*(uc+3*dist3d+i0*LLx*LLy)));
        ttg3=0.5*(tg3+am2*(*(uc+2*dist3d+in1*LLx*LLy)-*(uc+2*dist3d+i0*LLx*LLy)));
        ttg4=0.5*(tg4+am2*(*(uc+dist3d+in1*LLx*LLy)-*(uc+dist3d+i0*LLx*LLy)));
        ttg5=0.5*(tg5+am5*(*(uc+4*dist3d+in1*LLx*LLy)-*(uc+4*dist3d+i0*LLx*LLy)));
        ttg6=0.5*(tg6+am2*(*(uc+5*dist3d+in1*LLx*LLy)-*(uc+5*dist3d+i0*LLx*LLy)));*/
		ttg1=0.5*(tg1+am1*(fetch(texro,texi,texj+in1*LLy)-fetch(texro,texi,texj+i0*LLy)));
		ttg2=0.5*(tg2+am2*(fetch(texvz,texi,texj+in1*LLy)-fetch(texvz,texi,texj+i0*LLy)));
		ttg3=0.5*(tg3+am2*(fetch(texvy,texi,texj+in1*LLy)-fetch(texvy,texi,texj+i0*LLy)));
		ttg4=0.5*(tg4+am2*(fetch(texvx,texi,texj+in1*LLy)-fetch(texvx,texi,texj+i0*LLy)));
		ttg5=0.5*(tg5+am5*(fetch(texeng,texi,texj+in1*LLy)-fetch(texeng,texi,texj+i0*LLy)));
		ttg6=0.5*(tg6+am2*(fetch(textracer,texi,texj+in1*LLy)-fetch(textracer,texi,texj+i0*LLy)));
		*(gg1+(i0-2)*offset)=ttg1;
		*(gg1+dist3dgg+(i0-2)*offset)=ttg2;
		*(gg1+2*dist3dgg+(i0-2)*offset)=ttg3;
		*(gg1+3*dist3dgg+(i0-2)*offset)=ttg4;
		*(gg1+4*dist3dgg+(i0-2)*offset)=ttg5;
		*(gg1+5*dist3dgg+(i0-2)*offset)=ttg6;

		*(gg2+(i0-2)*offset)=ttg1-tg1;
		*(gg2+dist3dgg+(i0-2)*offset)=ttg2-tg2;
		*(gg2+2*dist3dgg+(i0-2)*offset)=ttg3-tg3;
		*(gg2+3*dist3dgg+(i0-2)*offset)=ttg4-tg4;
		*(gg2+4*dist3dgg+(i0-2)*offset)=ttg5-tg5;
		*(gg2+5*dist3dgg+(i0-2)*offset)=ttg6-tg6;

    }//5~8

    for(i=9;i<LL;i++)
	{
		ip4=i-7;
		ip3=i-6;
        ip2=i-5;
		ip1=i-4;
		i0=i-3;
		in1=i-2;
		in2=i-1;

/*		den=*(uc+i*LLx*LLy);
		xmt=*(uc+3*dist3d+i*LLx*LLy);
		ymt=*(uc+2*dist3d+i*LLx*LLy);
		zmt=*(uc+dist3d+i*LLx*LLy);
		eng=*(uc+4*dist3d+i*LLx*LLy);
		tracer=*(uc+5*dist3d+i*LLx*LLy);*/
		den=fetch(texro,texi,texj+i*LLy);
		xmt=fetch(texvz,texi,texj+i*LLy);
		ymt=fetch(texvy,texi,texj+i*LLy);
		zmt=fetch(texvx,texi,texj+i*LLy);
		eng=fetch(texeng,texi,texj+i*LLy);
		tracer=fetch(textracer,texi,texj+i*LLy);

		aden=1.0/den;   ///////////////////////////cuda dividing is slow?????????????????????
		vxm=xmt*aden;
		vym=ymt*aden;
		vzm=zmt*aden;
        vtemp=0.5*(vxm*vxm+vym*vym+vzm*vzm);
		pre=gm1*(eng-vtemp*den);
         
		*(fd+i*offset)=xmt;
		*(fd+dist3dfd+i*offset)=xmt*vxm+pre;
		*(fd+2*dist3dfd+i*offset)=ymt*vxm;
		*(fd+3*dist3dfd+i*offset)=zmt*vxm;
		*(fd+4*dist3dfd+i*offset)=vxm*(pre+eng);
		*(fd+5*dist3dfd+i*offset)=tracer*vxm;

        *(cs+i*offset)=sqrt(fabs(gamma*pre*aden)); 
        *(vx+i*offset)=vxm;
        *(vy+i*offset)=vym;
        *(vz+i*offset)=vzm;
        *(w+i*offset)=sqrt(fabs(den));
        *(h+i*offset)=fabs(gamma*pre*aden)/gm1+vtemp;

        ta1=fmax(fabs(*(vx+ip2*offset)-*(cs+ip2*offset)),fabs(*(vx+ip1*offset)-*(cs+ip1*offset)));
        ta2=fmax(fabs(*(vx+ip2*offset)),fabs(*(vx+ip1*offset)));
        ta3=fmax(fabs(*(vx+ip2*offset)+*(cs+ip2*offset)),fabs(*(vx+ip1*offset)+*(cs+ip1*offset)));
        ta4=fmax(fabs(*(vx+i0*offset)-*(cs+i0*offset)),fabs(*(vx+in1*offset)-*(cs+in1*offset)));
        ta5=fmax(fabs(*(vx+i0*offset)),fabs(*(vx+in1*offset)));
        ta6=fmax(fabs(*(vx+i0*offset)+*(cs+i0*offset)),fabs(*(vx+in1*offset)+*(cs+in1*offset)));
        ta1=fmax(fabs(*(vx+in2*offset)-*(cs+in2*offset)),ta1);
        ta2=fmax(fabs(*(vx+in2*offset)),ta2);
        ta3=fmax(fabs(*(vx+in2*offset)+*(cs+in2*offset)),ta3);
        ta4=fmax(fabs(*(vx+i*offset)-*(cs+i*offset)),ta4);
        ta5=fmax(fabs(*(vx+i*offset)),ta5);
        ta6=fmax(fabs(*(vx+i*offset)+*(cs+i*offset)),ta6);
		am1=fmax(ta1,ta4)*ama;
		am2=fmax(ta2,ta5)*ama;
		am5=fmax(ta3,ta6)*ama;///////////////////////////////////////////////////////////!!!!!!!!!!!!!!!!!!!!!!
        
		tg1=*(fd+in1*offset)-*(fd+i0*offset);
		tg2=*(fd+dist3dfd+in1*offset)-*(fd+dist3dfd+i0*offset);
		tg3=*(fd+2*dist3dfd+in1*offset)-*(fd+2*dist3dfd+i0*offset);
		tg4=*(fd+3*dist3dfd+in1*offset)-*(fd+3*dist3dfd+i0*offset);
		tg5=*(fd+4*dist3dfd+in1*offset)-*(fd+4*dist3dfd+i0*offset);
		tg6=*(fd+5*dist3dfd+in1*offset)-*(fd+5*dist3dfd+i0*offset);

/*		ttg1=0.5*(tg1+am1*(*(uc+in1*LLx*LLy)-*(uc+i0*LLx*LLy)));
        ttg2=0.5*(tg2+am2*(*(uc+3*dist3d+in1*LLx*LLy)-*(uc+3*dist3d+i0*LLx*LLy)));
        ttg3=0.5*(tg3+am2*(*(uc+2*dist3d+in1*LLx*LLy)-*(uc+2*dist3d+i0*LLx*LLy)));
        ttg4=0.5*(tg4+am2*(*(uc+dist3d+in1*LLx*LLy)-*(uc+dist3d+i0*LLx*LLy)));
        ttg5=0.5*(tg5+am5*(*(uc+4*dist3d+in1*LLx*LLy)-*(uc+4*dist3d+i0*LLx*LLy)));
        ttg6=0.5*(tg6+am2*(*(uc+5*dist3d+in1*LLx*LLy)-*(uc+5*dist3d+i0*LLx*LLy)));*/
		ttg1=0.5*(tg1+am1*(fetch(texro,texi,texj+in1*LLy)-fetch(texro,texi,texj+i0*LLy)));
		ttg2=0.5*(tg2+am2*(fetch(texvz,texi,texj+in1*LLy)-fetch(texvz,texi,texj+i0*LLy)));
		ttg3=0.5*(tg3+am2*(fetch(texvy,texi,texj+in1*LLy)-fetch(texvy,texi,texj+i0*LLy)));
		ttg4=0.5*(tg4+am2*(fetch(texvx,texi,texj+in1*LLy)-fetch(texvx,texi,texj+i0*LLy)));
		ttg5=0.5*(tg5+am5*(fetch(texeng,texi,texj+in1*LLy)-fetch(texeng,texi,texj+i0*LLy)));
		ttg6=0.5*(tg6+am2*(fetch(textracer,texi,texj+in1*LLy)-fetch(textracer,texi,texj+i0*LLy)));
		*(gg1+(i0-2)*offset)=ttg1;
		*(gg1+dist3dgg+(i0-2)*offset)=ttg2;
		*(gg1+2*dist3dgg+(i0-2)*offset)=ttg3;
		*(gg1+3*dist3dgg+(i0-2)*offset)=ttg4;
		*(gg1+4*dist3dgg+(i0-2)*offset)=ttg5;
		*(gg1+5*dist3dgg+(i0-2)*offset)=ttg6;

		*(gg2+(i0-2)*offset)=ttg1-tg1;
		*(gg2+dist3dgg+(i0-2)*offset)=ttg2-tg2;
		*(gg2+2*dist3dgg+(i0-2)*offset)=ttg3-tg3;
		*(gg2+3*dist3dgg+(i0-2)*offset)=ttg4-tg4;
		*(gg2+4*dist3dgg+(i0-2)*offset)=ttg5-tg5;
		*(gg2+5*dist3dgg+(i0-2)*offset)=ttg6-tg6;
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
        te0=*(w+ip2*offset)/(*(w+ip2*offset)+*(w+ip1*offset));

		te1=1-te0;
		te2=*(vx+ip2*offset);
		te3=*(vy+ip2*offset);
		te4=*(vz+ip2*offset);
		te5=*(vx+ip1*offset);
		te6=*(vy+ip1*offset);
		te7=*(vz+ip1*offset);
		tea=te0*te2+te1*te5;
		teb=te0*te3+te1*te6;
		tec=te0*te4+te1*te7;
        
		hm=te0*(*(h+ip2*offset))+te1*(*(h+ip1*offset));
		qm=0.5*(tea*tea+teb*teb+tec*tec);
		te8=gm1*qm;
        cm2=te0*pow(*(cs+ip2*offset),2)+te1*(pow(*(cs+ip1*offset),2))+0.5*te0*te1*gm1*((te2-te5)*(te2-te5)+(te3-te6)*(te3-te6)+(te4-te7)*(te4-te7));
		cm=sqrt(cm2);
		cmm=1.0/(3.0*cm2);
		te9=tea*cm;
		tea*=gm1;
		teb*=gm1;
		tec*=gm1;

		evr15=hm-te9;
		evr55=hm+te9;

		evl11=te8+te9;
		evl21=-(cm+tea);
		evl12=-teb/gm1*cm2;
		evl13=-tec/gm1*cm2;
		evl14=cm2-te8;
		evl15=te8-te9;
		evl25=cm-tea;

		tg1=*(gg1+(ip4-2)*offset);
		tg2=*(gg1+dist3dgg+(ip4-2)*offset);
		tg3=*(gg1+2*dist3dgg+(ip4-2)*offset);
		tg4=*(gg1+3*dist3dgg+(ip4-2)*offset);
		tg5=*(gg1+4*dist3dgg+(ip4-2)*offset);
		tg6=*(gg1+5*dist3dgg+(ip4-2)*offset);

	    h11=(evl11*tg1+evl21*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h21=evl12*tg1+cm2*tg3;
		h31=evl13*tg1+cm2*tg4;
		h41=evl14*tg1+tea*tg2+teb*tg3+tec*tg4-gm1*tg5;
	    h51=(evl15*tg1+evl25*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h61=tg6;

		tg1=*(gg1+(ip3-2)*offset);
		tg2=*(gg1+dist3dgg+(ip3-2)*offset);
		tg3=*(gg1+2*dist3dgg+(ip3-2)*offset);
		tg4=*(gg1+3*dist3dgg+(ip3-2)*offset);
		tg5=*(gg1+4*dist3dgg+(ip3-2)*offset);
		tg6=*(gg1+5*dist3dgg+(ip3-2)*offset);
	  
	    h12=(evl11*tg1+evl21*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h22=evl12*tg1+cm2*tg3;
		h32=evl13*tg1+cm2*tg4;
		h42=evl14*tg1+tea*tg2+teb*tg3+tec*tg4-gm1*tg5;
	    h52=(evl15*tg1+evl25*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h62=tg6;
		
		tg1=*(gg1+(ip2-2)*offset);
		tg2=*(gg1+dist3dgg+(ip2-2)*offset);
		tg3=*(gg1+2*dist3dgg+(ip2-2)*offset);
		tg4=*(gg1+3*dist3dgg+(ip2-2)*offset);
		tg5=*(gg1+4*dist3dgg+(ip2-2)*offset);
		tg6=*(gg1+5*dist3dgg+(ip2-2)*offset);

	    h13=(evl11*tg1+evl21*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h23=evl12*tg1+cm2*tg3;
		h33=evl13*tg1+cm2*tg4;
		h43=evl14*tg1+tea*tg2+teb*tg3+tec*tg4-gm1*tg5;
	    h53=(evl15*tg1+evl25*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h63=tg6;

		tg1=*(gg1+(ip1-2)*offset);
		tg2=*(gg1+dist3dgg+(ip1-2)*offset);
		tg3=*(gg1+2*dist3dgg+(ip1-2)*offset);
		tg4=*(gg1+3*dist3dgg+(ip1-2)*offset);
		tg5=*(gg1+4*dist3dgg+(ip1-2)*offset);
		tg6=*(gg1+5*dist3dgg+(ip1-2)*offset);

	    h14=(evl11*tg1+evl21*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h24=evl12*tg1+cm2*tg3;
		h34=evl13*tg1+cm2*tg4;
		h44=evl14*tg1+tea*tg2+teb*tg3+tec*tg4-gm1*tg5;
	    h54=(evl15*tg1+evl25*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h64=tg6;

		tg1=*(gg2+(i0-2)*offset);
		tg2=*(gg2+dist3dgg+(i0-2)*offset);
		tg3=*(gg2+2*dist3dgg+(i0-2)*offset);
		tg4=*(gg2+3*dist3dgg+(i0-2)*offset);
		tg5=*(gg2+4*dist3dgg+(i0-2)*offset);
		tg6=*(gg2+5*dist3dgg+(i0-2)*offset);

	    h15=(evl11*tg1+evl21*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h25=evl12*tg1+cm2*tg3;
		h35=evl13*tg1+cm2*tg4;
		h45=evl14*tg1+tea*tg2+teb*tg3+tec*tg4-gm1*tg5;
	    h55=(evl15*tg1+evl25*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h65=tg6;

		tg1=*(gg2+(ip1-2)*offset);
		tg2=*(gg2+dist3dgg+(ip1-2)*offset);
		tg3=*(gg2+2*dist3dgg+(ip1-2)*offset);
		tg4=*(gg2+3*dist3dgg+(ip1-2)*offset);
		tg5=*(gg2+4*dist3dgg+(ip1-2)*offset);
		tg6=*(gg2+5*dist3dgg+(ip1-2)*offset);

	    h16=(evl11*tg1+evl21*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h26=evl12*tg1+cm2*tg3;
		h36=evl13*tg1+cm2*tg4;
		h46=evl14*tg1+tea*tg2+teb*tg3+tec*tg4-gm1*tg5;
	    h56=(evl15*tg1+evl25*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h66=tg6;

		tg1=*(gg2+(ip2-2)*offset);
		tg2=*(gg2+dist3dgg+(ip2-2)*offset);
		tg3=*(gg2+2*dist3dgg+(ip2-2)*offset);
		tg4=*(gg2+3*dist3dgg+(ip2-2)*offset);
		tg5=*(gg2+4*dist3dgg+(ip2-2)*offset);
		tg6=*(gg2+5*dist3dgg+(ip2-2)*offset);

	    h17=(evl11*tg1+evl21*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h27=evl12*tg1+cm2*tg3;
		h37=evl13*tg1+cm2*tg4;
		h47=evl14*tg1+tea*tg2+teb*tg3+tec*tg4-gm1*tg5;
	    h57=(evl15*tg1+evl25*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h67=tg6;

		tg1=*(gg2+(ip3-2)*offset);
		tg2=*(gg2+dist3dgg+(ip3-2)*offset);
		tg3=*(gg2+2*dist3dgg+(ip3-2)*offset);
		tg4=*(gg2+3*dist3dgg+(ip3-2)*offset);
		tg5=*(gg2+4*dist3dgg+(ip3-2)*offset);
		tg6=*(gg2+5*dist3dgg+(ip3-2)*offset);

	    h18=(evl11*tg1+evl21*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h28=evl12*tg1+cm2*tg3;
		h38=evl13*tg1+cm2*tg4;
		h48=evl14*tg1+tea*tg2+teb*tg3+tec*tg4-gm1*tg5;
	    h58=(evl15*tg1+evl25*tg2-teb*tg3-tec*tg4+gm1*tg5)*0.5;
		h68=tg6;

		ts11=h11-h12;
		ts12=h12-h13;
		ts13=h13-h14;
		ts14=h15-h16;
		ts15=h16-h17;
		ts16=h17-h18;

		tt11=13.0*ts11*ts11+3.0*(h11-3.0*h12)*(h11-3.0*h12)+epweno;
		tt12=13.0*ts12*ts12+3.0*(h12+h13)*(h12+h13)+epweno;
		tt13=13.0*ts13*ts13+3.0*(3.0*h13-h14)*(3.0*h13-h14)+epweno;
		tt14=13.0*ts14*ts14+3.0*(h15-3.0*h16)*(h15-3.0*h16)+epweno;
		tt15=13.0*ts15*ts15+3.0*(h16+h17)*(h16+h17)+epweno;
		tt16=13.0*ts16*ts16+3.0*(3.0*h17-h18)*(3.0*h17-h18)+epweno;

		h11=(ts12-ts11)/(1.0+6.0*(tt11/tt12)*(tt11/tt12)+3.0*(tt11/tt13)*(tt11/tt13));
		h12=(ts13-ts12)/(2.0+(tt13/tt11)*(tt13/tt11)/1.5+4.0*(tt13/tt12)*(tt13/tt12))-0.25*(ts13-ts12);
		h13=(ts15-ts14)/(1.0+6.0*(tt14/tt15)*(tt14/tt15)+3.0*(tt14/tt16)*(tt14/tt16));
		h14=(ts16-ts15)/(2.0+(tt16/tt14)*(tt16/tt14)/1.5+4.0*(tt16/tt15)*(tt16/tt15))-0.25*(ts16-ts15);
		h15=(h11+h12+h13+h14)*cmm;

		ts11=h21-h22;
		ts12=h22-h23;
		ts13=h23-h24;
		ts14=h25-h26;
		ts15=h26-h27;
		ts16=h27-h28;

		tt11=13.0*ts11*ts11+3.0*(h21-3.0*h22)*(h21-3.0*h22)+epweno;
		tt12=13.0*ts12*ts12+3.0*(h22+h23)*(h22+h23)+epweno;
		tt13=13.0*ts13*ts13+3.0*(3.0*h23-h24)*(3.0*h23-h24)+epweno;
		tt14=13.0*ts14*ts14+3.0*(h25-3.0*h26)*(h25-3.0*h26)+epweno;
		tt15=13.0*ts15*ts15+3.0*(h26+h27)*(h26+h27)+epweno;
		tt16=13.0*ts16*ts16+3.0*(3.0*h27-h28)*(3.0*h27-h28)+epweno;

		h21=(ts12-ts11)/(1.0+6.0*(tt11/tt12)*(tt11/tt12)+3.0*(tt11/tt13)*(tt11/tt13));
		h22=(ts13-ts12)/(2.0+(tt13/tt11)*(tt13/tt11)/1.5+4.0*(tt13/tt12)*(tt13/tt12))-0.25*(ts13-ts12);
		h23=(ts15-ts14)/(1.0+6.0*(tt14/tt15)*(tt14/tt15)+3.0*(tt14/tt16)*(tt14/tt16));
		h24=(ts16-ts15)/(2.0+(tt16/tt14)*(tt16/tt14)/1.5+4.0*(tt16/tt15)*(tt16/tt15))-0.25*(ts16-ts15);
		h25=(h21+h22+h23+h24)*cmm;

		ts11=h31-h32;
		ts12=h32-h33;
		ts13=h33-h34;
		ts14=h35-h36;
		ts15=h36-h37;
		ts16=h37-h38;

		tt11=13.0*ts11*ts11+3.0*(h31-3.0*h32)*(h31-3.0*h32)+epweno;
		tt12=13.0*ts12*ts12+3.0*(h32+h33)*(h32+h33)+epweno;
		tt13=13.0*ts13*ts13+3.0*(3.0*h33-h34)*(3.0*h33-h34)+epweno;
		tt14=13.0*ts14*ts14+3.0*(h35-3.0*h36)*(h35-3.0*h36)+epweno;
		tt15=13.0*ts15*ts15+3.0*(h36+h37)*(h36+h37)+epweno;
		tt16=13.0*ts16*ts16+3.0*(3.0*h37-h38)*(3.0*h37-h38)+epweno;

		h31=(ts12-ts11)/(1.0+6.0*(tt11/tt12)*(tt11/tt12)+3.0*(tt11/tt13)*(tt11/tt13));
		h32=(ts13-ts12)/(2.0+(tt13/tt11)*(tt13/tt11)/1.5+4.0*(tt13/tt12)*(tt13/tt12))-0.25*(ts13-ts12);
		h33=(ts15-ts14)/(1.0+6.0*(tt14/tt15)*(tt14/tt15)+3.0*(tt14/tt16)*(tt14/tt16));
		h34=(ts16-ts15)/(2.0+(tt16/tt14)*(tt16/tt14)/1.5+4.0*(tt16/tt15)*(tt16/tt15))-0.25*(ts16-ts15);
		h35=(h31+h32+h33+h34)*cmm;

		ts11=h41-h42;
		ts12=h42-h43;
		ts13=h43-h44;
		ts14=h45-h46;
		ts15=h46-h47;
		ts16=h47-h48;

		tt11=13.0*ts11*ts11+3.0*(h41-3.0*h42)*(h41-3.0*h42)+epweno;
		tt12=13.0*ts12*ts12+3.0*(h42+h43)*(h42+h43)+epweno;
		tt13=13.0*ts13*ts13+3.0*(3.0*h43-h44)*(3.0*h43-h44)+epweno;
		tt14=13.0*ts14*ts14+3.0*(h45-3.0*h46)*(h45-3.0*h46)+epweno;
		tt15=13.0*ts15*ts15+3.0*(h46+h47)*(h46+h47)+epweno;
		tt16=13.0*ts16*ts16+3.0*(3.0*h47-h48)*(3.0*h47-h48)+epweno;

		h41=(ts12-ts11)/(1.0+6.0*(tt11/tt12)*(tt11/tt12)+3.0*(tt11/tt13)*(tt11/tt13));
		h42=(ts13-ts12)/(2.0+(tt13/tt11)*(tt13/tt11)/1.5+4.0*(tt13/tt12)*(tt13/tt12))-0.25*(ts13-ts12);
		h43=(ts15-ts14)/(1.0+6.0*(tt14/tt15)*(tt14/tt15)+3.0*(tt14/tt16)*(tt14/tt16));
		h44=(ts16-ts15)/(2.0+(tt16/tt14)*(tt16/tt14)/1.5+4.0*(tt16/tt15)*(tt16/tt15))-0.25*(ts16-ts15);
		h45=(h41+h42+h43+h44)*cmm;

		ts11=h51-h52;
		ts12=h52-h53;
		ts13=h53-h54;
		ts14=h55-h56;
		ts15=h56-h57;
		ts16=h57-h58;

		tt11=13.0*ts11*ts11+3.0*(h51-3.0*h52)*(h51-3.0*h52)+epweno;
		tt12=13.0*ts12*ts12+3.0*(h52+h53)*(h52+h53)+epweno;
		tt13=13.0*ts13*ts13+3.0*(3.0*h53-h54)*(3.0*h53-h54)+epweno;
		tt14=13.0*ts14*ts14+3.0*(h55-3.0*h56)*(h55-3.0*h56)+epweno;
		tt15=13.0*ts15*ts15+3.0*(h56+h57)*(h56+h57)+epweno;
		tt16=13.0*ts16*ts16+3.0*(3.0*h57-h58)*(3.0*h57-h58)+epweno;

		h51=(ts12-ts11)/(1.0+6.0*(tt11/tt12)*(tt11/tt12)+3.0*(tt11/tt13)*(tt11/tt13));
		h52=(ts13-ts12)/(2.0+(tt13/tt11)*(tt13/tt11)/1.5+4.0*(tt13/tt12)*(tt13/tt12))-0.25*(ts13-ts12);
		h53=(ts15-ts14)/(1.0+6.0*(tt14/tt15)*(tt14/tt15)+3.0*(tt14/tt16)*(tt14/tt16));
		h54=(ts16-ts15)/(2.0+(tt16/tt14)*(tt16/tt14)/1.5+4.0*(tt16/tt15)*(tt16/tt15))-0.25*(ts16-ts15);
		h55=(h51+h52+h53+h54)*cmm;

		ts11=h61-h62;
		ts12=h62-h63;
		ts13=h63-h64;
		ts14=h65-h66;
		ts15=h66-h67;
		ts16=h67-h68;

		tt11=13.0*ts11*ts11+3.0*(h61-3.0*h62)*(h61-3.0*h62)+epweno;
		tt12=13.0*ts12*ts12+3.0*(h62+h63)*(h62+h63)+epweno;
		tt13=13.0*ts13*ts13+3.0*(3.0*h63-h64)*(3.0*h63-h64)+epweno;
		tt14=13.0*ts14*ts14+3.0*(h65-3.0*h66)*(h65-3.0*h66)+epweno;
		tt15=13.0*ts15*ts15+3.0*(h66+h67)*(h66+h67)+epweno;
		tt16=13.0*ts16*ts16+3.0*(3.0*h67-h68)*(3.0*h67-h68)+epweno;

		h61=(ts12-ts11)/(3.0+18.0*(tt11/tt12)*(tt11/tt12)+9.0*(tt11/tt13)*(tt11/tt13));
		h62=(ts13-ts12)/(6.0+(tt13/tt11)*(tt13/tt11)*2.0+12.0*(tt13/tt12)*(tt13/tt12))-0.25*(ts13-ts12)/3.0;
		h63=(ts15-ts14)/(3.0+18.0*(tt14/tt15)*(tt14/tt15)+9.0*(tt14/tt16)*(tt14/tt16));
		h64=(ts16-ts15)/(6.0+(tt16/tt14)*(tt16/tt14)*2.0+12.0*(tt16/tt15)*(tt16/tt15))-0.25*(ts16-ts15)/3.0;
		h65=h61+h62+h63+h64;

		tfh1=h15+h45+h55;
		*(fh+(ip2-4)*offset)=(-*(fd+ip3*offset)+7.0*(*(fd+ip2*offset)+*(fd+ip1*offset))-*(fd+i0*offset))/12.0+tfh1;
		*(fh+dist3dfh+(ip2-4)*offset)=(-*(fd+dist3dfd+ip3*offset)+7.0*(*(fd+dist3dfd+ip2*offset)+*(fd+dist3dfd+ip1*offset))-*(fd+dist3dfd+i0*offset))/12.0+tfh1*tea/gm1-(h15-h55)*cm;
		*(fh+2*dist3dfh+(ip2-4)*offset)=(-*(fd+2*dist3dfd+ip3*offset)+7.0*(*(fd+2*dist3dfd+ip2*offset)+*(fd+2*dist3dfd+ip1*offset))-*(fd+2*dist3dfd+i0*offset))/12.0+tfh1*teb/gm1+h25;
		*(fh+3*dist3dfh+(ip2-4)*offset)=(-*(fd+3*dist3dfd+ip3*offset)+7.0*(*(fd+3*dist3dfd+ip2*offset)+*(fd+3*dist3dfd+ip1*offset))-*(fd+3*dist3dfd+i0*offset))/12.0+tfh1*tec/gm1+h35;
		*(fh+4*dist3dfh+(ip2-4)*offset)=(-*(fd+4*dist3dfd+ip3*offset)+7.0*(*(fd+4*dist3dfd+ip2*offset)+*(fd+4*dist3dfd+ip1*offset))-*(fd+4*dist3dfd+i0*offset))/12.0+evr15*h15+teb/gm1*h25+tec/gm1*h35+qm*h45+evr55*h55;
		*(fh+5*dist3dfh+(ip2-4)*offset)=(-*(fd+5*dist3dfd+ip3*offset)+7.0*(*(fd+5*dist3dfd+ip2*offset)+*(fd+5*dist3dfd+ip1*offset))-*(fd+5*dist3dfd+i0*offset))/12.0+h65;

    }//9~nLx+9
}
__global__ void hzleft_gpu(double* rhs,int nLx,int nLy,int nLz,double*fh)
{
	int LL2=nLz+1;
	int offset=nLy*nLx;
	int dist3dfh=LL2*nLx*nLy,dist3drhs=nLx*nLy*nLz;
	unsigned int tx=threadIdx.x;
	unsigned int ty=threadIdx.y;
	unsigned int bx=blockIdx.x;
	unsigned int by=blockIdx.y;
	unsigned int xid=bx*BLOCK_SIZE_x+tx;
	unsigned int yid=by*BLOCK_SIZE_yh+ty;
	rhs+=yid*nLx+xid;
	fh+=yid*nLx+xid;
    		for(unsigned int i=0;i<nLz;i++){
					*(rhs+i*offset)+=*(fh+i*offset)-*(fh+(i+1)*offset);
					*(rhs+dist3drhs+i*offset)+=*(fh+3*dist3dfh+i*offset)-*(fh+3*dist3dfh+(i+1)*offset);
					*(rhs+2*dist3drhs+i*offset)+=*(fh+2*dist3dfh+i*offset)-*(fh+2*dist3dfh+(i+1)*offset);
					*(rhs+3*dist3drhs+i*offset)+=*(fh+dist3dfh+i*offset)-*(fh+dist3dfh+(i+1)*offset);
					*(rhs+4*dist3drhs+i*offset)+=*(fh+4*dist3dfh+i*offset)-*(fh+4*dist3dfh+(i+1)*offset);
					*(rhs+5*dist3drhs+i*offset)+=*(fh+5*dist3dfh+i*offset)-*(fh+5*dist3dfh+(i+1)*offset);
			}
}
__global__ void fxleft_gpu(double* rhs,int nLx,int nLy,int nLz,double*fh)
{
	int dist3dfh=(nLx+1)*nLy*nLz,dist3drhs=nLx*nLy*nLz;
	unsigned int tx=threadIdx.x;
	unsigned int ty=threadIdx.y;
    fh+=(blockIdx.y*BLOCK_SIZE_yf+ty)*nLy*(nLx+1)+blockIdx.x*BLOCK_SIZE_xf+tx;
	rhs+=((blockIdx.y*BLOCK_SIZE_yf+ty)*nLy+blockIdx.x*BLOCK_SIZE_xf)*nLx+tx;
    __shared__ double shm[BLOCK_SIZE_xf*BLOCK_SIZE_yf][17];
	unsigned int lid=ty*BLOCK_SIZE_xf+tx;
			unsigned int i,k;
			double* lrhs=rhs;
    		for(i=0;i<nLx;i+=16){
					for(k=0;k<16;k++)
						shm[lid][k]=*(fh+(i+k)*nLy)-*(fh+(i+k+1)*nLy);	
					__syncthreads();
					for(k=0;k<16;k++)
						*(lrhs+k*nLx)+=shm[ty*16+k][tx];	
					__syncthreads();
					lrhs+=BLOCK_SIZE_xf;
			}
			lrhs=rhs+dist3drhs;
    		for(i=0;i<nLx;i+=16){
					for(k=0;k<16;k++){
						shm[lid][k]=*(fh+dist3dfh+(i+k)*nLy)-*(fh+dist3dfh+(i+k+1)*nLy);	
					}
					__syncthreads();
					for(k=0;k<16;k++){
						*(lrhs+k*nLx)+=shm[ty*16+k][tx];	
					}
					__syncthreads();
					lrhs+=BLOCK_SIZE_xf;
			}
			lrhs=rhs+2*dist3drhs;
    		for(i=0;i<nLx;i+=16){
					for(k=0;k<16;k++){
						shm[lid][k]=*(fh+2*dist3dfh+(i+k)*nLy)-*(fh+2*dist3dfh+(i+k+1)*nLy);	
					}
					__syncthreads();
					for(k=0;k<16;k++){
						*(lrhs+k*nLx)+=shm[ty*16+k][tx];	
					}
					__syncthreads();
					lrhs+=BLOCK_SIZE_xf;
			}
			lrhs=rhs+3*dist3drhs;
    		for(i=0;i<nLx;i+=16){
					for(k=0;k<16;k++){
						shm[lid][k]=*(fh+3*dist3dfh+(i+k)*nLy)-*(fh+3*dist3dfh+(i+k+1)*nLy);	
					}
					__syncthreads();
					for(k=0;k<16;k++){
						*(lrhs+k*nLx)+=shm[ty*16+k][tx];	
					}
					__syncthreads();
					lrhs+=BLOCK_SIZE_xf;
			}
			lrhs=rhs+4*dist3drhs;
    		for(i=0;i<nLx;i+=16){
					for(k=0;k<16;k++){
						shm[lid][k]=*(fh+4*dist3dfh+(i+k)*nLy)-*(fh+4*dist3dfh+(i+k+1)*nLy);	
					}
					__syncthreads();
					for(k=0;k<16;k++){
						*(lrhs+k*nLx)+=shm[ty*16+k][tx];	
					}
					__syncthreads();
					lrhs+=BLOCK_SIZE_xf;
			}
			lrhs=rhs+5*dist3drhs;
    		for(i=0;i<nLx;i+=16){
					for(k=0;k<16;k++){
						shm[lid][k]=*(fh+5*dist3dfh+(i+k)*nLy)-*(fh+5*dist3dfh+(i+k+1)*nLy);	
					}
					__syncthreads();
					for(k=0;k<16;k++){
						*(lrhs+k*nLx)+=shm[ty*16+k][tx];	
					}
					__syncthreads();
					lrhs+=BLOCK_SIZE_xf;
			}
	/*		rhs+=(yid*nLy+xid)*nLx;
    		for(unsigned int i=0;i<nLx;i++){
					*(rhs+i)+=*(fh+i*nLy)-*(fh+(i+1)*nLy);
					*(rhs+dist3drhs+i)+=*(fh+dist3dfh+i*nLy)-*(fh+dist3dfh+(i+1)*nLy);
					*(rhs+2*dist3drhs+i)+=*(fh+2*dist3dfh+i*nLy)-*(fh+2*dist3dfh+(i+1)*nLy);
					*(rhs+3*dist3drhs+i)+=*(fh+3*dist3dfh+i*nLy)-*(fh+3*dist3dfh+(i+1)*nLy);
					*(rhs+4*dist3drhs+i)+=*(fh+4*dist3dfh+i*nLy)-*(fh+4*dist3dfh+(i+1)*nLy);
					*(rhs+5*dist3drhs+i)+=*(fh+5*dist3dfh+i*nLy)-*(fh+5*dist3dfh+(i+1)*nLy);
			}*/
}
