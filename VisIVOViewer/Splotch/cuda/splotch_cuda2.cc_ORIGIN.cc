// each mpi_task manage 1 GPU
#include "cuda/splotch_cuda2.h"
#include "cuda/splotch_cuda.h"
#include "cxxsupport/string_utils.h"
#include "cuda/CuPolicy.h"
#include "cuda/CuRender.h"

using namespace std;
#ifdef SPLVISIVO
void cuda_rendering(int mydevID, arr2<COLOUR> &pic, vector<particle_sim> &particle, const vec3 &campos, const vec3 &lookat, vec3 &sky, vector<COLOURMAP> &amap, float b_brightness, paramfile &g_params,VisIVOServerOptions opt)
#else
void cuda_rendering(int mydevID, arr2<COLOUR> &pic, vector<particle_sim> &particle, const vec3 &campos, const vec3 &lookat, vec3 &sky, vector<COLOURMAP> &amap, float b_brightness, paramfile &g_params)
#endif
{
    tstack_push("CUDA");
    tstack_push("Device setup");
    long int nP = particle.size();
    pic.fill(COLOUR(0.0, 0.0, 0.0));
    int xres = pic.size1();
    int yres = pic.size2();
    cout << "resolution = " << xres << " x " << yres << endl;
    arr2<COLOUR> Pic_host(xres,yres);
    int ptypes = g_params.find<int>("ptypes",1);
    
    // Initialize policy class
    CuPolicy *policy = new CuPolicy(xres, yres, g_params);
    int ntiles = policy->GetNumTiles();
    
    // num particles to manage at once
    float factor = g_params.find<float>("particle_mem_factor", 3);
    long int len = cu_get_chunk_particle_count(policy, sizeof(cu_particle_sim), ntiles, factor);
    
    if (len <= 0)
    {
        cout << "Graphics memory setting error" << endl;
        mpiMgr.abort();
    }
    
    //CUDA Init
    cu_gpu_vars gv;
    memset(&gv, 0, sizeof(cu_gpu_vars));
    gv.policy = policy;
    // enable device and allocate arrays
    bool doLogs = true;
    int error = cu_init(mydevID, len, ntiles, &gv, g_params, campos, lookat, sky, b_brightness, doLogs);
    tstack_pop("Device setup");
    if (!error)
    {
        //a new linear pic object that will carry the result
        setup_colormap(ptypes, amap, &gv);
        float64 grayabsorb = g_params.find<float>("gray_absorption",0.2);
        bool a_eq_e = g_params.find<bool>("a_eq_e",true);
        
        int endP = 0;
        int startP = 0;
        int nPR = 0;
        
        while(endP < nP)
        {
            endP = startP + len;   //set range
            if (endP > nP) endP = nP;
            nPR += cu_draw_chunk(mydevID, (cu_particle_sim *) &(particle[startP]), endP-startP, Pic_host, &gv, a_eq_e, grayabsorb, xres, yres, doLogs);
            // combine host results of chunks
            tstack_push("combine images");
            for (int x=0; x<xres; x++)
                for (int y=0; y<yres; y++)
                    pic[x][y] += Pic_host[x][y];
            tstack_pop("combine images");
            cout << "Rank " << mpiMgr.rank() << ": Rendered " << nPR << "/" << nP << " particles" << endl << endl;
            startP = endP;
        }
        add_device_image(pic, &gv, xres, yres);
        tstack_pop("CUDA");
        cu_end(&gv);
    }
    
}


void setup_colormap(int ptypes, vector<COLOURMAP> &amap, cu_gpu_vars* gv)
{
    //init C style colormap
    cu_color_map_entry *amapD;//amap for Device
    int *amapDTypeStartPos; //begin indexes of ptypes in the linear amapD[]
    amapDTypeStartPos =new int[ptypes];
    int curPtypeStartPos =0;
    int size =0;
    //first we need to count all the entries to get colormap size
    for (int i=0; i<amap.size(); i++)
        size += amap[i].size();
    
    //then fill up the colormap amapD
    amapD =new cu_color_map_entry[size];
    int j,index =0;
    for(int i=0; i<amap.size(); i++)
    {
        for (j=0; j<amap[i].size(); j++)
        {
            amapD[index].val = amap[i].getX(j);
            COLOUR c (amap[i].getY(j));
            amapD[index].color.r = c.r;
            amapD[index].color.g = c.g;
            amapD[index].color.b = c.b;
            index++;
        }
        amapDTypeStartPos[i] = curPtypeStartPos;
        curPtypeStartPos += j;
    }
    //now let cuda init colormap on device
    cu_colormap_info clmp_info;
    clmp_info.map =amapD;
    clmp_info.mapSize =size;
    clmp_info.ptype_points =amapDTypeStartPos;
    clmp_info.ptypes =ptypes;
    cu_init_colormap(clmp_info, gv);
    
    delete []amapD;
    delete []amapDTypeStartPos;
}








/*









#ifndef NO_WIN_THREAD
#include <pthread.h>
#endif

#include "cuda/splotch_cuda2.h"
#include "cxxsupport/string_utils.h"
#include "cuda/CuPolicy.h"

using namespace std;

paramfile *g_params;
int ptypes = 0;
vector<particle_sim> particle_data;   //raw data from file
vec3 campos, lookat, sky;
vector<COLOURMAP> amap;
wallTimerSet cuWallTimers;
#ifdef SPLVISIVO
VisIVOServerOptions opt;
#endif

#ifdef SPLVISIVO
void cuda_rendering(int mydevID, int nDev, arr2<COLOUR> &pic,VisIVOServerOptions &opt2)
#else
void cuda_rendering(int mydevID, int nDev, arr2<COLOUR> &pic)
#endif
  {
  #ifdef SPLVISIVO
	opt=opt2;	
  #endif
  //see if host must be a working thread
  bool bHostThread = g_params->find<bool>("use_host_as_thread", false);
  int nThread = bHostThread? nDev+1: nDev;
  //init array info for threads control
  thread_info *tInfo = new thread_info[nThread];
  tInfo[0].pPic = &pic;      //local var pic is assigned to the first thread
  tInfo[0].devID = mydevID;
//  tInfo[0].npart_all = npart_all;
  for (int i=1; i<nDev; i++)
    {
    tInfo[i].devID = mydevID+i;
//    tInfo[i].npart_all = npart_all;
    tInfo[i].pPic = new arr2<COLOUR>(pic.size1(), pic.size2());
    }
  //make the last one work for host thread
  if (bHostThread)
    {
    tInfo[nThread-1].devID =-1;
//    tInfo[nThread-1].npart_all = npart_all;
    if (nThread-1 != 0)
      tInfo[nThread-1].pPic = new arr2<COLOUR>(pic.size1(), pic.size2());
    }
  //decide how to divide particles range by another function
  DevideThreadsTasks(tInfo, nThread, bHostThread);

#ifndef NO_WIN_THREAD // create cuda threads on Windows using CreateThread function
  HANDLE *tHandle = new HANDLE[nThread];
  //issue the threads
  for (int i=0; i<nDev; i++)
    tHandle[i] = CreateThread( NULL, 0,
      (LPTHREAD_START_ROUTINE)cu_thread_func,&(tInfo[i]), 0, NULL );
  //issue the host thread too
  if (bHostThread)
    tHandle[nDev] = CreateThread( NULL, 0,
      (LPTHREAD_START_ROUTINE)host_thread_func,&(tInfo[nDev]), 0, NULL );
  WaitForMultipleObjects(nThread, tHandle, true, INFINITE);

#else // create cuda threads on Linux using pthread_create function

//  planck_assert(nDev <= 1, "can't have multiple cuda threads on Linux (yet), so 'gpu_number' must be 1");
  pthread_t *tHandle = new pthread_t[nThread];
  for (int i=0; i<nDev; i++)
     pthread_create(&(tHandle[i]), NULL, cu_thread_func, (void *) &(tInfo[i]) );
  if (bHostThread)
     pthread_create(&(tHandle[nDev]), NULL, host_thread_func, (void *) &(tInfo[nDev]) );
  void *status[nThread];
  for (int i=0; i <nThread; ++i) pthread_join(tHandle[i], &status[i]);
//  cu_thread_func (&(tInfo[0])); //just call it as normal function
#endif  //if not NO_WIN_THREAD

  // combine the results of multiple threads(devices + host) to pic
  for (int i=1; i<nThread; i++)
      for (int x=0; x<pic.size1(); x++)
        for (int y=0; y<pic.size2(); y++)
              pic[x][y] = pic[x][y] + (*tInfo[i].pPic)[x][y];

  if (g_params->getVerbosity())
   for (int i=0; i<nThread; i++)
    {
    if (tInfo[i].devID!=-1)
      {
      cout<< endl <<"Rank " << mpiMgr.rank() << ": Times of GPU" << i << ":" <<endl;
      GPUReport(tInfo[i].times);
      cout<<endl;
      }
    }
  
  if (mpiMgr.master()) cuWallTimers = tInfo[0].times;

  for (int i=1; i<nThread; i++)
    delete tInfo[i].pPic;
  delete [] tInfo;
  delete [] tHandle;
  }


void DevideThreadsTasks(thread_info *tInfo, int nThread, bool bHostThread)
  {
  bool bTestLoadBalancing = g_params->find<bool>("test_load_balancing", false);
  unsigned int curStart = 0;
  int hostLoad = bHostThread? g_params->find<int>("host_load",0): 0;
  int nDev = bHostThread? nThread-1: nThread;
  int onePercent = particle_data.size()/100;
  int averageDevLen = (nDev!=0)? onePercent *(100-hostLoad)/nDev : 0;

  for (int i=0; i<nThread; i++)
    {
    tInfo[i].startP = curStart;
    if (tInfo[i].devID != -1) //not a host
      {
      if (bTestLoadBalancing)
        {
        int gpuLoad = g_params->find<int>("gpu_load"+dataToString(i),0);
        tInfo[i].endP = curStart + gpuLoad * onePercent - 1;
        }
      else
        tInfo[i].endP = curStart + averageDevLen - 1;
      }
    else //if this is a host
      {
      tInfo[i].endP = curStart + hostLoad * onePercent - 1;
      }
    curStart = tInfo[i].endP + 1;
    }

  tInfo[nThread-1].endP = particle_data.size()-1;
  }



THREADFUNC cu_thread_func(void *pinfo)
 {
  //a new thread info object that will carry each chunk's drawing
  thread_info *pInfoOutput = (thread_info*) pinfo;
  thread_info ti = *pInfoOutput;

  //do some cleaning for final thread_info
  pInfoOutput->pPic->fill(COLOUR(0.0, 0.0, 0.0));

  // Initialize policy class
  CuPolicy *policy = new CuPolicy(*g_params); 

  //a new pic object residing in ti that will carry the result
  arr2<COLOUR> pic(pInfoOutput->pPic->size1(), pInfoOutput->pPic->size2());
  ti.pPic = &pic;

  // num particles to manage at once
  float factor = g_params->find<float>("particle_mem_factor", 3);
  int len = cu_get_chunk_particle_count(*g_params, policy, sizeof(cu_particle_sim), factor); 
  if (len == 0)
    {
    printf("\nGraphics memory setting error\n");
    mpiMgr.abort();
    }

  //CUDA Init
  cu_gpu_vars gv; //for each gpu/thread a variable pack is needed
  memset(&gv, 0, sizeof(cu_gpu_vars));
  gv.policy = policy;
  // enable device and allocate arrays
  cu_init(pInfoOutput->devID, len, &gv, *g_params, campos, lookat, sky);

  //CUDA Coloring
  setup_colormap(ptypes, &gv);

  int endP = ti.endP;
  ti.endP = ti.startP;
  while(ti.endP < endP)
    {
    ti.endP =ti.startP + len - 1;   //set range
    if (ti.endP > endP) ti.endP = endP; 
    cu_draw_chunk(&ti, &gv);
    // combine results of chunks
    pInfoOutput->times.start("gcombine");
    for (int x=0; x<pic.size1(); x++)
      for (int y=0; y<pic.size2(); y++)
        (*(pInfoOutput->pPic))[x][y] += pic[x][y];
    pInfoOutput->times.stop("gcombine");
    ti.startP = ti.endP + 1;
    }
  pInfoOutput->times = ti.times;

  cu_end(&gv);
 }


void cu_draw_chunk(void *pinfo, cu_gpu_vars* gv)
  {

  //get the input info
  thread_info *tInfo = (thread_info*)pinfo;
  tInfo->times.start("gpu_thread");

  int nParticle = tInfo->endP - tInfo->startP + 1;
  printf("Rank %d - GPU %d : Processing %d particles\n", mpiMgr.rank(), tInfo->devID, nParticle); fflush(stdout);

  paramfile &params(*g_params);

  //copy data address to local C-like array pointer d_particle_data
  tInfo->times.start("gcopy");
  cu_particle_sim *d_particle_data = &(particle_data[tInfo->startP]);
  cu_copy_particles_to_device(d_particle_data, nParticle, gv);
  tInfo->times.stop("gcopy");

  //init cu_particle_splotch array memory
  cu_particle_splotch *cu_ps;
  cu_ps = new cu_particle_splotch[nParticle];
  memset(cu_ps, 0, nParticle);

  //CUDA Transformation
  tInfo->times.start("gtransform");
  cu_transform(nParticle, cu_ps, gv);
  tInfo->times.stop("gtransform");

/* temporarily ignore sorting 191109.
   it becomes complicated when using multiple threads with sorting
   //then copy particles back to host for sorting
   for (int i=0; i<nParticle; i++)
   memcpy( &(particle_data[iWRONG]),&(d_particle_data[i]), sizeof(cu_particle_sim));
PROBLEM HERE!

// --------------------------------
// ----------- Sorting ------------
// --------------------------------
// cout << endl << "applying sort (" << npart << ") ..." << endl;

   int sort_type = params.find<int>("sort_type",1);
   particle_sort(particle_data,sort_type,true);
     //we can sort by size(r) for balancing with cuda
     //particle_sort(particle_data,4,true);

     //copy sorted data back to device, not a must!
     //first to C-style object
     for(int i=0; i<particle_data.size(); i++)
     memcpy( &(d_particle_data[i]), &(particle_data[i]), sizeof(cu_particle_sim));
    cu_copy_particle_sim_to_device(d_particle_data, particle_data.size());
*/
 

// ----------------------------------
// ----------- Rendering ------------
// ----------------------------------
/*
  //get parameters for rendering
  int xres = params.find<int>("xres",800),
      yres = params.find<int>("yres",xres);
//  long nsplotch=pFiltered;
//  long nsplotch_all=nsplotch;
//  mpiMgr.allreduce(nsplotch_all,MPI_Manager::Sum);
  float64 grayabsorb = params.find<float>("gray_absorption",0.2);
  bool a_eq_e = params.find<bool>("a_eq_e",true);

  //prepare fragment buffer memory space first
  void *fragBuf;
  size_t nFBufInByte = gv->policy->GetFBufSize() <<20;

  int nFBufInCell;
  if (a_eq_e)
    {
    nFBufInCell = nFBufInByte/sizeof(cu_fragment_AeqE);
    fragBuf = new cu_fragment_AeqE[nFBufInCell];
    }
  else
    {
    nFBufInCell = nFBufInByte/sizeof(cu_fragment_AneqE);
    fragBuf = new cu_fragment_AneqE[nFBufInCell];
    }

  int maxRegion = gv->policy->GetMaxRegion();
  int chunk_dim = nParticle;
  // new array of particles produced after filter and splitting
  cu_particle_splotch *cu_ps_filtered;
  cu_ps_filtered = new cu_particle_splotch[chunk_dim];

  //clear the output pic
  tInfo->pPic ->fill(COLOUR(0.0, 0.0, 0.0));
 
  int End_cu_ps = 0, Start_cu_ps=0;
  while (End_cu_ps < nParticle)
  {
   int nFragments2Render = 0;
   //filter and split particles to a cu_ps_filtered
   tInfo->times.start("gfilter");
   int pFiltered = filter_chunk(Start_cu_ps, chunk_dim, nParticle, maxRegion, 
                                nFBufInCell, cu_ps, cu_ps_filtered, &End_cu_ps,
                                &nFragments2Render);
   tInfo->times.stop("gfilter");

   tInfo->times.start("gcopy");
   cu_copy_particles_to_render(cu_ps_filtered, pFiltered, gv);
   tInfo->times.start("gcopy");

   // render chunks of pFiltered particles
//   tInfo->times.start("gcoloring");  moved inside the render kernel
//   cu_colorize(pFiltered, gv);
//   tInfo->times.stop("gcoloring");

   tInfo->times.start("grender");
   cu_render1(pFiltered, a_eq_e, (float) grayabsorb, gv);
   tInfo->times.stop("grender");

    //collect result
    tInfo->times.start("gcopy-fbuf");
    cu_get_fbuf(fragBuf, a_eq_e, nFragments2Render, gv);
    tInfo->times.stop("gcopy-fbuf");
 
    //combine chunks  
    //cu_ps_filtered:       the particle array
    tInfo->times.start("gcombine");
    combine_chunk(0, pFiltered-1, cu_ps_filtered, fragBuf, a_eq_e, grayabsorb, *(tInfo->pPic));
    tInfo->times.stop("gcombine");
 /*  render_chunk(pFiltered, nFBufInCell, cu_ps_filtered,
                fragBuf, gv, a_eq_e, grayabsorb, *(tInfo->pPic), tInfo->times);*/
7
printf("Rank %d - GPU %d : Rendered %d/%d particles \n",  mpiMgr.rank(), tInfo->devID, End_cu_ps, nParticle);
   
   Start_cu_ps = End_cu_ps;
  }

  delete []cu_ps;
  delete []cu_ps_filtered;
  if (a_eq_e)
    delete [] ((cu_fragment_AeqE *) fragBuf);
  else
    delete [] ((cu_fragment_AneqE *) fragBuf);

  tInfo->times.stop("gpu_thread");
  }


//filter and split particles to a cu_ps_filtered of size nParticles
int filter_chunk(int StartP, int chunk_dim, int nParticle, int maxRegion, 
                 int nFBufInCell, cu_particle_splotch *cu_ps, 
                 cu_particle_splotch *cu_ps_filtered, int *End_cu_ps, 
                 int *nFragments2Render)
{
  cu_particle_splotch p, pNew;
  int region, nsplit;
  bool finished = false;

  unsigned long posInFragBuf = 0; 
  int pFiltered = 0;
  int i=StartP;  // start chunk position in cu_ps

  // filter particles until cu_ps is finished or cu_ps_filtered array is full
  while(!finished && (i < nParticle))
   {
     //select valid ones
     p = cu_ps[i];
     if (p.isValid)
     {
       int h = p.maxy - p.miny;
       int w = p.maxx - p.minx;
       region = h*w;

       if (region <= maxRegion)
//       if(p.r <= 32.0)
       {
         // set the start position of the particle in fragment buffer
         if ((pFiltered + 1 <= chunk_dim) && (posInFragBuf+region < nFBufInCell)) 
         {
           p.posInFragBuf = posInFragBuf;
           cu_ps_filtered[pFiltered] = p;
           pFiltered++;
           posInFragBuf += region;
         }
         else finished = true; 
       }
       else
       { //particle too big -> split along y direction
         pNew = p;
         int w1 = (maxRegion%h == 0) ? (maxRegion/h):(maxRegion/h + 1);
         nsplit = w/w1 + 1;
         if ((pFiltered + nsplit <= chunk_dim) && (posInFragBuf+region < nFBufInCell))
         {
           for (int minx = p.minx; minx < p.maxx; minx += w1)
           {
             pNew.minx = minx;  //minx,maxx of pNew need to be set
             pNew.maxx = (minx+w1 >= p.maxx) ? p.maxx : minx+w1; 
             // set the start position of the particle in fragment buffer
             pNew.posInFragBuf = posInFragBuf; 
             cu_ps_filtered[pFiltered] = pNew;

             pFiltered++;
             int newregion = (pNew.maxx - pNew.minx) * (pNew.maxy - pNew.miny);
             posInFragBuf += newregion;
           }
         }
         else finished = true; 
       }
      }
      i++;
    }

   *End_cu_ps = i;
   *nFragments2Render = posInFragBuf;
   return pFiltered;  // return chunk position reached in cu_ps
}


void combine_chunk(int StartP, int EndP, cu_particle_splotch *cu_ps_filtered, 
                  void *fragBuf, bool a_eq_e, float64 grayabsorb, arr2<COLOUR> &pPic)
{

    if (a_eq_e)
    {
      cu_fragment_AeqE *fragBufAeqE = (cu_fragment_AeqE *)fragBuf;
      for (int pPos=StartP, fPos=0; pPos<EndP; pPos++)
      {
        for (int x =cu_ps_filtered[pPos].minx; x <cu_ps_filtered[pPos].maxx; x++)
        {
          for (int y =cu_ps_filtered[pPos].miny; y <cu_ps_filtered[pPos].maxy; y++)
          {
            pPic[x][y].r += fragBufAeqE[fPos].aR;
            pPic[x][y].g += fragBufAeqE[fPos].aG;
            pPic[x][y].b += fragBufAeqE[fPos].aB;
            fPos ++;
          }
        }
      }
    }
    else
    {
      cu_fragment_AneqE *fragBufAneqE = (cu_fragment_AneqE *)fragBuf;
      for (int pPos=StartP, fPos=0; pPos<EndP; pPos++)
      {
        for (int x =cu_ps_filtered[pPos].minx; x <cu_ps_filtered[pPos].maxx; x++)
        {
          for (int y =cu_ps_filtered[pPos].miny; y <cu_ps_filtered[pPos].maxy; y++) 
          {
            pPic[x][y].r += fragBufAneqE[fPos].aR *
				   (pPic[x][y].r - fragBufAneqE[fPos].qR);
            pPic[x][y].g += fragBufAneqE[fPos].aG *
				   (pPic[x][y].g - fragBufAneqE[fPos].qG); 
            pPic[x][y].b += fragBufAneqE[fPos].aB *
                                   (pPic[x][y].b - fragBufAneqE[fPos].qB);
            fPos ++;
          }
        }
      }
    }
}

void setup_colormap(int ptypes, cu_gpu_vars* gv)
{
//init C style colormap
  cu_color_map_entry *amapD;//amap for Device
  int *amapDTypeStartPos; //begin indexes of ptypes in the linear amapD[]
  amapDTypeStartPos =new int[ptypes];
  int curPtypeStartPos =0;
  int size =0;
  //first we need to count all the entries to get colormap size
  for (int i=0; i<amap.size(); i++)
    size += amap[i].size();

  //then fill up the colormap amapD
  amapD =new cu_color_map_entry[size];
  int j,index =0;
  for(int i=0; i<amap.size(); i++)
    {
    for (j=0; j<amap[i].size(); j++)
      {
      amapD[index].val = amap[i].getX(j);
      COLOUR c (amap[i].getY(j));
      amapD[index].color.r = c.r;
      amapD[index].color.g = c.g;
      amapD[index].color.b = c.b;
      index++;
      }
    amapDTypeStartPos[i] = curPtypeStartPos;
    curPtypeStartPos += j;
    }
  //now let cuda init colormap on device
  cu_colormap_info clmp_info;
  clmp_info.map =amapD;
  clmp_info.mapSize =size;
  clmp_info.ptype_points =amapDTypeStartPos;
  clmp_info.ptypes =ptypes;
  cu_init_colormap(clmp_info, gv);

  delete []amapD;
  delete []amapDTypeStartPos;
}



THREADFUNC host_thread_func(void *p)
  {
  thread_info *tInfo = (thread_info*)p;

  vector<particle_sim>::iterator i1,i2;
  i1 =particle_data.begin() + tInfo->startP;
  i2 =particle_data.begin() + tInfo->endP + 1;
  vector<particle_sim> particles(i1,i2);
#ifdef SPLVISIVO
  host_rendering(*g_params, particles, *(tInfo->pPic), campos, lookat, sky, amap,0.0,opt);
#else
  host_rendering(*g_params, particles, *(tInfo->pPic), campos, lookat, sky, amap,0.0);
#endif

  }


void GPUReport(wallTimerSet &cuTimers)
  {
    cout << "Copy  (secs)               : " << cuTimers.acc("gcopy") << endl;
    cout << "Copy-fbuf  (secs)          : " << cuTimers.acc("gcopy-fbuf") << endl;
    cout << "Transforming Data (secs)   : " << cuTimers.acc("gtransform") << endl;
//    cout << "Sorting Data (secs)        : " << cuTimers.acc("gsort") << endl;
//    cout << "Coloring Sub-Data (secs)   : " << cuTimers.acc("gcoloring") << endl;
    cout << "Filter Sub-Data (secs)     : " << cuTimers.acc("gfilter") << endl;
    cout << "Rendering Sub-Data (secs)  : " << cuTimers.acc("grender") << endl;
    cout << "Combine Sub-image (secs)   : " << cuTimers.acc("gcombine") << endl;
    cout << "Cuda thread (secs)         : " << cuTimers.acc("gpu_thread") << endl;
  }

void cuda_timeReport(paramfile &params)
  {
  if (mpiMgr.master())
    {
    cout << endl << "--------------------------------------------" << endl;
    cout << "Summary of timings" << endl;
    cout << "--------------------------------------------" << endl;
    cout<< endl <<"Times of GPU:" <<endl;
    GPUReport (cuWallTimers);
    cout <<  "--------------------------------------------" << endl;

    if (params.find<bool>("use_host_as_thread", false))
      {
      cout<< endl <<"Times of CPU HOST as threads:" <<endl;
      hostTimeReport(wallTimers);
      cout << "Host thread (secs)         : " << wallTimers.acc("host_thread") << endl;
      cout << "--------------------------------------------" << endl;
      }
    cout << "Setup Data (secs)          : " << wallTimers.acc("setup") << endl;
    cout << "Read Data (secs)           : " << wallTimers.acc("read") << endl;
    cout << "Ranging Data (secs)        : " << wallTimers.acc("range") << endl;
    cout << "Postprocessing (secs)      : " << wallTimers.acc("postproc") << endl;
    cout << "Write Data (secs)          : " << wallTimers.acc("write") << endl;
    cout << "Total (secs)               : " << wallTimers.acc("full") << endl;
    }
  }
*/
