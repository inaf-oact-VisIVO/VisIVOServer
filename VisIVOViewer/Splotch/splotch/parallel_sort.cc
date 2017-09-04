#include "cxxsupport/datatypes.h"
#include "cxxsupport/mpi_support.h"

using namespace std;

template<typename T> void parallel_sort(arr<T> &data)
  {
  data.sort();
  arr<size_t> nlist;
  mpiMgr.allgather(data.size(),nlist);
  tsize nmax=0;
  for (tsize i=0; i<mpiMgr.num_ranks(); ++i)
    nmax = max(nmax,nlist[i]);

  for(tsize ncpu_in_group=2; ncpu_in_group<=mpiMgr.num_ranks(); ncpu_in_group*=2)
    {
    tsize groupnr = mpiMgr.rank()/ncpu_in_group;
    tsize master = ncpu_in_group*groupnr;
    parallel_merge(master, ncpu_in_group, nlist, nmax, data);
    }
  }

template void parallel_sort(arr<double> &data);

template<typename T> void parallel_merge (tsize master, tsize ncpu, const arr<size_t> &nlist,
  size_t nmax, arr<T> &data)
  {
  tsize ndata = data.size();

  if(master + ncpu/2 >= mpiMgr.num_ranks())  /* nothing to do */
    return;

  if (mpiMgr.rank() != master)
    if (ndata>0)
      mpiMgr.sendrecv_replaceRawVoid(&data[0], NAT_CHAR, ndata*sizeof(T),
        master, master);
  else
    {
    arr<T> list_a(nmax), list_b(nmax), list_r(nmax);

    tsize cpua = master,
          cpub = master + ncpu/2,
          cpur = master;

    tsize na=0, nb=0, nr=0;

    memcpy(&list_a[0], &data[0], ndata*sizeof(T));
    if (nlist[cpub])
      mpiMgr.recvRawVoid(&list_b[0], NAT_CHAR, nlist[cpub]*sizeof(T), cpub);

    while (cpur < min(master+ncpu,tsize(mpiMgr.num_ranks())))
      {
      while(na>=nlist[cpua] && cpua<master+ncpu/2-1)
        {
        ++cpua;
        if (nlist[cpua])
          mpiMgr.recvRawVoid(&list_a[0], NAT_CHAR, nlist[cpua]*sizeof(T), cpua);
        na=0;
        }
      while (nb>=nlist[cpub] && cpub<master+ncpu-1 && cpub<mpiMgr.num_ranks()-1)
        {
        ++cpub;
        if (nlist[cpub])
          mpiMgr.recvRawVoid(&list_b[0], NAT_CHAR, nlist[cpub]*sizeof(T), cpub);
        nb=0;
        }

      while (nr>=nlist[cpur])
        {
        if (cpur==master)
          memcpy(&data[0],&list_r[0],nr*sizeof(T));
        else
          if (nlist[cpur])
            mpiMgr.sendRawVoid(&list_r[0], NAT_CHAR, nlist[cpur]*sizeof(T), cpur);
        nr=0;
        ++cpur;
        }

      if (na<nlist[cpua] && nb<nlist[cpub])
        list_r[nr++] = (list_a[na]<list_b[nb]) ? list_a[na++] : list_b[nb++];
      else if (na<nlist[cpua])
        list_r[nr++] = list_a[na++];
      else if (nb<nlist[cpub])
        list_r[nr++] = list_b[nb++];
      }
    }
  }
