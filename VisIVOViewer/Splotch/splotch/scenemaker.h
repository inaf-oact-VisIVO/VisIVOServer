#ifndef SPLOTCH_SCENE_MAKER_H
#define SPLOTCH_SCENE_MAKER_H

#include <vector>
#include <string>
#include <fstream>

#include "cxxsupport/paramfile.h"
#include "splotch/splotchutils.h"

#ifdef SPLVISIVO
#include "optionssetter.h"
#endif

class sceneMaker
  {
  private:
    struct scene
      {
      vec3 campos, lookat, sky;
      double fidx;
      std::string outname;
      bool keep_particles, reuse_particles;
      scene (const vec3 &c, const vec3 &l, const vec3 &s, double fdx,
             const std::string &oname, bool keep, bool reuse)
        : campos(c), lookat(l), sky(s), fidx(fdx), outname(oname),
          keep_particles(keep), reuse_particles(reuse) {}
	
      };

    std::vector<scene> scenes;
    int cur_scene;

    paramfile &params;

    int interpol_mode;

    double boxsize;
// only used if interpol_mode>0
    std::vector<particle_sim> p1,p2;
    std::vector<MyIDType> id1,id2;
    std::vector<uint32> idx1,idx2;
    int snr1_now,snr2_now;
    double time1,time2;
// only used if interpol_mode>1
    std::vector<vec3f> vel1,vel2;

// only used if interpol_mode>0
    void particle_interpolate(std::vector<particle_sim> &p, double frac);

// only used if the same particles are used for more than one scene
    std::vector<particle_sim> p_orig;

#ifdef SPLVISIVO
   void fetchFiles(std::vector<particle_sim> &particle_data, double fidx,VisIVOServerOptions &opt);
#else
    void fetchFiles(std::vector<particle_sim> &particle_data, double fidx);
#endif
  public:
#ifdef SPLVISIVO
    sceneMaker (paramfile &par, VisIVOServerOptions &opt);
    bool getNextScene (std::vector<particle_sim> &particle_data, 
                       std::vector<particle_sim> &r_points, vec3 &campos,
                       vec3 &lookat, vec3 &sky, std::string &outfile, VisIVOServerOptions &opt);
#else
    sceneMaker (paramfile &par);
    bool getNextScene (std::vector<particle_sim> &particle_data, 
                       std::vector<particle_sim> &r_points, vec3 &campos,
                       vec3 &lookat, vec3 &sky, std::string &outfile);
#endif
  };

#endif
