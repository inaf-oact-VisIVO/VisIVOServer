/*
 * Copyright (c) 2004-2010
 *              Martin Reinecke (1), Klaus Dolag (1)
 *               (1) Max-Planck-Institute for Astrophysics
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */
#include <iostream>
#include <cmath>
#include <algorithm>
#include <unistd.h>

#ifdef SPLVISIVO
#include "optionssetter.h"
#endif

#include "splotch/scenemaker.h"
#include "splotch/splotch_host.h"
#include "splotch/splotchutils.h"
#include "cxxsupport/walltimer.h"
#include "cxxsupport/ls_image.h"
#include "cxxsupport/announce.h"

#ifdef CUDA
#include "cuda/splotch_cuda2.h"
#endif
#ifdef OPENCL
#include "opencl/splotch_cuda2.h"
#endif

#ifdef PREVIEWER
#include "previewer/simple_gui/SimpleGUI.h"
#include "previewer/libs/core/FileLib.h"
#include <string>
#include <iostream>
#endif

using namespace std;
#ifdef SPLVISIVO
int splotchMain (VisIVOServerOptions opt)
#else
int main (int argc, const char **argv)
#endif
{
    
#ifdef PREVIEWER
    bool pvMode = false;
    int paramFileArg = -1;
    for(int i = 0; i < argc; i++)
    {
        // Look for command line switch (-pv launches previewer)
        if(std::string(argv[i]) == std::string("-pv"))
        {
            pvMode = true;
        }
        // Look for the parameter file
        if(previewer::FileLib::FileExists(argv[i]))
        {
            paramFileArg = i;
        }
    }
    // Preview mode is enabled
    if(pvMode)
    {
        // If param file exists launch app
        if(paramFileArg > 0)
        {
            // Launch app with simple GUI
            previewer::simple_gui::SimpleGUI program;
            program.Load(std::string(argv[paramFileArg]));
        }
        else
        {
            // Output a message as input param file does not exist
            std::cout << "Invalid input parameter file." << std::endl;
            return 0;
        }
    }
#endif
    
    
    tstack_push("Splotch");
    tstack_push("Setup");
    bool master = mpiMgr.master();
    
    // Prevent the rendering to quit if the "stop" file
    // has been forgotten before Splotch was started
    if (master)
    {
        if (fopen("stop", "r"))
            unlink("stop");
    }
    
#ifdef SPLVISIVO
    planck_assert(!opt.splotchpar.empty(),"usage: --splotch <parameter file>");
    paramfile params (opt.splotchpar.c_str(),false);
    
    params.setParam("camera_x",opt.spPosition[0]);
    params.setParam("camera_y",opt.spPosition[1]);
    params.setParam("camera_z",opt.spPosition[2]);
    
    params.setParam("lookat_x",opt.spLookat[0]);
    params.setParam("lookat_y",opt.spLookat[1]);
    params.setParam("lookat_z",opt.spLookat[2]);
    
    params.setParam("sky_x",opt.spSky[0]);
    params.setParam("sky_y",-opt.spSky[1]);
    params.setParam("sky_z",opt.spSky[2]);
        
    params.setParam("colorbar",opt.showLut);
    //params.setParam("palette0",opt.colorTable);

    
    params.setParam("autocalc_cam",opt.autocalcCam);

    
    params.setParam("outfile",opt.imageName);
    params.setParam("interpolation_mode",0);
    params.setParam("simtype",10);
    params.setParam("fov",opt.fov);
    params.setParam("ptypes",1);
#else
    module_startup ("splotch",argc,argv,2,"<parameter file>",master);
    paramfile params (argv[1],false);
    
    planck_assert(!params.param_present("geometry_file"),
                  "Creating of animations has been changed, use \"scene_file\" instead of "
                  "\"geometry_file\" ... ");
#endif
    
    vector<particle_sim> particle_data; //raw data from file
    vector<particle_sim> r_points;
    
#if (!defined(OPENCL))
    vec3 campos, lookat, sky;
    vector<COLOURMAP> amap;
#else
    ptypes = params.find<int>("ptypes",1);
    g_params =&params;
#endif
#if (defined(CUDA) || defined(OPENCL))
    int myID = mpiMgr.rank();
    int nDevNode = check_device(myID);     // number of GPUs available per node
#ifdef CUDA
    int nDevProc = 1;   // number of GPU required per process
    cout << "Configuration supported is 1 gpu for each mpi process" << endl;
#else
    int nDevProc = params.find<int>("gpu_number",1);
#endif
    int nCoreNode = params.find<int>("cores_number",12); //number of cores per node
    int mydevID = -1;
    // We assume a geometry where processes use only one gpu if available
    if (nDevNode > 0 && nDevProc == 1)
    {
        mydevID = myID%nCoreNode; //ID within the node
        if (mydevID>=nDevNode)
        {
            cout << "There isn't a gpu available for process = " << myID << " computation will be performed on the host" << endl;
            mydevID = -1;
        }
    }
#ifdef OPENCL
    // b) or all processes use a number of GPUs > 1 and <= nDevNode
    else
        planck_assert(nDevNode>=nDevProc, string("Number of GPUs available = ")
                      +dataToString(nDevNode) + " is lower than the number of GPUs required = "
                      +dataToString(nDevProc));
#endif
    bool gpu_info = params.find<bool>("gpu_info",true);
    if (gpu_info)
        if (mydevID >= 0) print_device_info(myID, mydevID);
#endif // CUDA
    
    
#ifdef SPLVISIVO
    get_colourmaps(params,amap,opt);
#else
    get_colourmaps(params,amap);
#endif
    tstack_pop("Setup");
    string outfile;
    
    
#ifdef SPLVISIVO
    sceneMaker sMaker(params,opt);
    while (sMaker.getNextScene (particle_data, r_points, campos, lookat, sky, outfile,opt))
#else
        sceneMaker sMaker(params);
    while (sMaker.getNextScene (particle_data, r_points, campos, lookat, sky, outfile))
#endif
    {
        bool a_eq_e = params.find<bool>("a_eq_e",true);
        int xres = params.find<int>("xres",800),
        yres = params.find<int>("yres",xres);
        arr2<COLOUR> pic(xres,yres);
        
        // calculate boost factor for brightness
        bool boost = params.find<bool>("boost",false);
        float b_brightness = boost ?
        float(particle_data.size())/float(r_points.size()) : 1.0;
        if(particle_data.size()>0)
        {
#if (!defined(CUDA) && !defined(OPENCL))
            if(boost)
            {
#ifdef SPLVISIVO
                host_rendering(params, r_points, pic, campos, lookat, sky, amap, b_brightness, opt);
#else
            host_rendering(params, r_points, pic, campos, lookat, sky, amap, b_brightness);
#endif
        }else
        {
#ifdef SPLVISIVO
            host_rendering(params, particle_data, pic, campos, lookat, sky, amap, b_brightness, opt);
#else
            host_rendering(params, particle_data, pic, campos, lookat, sky, amap, b_brightness);
#endif
        }
#else
            if (mydevID >= 0)
            {
#ifdef CUDA
                if (!a_eq_e) planck_fail("CUDA only supported for A==E so far");
                tstack_push("CUDA");
                if(boost)
                {
#ifdef SPLVISIVO
                    cuda_rendering(mydevID, pic, r_points, campos, lookat, sky, amap, b_brightness, params,opt);
#else
                    cuda_rendering(mydevID, pic, r_points, campos, lookat, sky, amap, b_brightness, params);

#endif
                }else
                {
#ifdef SPLVISIVO
                    cuda_rendering(mydevID, pic, particle_data, campos, lookat, sky, amap, b_brightness, params,opt);
#else
                    cuda_rendering(mydevID, pic, particle_data, campos, lookat, sky, amap, b_brightness, params);

#endif
                }
                tstack_pop("CUDA");
#endif
#ifdef OPENCL
                tstack_push("OPENCL");
#ifdef SPLVISIVO
                opencl_rendering(mydevID, particle_data, nDevProc, pic,opt);
#else
                opencl_rendering(mydevID, particle_data, nDevProc, pic);
#endif
                tstack_pop("OPENCL");
#endif
            }
            else
            {
                if(boost)
                {
#ifdef SPLVISIVO

                    host_rendering(params, r_points, pic, campos, lookat, sky, amap, b_brightness,opt);
#else
                    host_rendering(params, r_points, pic, campos, lookat, sky, amap, b_brightness);
#endif
                }
                
                else
                {
#ifdef SPLVISIVO            
                    host_rendering(params, particle_data, pic, campos, lookat, sky, amap, b_brightness,opt);
#else
                    host_rendering(params, particle_data, pic, campos, lookat, sky, amap, b_brightness);
#endif
                    
            
                }
            }
#endif
        }
        
        
        tstack_push("Post-processing");
        mpiMgr.allreduceRaw
        (reinterpret_cast<float *>(&pic[0][0]),3*xres*yres,MPI_Manager::Sum);
        
        exptable<float32> xexp(-20.0);
        if (mpiMgr.master() && a_eq_e)
            for (int ix=0;ix<xres;ix++)
                for (int iy=0;iy<yres;iy++)
                {
                    pic[ix][iy].r=-xexp.expm1(pic[ix][iy].r);
                    pic[ix][iy].g=-xexp.expm1(pic[ix][iy].g);
                    pic[ix][iy].b=-xexp.expm1(pic[ix][iy].b);
                }
        
        tstack_replace("Post-processing","Output");
        
        if (master && params.find<bool>("colorbar",false))
        {
            cout << endl << "creating color bar ..." << endl;
            add_colorbar(params,pic,amap);
        }
        
        if(!params.find<bool>("AnalyzeSimulationOnly"))
        {
            if (master)
            {
                cout << endl << "saving file ..." << endl;
                
                LS_Image img(pic.size1(),pic.size2());
                
                for (tsize i=0; i<pic.size1(); ++i)
                    for (tsize j=0; j<pic.size2(); ++j)
                        img.put_pixel(i,j,Colour(pic[i][j].r,pic[i][j].g,pic[i][j].b));
                int pictype = params.find<int>("pictype",0);
                switch(pictype)
                {
                    case 0:
                        img.write_TGA(outfile+".tga");
                        break;
                    case 1:
                        planck_fail("ASCII PPM no longer supported");
                        break;
                    case 2:
                        img.write_PPM(outfile+".ppm");
                        break;
                    case 3:
                        img.write_TGA_rle(outfile+".tga");
                        break;
                    default:
                        planck_fail("No valid image file type given ...");
                        break;
                }
            }
        }
        
        tstack_pop("Output");

        // Also meant to happen when using CUDA - unimplemented.
#if (defined(OPENCL))
        cuda_timeReport();
#endif
        timeReport();
        
        mpiMgr.barrier();
        // Abandon ship if a file named "stop" is found in the working directory.
        // ==>  Allows to stop rendering conveniently using a simple "touch stop".
        if (fopen("stop", "r"))
        {
            mpiMgr.~MPI_Manager();
            break;
        }
        
    }

#ifdef VS
    //Just to hold the screen to read the messages when debugging
    cout << endl << "Press any key to end..." ;
    getchar();
#endif
}
