/*
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
 * 		Tim Dykes and Ian Cant
 * 		University of Portsmouth
 *
 */

#include "ParticleData.h"
// Need previewer include, in cpp to avoid circular deps
#include "previewer/Previewer.h"

namespace previewer
{

	void ParticleData::Load()
	{
		
		paramfile* splotchParams = Previewer::parameterInfo.GetParamFileReference();
		sceneMaker sMaker(*splotchParams);

		vec3 dummyCam = vec3(0,0,0);
		vec3 dummyAt = vec3(0,0,0);
		vec3 dummySky = vec3(0,0,0);
		std::string dummyOutfile = " ";
	
		// Splotch param find
		bool boost = splotchParams->find<bool>("boost",false);

		// Previewer param find
		//bool boost = param.GetParameter<bool>("boost");

		if(boost)
		{
			//if boost particledata is second argument(boosted data)
			bool ret = sMaker.getNextScene(dummySplotchParticleData, particleList, dummyCam, dummyAt, dummySky, dummyOutfile);
			if(!ret) std::cout << "sMaker.getNextScene failed." << std::endl;
		}
		else 
		{
			//else use full data set (first arg)
			bool ret = sMaker.getNextScene(particleList, dummySplotchParticleData, dummyCam, dummyAt, dummySky, dummyOutfile);
			if(!ret) std::cout << "sMaker.getNextScene failed." << std::endl;		
		}

		// Get colourmap
		int ret = FileLib::LoadColourPalette(Previewer::parameterInfo, colourMaps);
		if(!ret)
			std::cout << "Colour palette in param file is invalid" << std::endl;

		// Get brightness (currently not being used)
		const int numTypes = splotchParams->find<int>("ptypes",1);
		std::vector<float> brightness;

		for(int i = 0; i < numTypes; i++)
		{
			brightness.push_back(splotchParams->find<float>("brightness"+dataToString(i),1.f));
		}
		

		// Copy data into our own structure
		OriginalRGBData.resize(particleList.size());

		for(uint i = 0; i < particleList.size(); i++)
		{
		
			//Hide unsightly large red particles and enlarge all others for smoother blending
			if(particleList[i].type == 1)
				particleList[i].r *= 0.25;
			else particleList[i].r *= splotchParams->find<double>("preview_radial_mod",2);

			// Store orginal colour, to use in regeneration of particle colour with new colourmap
			OriginalRGBData[i].x = particleList[i].e.r;
			OriginalRGBData[i].y = particleList[i].e.g;
			OriginalRGBData[i].z = particleList[i].e.b;

			// Generate colour in same way splotch does (Add brightness here):
			particleList[i].e = colourMaps[particleList[i].type].getVal_const(particleList[i].e.r) * particleList[i].I;
		}

		// Compute and store bounding box of data
		BBox.Compute(particleList);

		DebugPrint("BBox Data");
		DebugPrint("minX: ", BBox.minX);
		DebugPrint("maxX: ", BBox.maxX);
		DebugPrint("miny: ", BBox.minY);
		DebugPrint("maxY: ", BBox.maxY);
		DebugPrint("minZ: ", BBox.minZ);
		DebugPrint("maxZ: ", BBox.maxZ);

		DebugPrint("Particle Data has been loaded");
		DebugPrint("Number of particles to be rendered: ", (int)particleList.size());

	}	

	void ParticleData::ReloadColourData()
	{
		DebugPrint("Reloading colour data");

		// Get new colourmaps
		int ret = FileLib::LoadColourPalette(Previewer::parameterInfo, colourMaps);
		if(!ret)
		{
			std::cout << "Colour palette is invalid" << std::endl;
			return;
		}

		//Reload colour data with new palette
		for(uint i = 0; i < particleList.size(); i++)
		{
			// Use original RGB data as particleList's colour data will have been overwritten
			particleList[i].e = (colourMaps[particleList[i].type].getVal_const(OriginalRGBData[i].x) * particleList[i].I) ;
		}

		DebugPrint("Loaded new colours");

	}

	void ParticleData::SetPalette(std::string paletteFilename, int particleType)
	{
		DebugPrint("Setting new palette");

		paramfile* splotchParams = Previewer::parameterInfo.GetParamFileReference();

		DebugPrint("Got file ref");

		// Check if particle type has vector colour
		if (splotchParams->find<bool>("color_is_vector"+dataToString(particleType),false))
		{
			std::cout << " color of ptype " << particleType << " is vector, so cannot set new colourMap ..." << std::endl;
			return;
		}

		DebugPrint("Checked for vector colour");
 
		// Return if there is no palette parameter available for specified particle type
		// E.G. particle type is invalid
		if(!splotchParams->param_present("palette"+dataToString(particleType)))
		{
			std::cout << "Cannot set new palette for particle type " << particleType << std::endl;
			return;
		}

		DebugPrint("Checked for present parameter");

		// Write new palette filename to parameter file
		splotchParams->setParam<std::string>("palette"+dataToString(particleType),paletteFilename);

		DebugPrint("Set new parameter");
	}

	std::string ParticleData::GetPalette(int particleType)
	{
		paramfile* splotchParams = Previewer::parameterInfo.GetParamFileReference();

		DebugPrint("Got file ref");

		if(!splotchParams->param_present("palette"+dataToString(particleType)))
		{
			std::cout << "No palette in parameter file for particle type: " << particleType << std::endl;
			std::string ret = "NULL";
			return ret;
		}
		else
		{
			return splotchParams->find<std::string>("palette"+dataToString(particleType),"NULL");
		}

	}

	BoundingBox ParticleData::GetBoundingBox() const
	{
		return BBox;
	}

	ParticleList ParticleData::GetParticleList() const
	{
		return particleList;
	}


}