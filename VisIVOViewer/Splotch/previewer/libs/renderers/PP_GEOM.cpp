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
 
#include "PP_GEOM.h"
 #include "previewer/Previewer.h"

namespace previewer
{
	void PP_GEOM::Load(const ParticleData& pData)
	{
		//Store passed in data and bounding box
		particleList = pData.GetParticleList();
		dataBBox = pData.GetBoundingBox();

		float maxI = 0;
		float minI = 100;
		float maxR = 0; 
		float minR = 100;
		// Check max/mins
		for(unsigned int i = 0; i < particleList.size(); i++)
		{
			// Check for maximuma nd minimum radii/intensities
			maxI = ((particleList[i].I > maxI) ? particleList[i].I : maxI);
			minI = ((particleList[i].I < minI) ? particleList[i].I : minI);
			maxR = ((particleList[i].r > maxR) ? particleList[i].r : maxR);
			minR = ((particleList[i].r < minR) ? particleList[i].r : minR);
			// Compress colours to avoid saturation when drawing large numbers of particles
			particleList[i].e *= 0.075;
		}

		std::cout << "Min I: " << minI << " Max I: " << maxI << std::endl;
		std::cout << "Min R: " << minR << " Max R: " << maxR << std::endl;

		// Create the VBO for particle drawing
		genVBO();

		PrintOpenGLError();

		// Set up GUI here if necessary...

		// Create material
		material = new PP_ParticleMaterial();

		// Load shader program, with geometry shader
		material->Load("PP_GEOM", true);

		// Set up rest of material
		material->SetBlend(true);
		material->SetBlendSrc(GL_SRC_ALPHA);
		material->SetBlendDst(GL_ONE);
		material->SetTexture(true);
		material->LoadTexture("previewer/data/textures/particle.tga", GL_TEXTURE_2D);

		// Load and position camera
		camera.SetPerspectiveProjection(ParticleSimulation::GetFOV(), ParticleSimulation::GetAspectRatio(), 1, 200000);

		// Check if recalculation is required
		paramfile* splotchParams = Previewer::parameterInfo.GetParamFileReference();
		bool recalc = splotchParams->find<bool>("recalc_preview_cam",true);
		if(recalc)
			camera.Create(dataBBox);
		else
		{
			vec3f lookat(splotchParams->find<double>("lookat_x"),splotchParams->find<double>("lookat_y"),splotchParams->find<double>("lookat_z"));
			vec3f sky(splotchParams->find<double>("sky_x",0),splotchParams->find<double>("sky_y",0),splotchParams->find<double>("sky_z",1));
			vec3f campos(splotchParams->find<double>("camera_x",0),splotchParams->find<double>("camera_y",0),splotchParams->find<double>("camera_z",1));

			camera.Create(campos,lookat,(sky*=-1));
		}

		camera.SetMainCameraStatus(true);

		std::cout <<  " Renderer loaded " << std::endl;
		
	}

	void PP_GEOM::Draw()
	{
		glEnable(GL_SCISSOR_TEST);

		glViewport(0,0, ParticleSimulation::GetSimWindowWidth(),ParticleSimulation::GetSimWindowHeight());
		glScissor(0,0, ParticleSimulation::GetSimWindowWidth(),ParticleSimulation::GetSimWindowHeight());
		Clear(0.2,0.2,0.2,1.0);

		// Set 3d rendering viewport size 
		float viewXmin = ParticleSimulation::GetRenderXMin();
		float viewYmin = ParticleSimulation::GetRenderYMin();

		// Set 3d render viewport to scissor and clear
		glViewport(viewXmin,viewYmin, ParticleSimulation::ParticleSimulation::GetRenderWidth(), ParticleSimulation::GetRenderHeight());
		glScissor(viewXmin,viewYmin, ParticleSimulation::GetRenderWidth(),ParticleSimulation::GetRenderHeight());

		glClearColor(0.0, 0.0, 0.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// Bind material
		material->Bind(camera.GetMVPMatrix());


		// Bind VBO
		glBindBuffer(GL_ARRAY_BUFFER, VBO);

		// Set vertex/colour pointers (using offsets for interleaved data)
  		glVertexPointer( 3, GL_FLOAT, byteInterval, (void*)sizeof(particleList[0].e) );
 		glColorPointer(  3, GL_FLOAT, byteInterval, 0 );
 		// Set normal pointer (used to pass extra particle data for use in shader)
 		glNormalPointer(GL_FLOAT, byteInterval, (void*)(sizeof(particleList[0].e) +  sizeof(particleList[0].x)*3));

 		// Enable client states
	    glEnableClientState(GL_VERTEX_ARRAY);
	    glEnableClientState(GL_COLOR_ARRAY);
	    glEnableClientState(GL_NORMAL_ARRAY);

	   	// Draw
	    glDrawArrays(GL_POINTS, 0, particleList.size() );

	    //Disable client states
	    glDisableClientState(GL_VERTEX_ARRAY);
	    glDisableClientState(GL_COLOR_ARRAY);
	    glDisableClientState(GL_NORMAL_ARRAY);

		// Unbind vbo
		glBindBuffer(GL_ARRAY_BUFFER, 0);		

		// Unbind material
		material->Unbind();

		//Check for error
		//PrintOpenGLError();

		// Draw GUI

		// Draw labels
		glDisable(GL_SCISSOR_TEST);
	}

	void PP_GEOM::Unload()
	{

	}

	void PP_GEOM::Update()
	{
		camera.SetPerspectiveProjection(ParticleSimulation::GetFOV(), ParticleSimulation::GetAspectRatio(), 1, 200000);
	}

	void PP_GEOM::OnKeyPress(Event ev)
	{
		// Movement * seconds per frame for uniform movement regardless of framerate 
		// (* 100 / *0.1 to make scale similar for setters move and rotation)
		if(ev.keyID == "w")
				camera.MoveForward(ParticleSimulation::GetMoveSpeed()*100*WindowManager::GetSPF());

		else if(ev.keyID == "a")
				camera.MoveRight(-ParticleSimulation::GetMoveSpeed()*100*WindowManager::GetSPF());

		else if(ev.keyID == "s")
				camera.MoveForward(-ParticleSimulation::GetMoveSpeed()*100*WindowManager::GetSPF());

		else if(ev.keyID == "d")
				camera.MoveRight(ParticleSimulation::GetMoveSpeed()*100*WindowManager::GetSPF());

		else if(ev.keyID == "q")
				camera.MoveUpward(ParticleSimulation::GetMoveSpeed()*50*WindowManager::GetSPF());

		else if(ev.keyID == "e")
				camera.MoveUpward(-ParticleSimulation::GetMoveSpeed()*50*WindowManager::GetSPF());

		else if(ev.keyID == "i")
				camera.RotateAroundTarget(dataBBox.centerPoint, 0.0,ParticleSimulation::GetRotationSpeed()*0.1*WindowManager::GetSPF(),0.0);

		else if(ev.keyID == "j")
				camera.RotateAroundTarget(dataBBox.centerPoint, ParticleSimulation::GetRotationSpeed()*0.1*WindowManager::GetSPF(),0.0,0.0);

		else if(ev.keyID == "k")
				camera.RotateAroundTarget(dataBBox.centerPoint, 0.0,-ParticleSimulation::GetRotationSpeed()*0.1*WindowManager::GetSPF(),0.0);

		else if(ev.keyID == "l")
				camera.RotateAroundTarget(dataBBox.centerPoint, -ParticleSimulation::GetRotationSpeed()*0.1*WindowManager::GetSPF(),0.0,0.0);
	}

	void PP_GEOM::OnMotion(Event ev)
	{
		//account for mouse moving around screen between clicks or going off the render screen
		if( (ev.mouseX > (mouseMotionX + 10.f)) || (ev.mouseX < (mouseMotionX - 10.f)) || 
			(ev.mouseX > ParticleSimulation::GetRenderXMax()) || (ev.mouseX < ParticleSimulation::GetRenderXMin()) ||
			 (ev.mouseY > ParticleSimulation::GetRenderYMax()) || (ev.mouseY < ParticleSimulation::GetRenderYMin())  )

			mouseMotionX = ev.mouseX;
		if( (ev.mouseY > (mouseMotionY + 10.f)) || (ev.mouseY < (mouseMotionY - 10.f)) ||
		 	(ev.mouseY > ParticleSimulation::GetRenderYMax()) || (ev.mouseY < ParticleSimulation::GetRenderYMin()) ||
		 	(ev.mouseX > ParticleSimulation::GetRenderXMax()) || (ev.mouseX < ParticleSimulation::GetRenderXMin()) )
			mouseMotionY = ev.mouseY;

		float xRot = (mouseMotionX - ev.mouseX)*0.1f; //*time elapsed in frame
		float yRot = (mouseMotionY - ev.mouseY)*0.1f;

		camera.Rotate(xRot, yRot, 0.f);

		mouseMotionX = ev.mouseX;
		mouseMotionY = ev.mouseY;
	} 

	void PP_GEOM::genVBO()
	{

		// Generate Vertex Buffer Object with space reserved for data, then insert interleaved data.
		int bufferSize = (particleList.size()*sizeof(particle_sim));

		glGenBuffers(1, &VBO);	
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferData(GL_ARRAY_BUFFER, bufferSize, 0, GL_STATIC_DRAW); //pass 0 to reserve space
		glBufferSubData(GL_ARRAY_BUFFER, 0, bufferSize, &particleList[0]); //enum,offset,size,start
		glBindBuffer(GL_ARRAY_BUFFER, 0);

		// Byte interval = size of particle_sim
		byteInterval = (int)sizeof(particle_sim);
	}

}

