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
 
#include "GUICommand.h"
#include "SimpleGUI.h"

namespace previewer
{
	namespace simple_gui
	{
		GUICommand::GUICommand()
		{
			isTerminating = false;

			// Set up confirmation list
			confirmationsList.push_back("PARAMETER FILE WRITTEN TO");
			confirmationsList.push_back("PALETTE FILE LOADED FOR");
			confirmationsList.push_back("COLOURS RELOADED");
			confirmationsList.push_back("FPS");
			confirmationsList.push_back("ANIMATION POINT SAVED AT");
			confirmationsList.push_back("LAST ANIMATION POINT REMOVED");
			confirmationsList.push_back("MOVEMENT SPEED SET TO");
			confirmationsList.push_back("ROTATION SPEED SET TO");
			confirmationsList.push_back("ANIMATION PREVIEW IN PROGRESS");
			confirmationsList.push_back("ANIMATION FILE SAVED TO");
			confirmationsList.push_back("ANIMATION FILE LOADED FROM");
			confirmationsList.push_back("X RESOLUTION SET TO");
			confirmationsList.push_back("Y RESOLUTION SET TO");
			confirmationsList.push_back("RESOLUTION SET TO");
			confirmationsList.push_back("FIELD OF VIEW SET TO");
			confirmationsList.push_back("INTERPOLATED - CHECK CLI FOR MORE DETAILS");
			confirmationsList.push_back("WRITTEN ANIMATION FILE - CHECK CLI FOR MORE DETAILS");
			confirmationsList.push_back("CAMERA INTERPOLATION TYPE SET TO");
			confirmationsList.push_back("LOOKAT INTERPOLATION TYPE SET TO");
			confirmationsList.push_back("PROGRAM RUNNING, CHECK CLI FOR MORE DETAILS");
			confirmationsList.push_back("VIEWING IMAGE, TYPE STOP VIEWING TO RETURN TO PREVIEW");
			confirmationsList.push_back("VIEWING STOPPED");
			confirmationsList.push_back("PARAMETER SET");	
		}

		void GUICommand::Undo()
		{
			//
		}
		void GUICommand::Redo()
		{
			//
		}

		std::string GUICommand::GetCurrentCommandLine()
		{
			return currentCommandLine;
		}
		void GUICommand::HandleKeyEvent(Event ev, Previewer& pv)
		{
			if(ev.keyID == "#")
			{
				SimpleGUI::guiActive = !SimpleGUI::guiActive;
				currentCommandLine = "";
			}

			else if(SimpleGUI::guiActive)
			{
				if(ev.keyID == "RETURN")
				{
					ProcessCommand(pv);
				}
				else if(ev.keyID == "BACKSPACE")
				{
					currentCommandLine = currentCommandLine.substr(0, currentCommandLine.size()-1);
				}
				else if(ev.keyID == "SPACE")
				{
					currentCommandLine += " ";
				}
				else
				{
					if(ev.keyID.length() == 1)
						currentCommandLine += ev.keyID;
				}
			}
		}

		void GUICommand::ProcessCommand(Previewer& pv)
		{
			std::string command;
			command = currentCommandLine;
			std::transform(command.begin(), command.end(), command.begin(), ::toupper);

			// Quit Event
			if(DoesStringBegin(command, "QUIT"))
			{
				// Quit application
				currentCommandLine.clear();
				isTerminating = true;
				pv.TriggerOnQuitApplicationEvent();
			}

			// Write Parameter File
			if(DoesStringBegin(command, "WRITE PARAM FILE"))
			{
				if(GetArgFromString(currentCommandLine, 4) != "ARG_NOT_FOUND")
				{
					// Writing parameter file
					pv.WriteParameterFile(GetArgFromString(currentCommandLine, 4));
					currentCommandLine = "Parameter file written to " + GetArgFromString(currentCommandLine, 4);
				}
			}

			// Set Palette
			if(DoesStringBegin(command, "SET PALETTE"))
			{
				// Chaging the palette
				pv.SetPalette(GetArgFromString(currentCommandLine, 3), atoi(GetArgFromString(currentCommandLine, 4).c_str()));
				currentCommandLine = "Palette file loaded for particle type " + GetArgFromString(currentCommandLine, 4);
			}

			// Reload Colours
			if(DoesStringBegin(command, "RELOAD COLOURS"))
			{
				// Reloading the colour data
				pv.ReloadColourData();
				currentCommandLine = "Colours reloaded";
			}

			// Get FPS
			if(DoesStringBegin(command, "GET FPS"))
			{	
				std::stringstream sstream;
				sstream << pv.GetFPS();
				std::string str;
				sstream >> str;
				// Output the fps
				currentCommandLine = "FPS " + str;
			}

			// Set Animation Point
			if(DoesStringBegin(command, "SET ANIMATION POINT"))
			{
				pv.SetAnimationPoint(atoi(GetArgFromString(currentCommandLine, 4).c_str()));
				pv.UnloadAnimationFile();
				currentCommandLine = "Animation point saved at time: " + GetArgFromString(currentCommandLine, 4);
			}

			// Remove Animation Point
			if(DoesStringBegin(command, "REMOVE ANIMATION POINT"))
			{
				pv.RemoveAnimationPoint();
				currentCommandLine = "Last animation point removed";
			}

			// Set Movement Speed
			if(DoesStringBegin(command, "SET MOVE SPEED"))
			{
				pv.SetMovementSpeed(atoi(GetArgFromString(currentCommandLine, 4).c_str()));
				currentCommandLine = "Movement speed set to " + GetArgFromString(currentCommandLine, 4);
			}

			// Set Rotation Speed
			if(DoesStringBegin(command, "SET ROTATION SPEED"))
			{
				pv.SetRotationSpeed(atoi(GetArgFromString(currentCommandLine, 4).c_str()));
				currentCommandLine = "Rotation speed set to " + GetArgFromString(currentCommandLine, 4);
			}

			// Preview Animation
			if(DoesStringBegin(command, "PREVIEW ANIMATION"))
			{
				//pv.UnloadAnimationFile();
				pv.Interpolate();
				pv.PreviewAnimation();
				currentCommandLine = "Animation preview in progress";
			}

			// Save animation file 
			if(DoesStringBegin(command, "SAVE ANIMATION"))
			{
				pv.Interpolate();
				pv.SaveAnimationFile(GetArgFromString(currentCommandLine, 3));
				currentCommandLine = "Animation file saved to " + GetArgFromString(currentCommandLine, 3);
			}

			// Load animation file
			if(DoesStringBegin(command, "LOAD ANIMATION"))
			{
				pv.LoadAnimationFile(GetArgFromString(currentCommandLine, 3));
				currentCommandLine = "Animation file loaded from " + GetArgFromString(currentCommandLine, 3);
			}

			// Set X Resolution
			if(DoesStringBegin(command, "SET XRES"))
			{
				pv.SetXRes(atoi(GetArgFromString(currentCommandLine, 3).c_str()));
				currentCommandLine = "X resolution set to " + GetArgFromString(currentCommandLine, 3);
			}

			// Set Y Resolution
			if(DoesStringBegin(command, "SET YRES"))
			{
				pv.SetYRes(atoi(GetArgFromString(currentCommandLine, 3).c_str()));
				currentCommandLine = "Y resolution set to " + GetArgFromString(currentCommandLine, 3);
			}

			// Set Resolution
			if(DoesStringBegin(command, "SET RESOLUTION"))
			{
				pv.SetXRes(atoi(GetArgFromString(currentCommandLine, 3).c_str()));
				pv.SetYRes(atoi(GetArgFromString(currentCommandLine, 4).c_str()));
				currentCommandLine = "Resolution set to X:" + GetArgFromString(currentCommandLine, 3) + " Y:" + GetArgFromString(currentCommandLine, 4);
			}

			// Field of view
			if(DoesStringBegin(command, "SET FOV"))
			{
				pv.SetFOV(atoi(GetArgFromString(currentCommandLine, 3).c_str()));
				currentCommandLine = "Field of view set to " + GetArgFromString(currentCommandLine, 3);
			}

			// Field of view
			if(DoesStringBegin(command, "INTERPOLATE"))
			{
				pv.Interpolate();
				currentCommandLine = "Interpolated - check cli for more details";
			}

			// Write splotch "scene file"
			if(DoesStringBegin(command, "WRITE SPLOTCH ANIMATION"))
			{
				pv.Interpolate();
				pv.WriteSplotchAnimationFile(GetArgFromString(currentCommandLine, 4));
				currentCommandLine = "Written animation file - check cli for more details";
			}

			// Allow you to set camera interpolation type			
			if(DoesStringBegin(command, "SET CAMERA INTERPOLATION"))
			{
				if(GetArgFromString(currentCommandLine, 4) == "linear")
				{
					pv.SetCameraInterpolation(LINEAR);
					currentCommandLine = "Camera interpolation type set to: LINEAR";
				}
				else if(GetArgFromString(currentCommandLine, 4) == "cubic")
				{
					pv.SetCameraInterpolation(CUBIC);
					currentCommandLine = "Camera interpolation type set to: CUBIC";
				}
				else 
				{
					currentCommandLine = "Camera interpolation type set to: UNRECOGNISED";
				}

			}

			// Allow you to set lookat interpolation type
			if(DoesStringBegin(command, "SET LOOKAT INTERPOLATION"))
			{
				if(GetArgFromString(currentCommandLine, 4) == "linear")
				{
					pv.SetLookatInterpolation(LINEAR);
					currentCommandLine = "Lookat interpolation type set to: LINEAR";
				}
				else if(GetArgFromString(currentCommandLine, 4) == "cubic")
				{
					pv.SetLookatInterpolation(CUBIC);
					currentCommandLine = "Lookat interpolation type set to: CUBIC";
				}
				else 
				{
					currentCommandLine = "Lookat interpolation type set to: UNRECOGNISED";
				}

			}

			// Run splotch as a new process
			if(DoesStringBegin(command, "RUN"))
			{
				std::string arg = "./"+GetArgFromString(currentCommandLine, 2)+" "+GetArgFromString(currentCommandLine, 3);
				system(arg.c_str());
				currentCommandLine = "Program running, check cli for more details";
			}

			// View splotch image
			if(DoesStringBegin(command, "VIEW IMAGE"))
			{
				pv.ViewImage(GetArgFromString(currentCommandLine, 3));
				currentCommandLine = "Viewing image, type stop viewing to return to preview";
			}

			// Stop viewing splotch image
			if(DoesStringBegin(command, "STOP VIEWING"))
			{
				pv.StopViewingImage();
				currentCommandLine = "Viewing stopped";
			}

			// Set other parameters not already settable
			if(DoesStringBegin(command, "SET PARAM"))
			{
				pv.SetParameter(GetArgFromString(currentCommandLine, 3), GetArgFromString(currentCommandLine, 4));
				currentCommandLine = "Parameter set";
			}		

			// Get other parameters not already gettable
			if(DoesStringBegin(command, "GET PARAM"))
			{
				std::string param = pv.GetParameter(GetArgFromString(currentCommandLine, 3));
				currentCommandLine = GetArgFromString(currentCommandLine, 3)+" = "+param;
				std::string confirm = GetArgFromString(currentCommandLine, 3)+" = "+param;
				std::transform(confirm.begin(), confirm.end(), confirm.begin(), ::toupper);
				confirmationsList.push_back(confirm);	
			}
		


			// Handle response messages
			for(std::list<std::string>::iterator it = confirmationsList.begin(); it != confirmationsList.end(); it++)
			{
				if(DoesStringBegin(command, *it))
				{
					currentCommandLine.clear();
				}
			}
		}

		//void GUICommand::ProcessReverseCommand(std::string command)
		//{
			// Undo command feature not written yet
		//}

		bool GUICommand::IsTerminating()
		{
			return isTerminating;
		}

		std::string GUICommand::GetArgFromString(std::string str, int arg)
		{
			bool protectedMode = false; // not in " protected mode
			int currentArg = 1; // analysing first argument
			std::string currentArgStr;

			for(std::string::iterator it = str.begin(); it != str.end(); it++)
			{
				// Test for space
				if(*it == ' ')
				{
					if(protectedMode)
					{
						// Ignore the space
						currentArgStr += *it;
					}
					else
					{
						// End of argument
						if(currentArg == arg)
						{
							return currentArgStr;
						}
						else
						{
							currentArg++;
							currentArgStr.clear();
						}
					}
				}

				// Test for protected mode character
				else if( *it == '"')
				{
					if(protectedMode)
						protectedMode = false;
					else
						protectedMode = true;
				}

				// Any other character
				else
				{
					currentArgStr += *it;
				}
			}

			if(currentArg != arg)
			{
				return std::string("ARG_NOT_FOUND");
			}
			else
			{
				return currentArgStr;
			}
		}

		bool GUICommand::DoesStringBegin(std::string str1, std::string str2)
		{
			return (str1.compare(0, str2.length(), str2) == 0) ? true : false;
		}
	}
}