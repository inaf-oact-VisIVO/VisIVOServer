#version 120

uniform mat4 MVP; //ModelViewProjectionMatrix

void main()
{	

    //pass through colour and position after multiplying pos by matrices
    gl_FrontColor = vec4(gl_Color.x, gl_Color.y, gl_Color.z, 0.2);

    // clamp minimum size 
    if(gl_Normal.x < 200)
        gl_PointSize = 0.2;
    else
	  gl_PointSize = (gl_Normal.x/1000);


    gl_Position = MVP * gl_Vertex;
}


